#!/usr/bin/env python3
import argparse
import os
import re
import sys
import shlex
import subprocess
import csv
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime


def run_cmd(cmd, log_fh=None, env=None):
    """Run a shell command, writing full output to logs and a concise terminal summary."""
    cmd_str = " ".join(cmd)
    print(f"\n[Running] {cmd_str}")
    if log_fh is not None:
        log_fh.write(f"\n$ {cmd_str}\n")
        log_fh.flush()

    process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        env=env,
    )
    line_count = 0
    for line in process.stdout:
        line_count += 1
        if log_fh is not None:
            log_fh.write(line)
    process.wait()
    if log_fh is not None:
        log_fh.flush()

    if process.returncode != 0:
        print(f"[Failed] {cmd[0]} exited with code {process.returncode}. See logs for details.")
        raise RuntimeError(f"Command failed ({process.returncode}): {' '.join(cmd)}")
    print(f"[Done] {cmd[0]} ({line_count} log lines)")


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Run evo pipeline per chromosome, aggregate genome-level stats, "
            "and write run/chromosome logs."
        )
    )
    parser.add_argument("input_vcf_gz", help="Input VCF.gz")
    parser.add_argument("haplotagged_bam", help="Input haplotagged BAM")
    parser.add_argument("output_dir", help="Output directory")
    parser.add_argument("vcf_sample_name", nargs="?", default="", help="Optional VCF sample name")
    parser.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="Number of chromosomes to process in parallel (default: 1)",
    )
    args = parser.parse_args()
    if args.jobs <= 0:
        parser.error("--jobs must be >= 1")
    return args


def get_chromosomes_from_bam(bam_path):
    """Return (chrom, length) pairs for standard chromosomes only (chr or no-chr)."""
    cmd = ["samtools", "idxstats", bam_path]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)

    chromosomes = []
    for line in result.stdout.strip().split("\n"):
        fields = line.split("\t")
        if len(fields) < 2:
            continue

        chrom = fields[0]
        length = int(fields[1])

        # Normalize: drop optional "chr" prefix for checking
        c = chrom[3:] if chrom.startswith("chr") else chrom

        # Keep only 1-22, X, Y, M, MT
        if re.fullmatch(r"(?:[1-9]|1[0-9]|2[0-2]|X|Y|M|MT)", c) is None:
            print(f"Skipping non-standard contig {chrom}")
            continue

        chromosomes.append((chrom, length))

    return chromosomes


def merged_interval_coverage_bp(intervals):
    """Return covered bp after merging inclusive [start, end] intervals."""
    if not intervals:
        return 0

    intervals = sorted(intervals)
    merged_start, merged_end = intervals[0]
    covered = 0

    for start, end in intervals[1:]:
        if start <= merged_end + 1:
            merged_end = max(merged_end, end)
        else:
            covered += merged_end - merged_start + 1
            merged_start, merged_end = start, end

    covered += merged_end - merged_start + 1
    return covered


def component_span_stats_from_stats(component_stats_path):
    """
    Compute two connected-component span metrics from component_statistics.txt:
      1) summed_component_bp: sum of all component spans (counts overlaps multiple times)
      2) unique_component_bp: covered bp after merging overlapping component intervals
         (counts overlaps once)
    """
    if not os.path.exists(component_stats_path):
        return 0, 0

    intervals = []
    summed_component_bp = 0
    with open(component_stats_path, "r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            try:
                start = int(row["min_pos"])
                end = int(row["max_pos"])
            except (KeyError, TypeError, ValueError):
                continue
            if end < start:
                continue
            intervals.append((start, end))
            # Inclusive genomic coordinates: [start, end]
            summed_component_bp += (end - start + 1)

    unique_component_bp = merged_interval_coverage_bp(intervals)
    return summed_component_bp, unique_component_bp



def read_mutation_stats_tsv(tsv_path):
    """Read mutation_stats.tsv and return a dict of metric -> count."""
    stats = {}
    if not os.path.exists(tsv_path):
        return stats
    with open(tsv_path, "r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            stats[row["metric"]] = int(row["count"])
    # Backward compatibility with older metric name.
    if "orphaned" not in stats and "rejected" in stats:
        stats["orphaned"] = stats["rejected"]
    return stats


def extract_max_span_scan_from_chrom_log(chrom_log_path):
    """
    Parse a '[max_span_scan]' line written by ./main from a chromosome log.
    Returns a dict of parsed key/value fields plus the raw line.
    """
    if not os.path.exists(chrom_log_path):
        return {}
    last_line = ""
    with open(chrom_log_path, "r") as f:
        for line in f:
            if line.startswith("[max_span_scan]"):
                last_line = line.strip()
    if not last_line:
        return {}

    out = {"raw": last_line}
    for tok in last_line.split():
        if "=" not in tok:
            continue
        k, v = tok.split("=", 1)
        out[k] = v
    return out


def process_chromosome(
    chrom,
    length,
    bam_path,
    vcf_path,
    outdir,
    log_dir,
    main_program,
    vcf_sample_name,
    tool_threads,
    force_single_thread_tools,
):
    region = f"{chrom}:1-{length}"
    chrom_log_path = os.path.join(log_dir, f"{chrom}.log")
    pre_dir = os.path.join(outdir, f"{chrom}_pre")
    out_dir = os.path.join(outdir, f"{chrom}_out")
    os.makedirs(pre_dir, exist_ok=True)
    os.makedirs(out_dir, exist_ok=True)

    chrom_t0 = time.monotonic()
    cmd_env = os.environ.copy()
    if force_single_thread_tools:
        cmd_env["OMP_NUM_THREADS"] = "1"

    print(f"\n=== Processing {chrom} ({region}) ===")
    print(f"Preprocessing directory: {pre_dir}")
    print(f"Main program output directory: {out_dir}")
    print(f"Chromosome log: {chrom_log_path}")

    with open(chrom_log_path, "w") as chrom_log:
        chrom_log.write(f"Chromosome: {chrom}\n")
        chrom_log.write(f"Region: {region}\n")
        chrom_log.write(f"Preprocessing directory: {pre_dir}\n")
        chrom_log.write(f"Main output directory: {out_dir}\n")
        chrom_log.write(f"tool_threads_per_subprocess={tool_threads}\n")
        chrom_log.write(f"force_single_thread_tools={int(force_single_thread_tools)}\n")
        chrom_log.flush()

        # --- Step 1: Split BAM for this chromosome with filters ---
        chrom_bam = os.path.join(pre_dir, f"{chrom}.bam")
        # Exclude unmapped, secondary, supplementary reads; MAPQ < 20
        run_cmd([
            "samtools", "view", "-b", "-@", str(tool_threads), "-F", "2308", "-q", "20",
            bam_path, region, "-o", chrom_bam
        ], log_fh=chrom_log, env=cmd_env)
        run_cmd(["samtools", "index", "-@", str(tool_threads), chrom_bam], log_fh=chrom_log, env=cmd_env)

        # --- Step 2: Split VCF for this chromosome with PASS + SNP-only filter ---
        chrom_vcf = os.path.join(pre_dir, f"{chrom}.vcf.gz")
        chrom_vcf_plain = os.path.join(pre_dir, f"{chrom}.vcf")
        run_cmd([
            "bcftools", "view", "--threads", str(tool_threads), "-Oz", "-r", chrom, "-f", "PASS", "-v", "snps",
            vcf_path, "-o", chrom_vcf
        ], log_fh=chrom_log, env=cmd_env)
        run_cmd(["bcftools", "index", "--threads", str(tool_threads), chrom_vcf], log_fh=chrom_log, env=cmd_env)
        run_cmd(
            ["bcftools", "view", "--threads", str(tool_threads), "-Ov", chrom_vcf, "-o", chrom_vcf_plain],
            log_fh=chrom_log,
            env=cmd_env,
        )

        # --- Step 3: Run main program ---
        main_cmd = [
            main_program,
            "--vcf", chrom_vcf,
            "--bam", chrom_bam,
            "--output-dir", out_dir,
        ]
        if vcf_sample_name:
            main_cmd.extend(["--vcf-sample-name", vcf_sample_name])
        run_cmd(main_cmd, log_fh=chrom_log, env=cmd_env)

    component_stats_path = os.path.join(out_dir, "component_statistics.txt")
    chrom_component_bp_summed, chrom_component_bp_unique = component_span_stats_from_stats(
        component_stats_path
    )
    mutation_stats_tsv = os.path.join(out_dir, "mutation_stats.tsv")
    chrom_mstats = read_mutation_stats_tsv(mutation_stats_tsv)
    chrom_pct_summed = (100.0 * chrom_component_bp_summed / length) if length > 0 else 0.0
    chrom_pct_unique = (100.0 * chrom_component_bp_unique / length) if length > 0 else 0.0
    elapsed_seconds = time.monotonic() - chrom_t0

    print(
        f"[Connected components] {chrom}: "
        f"summed={chrom_component_bp_summed} bp ({chrom_pct_summed:.4f}%), "
        f"unique={chrom_component_bp_unique} bp ({chrom_pct_unique:.4f}%)"
    )
    print(f"=== Done {chrom} ({elapsed_seconds:.2f}s) ===")

    return {
        "chrom": chrom,
        "length": length,
        "chrom_log_path": chrom_log_path,
        "out_dir": out_dir,
        "component_bp_summed": chrom_component_bp_summed,
        "component_bp_unique": chrom_component_bp_unique,
        "chrom_pct_summed": chrom_pct_summed,
        "chrom_pct_unique": chrom_pct_unique,
        "mutation_stats": chrom_mstats,
        "elapsed_seconds": elapsed_seconds,
    }


def main():
    args = parse_args()
    vcf_path = args.input_vcf_gz
    bam_path = args.haplotagged_bam
    outdir = os.path.abspath(args.output_dir)
    vcf_sample_name = args.vcf_sample_name
    jobs = args.jobs

    main_program = "./main"

    # Validate inputs
    for f in [vcf_path, bam_path, main_program]:
        if not os.path.exists(f):
            sys.exit(f"Error: required file not found: {f}")

    os.makedirs(outdir, exist_ok=True)
    log_dir = os.path.join(outdir, "logs")
    os.makedirs(log_dir, exist_ok=True)
    run_stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_log_path = os.path.join(log_dir, f"run_{run_stamp}.log")

    # Derive genome name
    genome_name = os.path.basename(bam_path).replace(".bam", "")
    print(f"Genome name: {genome_name}")
    print(f"Logs directory: {log_dir}")
    print(f"Run log: {run_log_path}")

    # Get chromosomes
    print("Extracting chromosome list from BAM header using samtools...")
    chromosomes = get_chromosomes_from_bam(bam_path)
    print(f"Detected {len(chromosomes)} chromosomes.")
    genome_bp = sum(length for _, length in chromosomes)
    total_component_bp_summed = 0
    total_component_bp_unique = 0
    genome_mutation_stats = {
        "total_snvs": 0,
        "in_relations": 0,
        "in_graph": 0,
        "with_timing": 0,
        "orphaned": 0,
        "singleton": 0,
        "accepted_edges": 0,
        "timing_edges": 0,
    }

    run_start_wall = datetime.now()
    run_start_monotonic = time.monotonic()
    invocation_cmd = " ".join(shlex.quote(x) for x in [sys.executable, *sys.argv])
    tool_threads = 1 if jobs > 1 else 4
    force_single_thread_tools = jobs > 1

    with open(run_log_path, "w") as run_log:
        run_log.write(
            f"Run started (local): {run_start_wall.strftime('%Y-%m-%d %H:%M:%S')}\n"
        )
        run_log.write(f"Run started: {run_start_wall.isoformat(timespec='seconds')}\n")
        run_log.write(f"Invocation command: {invocation_cmd}\n")
        run_log.write(f"Genome: {genome_name}\n")
        run_log.write(f"VCF: {vcf_path}\n")
        run_log.write(f"BAM: {bam_path}\n")
        run_log.write(f"Output: {outdir}\n")
        if vcf_sample_name:
            run_log.write(f"VCF sample name: {vcf_sample_name}\n")
        run_log.write("\n=== Default/auto-selected parameters ===\n")
        run_log.write(f"main_default_min_reads=2\n")
        run_log.write("max_pair_distance_source=computed_by_main_per_chrom_bam\n")
        run_log.write(f"jobs={jobs}\n")
        run_log.write(f"tool_threads_per_subprocess={tool_threads}\n")
        run_log.write("parallel_mode=chromosome_thread_pool\n")
        run_log.write(f"Detected chromosomes: {len(chromosomes)}\n")
        run_log.write(f"Genome bp (selected chromosomes): {genome_bp}\n")
        run_log.write(
            "per_chrom_max_span_summary_source=parsed_from_chrom_log_[max_span_scan]_line\n"
        )
        run_log.flush()

        chrom_tasks = sorted(chromosomes, key=lambda x: x[1], reverse=True)
        if jobs == 1:
            for chrom, length in chrom_tasks:
                result = process_chromosome(
                    chrom=chrom,
                    length=length,
                    bam_path=bam_path,
                    vcf_path=vcf_path,
                    outdir=outdir,
                    log_dir=log_dir,
                    main_program=main_program,
                    vcf_sample_name=vcf_sample_name,
                    tool_threads=tool_threads,
                    force_single_thread_tools=force_single_thread_tools,
                )
                total_component_bp_summed += result["component_bp_summed"]
                total_component_bp_unique += result["component_bp_unique"]
                for key in genome_mutation_stats:
                    genome_mutation_stats[key] += result["mutation_stats"].get(key, 0)
                run_log.write(
                    f"{result['chrom']}\tOK\t{result['chrom_log_path']}\t"
                    f"elapsed_seconds={result['elapsed_seconds']:.3f}\n"
                )
                run_log.write(
                    f"{result['chrom']}\tconnected_component_bp_summed={result['component_bp_summed']}\t"
                    f"connected_component_bp_unique={result['component_bp_unique']}\t"
                    f"chrom_bp={result['length']}\t"
                    f"chrom_pct_summed={result['chrom_pct_summed']:.6f}\t"
                    f"chrom_pct_unique={result['chrom_pct_unique']:.6f}\n"
                )
                max_span = extract_max_span_scan_from_chrom_log(result["chrom_log_path"])
                if max_span:
                    run_log.write(f"{result['chrom']}\t{max_span['raw']}\n")
                else:
                    run_log.write(f"{result['chrom']}\tmax_span_scan=not_found_in_chrom_log\n")
                run_log.flush()
        else:
            with ThreadPoolExecutor(max_workers=jobs) as executor:
                futures = {
                    executor.submit(
                        process_chromosome,
                        chrom,
                        length,
                        bam_path,
                        vcf_path,
                        outdir,
                        log_dir,
                        main_program,
                        vcf_sample_name,
                        tool_threads,
                        force_single_thread_tools,
                    ): (chrom, length)
                    for chrom, length in chrom_tasks
                }
                try:
                    for future in as_completed(futures):
                        result = future.result()
                        total_component_bp_summed += result["component_bp_summed"]
                        total_component_bp_unique += result["component_bp_unique"]
                        for key in genome_mutation_stats:
                            genome_mutation_stats[key] += result["mutation_stats"].get(key, 0)
                        run_log.write(
                            f"{result['chrom']}\tOK\t{result['chrom_log_path']}\t"
                            f"elapsed_seconds={result['elapsed_seconds']:.3f}\n"
                        )
                        run_log.write(
                            f"{result['chrom']}\tconnected_component_bp_summed={result['component_bp_summed']}\t"
                            f"connected_component_bp_unique={result['component_bp_unique']}\t"
                            f"chrom_bp={result['length']}\t"
                            f"chrom_pct_summed={result['chrom_pct_summed']:.6f}\t"
                            f"chrom_pct_unique={result['chrom_pct_unique']:.6f}\n"
                        )
                        max_span = extract_max_span_scan_from_chrom_log(result["chrom_log_path"])
                        if max_span:
                            run_log.write(f"{result['chrom']}\t{max_span['raw']}\n")
                        else:
                            run_log.write(f"{result['chrom']}\tmax_span_scan=not_found_in_chrom_log\n")
                        run_log.flush()
                except Exception as exc:
                    for f in futures:
                        f.cancel()
                    run_log.write(f"ERROR\tparallel_chromosome_processing_failed\t{exc}\n")
                    run_log.flush()
                    raise

        genome_pct_summed = (100.0 * total_component_bp_summed / genome_bp) if genome_bp > 0 else 0.0
        genome_pct_unique = (100.0 * total_component_bp_unique / genome_bp) if genome_bp > 0 else 0.0
        print(
            f"\nConnected component total span (summed): {total_component_bp_summed} bp "
            f"({genome_pct_summed:.6f}% of selected genome)"
        )
        print(
            f"Connected component total span (unique): {total_component_bp_unique} bp "
            f"({genome_pct_unique:.6f}% of selected genome)"
        )
        run_log.write(
            f"TOTAL\tconnected_component_bp_summed={total_component_bp_summed}\t"
            f"connected_component_bp_unique={total_component_bp_unique}\t"
            f"genome_bp={genome_bp}\t"
            f"genome_pct_summed={genome_pct_summed:.6f}\t"
            f"genome_pct_unique={genome_pct_unique:.6f}\n"
        )
        run_log.flush()

        # ---- Genome-wide mutation connectivity statistics ----
        g = genome_mutation_stats
        def _pct(num, denom):
            return f"{100.0 * num / denom:.1f}%" if denom > 0 else "N/A"

        mut_lines = [
            "\n=== Genome-wide mutation connectivity statistics ===",
        ]
        if g["total_snvs"] > 0:
            mut_lines.append(f"Total mutations (filtered VCF):        {g['total_snvs']}")
        mut_lines.append(
            f"  In pairwise relations:               {g['in_relations']}"
            + (f"  ({_pct(g['in_relations'], g['total_snvs'])} of total)" if g["total_snvs"] > 0 else "")
        )
        mut_lines.append(
            f"    In graph (accepted edges):          {g['in_graph']}"
            + (f"  ({_pct(g['in_graph'], g['in_relations'])} of connected)" if g["in_relations"] > 0 else "")
        )
        mut_lines.append(
            f"      - with timing edge:              {g['with_timing']}"
            + (f"  ({_pct(g['with_timing'], g['in_graph'])} of in-graph)" if g["in_graph"] > 0 else "")
        )
        mut_lines.append(
            f"    Orphaned (all edges rejected):      {g['orphaned']}"
            + (f"  ({_pct(g['orphaned'], g['in_relations'])} of connected)" if g["in_relations"] > 0 else "")
        )
        if g["total_snvs"] > 0:
            mut_lines.append(
                f"  Singleton (no pairwise relations):    {g['singleton']}"
                f"  ({_pct(g['singleton'], g['total_snvs'])} of total)"
            )
        mut_lines.append(f"Total accepted edges:                  {g['accepted_edges']}")
        mut_lines.append(f"  - timing edges:                      {g['timing_edges']}")

        mut_summary = "\n".join(mut_lines)
        print(mut_summary)
        run_log.write(mut_summary + "\n")

        run_end_wall = datetime.now()
        elapsed_seconds = time.monotonic() - run_start_monotonic
        elapsed_h = int(elapsed_seconds // 3600)
        elapsed_m = int((elapsed_seconds % 3600) // 60)
        elapsed_s = elapsed_seconds % 60
        elapsed_hms = f"{elapsed_h:02d}:{elapsed_m:02d}:{elapsed_s:06.3f}"
        runtime_line = (
            f"TOTAL_RUNTIME\tseconds={elapsed_seconds:.3f}\thms={elapsed_hms}\t"
            f"started={run_start_wall.isoformat(timespec='seconds')}\t"
            f"ended={run_end_wall.isoformat(timespec='seconds')}\n"
        )
        run_log.write(runtime_line)
        run_log.flush()

    print(f"\nAll chromosomes processed for {genome_name}. Results saved in {outdir}")
    print(f"Total runtime: {elapsed_hms} ({elapsed_seconds:.3f}s)")

if __name__ == "__main__":
    main()
