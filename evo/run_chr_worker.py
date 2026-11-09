#!/usr/bin/env python3

import os
import subprocess
import time

import pysam

from run_chr_log_stats import (
    component_span_stats_from_stats,
    read_chromosome_log_stats,
    parse_haplotag_detected_value,
)
from run_chr_runlog import write_chromosome_edge_summary


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


def assign_unique_alignment_names(bam_path, prefix):
    """Give every alignment record a unique name for independent counting."""
    output_path = f"{bam_path}.alignment_records.bam"
    count = 0
    try:
        with pysam.AlignmentFile(bam_path, "rb") as source:
            with pysam.AlignmentFile(output_path, "wb", template=source) as output:
                for alignment in source.fetch(until_eof=True):
                    count += 1
                    alignment.query_name = f"{prefix}_aln{count}"
                    output.write(alignment)
        os.replace(output_path, bam_path)
    finally:
        if os.path.exists(output_path):
            os.unlink(output_path)
    return count


def scan_vcf_filter_pass_presence(vcf_path):
    """
    Scan all VCF entries and return:
      - has_pass: True if at least one entry has FILTER=PASS (case-insensitive)
      - total_records: number of scanned variant records
      - pass_records: number of records with PASS in FILTER
    """
    cmd = ["bcftools", "query", "-f", "%FILTER\\n", vcf_path]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    total_records = 0
    pass_records = 0
    for line in result.stdout.splitlines():
        total_records += 1
        filters = [tok.strip().lower() for tok in line.split(";") if tok.strip()]
        if "pass" in filters:
            pass_records += 1
    return pass_records > 0, total_records, pass_records


def prepare_chromosome_bams(chrom, length, pre_dir, args, chrom_log, cmd_env):
    """Extract, rename, and index one chromosome BAM per input."""
    started = time.monotonic()
    chrom_bams = []
    for bam_index, bam_path in enumerate(args.bams, start=1):
        filename = f"{chrom}.bam" if len(args.bams) == 1 else f"{chrom}.input{bam_index}.bam"
        chrom_bam = os.path.join(pre_dir, filename)
        chrom_bams.append(chrom_bam)
        # Exclude unmapped/secondary and MAPQ <20; retain supplementary alignments.
        run_cmd([
            "samtools", "view", "-b", "-@", str(args.tool_threads),
            "-F", "260", "-q", "20", bam_path, f"{chrom}:1-{length}",
            "-o", chrom_bam,
        ], log_fh=chrom_log, env=cmd_env)
        alignment_count = assign_unique_alignment_names(chrom_bam, f"input{bam_index}")
        chrom_log.write(f"alignment_records_input{bam_index}={alignment_count}\n")
        run_cmd(
            ["samtools", "index", "-@", str(args.tool_threads), chrom_bam],
            log_fh=chrom_log, env=cmd_env,
        )
    chrom_log.write(f"[stage_timing] stage=bam_preprocess seconds={time.monotonic() - started:.6f}\n")
    chrom_log.flush()
    return chrom_bams


def prepare_germline_vcf(source_path, output_path, chrom, sample_name, require_pass):
    """Select the requested sample's biallelic, phased heterozygous SNPs."""
    total = eligible = 0
    with pysam.VariantFile(source_path) as source:
        source.subset_samples([sample_name])
        with pysam.VariantFile(output_path, "wz", header=source.header) as output:
            records = source.fetch(chrom) if chrom in source.header.contigs else ()
            for record in records:
                total += 1
                if require_pass and "PASS" not in record.filter:
                    continue
                if (len(record.alleles) != 2 or
                        any(len(a) != 1 for a in record.alleles) or
                        record.alts[0] in ("*", ".", "X") or record.ref == record.alts[0]):
                    continue
                sample = record.samples[sample_name]
                genotype = sample.get("GT")
                if not genotype or genotype[:2] not in ((0, 1), (1, 0)):
                    continue
                # Match the previous reader: the first two alleles must be phased.
                # For polyploid calls, later separators do not affect eligibility.
                phased = sample.phased
                if len(genotype) > 2:
                    gt_index = list(record.format).index("GT")
                    gt_text = str(record).split("\t")[9].split(":")[gt_index]
                    phased = gt_text.startswith(("0|1", "1|0"))
                if phased:
                    output.write(record)
                    eligible += 1
    return total, eligible


def prepare_chromosome_vcfs(chrom, pre_dir, args, chrom_log, cmd_env):
    """Extract and index somatic and optional germline SNP VCFs."""
    inputs = [
        ("somatic", args.input_vcf_gz, args.require_pass_filter),
        ("germline", args.germline_vcf_gz, args.germline_require_pass_filter),
    ]
    prepared = {}
    for kind, source_vcf, require_pass in inputs:
        if not source_vcf:
            continue
        started = time.monotonic()
        suffix = "" if kind == "somatic" else ".germline"
        chrom_vcf = os.path.join(pre_dir, f"{chrom}{suffix}.vcf.gz")
        filter_label = "vcf" if kind == "somatic" else "germline_vcf"
        chrom_log.write(f"{filter_label}_require_filter_PASS={int(require_pass)}\n")
        if kind == "germline":
            total, eligible = prepare_germline_vcf(
                source_vcf, chrom_vcf, chrom, args.germline_vcf_sample_name, require_pass
            )
            chrom_log.write(f"[germline_markers] input_records={total} eligible_records={eligible}\n")
        else:
            command = ["bcftools", "view", "--threads", str(args.tool_threads), "-Oz", "-r", chrom]
            if require_pass:
                command.extend(["-f", "PASS"])
            command.extend(["-v", "snps", source_vcf, "-o", chrom_vcf])
            run_cmd(command, log_fh=chrom_log, env=cmd_env)
        run_cmd(
            ["bcftools", "index", "--threads", str(args.tool_threads), chrom_vcf],
            log_fh=chrom_log, env=cmd_env,
        )
        if kind == "somatic":
            plain_vcf = os.path.join(pre_dir, f"{chrom}.vcf")
            run_cmd(
                ["bcftools", "view", "--threads", str(args.tool_threads), "-Ov", chrom_vcf, "-o", plain_vcf],
                log_fh=chrom_log, env=cmd_env,
            )
        chrom_log.write(
            f"[stage_timing] stage={kind}_vcf_preprocess seconds={time.monotonic() - started:.6f}\n"
        )
        prepared[kind] = chrom_vcf
    chrom_log.flush()
    return prepared["somatic"], prepared.get("germline")


def run_chromosome_main(chrom_bams, chrom_vcfs, out_dir, args, chrom_log, cmd_env):
    """Run inference and graph construction, then remove temporary BAMs on success."""
    somatic_vcf, germline_vcf = chrom_vcfs
    command = [args.main_program, "--somatic-vcf", somatic_vcf]
    for bam in chrom_bams:
        command.extend(["--bam", bam])
    command.extend([
        "--output-dir", out_dir,
        "--window-size", str(args.window_size),
        "--window-stride", str(args.window_stride),
    ])
    if germline_vcf:
        command.extend(["--germline-vcf", germline_vcf])
    for bed_file in args.exclude_regions_bed or []:
        command.extend(["--exclude-regions-bed", bed_file])
    if args.divergent_same_hp:
        command.append("--divergent-same-hp")

    started = time.monotonic()
    run_cmd(command, log_fh=chrom_log, env=cmd_env)
    chrom_log.write(f"[stage_timing] stage=main seconds={time.monotonic() - started:.6f}\n")
    for bam in chrom_bams:
        for temporary_file in (bam, bam + ".bai"):
            if os.path.exists(temporary_file):
                os.remove(temporary_file)
                chrom_log.write(f"removed_temporary_file={temporary_file}\n")
    chrom_log.flush()


def collect_chromosome_results(chrom, length, out_dir, chrom_log, started):
    """Write the chromosome summary and return statistics to the genome runner."""
    chrom_log.flush()
    chrom_log_path = chrom_log.name
    mutation_stats, edge_stats, haplotag_status = read_chromosome_log_stats(chrom_log_path)
    haplotag_detected = parse_haplotag_detected_value(haplotag_status)
    write_chromosome_edge_summary(chrom_log, haplotag_detected, edge_stats)
    chrom_log.flush()

    stats_path = os.path.join(out_dir, "component_statistics.txt")
    summed, unique = component_span_stats_from_stats(stats_path)
    pct_summed = 100.0 * summed / length if length > 0 else 0.0
    pct_unique = 100.0 * unique / length if length > 0 else 0.0
    elapsed_seconds = time.monotonic() - started
    print(
        f"[Connected components] {chrom}: "
        f"summed={summed} bp ({pct_summed:.4f}%), "
        f"unique={unique} bp ({pct_unique:.4f}%)"
    )
    print(f"=== Done {chrom} ({elapsed_seconds:.2f}s) ===")
    return {
        "chrom": chrom,
        "length": length,
        "chrom_log_path": chrom_log_path,
        "out_dir": out_dir,
        "component_bp_summed": summed,
        "component_bp_unique": unique,
        "chrom_pct_summed": pct_summed,
        "chrom_pct_unique": pct_unique,
        "mutation_stats": mutation_stats,
        "edge_haplotype_stats": edge_stats,
        "haplotag_status": haplotag_status,
        "elapsed_seconds": elapsed_seconds,
    }


def process_chromosome(chrom, length, args):
    """Coordinate the four steps for one chromosome."""
    pre_dir = os.path.join(args.output_dir, f"{chrom}_pre")
    out_dir = os.path.join(args.output_dir, f"{chrom}_out")
    chrom_log_path = os.path.join(args.log_dir, f"{chrom}.log")
    os.makedirs(pre_dir, exist_ok=True)
    os.makedirs(out_dir, exist_ok=True)
    started = time.monotonic()
    cmd_env = os.environ.copy()
    if args.force_single_thread_tools:
        cmd_env["OMP_NUM_THREADS"] = "1"

    print(f"\n=== Processing {chrom} ({chrom}:1-{length}) ===")
    print(f"Preprocessing directory: {pre_dir}")
    print(f"Main program output directory: {out_dir}")
    print(f"Chromosome log: {chrom_log_path}")
    with open(chrom_log_path, "w") as chrom_log:
        chrom_log.write(
            f"Chromosome: {chrom}\n"
            f"Region: {chrom}:1-{length}\n"
            f"Preprocessing directory: {pre_dir}\n"
            f"Main output directory: {out_dir}\n"
            f"tool_threads_per_subprocess={args.tool_threads}\n"
            f"input_bam_count={len(args.bams)}\n"
            f"force_single_thread_tools={int(args.force_single_thread_tools)}\n"
            f"haplotag_mode={'germline_vcf' if args.germline_vcf_gz else 'none'}\n"
            f"divergent_same_hp={int(args.divergent_same_hp)}\n"
            f"maximum_mutations={args.window_size}\n"
            f"group_stride={args.window_stride}\n"
        )
        for bed_file in args.exclude_regions_bed or []:
            chrom_log.write(f"exclude_regions_bed={bed_file}\n")
        chrom_log.flush()

        prepared_bams = prepare_chromosome_bams(chrom, length, pre_dir, args, chrom_log, cmd_env)
        chrom_vcfs = prepare_chromosome_vcfs(chrom, pre_dir, args, chrom_log, cmd_env)
        run_chromosome_main(prepared_bams, chrom_vcfs, out_dir, args, chrom_log, cmd_env)
        return collect_chromosome_results(chrom, length, out_dir, chrom_log, started)
