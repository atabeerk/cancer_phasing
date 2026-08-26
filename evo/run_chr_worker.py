#!/usr/bin/env python3

import os
import re
import subprocess
import time

from run_chr_log_stats import (
    RELATIONS,
    component_span_stats_from_stats,
    extract_edge_haplotype_summary_from_chrom_log,
    extract_haplotag_status_from_chrom_log,
    extract_mutation_stats_from_chrom_log,
    parse_haplotag_detected_value,
)


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


def get_vcf_samples(vcf_path):
    """Return sample names from a VCF header in declared order."""
    result = subprocess.run(
        ["bcftools", "query", "-l", vcf_path],
        capture_output=True,
        text=True,
        check=True,
    )
    return [line.strip() for line in result.stdout.splitlines() if line.strip()]


def process_chromosome(
    chrom,
    length,
    bam_paths,
    vcf_path,
    germline_vcf_path,
    outdir,
    log_dir,
    main_program,
    vcf_sample_name,
    germline_vcf_sample_name,
    exclude_regions_bed,
    require_pass_filter,
    germline_require_pass_filter,
    divergent_same_hp,
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
        chrom_log.write(f"input_bam_count={len(bam_paths)}\n")
        chrom_log.write(f"force_single_thread_tools={int(force_single_thread_tools)}\n")
        chrom_log.write(
            "haplotag_mode="
            + ("germline_vcf" if germline_vcf_path else "bam_hp_tags")
            + "\n"
        )
        chrom_log.write(
            f"divergent_same_hp={int(divergent_same_hp)}\n"
        )
        if exclude_regions_bed:
            chrom_log.write(
                f"exclude_regions_bed={exclude_regions_bed}\n"
            )
        chrom_log.flush()

        # --- Step 1: Split BAM for this chromosome with filters ---
        bam_preprocess_t0 = time.monotonic()
        chrom_bams = []
        for bam_index, bam_path in enumerate(bam_paths, start=1):
            if len(bam_paths) == 1:
                chrom_bam = os.path.join(pre_dir, f"{chrom}.bam")
            else:
                chrom_bam = os.path.join(
                    pre_dir, f"{chrom}.input{bam_index}.bam"
                )
            chrom_bams.append(chrom_bam)
            # Exclude unmapped, secondary, supplementary reads; MAPQ < 20
            run_cmd([
                "samtools", "view", "-b", "-@", str(tool_threads),
                "-F", "2308", "-q", "20", bam_path, region,
                "-o", chrom_bam,
            ], log_fh=chrom_log, env=cmd_env)
            run_cmd(
                ["samtools", "index", "-@", str(tool_threads), chrom_bam],
                log_fh=chrom_log,
                env=cmd_env,
            )
        chrom_log.write(
            f"[stage_timing] stage=bam_preprocess seconds="
            f"{time.monotonic() - bam_preprocess_t0:.6f}\n"
        )
        chrom_log.flush()

        # --- Step 2: Split VCF for this chromosome with SNP-only filter.
        # Require FILTER=PASS only when PASS is defined in VCF header metadata.
        chrom_vcf = os.path.join(pre_dir, f"{chrom}.vcf.gz")
        chrom_vcf_plain = os.path.join(pre_dir, f"{chrom}.vcf")
        somatic_vcf_t0 = time.monotonic()
        chrom_log.write(f"vcf_require_filter_PASS={int(require_pass_filter)}\n")
        view_cmd = [
            "bcftools",
            "view",
            "--threads",
            str(tool_threads),
            "-Oz",
            "-r",
            chrom,
        ]
        if require_pass_filter:
            view_cmd.extend(["-f", "PASS"])
        view_cmd.extend(["-v", "snps", vcf_path, "-o", chrom_vcf])
        run_cmd(view_cmd, log_fh=chrom_log, env=cmd_env)
        run_cmd(["bcftools", "index", "--threads", str(tool_threads), chrom_vcf], log_fh=chrom_log, env=cmd_env)
        run_cmd(
            ["bcftools", "view", "--threads", str(tool_threads), "-Ov", chrom_vcf, "-o", chrom_vcf_plain],
            log_fh=chrom_log,
            env=cmd_env,
        )
        chrom_log.write(
            f"[stage_timing] stage=somatic_vcf_preprocess seconds="
            f"{time.monotonic() - somatic_vcf_t0:.6f}\n"
        )

        chrom_germline_vcf = None
        if germline_vcf_path:
            germline_vcf_t0 = time.monotonic()
            chrom_germline_vcf = os.path.join(
                pre_dir, f"{chrom}.germline.vcf.gz"
            )
            chrom_log.write(
                "germline_vcf_require_filter_PASS="
                f"{int(germline_require_pass_filter)}\n"
            )
            germline_view_cmd = [
                "bcftools",
                "view",
                "--threads",
                str(tool_threads),
                "-Oz",
                "-r",
                chrom,
            ]
            if germline_require_pass_filter:
                germline_view_cmd.extend(["-f", "PASS"])
            germline_view_cmd.extend([
                "-v",
                "snps",
                germline_vcf_path,
                "-o",
                chrom_germline_vcf,
            ])
            run_cmd(germline_view_cmd, log_fh=chrom_log, env=cmd_env)
            run_cmd(
                [
                    "bcftools",
                    "index",
                    "--threads",
                    str(tool_threads),
                    chrom_germline_vcf,
                ],
                log_fh=chrom_log,
                env=cmd_env,
            )
            chrom_log.write(
                f"[stage_timing] stage=germline_vcf_preprocess seconds="
                f"{time.monotonic() - germline_vcf_t0:.6f}\n"
            )
        chrom_log.flush()

        # --- Step 3: Run main program ---
        main_cmd = [main_program, "--somatic-vcf", chrom_vcf]
        for chrom_bam in chrom_bams:
            main_cmd.extend(["--bam", chrom_bam])
        main_cmd.extend(["--output-dir", out_dir])
        if vcf_sample_name:
            main_cmd.extend(["--somatic-sample-name", vcf_sample_name])
        if chrom_germline_vcf:
            main_cmd.extend(["--germline-vcf", chrom_germline_vcf])
            main_cmd.extend([
                "--germline-vcf-sample-name",
                germline_vcf_sample_name,
            ])
        if exclude_regions_bed:
            main_cmd.extend([
                "--exclude-regions-bed",
                exclude_regions_bed,
            ])
        if divergent_same_hp:
            main_cmd.append("--divergent-same-hp")

        main_t0 = time.monotonic()
        run_cmd(main_cmd, log_fh=chrom_log, env=cmd_env)
        chrom_log.write(
            f"[stage_timing] stage=main seconds="
            f"{time.monotonic() - main_t0:.6f}\n"
        )

        # The chromosome BAM is no longer needed after ./main succeeds.
        for chrom_bam in chrom_bams:
            for temporary_bam_file in (chrom_bam, chrom_bam + ".bai"):
                if os.path.exists(temporary_bam_file):
                    os.remove(temporary_bam_file)
                    chrom_log.write(
                        f"removed_temporary_file={temporary_bam_file}\n"
                    )
        chrom_log.flush()

        # --- Step 4: Chromosome-wide edge summary ---
        chrom_log.flush()
        haplotag_status = extract_haplotag_status_from_chrom_log(chrom_log_path)
        haplotag_detected = parse_haplotag_detected_value(haplotag_status)
        chrom_edge_hap = extract_edge_haplotype_summary_from_chrom_log(chrom_log_path)
        chrom_mstats = extract_mutation_stats_from_chrom_log(chrom_log_path)

        def _pct(num, denom):
            return f"{100.0 * num / denom:.1f}%" if denom > 0 else "N/A"

        relation_label = {
            "cooccurring": "Co-occurring",
            "timing": "Timing",
            "divergent": "Divergent",
        }
        chrom_log.write("\n=== Chromosome-wide edge summary ===\n")
        if haplotag_detected in (0, 1):
            chrom_log.write(f"Haplotagging detected for this chromosome: {haplotag_detected}\n")
        else:
            chrom_log.write("Haplotagging detected for this chromosome: unknown\n")
        if haplotag_detected == 0:
            chrom_log.write(
                "Haplotype-splitting entries omitted: haplotagging was not detected for this chromosome.\n"
            )

        for relation in RELATIONS:
            rec = chrom_edge_hap[relation]
            chrom_log.write(f"\n{relation_label[relation]} edges (accepted):        {rec['total']}\n")
            if haplotag_detected != 0:
                chrom_log.write(
                    f"  Same haplotype:                    {rec['same_haplotype']}"
                    + (f"  ({_pct(rec['same_haplotype'], rec['total'])} of relation)" if rec["total"] > 0 else "")
                    + "\n"
                )
                chrom_log.write(
                    f"    - HP1-HP1:                       {rec['hp1_hp1']}"
                    + (f"  ({_pct(rec['hp1_hp1'], rec['same_haplotype'])} of same-haplotype)" if rec["same_haplotype"] > 0 else "")
                    + "\n"
                )
                chrom_log.write(
                    f"    - HP2-HP2:                       {rec['hp2_hp2']}"
                    + (f"  ({_pct(rec['hp2_hp2'], rec['same_haplotype'])} of same-haplotype)" if rec["same_haplotype"] > 0 else "")
                    + "\n"
                )
                chrom_log.write(
                    f"  Different haplotypes (HP1-HP2):    {rec['different_haplotype']}"
                    + (f"  ({_pct(rec['different_haplotype'], rec['total'])} of relation)" if rec["total"] > 0 else "")
                    + "\n"
                )
                chrom_log.write(
                    f"  At least one unknown or mixed:     {rec['at_least_one_unknown']}"
                    + (f"  ({_pct(rec['at_least_one_unknown'], rec['total'])} of relation)" if rec["total"] > 0 else "")
                    + "\n"
                )
        chrom_log.flush()

    component_stats_path = os.path.join(out_dir, "component_statistics.txt")
    chrom_component_bp_summed, chrom_component_bp_unique = component_span_stats_from_stats(
        component_stats_path
    )
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
        "edge_haplotype_stats": chrom_edge_hap,
        "elapsed_seconds": elapsed_seconds,
    }
