#!/usr/bin/env python3
import os
import sys
import shlex
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime

from run_chr_cli import (
    chrom_is_excluded,
    parse_args,
    parse_exclude_chroms,
    preflight_requirements,
)
from run_chr_log_stats import (
    add_edge_haplotype_summary,
    new_edge_haplotype_summary,
)
from run_chr_runlog import (
    append_chromosome_result,
    build_mutation_summary,
    format_elapsed_hms,
    write_genome_edge_summary,
    write_parallel_failure,
    write_postprocessing_header,
    write_run_header,
    write_runtime_line,
    write_total_component_line,
)
from run_chr_worker import (
    get_chromosomes_from_bam,
    get_vcf_samples,
    process_chromosome,
    run_cmd,
    scan_vcf_filter_pass_presence,
)


def merge_haplotagged_snv_tsvs(outdir, chromosomes):
    """
    Merge per-chromosome haplotagged SNV tables into one top-level TSV.
    Inputs are expected at:
      <outdir>/<chrom>_out/haplotagged_snvs.tsv
    Output is written to:
      <outdir>/haplotagged_snvs.tsv
    """
    merged_path = os.path.join(outdir, "haplotagged_snvs.tsv")
    wrote_header = False
    merged_files = 0
    rows_written = 0

    with open(merged_path, "w", encoding="utf-8") as out:
        for chrom, _length in chromosomes:
            per_chrom = os.path.join(outdir, f"{chrom}_out", "haplotagged_snvs.tsv")
            if not os.path.exists(per_chrom):
                continue

            with open(per_chrom, "r", encoding="utf-8") as inp:
                header = inp.readline()
                if not header:
                    continue
                if not wrote_header:
                    out.write(header.rstrip("\n") + "\n")
                    wrote_header = True

                for line in inp:
                    if not line.strip():
                        continue
                    out.write(line.rstrip("\n") + "\n")
                    rows_written += 1
            merged_files += 1

    # Ensure a valid header-only file exists even if no per-chrom files were found.
    if not wrote_header:
        with open(merged_path, "w", encoding="utf-8") as out:
            out.write(
                "CHR\tPOSITION\tTOTAL_COVERAGE\tHP_COUNTS"
                "\tHP1_HVAF\tHP2_HVAF\tHP_ASSIGNMENT\n"
            )

    return merged_path, merged_files, rows_written


def main():
    # Parse and normalize CLI inputs used by the orchestration layer.
    args = parse_args()
    vcf_path = args.input_vcf_gz
    bam_paths = args.haplotagged_bams
    outdir = os.path.abspath(args.output_dir)
    vcf_sample_name = args.vcf_sample_name
    germline_vcf_path = args.germline_vcf_gz
    germline_vcf_sample_name = args.germline_vcf_sample_name
    exclude_regions_bed = args.exclude_regions_bed
    divergent_same_hp = args.divergent_same_hp
    jobs = args.jobs
    post_vcfs = args.vcfs
    post_cn_bed_hp1 = args.cn_bed_hp1
    post_cn_bed_hp2 = args.cn_bed_hp2
    post_tree = args.tree
    excluded_chroms = parse_exclude_chroms(args.exclude_chrom)

    main_program = "./main"

    try:
        preflight_requirements(main_program)
    except RuntimeError as exc:
        sys.exit(f"Error: {exc}")

    # Validate all required on-disk inputs before launching any work.
    required_files = [vcf_path, main_program, *bam_paths]
    if germline_vcf_path:
        required_files.append(germline_vcf_path)
    if exclude_regions_bed:
        required_files.append(exclude_regions_bed)
    for f in required_files:
        if not os.path.exists(f):
            sys.exit(f"Error: required file not found: {f}")

    somatic_samples = get_vcf_samples(vcf_path)
    if vcf_sample_name:
        if vcf_sample_name not in somatic_samples:
            sys.exit(
                f"Error: somatic VCF sample '{vcf_sample_name}' was not found."
            )
    elif len(somatic_samples) == 1:
        vcf_sample_name = somatic_samples[0]
        print(f"Using sole somatic VCF sample: {vcf_sample_name}")
    elif len(somatic_samples) > 1:
        sys.exit(
            "Error: --somatic-sample-name is required because the somatic VCF "
            f"contains {len(somatic_samples)} samples."
        )

    if germline_vcf_path:
        germline_samples = get_vcf_samples(germline_vcf_path)
        if germline_vcf_sample_name not in germline_samples:
            sys.exit(
                "Error: germline VCF sample "
                f"'{germline_vcf_sample_name}' was not found."
            )

    os.makedirs(outdir, exist_ok=True)
    log_dir = os.path.join(outdir, "logs")
    os.makedirs(log_dir, exist_ok=True)
    run_stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_log_path = os.path.join(log_dir, f"run_{run_stamp}.log")

    # Derive genome name
    genome_name = os.path.basename(bam_paths[0]).replace(".bam", "")
    print(f"Genome name: {genome_name}")
    print(f"Input BAMs: {len(bam_paths)}")
    print(f"Logs directory: {log_dir}")
    print(f"Run log: {run_log_path}")

    # Build the chromosome task list and initialize genome-level accumulators.
    print("Extracting chromosome list from BAM header using samtools...")
    detected_chromosomes = get_chromosomes_from_bam(bam_paths[0])
    for bam_path in bam_paths[1:]:
        bam_chromosomes = get_chromosomes_from_bam(bam_path)
        if bam_chromosomes != detected_chromosomes:
            sys.exit(
                "Error: BAM sequence names and lengths do not match the first "
                f"input BAM: {bam_path}"
            )
    chromosomes = [
        (chrom, length)
        for chrom, length in detected_chromosomes
        if not chrom_is_excluded(chrom, excluded_chroms)
    ]
    excluded_detected_chromosomes = [
        (chrom, length)
        for chrom, length in detected_chromosomes
        if chrom_is_excluded(chrom, excluded_chroms)
    ]
    print(f"Detected {len(detected_chromosomes)} chromosomes.")
    if excluded_chroms:
        requested = ", ".join(sorted(excluded_chroms))
        matched = ", ".join(chrom for chrom, _ in excluded_detected_chromosomes) or "none"
        print(f"Excluded chromosomes requested: {requested}")
        print(f"Excluded chromosomes found in BAM: {matched}")
    if not chromosomes:
        sys.exit("Error: no chromosomes left after applying --exclude-chrom.")
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
        "hp1": 0,
        "hp2": 0,
        "mixed": 0,
        "unknown": 0,
    }
    genome_edge_haplotype_stats = new_edge_haplotype_summary()
    chrom_haplotag_detected = {}

    run_start_wall = datetime.now()
    run_start_monotonic = time.monotonic()
    invocation_cmd = " ".join(shlex.quote(x) for x in [sys.executable, *sys.argv])
    tool_threads = 1 if jobs > 1 else 4
    force_single_thread_tools = jobs > 1
    require_pass_filter, total_vcf_records, pass_vcf_records = scan_vcf_filter_pass_presence(vcf_path)
    germline_require_pass_filter = False
    germline_total_vcf_records = 0
    germline_pass_vcf_records = 0
    if germline_vcf_path:
        (
            germline_require_pass_filter,
            germline_total_vcf_records,
            germline_pass_vcf_records,
        ) = scan_vcf_filter_pass_presence(germline_vcf_path)
    print(
        "VCF filter mode: "
        + (
            f"require FILTER=PASS (found {pass_vcf_records} PASS records out of {total_vcf_records})"
            if require_pass_filter
            else f"ignore FILTER=PASS (found 0 PASS records out of {total_vcf_records})"
        )
    )
    if germline_vcf_path:
        print(
            "Germline VCF filter mode: "
            + (
                "require FILTER=PASS "
                f"(found {germline_pass_vcf_records} PASS records out of "
                f"{germline_total_vcf_records})"
                if germline_require_pass_filter
                else "ignore FILTER=PASS "
                f"(found 0 PASS records out of {germline_total_vcf_records})"
            )
        )
    print(
        "Haplotag mode: "
        + ("germline_vcf" if germline_vcf_path else "bam_hp_tags")
    )

    with open(run_log_path, "w") as run_log:
        write_run_header(
            run_log=run_log,
            run_start_wall=run_start_wall,
            invocation_cmd=invocation_cmd,
            genome_name=genome_name,
            vcf_path=vcf_path,
            bam_paths=bam_paths,
            outdir=outdir,
            vcf_sample_name=vcf_sample_name,
            germline_vcf_path=germline_vcf_path,
            germline_vcf_sample_name=germline_vcf_sample_name,
            exclude_regions_bed=exclude_regions_bed,
            germline_total_vcf_records=germline_total_vcf_records,
            germline_pass_vcf_records=germline_pass_vcf_records,
            germline_require_pass_filter=germline_require_pass_filter,
            divergent_same_hp=divergent_same_hp,
            jobs=jobs,
            tool_threads=tool_threads,
            total_vcf_records=total_vcf_records,
            pass_vcf_records=pass_vcf_records,
            require_pass_filter=require_pass_filter,
            chromosomes=chromosomes,
            genome_bp=genome_bp,
            detected_chromosome_count=len(detected_chromosomes),
            excluded_chroms=excluded_chroms,
            excluded_detected_chromosomes=excluded_detected_chromosomes,
        )
        run_log.flush()

        # Execute per-chromosome jobs (serial or thread pool) and fold results into genome totals.
        chrom_tasks = sorted(chromosomes, key=lambda x: x[1], reverse=True)
        if jobs == 1:
            for chrom, length in chrom_tasks:
                result = process_chromosome(
                    chrom=chrom,
                    length=length,
                    bam_paths=bam_paths,
                    vcf_path=vcf_path,
                    germline_vcf_path=germline_vcf_path,
                    outdir=outdir,
                    log_dir=log_dir,
                    main_program=main_program,
                    vcf_sample_name=vcf_sample_name,
                    germline_vcf_sample_name=germline_vcf_sample_name,
                    exclude_regions_bed=exclude_regions_bed,
                    require_pass_filter=require_pass_filter,
                    germline_require_pass_filter=germline_require_pass_filter,
                    divergent_same_hp=divergent_same_hp,
                    tool_threads=tool_threads,
                    force_single_thread_tools=force_single_thread_tools,
                )
                total_component_bp_summed += result["component_bp_summed"]
                total_component_bp_unique += result["component_bp_unique"]
                for key in genome_mutation_stats:
                    genome_mutation_stats[key] += result["mutation_stats"].get(key, 0)
                chrom_edge_hap = result["edge_haplotype_stats"]
                add_edge_haplotype_summary(genome_edge_haplotype_stats, chrom_edge_hap)
                chrom_haplotag_detected[result["chrom"]] = append_chromosome_result(run_log, result)
                run_log.flush()
        else:
            with ThreadPoolExecutor(max_workers=jobs) as executor:
                futures = {
                    executor.submit(
                        process_chromosome,
                        chrom=chrom,
                        length=length,
                        bam_paths=bam_paths,
                        vcf_path=vcf_path,
                        germline_vcf_path=germline_vcf_path,
                        outdir=outdir,
                        log_dir=log_dir,
                        main_program=main_program,
                        vcf_sample_name=vcf_sample_name,
                        germline_vcf_sample_name=germline_vcf_sample_name,
                        exclude_regions_bed=exclude_regions_bed,
                        require_pass_filter=require_pass_filter,
                        germline_require_pass_filter=germline_require_pass_filter,
                        divergent_same_hp=divergent_same_hp,
                        tool_threads=tool_threads,
                        force_single_thread_tools=force_single_thread_tools,
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
                        chrom_edge_hap = result["edge_haplotype_stats"]
                        add_edge_haplotype_summary(genome_edge_haplotype_stats, chrom_edge_hap)
                        chrom_haplotag_detected[result["chrom"]] = append_chromosome_result(run_log, result)
                        run_log.flush()
                except Exception as exc:
                    for f in futures:
                        f.cancel()
                    write_parallel_failure(run_log, exc)
                    run_log.flush()
                    raise

        # Emit final genome-level summaries after all chromosomes complete.
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
        write_total_component_line(
            run_log=run_log,
            total_component_bp_summed=total_component_bp_summed,
            total_component_bp_unique=total_component_bp_unique,
            genome_bp=genome_bp,
            genome_pct_summed=genome_pct_summed,
            genome_pct_unique=genome_pct_unique,
        )
        run_log.flush()

        mut_summary = build_mutation_summary(genome_mutation_stats)
        print(mut_summary)
        run_log.write(mut_summary + "\n")

        write_genome_edge_summary(
            run_log=run_log,
            genome_edge_haplotype_stats=genome_edge_haplotype_stats,
            chrom_haplotag_detected=chrom_haplotag_detected,
        )
        run_log.flush()

        run_end_wall = datetime.now()
        elapsed_seconds = time.monotonic() - run_start_monotonic
        elapsed_hms = format_elapsed_hms(elapsed_seconds)
        write_runtime_line(run_log, elapsed_seconds, elapsed_hms, run_start_wall, run_end_wall)
        run_log.flush()

    # Run optional/standard postprocessing as the final pipeline stage.
    postprocess_script = os.path.join(os.path.dirname(os.path.abspath(__file__)), "postprocess", "run_postprocessing.py")
    if not os.path.exists(postprocess_script):
        raise RuntimeError(f"Default postprocess script not found: {postprocess_script}")

    post_cmd = [sys.executable, postprocess_script, "--outdir", outdir]
    if post_vcfs:
        post_cmd.extend(["--vcfs", os.path.abspath(post_vcfs)])
    if post_cn_bed_hp1 and post_cn_bed_hp2:
        post_cmd.extend(["--cn-bed-hp1", os.path.abspath(post_cn_bed_hp1)])
        post_cmd.extend(["--cn-bed-hp2", os.path.abspath(post_cn_bed_hp2)])
    if post_tree:
        post_cmd.extend(["--tree", os.path.abspath(post_tree)])

    with open(run_log_path, "a", encoding="utf-8") as run_log:
        merged_tsv, merged_files, merged_rows = merge_haplotagged_snv_tsvs(outdir, chromosomes)
        run_log.write(
            f"haplotagged_snv_table={merged_tsv}\t"
            f"source_chrom_files={merged_files}\trows={merged_rows}\n"
        )
        run_log.flush()

        write_postprocessing_header(run_log)
        run_log.flush()
        run_cmd(post_cmd, log_fh=run_log)

    print(f"\nAll chromosomes processed for {genome_name}. Results saved in {outdir}")
    print(f"Total runtime: {elapsed_hms} ({elapsed_seconds:.3f}s)")


if __name__ == "__main__":
    main()
