#!/usr/bin/env python3

from run_chr_log_stats import (
    RELATIONS,
    extract_haplotag_status_from_chrom_log,
    extract_max_span_scan_from_chrom_log,
    parse_haplotag_detected_value,
)


def _pct(num, denom):
    return f"{100.0 * num / denom:.1f}%" if denom > 0 else "N/A"


def write_run_header(
    run_log,
    run_start_wall,
    invocation_cmd,
    genome_name,
    vcf_path,
    bam_path,
    outdir,
    vcf_sample_name,
    jobs,
    tool_threads,
    total_vcf_records,
    pass_vcf_records,
    require_pass_filter,
    chromosomes,
    genome_bp,
):
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
    run_log.write("main_default_min_reads=2\n")
    run_log.write("max_pair_distance_source=computed_by_main_per_chrom_bam\n")
    run_log.write("haplotag_detect_source=computed_by_main_per_chrom_bam\n")
    run_log.write(f"jobs={jobs}\n")
    run_log.write(f"tool_threads_per_subprocess={tool_threads}\n")
    run_log.write(f"vcf_total_records_scanned={total_vcf_records}\n")
    run_log.write(f"vcf_pass_records_found={pass_vcf_records}\n")
    run_log.write(f"vcf_require_filter_PASS={int(require_pass_filter)}\n")
    run_log.write("parallel_mode=chromosome_thread_pool\n")
    run_log.write(f"Detected chromosomes: {len(chromosomes)}\n")
    run_log.write(f"Genome bp (selected chromosomes): {genome_bp}\n")
    run_log.write(
        "per_chrom_max_span_summary_source=parsed_from_chrom_log_[max_span_scan]_line\n"
    )


def append_chromosome_result(run_log, result):
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

    haplotag_status = extract_haplotag_status_from_chrom_log(result["chrom_log_path"])
    if "haplotagged" in haplotag_status:
        source = haplotag_status.get("source", "unknown")
        run_log.write(
            f"{result['chrom']}\thaplotag_detected={haplotag_status['haplotagged']}\t"
            f"haplotag_source={source}\n"
        )
    else:
        run_log.write(f"{result['chrom']}\thaplotag_detected=not_found_in_chrom_log\n")

    return parse_haplotag_detected_value(haplotag_status)


def write_parallel_failure(run_log, exc):
    run_log.write(f"ERROR\tparallel_chromosome_processing_failed\t{exc}\n")


def write_total_component_line(
    run_log,
    total_component_bp_summed,
    total_component_bp_unique,
    genome_bp,
    genome_pct_summed,
    genome_pct_unique,
):
    run_log.write(
        f"TOTAL\tconnected_component_bp_summed={total_component_bp_summed}\t"
        f"connected_component_bp_unique={total_component_bp_unique}\t"
        f"genome_bp={genome_bp}\t"
        f"genome_pct_summed={genome_pct_summed:.6f}\t"
        f"genome_pct_unique={genome_pct_unique:.6f}\n"
    )


def build_mutation_summary(genome_mutation_stats):
    g = genome_mutation_stats
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
        f"    In graph (has at least one accepted edge): {g['in_graph']}"
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
            f"  Singleton (has no pairwise relations): {g['singleton']}"
            f"  ({_pct(g['singleton'], g['total_snvs'])} of total)"
        )
    mut_lines.append(
        f"  HP1 assigned:                         {g['hp1']}"
        + (f"  ({_pct(g['hp1'], g['total_snvs'])} of total)" if g["total_snvs"] > 0 else "")
    )
    mut_lines.append(
        f"  HP2 assigned:                         {g['hp2']}"
        + (f"  ({_pct(g['hp2'], g['total_snvs'])} of total)" if g["total_snvs"] > 0 else "")
    )
    mut_lines.append(
        f"  MIXED assigned:                       {g['mixed']}"
        + (f"  ({_pct(g['mixed'], g['total_snvs'])} of total)" if g["total_snvs"] > 0 else "")
    )
    mut_lines.append(
        f"  UNKNOWN assigned:                     {g['unknown']}"
        + (f"  ({_pct(g['unknown'], g['total_snvs'])} of total)" if g["total_snvs"] > 0 else "")
    )
    mut_lines.append(f"Total accepted edges:                  {g['accepted_edges']}")
    mut_lines.append(f"  - timing edges:                      {g['timing_edges']}")
    return "\n".join(mut_lines)


def write_genome_edge_summary(run_log, genome_edge_haplotype_stats, chrom_haplotag_detected):
    run_log.write("\n=== Genome-wide edge summary ===\n")
    without_haplotag = sorted(
        chrom for chrom, detected in chrom_haplotag_detected.items() if detected == 0
    )
    any_haplotag_detected = any(
        detected == 1 for detected in chrom_haplotag_detected.values()
    )
    if without_haplotag:
        run_log.write(
            "Chromosomes without haplotagging detected: "
            + ", ".join(without_haplotag)
            + "\n"
        )
    if not any_haplotag_detected:
        run_log.write(
            "Haplotype-splitting entries omitted: no chromosome had haplotagging detected.\n"
        )

    relation_label = {
        "cooccurring": "Co-occurring",
        "timing": "Timing",
        "divergent": "Divergent",
    }
    for relation in RELATIONS:
        rec = genome_edge_haplotype_stats[relation]

        run_log.write(f"\n{relation_label[relation]} edges (accepted):        {rec['total']}\n")
        run_log.write(
            f"  Without loss:                      {rec['without_loss']}"
            + (f"  ({_pct(rec['without_loss'], rec['total'])} of relation)" if rec["total"] > 0 else "")
            + "\n"
        )
        run_log.write(
            f"  With loss:                         {rec['with_loss']}"
            + (f"  ({_pct(rec['with_loss'], rec['total'])} of relation)" if rec["total"] > 0 else "")
            + "\n"
        )

        if any_haplotag_detected:
            run_log.write(
                f"  Same haplotype:                    {rec['same_haplotype']}"
                + (f"  ({_pct(rec['same_haplotype'], rec['total'])} of relation)" if rec["total"] > 0 else "")
                + "\n"
            )
            run_log.write(
                f"    - HP1-HP1:                       {rec['hp1_hp1']}"
                + (f"  ({_pct(rec['hp1_hp1'], rec['same_haplotype'])} of same-haplotype)" if rec["same_haplotype"] > 0 else "")
                + "\n"
            )
            run_log.write(
                f"    - HP2-HP2:                       {rec['hp2_hp2']}"
                + (f"  ({_pct(rec['hp2_hp2'], rec['same_haplotype'])} of same-haplotype)" if rec["same_haplotype"] > 0 else "")
                + "\n"
            )
            run_log.write(
                f"  Different haplotypes (HP1-HP2):    {rec['different_haplotype']}"
                + (f"  ({_pct(rec['different_haplotype'], rec['total'])} of relation)" if rec["total"] > 0 else "")
                + "\n"
            )
            run_log.write(
                f"  At least one unknown or mixed:     {rec['at_least_one_unknown']}"
                + (f"  ({_pct(rec['at_least_one_unknown'], rec['total'])} of relation)" if rec["total"] > 0 else "")
                + "\n"
            )


def format_elapsed_hms(elapsed_seconds):
    elapsed_h = int(elapsed_seconds // 3600)
    elapsed_m = int((elapsed_seconds % 3600) // 60)
    elapsed_s = elapsed_seconds % 60
    return f"{elapsed_h:02d}:{elapsed_m:02d}:{elapsed_s:06.3f}"


def write_runtime_line(run_log, elapsed_seconds, elapsed_hms, run_start_wall, run_end_wall):
    runtime_line = (
        f"TOTAL_RUNTIME\tseconds={elapsed_seconds:.3f}\thms={elapsed_hms}\t"
        f"started={run_start_wall.isoformat(timespec='seconds')}\t"
        f"ended={run_end_wall.isoformat(timespec='seconds')}\n"
    )
    run_log.write(runtime_line)


def write_postprocessing_header(run_log):
    run_log.write("\n=== Postprocessing ===\n")
