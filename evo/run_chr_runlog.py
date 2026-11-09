#!/usr/bin/env python3

from run_chr_log_stats import (
    RELATIONS,
    parse_haplotag_detected_value,
)


def _pct(num, denom):
    return f"{100.0 * num / denom:.1f}%" if denom > 0 else "N/A"


def write_run_header(run_log, args):
    run_log.write(f"Run started: {args.run_start.isoformat(timespec='seconds')}\n")
    run_log.write(f"Invocation command: {args.invocation_cmd}\n")
    run_log.write(f"Genome: {args.genome_name}\n")
    run_log.write(f"Somatic VCF: {args.input_vcf_gz}\n")
    for bam_path in args.bams:
        run_log.write(f"BAM: {bam_path}\n")
    run_log.write(f"Output: {args.output_dir}\n")
    if args.germline_vcf_gz:
        run_log.write(f"Germline VCF: {args.germline_vcf_gz}\n")
        run_log.write(f"Germline VCF sample name: {args.germline_vcf_sample_name}\n")
    for bed_file in args.exclude_regions_bed:
        run_log.write(f"Excluded-regions BED: {bed_file}\n")
    if args.cn_bed:
        run_log.write(f"Copy-number BED (annotation only): {args.cn_bed}\n")
    run_log.write("\n=== Run parameters ===\n")
    haplotag_mode = "germline_vcf" if args.germline_vcf_gz else "none"
    run_log.write(f"haplotag_mode={haplotag_mode}\n")
    run_log.write(f"divergent_same_hp={int(args.divergent_same_hp)}\n")
    run_log.write(f"maximum_mutations={args.window_size}\n")
    run_log.write(f"group_stride={args.window_stride}\n")
    run_log.write(
        "haplotag_detect_source="
        + (
            "germline_vcf_candidate_reads"
            if args.germline_vcf_gz
            else "none"
        )
        + "\n"
    )
    run_log.write(f"jobs={args.jobs}\n")
    run_log.write(f"tool_threads_per_subprocess={args.tool_threads}\n")
    run_log.write(f"vcf_total_records_scanned={args.total_vcf_records}\n")
    run_log.write(f"vcf_pass_records_found={args.pass_vcf_records}\n")
    run_log.write(f"vcf_require_filter_PASS={int(args.require_pass_filter)}\n")
    if args.germline_vcf_gz:
        run_log.write(f"germline_vcf_total_records_scanned={args.germline_total_vcf_records}\n")
        run_log.write(f"germline_vcf_pass_records_found={args.germline_pass_vcf_records}\n")
        run_log.write(f"germline_vcf_require_filter_PASS={int(args.germline_require_pass_filter)}\n")
    run_log.write("parallel_mode=chromosome_thread_pool\n")
    run_log.write(f"Detected chromosomes before exclusion: {len(args.detected_chromosomes)}\n")
    if args.excluded_chroms:
        requested = ",".join(sorted(args.excluded_chroms))
        matched = ",".join(chrom for chrom, _ in args.excluded_detected_chromosomes) or "none"
        run_log.write(f"Excluded chromosomes requested: {requested}\n")
        run_log.write(f"Excluded chromosomes found in BAM: {matched}\n")
    run_log.write(f"Detected chromosomes: {len(args.chromosomes)}\n")
    genome_bp = sum(length for _, length in args.chromosomes)
    run_log.write(f"Genome bp (selected chromosomes): {genome_bp}\n")


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

    haplotag_status = result["haplotag_status"]
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


def _write_relation_summary(run_log, edge_stats, include_haplotypes):
    relation_label = {
        "cooccurring": "Co-occurring",
        "timing": "Timing",
        "divergent": "Divergent",
    }
    for relation in RELATIONS:
        rec = edge_stats[relation]
        run_log.write(f"\n{relation_label[relation]} edges (accepted):        {rec['total']}\n")
        if include_haplotypes:
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


def write_chromosome_edge_summary(chrom_log, haplotag_detected, chrom_edge_hap):
    chrom_log.write("\n=== Chromosome-wide edge summary ===\n")
    if haplotag_detected in (0, 1):
        chrom_log.write(f"Haplotagging detected for this chromosome: {haplotag_detected}\n")
    else:
        chrom_log.write("Haplotagging detected for this chromosome: unknown\n")
    if haplotag_detected == 0:
        chrom_log.write(
            "Haplotype-splitting entries omitted: haplotagging was not detected for this chromosome.\n"
        )
    _write_relation_summary(chrom_log, chrom_edge_hap, haplotag_detected != 0)


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

    _write_relation_summary(run_log, genome_edge_haplotype_stats, any_haplotag_detected)


def format_elapsed_hms(elapsed_seconds):
    elapsed_h = int(elapsed_seconds // 3600)
    elapsed_m = int((elapsed_seconds % 3600) // 60)
    elapsed_s = elapsed_seconds % 60
    return f"{elapsed_h:02d}:{elapsed_m:02d}:{elapsed_s:06.3f}"


def write_runtime_line(run_log, run_start, run_end):
    elapsed_seconds = (run_end - run_start).total_seconds()
    elapsed_hms = format_elapsed_hms(elapsed_seconds)
    run_log.write(
        f"TOTAL_RUNTIME\tseconds={elapsed_seconds:.3f}\thms={elapsed_hms}\t"
        f"started={run_start.isoformat(timespec='seconds')}\t"
        f"ended={run_end.isoformat(timespec='seconds')}\n"
    )
    return elapsed_seconds, elapsed_hms


def write_postprocessing_header(run_log):
    run_log.write("\n=== Postprocessing ===\n")
