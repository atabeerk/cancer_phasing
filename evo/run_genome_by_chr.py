#!/usr/bin/env python3
import os
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime

from graph_ops.simplify_graph import write_topn_graphs
from run_chr_cli import parse_and_validate_args
from run_chr_log_stats import MUTATION_STAT_FIELDS, add_edge_haplotype_summary, new_edge_haplotype_summary
from run_chr_runlog import (
    append_chromosome_result,
    build_mutation_summary,
    write_genome_edge_summary,
    write_parallel_failure,
    write_postprocessing_header,
    write_run_header,
    write_runtime_line,
)
from run_chr_worker import process_chromosome, scan_vcf_filter_pass_presence
from postprocess_results import postprocess_results


def main():
    args = parse_and_validate_args()
    args.run_start = datetime.now()
    args.invocation_cmd = f"{sys.executable} {' '.join(sys.argv)}"
    args.genome_name = os.path.basename(args.bams[0]).replace(".bam", "")
    args.log_dir = os.path.join(args.output_dir, "logs")
    args.tool_threads = 1 if args.jobs > 1 else 4
    args.force_single_thread_tools = args.jobs > 1

    os.makedirs(args.log_dir, exist_ok=True)
    run_stamp = args.run_start.strftime("%Y%m%d_%H%M%S")
    run_log_path = os.path.join(args.log_dir, f"run_{run_stamp}.log")
    print(f"Genome name: {args.genome_name}")
    print(f"Input BAMs: {len(args.bams)}")
    print(f"Logs directory: {args.log_dir}")
    print(f"Run log: {run_log_path}")

    (
        args.require_pass_filter,
        args.total_vcf_records,
        args.pass_vcf_records,
    ) = scan_vcf_filter_pass_presence(args.input_vcf_gz)
    args.germline_require_pass_filter = False
    args.germline_total_vcf_records = 0
    args.germline_pass_vcf_records = 0
    if args.germline_vcf_gz:
        (
            args.germline_require_pass_filter,
            args.germline_total_vcf_records,
            args.germline_pass_vcf_records,
        ) = scan_vcf_filter_pass_presence(args.germline_vcf_gz)
    print(
        "VCF filter mode: "
        + (
            f"require FILTER=PASS (found {args.pass_vcf_records} PASS records out of {args.total_vcf_records})"
            if args.require_pass_filter
            else f"ignore FILTER=PASS (found 0 PASS records out of {args.total_vcf_records})"
        )
    )
    if args.germline_vcf_gz:
        print(
            "Germline VCF filter mode: "
            + (
                f"require FILTER=PASS (found {args.germline_pass_vcf_records} PASS records out of "
                f"{args.germline_total_vcf_records})"
                if args.germline_require_pass_filter
                else f"ignore FILTER=PASS (found 0 PASS records out of {args.germline_total_vcf_records})"
            )
        )
    print(f"Haplotag mode: {'germline_vcf' if args.germline_vcf_gz else 'none'}")

    genome_mutation_stats = dict.fromkeys(MUTATION_STAT_FIELDS, 0)
    genome_edge_haplotype_stats = new_edge_haplotype_summary()
    chrom_haplotag_detected = {}

    with open(run_log_path, "w") as run_log:
        write_run_header(run_log, args)
        run_log.flush()

        chrom_tasks = sorted(args.chromosomes, key=lambda item: item[1], reverse=True)
        if args.jobs == 1:
            for chrom, length in chrom_tasks:
                result = process_chromosome(chrom, length, args)
                for key in genome_mutation_stats:
                    genome_mutation_stats[key] += result["mutation_stats"].get(key, 0)
                add_edge_haplotype_summary(genome_edge_haplotype_stats, result["edge_haplotype_stats"])
                chrom_haplotag_detected[result["chrom"]] = append_chromosome_result(run_log, result)
                run_log.flush()
        else:
            with ThreadPoolExecutor(max_workers=args.jobs) as executor:
                futures = {
                    executor.submit(process_chromosome, chrom, length, args): (chrom, length)
                    for chrom, length in chrom_tasks
                }
                try:
                    for future in as_completed(futures):
                        result = future.result()
                        for key in genome_mutation_stats:
                            genome_mutation_stats[key] += result["mutation_stats"].get(key, 0)
                        add_edge_haplotype_summary(genome_edge_haplotype_stats, result["edge_haplotype_stats"])
                        chrom_haplotag_detected[result["chrom"]] = append_chromosome_result(run_log, result)
                        run_log.flush()
                except Exception as exc:
                    for future in futures:
                        future.cancel()
                    write_parallel_failure(run_log, exc)
                    run_log.flush()
                    raise

        mut_summary = build_mutation_summary(genome_mutation_stats)
        print(mut_summary)
        run_log.write(mut_summary + "\n")
        write_genome_edge_summary(run_log, genome_edge_haplotype_stats, chrom_haplotag_detected)

        topn_paths = write_topn_graphs(args.output_dir, args.top_components, args.chromosomes)
        topn_message = "Top-N graphs: " + ", ".join(topn_paths)
        print(topn_message)
        run_log.write(topn_message + "\n")
        run_log.flush()

        write_postprocessing_header(run_log)
        postprocess_results(args.output_dir, args.chromosomes, args.cn_bed, run_log,
                            has_germline_vcf=bool(args.germline_vcf_gz))

        elapsed_seconds, elapsed_hms = write_runtime_line(run_log, args.run_start, datetime.now())

    print(f"\nAll chromosomes processed for {args.genome_name}. Results saved in {args.output_dir}")
    print(f"Total runtime: {elapsed_hms} ({elapsed_seconds:.3f}s)")


if __name__ == "__main__":
    main()
