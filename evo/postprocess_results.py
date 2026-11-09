"""Merge final mutation tables, annotate CN, report chains, summarize and package."""

import csv
import os
import sys
import tarfile
from contextlib import redirect_stdout
from pathlib import Path
from statistics import mean, median

from postprocess.annotate_copy_number import annotate_copy_number
from postprocess.report_condensed_timing_chains import write_timing_chains


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
                "CHR\tPOSITION\tTOTAL_COVERAGE"
                "\tHP1_SUPPORTING\tHP1_UNSUPPORTING"
                "\tHP2_SUPPORTING\tHP2_UNSUPPORTING"
                "\tUNKNOWN_SUPPORTING\tUNKNOWN_UNSUPPORTING"
                "\tHP1_HVAF\tHP2_HVAF\tHP_ASSIGNMENT\n"
            )

    return merged_path, merged_files, rows_written



SUMMARY_FIELDS = [
    "chromosome", "avg_span_bp", "median_span_bp", "top10_max_span_bp",
    "avg_num_nodes", "top10_max_num_nodes", "num_components",
]


def summarize_component_stats(stats_file):
    """Summarize a chromosome, including chromosomes with no accepted edges."""
    rows = []
    if stats_file.exists():
        with stats_file.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
    spans = [int(row["span_bp"]) for row in rows]
    sizes = [int(row["num_nodes"]) for row in rows]
    return {
        "avg_span_bp": f"{mean(spans) if spans else 0:.1f}",
        "median_span_bp": f"{median(spans) if spans else 0:.1f}",
        "top10_max_span_bp": ",".join(str(value) for value in sorted(spans, reverse=True)[:10]),
        "avg_num_nodes": f"{mean(sizes) if sizes else 0:.1f}",
        "top10_max_num_nodes": ",".join(str(value) for value in sorted(sizes, reverse=True)[:10]),
        "num_components": len(rows),
    }


def chrom_sort_key(chrom):
    name = chrom.removeprefix("chr").removeprefix("Chr")
    if name.isdigit():
        return int(name), ""
    return {"X": 23, "Y": 24}.get(name.upper(), 25), name


def write_global_summary(outdir, chromosomes):
    output = outdir / "global_summary.txt"
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for chrom, _length in sorted(chromosomes, key=lambda item: chrom_sort_key(item[0])):
            stats = summarize_component_stats(outdir / f"{chrom}_out" / "component_statistics.txt")
            writer.writerow({"chromosome": chrom, **stats})
    return output


def compress_all_results(outdir, chromosomes, run_log_path=None, *, has_germline_vcf=False):
    """Package final outputs and logs; omit temporary BAMs and pileups."""
    archive = outdir / f"{outdir.name}.tar.gz"
    with tarfile.open(archive, "w:gz") as tar:
        names = ["global_summary.txt", "timing_chains.tsv", "merged_haplotype_cn_segments.tsv", "topN_graphs"]
        if has_germline_vcf:
            names.append("haplotagged_snvs.tsv")
        for name in names:
            path = outdir / name
            if path.exists():
                tar.add(path, arcname=name)
        logs = outdir / "logs"
        if logs.exists():
            names = {f"{chrom}.log" for chrom, _length in chromosomes}
            if run_log_path is not None:
                names.add(Path(run_log_path).name)
            for name in sorted(names):
                path = logs / name
                if path.is_file():
                    tar.add(path, arcname=f"logs/{name}")
        chromosome_outputs = ["component_statistics.txt", "graphs", "inference"]
        if has_germline_vcf:
            chromosome_outputs.append("haplotagged_snvs.tsv")
        for chrom, _length in chromosomes:
            for name in chromosome_outputs:
                path = outdir / f"{chrom}_out" / name
                if path.exists():
                    tar.add(path, arcname=f"{chrom}/{name}")
    return archive


def postprocess_results(outdir, chromosomes, cn_bed=None, log_fh=None, *, has_germline_vcf=False):
    """Finish one genome run using its selected chromosomes."""
    outdir = Path(outdir)

    def report(message):
        print(message)
        if log_fh is not None:
            log_fh.write(message + "\n")
            log_fh.flush()

    if has_germline_vcf:
        merged, files, rows = merge_haplotagged_snv_tsvs(outdir, chromosomes)
        report(f"haplotagged_snv_table={merged}\tsource_chrom_files={files}\trows={rows}")
    if cn_bed is not None:
        with redirect_stdout(log_fh if log_fh is not None else sys.stdout):
            nodes, components = annotate_copy_number(outdir, chromosomes, cn_bed)
        report(f"CN annotations: {nodes} graph nodes, {components} components")
    chains, count = write_timing_chains(outdir, chromosomes)
    report(f"Timing chains: {count} complete paths written to {chains}")
    summary = write_global_summary(outdir, chromosomes)
    report(f"Global summary written to: {summary}")
    report("Packaging final results and logs")
    archive = compress_all_results(outdir, chromosomes, getattr(log_fh, "name", None),
                                   has_germline_vcf=has_germline_vcf)
    report(f"Results archive: {archive}")
