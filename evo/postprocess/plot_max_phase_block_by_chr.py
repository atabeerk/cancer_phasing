#!/usr/bin/env python3
"""
Plot maximum phase block span by chromosome from global_summary.txt files.

For each input file, this script reads `top10_max_span_bp`, takes the first value
as the per-chromosome maximum span, and writes a bar plot to the same directory.
"""

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter


def chrom_sort_key(chrom: str):
    name = chrom.strip()
    if name.lower().startswith("chr"):
        name = name[3:]
    name_upper = name.upper()

    if name.isdigit():
        return (0, int(name), "")
    if name_upper == "X":
        return (1, 23, "")
    if name_upper == "Y":
        return (1, 24, "")
    if name_upper in {"M", "MT"}:
        return (1, 25, "")
    return (2, 999, name_upper)


def extract_max_spans(summary_path: Path):
    rows = []
    with summary_path.open("r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"chromosome", "top10_max_span_bp"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"{summary_path} is missing required columns: {', '.join(sorted(missing))}"
            )

        for row in reader:
            chrom = (row.get("chromosome") or "").strip()
            top10 = (row.get("top10_max_span_bp") or "").strip()
            if not chrom or not top10:
                continue

            first_value = top10.split(",")[0].strip()
            try:
                max_span_bp = int(float(first_value))
            except ValueError:
                continue

            # Keep autosomes only (1-22), excluding X/Y/M.
            name = chrom
            if name.lower().startswith("chr"):
                name = name[3:]
            if not name.isdigit():
                continue
            chrom_num = int(name)
            if chrom_num < 1 or chrom_num > 22:
                continue

            normalized_chrom = f"chr{chrom_num}"
            rows.append((normalized_chrom, max_span_bp))

    rows.sort(key=lambda x: chrom_sort_key(x[0]))
    return rows


def plot_one_summary(summary_path: Path, output_name: str, dpi: int):
    values = extract_max_spans(summary_path)
    if not values:
        raise ValueError(f"No plottable chromosome spans found in {summary_path}")

    chromosomes = [chrom for chrom, _ in values]
    max_spans_mb = [span / 1_000_000.0 for _, span in values]

    fig, ax = plt.subplots(figsize=(12, 5))
    # Warm color-blind-friendly styling (Okabe-Ito inspired orange tones).
    ax.bar(
        chromosomes,
        max_spans_mb,
        width=0.8,
        color="#E69F00",
        edgecolor="#D55E00",
        linewidth=0.7,
    )
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("Maximum phase block length")
    ax.set_title("Maximum phase block length by chromosome")
    ax.tick_params(axis="x", rotation=45)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g} Mb"))
    ax.grid(axis="y", alpha=0.25, linewidth=0.7)
    fig.tight_layout()

    output_path = summary_path.parent / output_name
    fig.savefig(output_path, dpi=dpi)
    plt.close(fig)
    return output_path


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Generate chromosome-vs-maximum-phase-block-length plots from one or more "
            "global_summary.txt files."
        )
    )
    parser.add_argument(
        "summary_files",
        nargs="+",
        help="One or more paths to global_summary.txt files.",
    )
    parser.add_argument(
        "--output-name",
        default="max_phase_block_by_chr.png",
        help="Output filename written next to each input file.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=150,
        help="Output PNG resolution (default: 150).",
    )
    args = parser.parse_args()

    for path_str in args.summary_files:
        summary_path = Path(path_str).resolve()
        if not summary_path.exists():
            raise FileNotFoundError(f"Input file not found: {summary_path}")
        output_path = plot_one_summary(summary_path, args.output_name, args.dpi)
        print(f"Wrote: {output_path}")


if __name__ == "__main__":
    main()
