#!/usr/bin/env python3
import os
import pandas as pd
from pathlib import Path
import tarfile
from pathlib import Path


def summarize_component_stats(stats_file):
    import pandas as pd

    df = pd.read_csv(stats_file, sep="\t")

    # Handle empty files
    if df.empty:
        print(f"Warning: {stats_file} is empty — skipping.")
        return {
            "avg_span_bp": 0,
            "top10_max_span_bp": "",
            "avg_haplotypes": 0,
            "top10_max_haplotypes": "",
            "avg_num_nodes": 0,
            "top10_max_num_nodes": "",
            "num_components": 0,
        }

    # Convert numeric columns safely
    for col in ["span_bp", "haplotypes", "num_nodes"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    # Drop rows where all numeric values are NaN
    df = df.dropna(subset=["span_bp", "haplotypes", "num_nodes"])
    if df.empty:
        print(f"Warning: {stats_file} has no valid numeric rows — skipping.")
        return {
            "avg_span_bp": 0,
            "top10_max_span_bp": "",
            "avg_haplotypes": 0,
            "top10_max_haplotypes": "",
            "avg_num_nodes": 0,
            "top10_max_num_nodes": "",
            "num_components": 0,
        }

    return {
        "avg_span_bp": df["span_bp"].mean(),
        "top10_max_span_bp": ",".join(map(str, df["span_bp"].nlargest(10).tolist())),
        "avg_haplotypes": df["haplotypes"].mean(),
        "top10_max_haplotypes": ",".join(map(str, df["haplotypes"].nlargest(10).tolist())),
        "avg_num_nodes": df["num_nodes"].mean(),
        "top10_max_num_nodes": ",".join(map(str, df["num_nodes"].nlargest(10).tolist())),
        "num_components": len(df),
    }



def chrom_sort_key(chrom):
    """Sort chromosomes numerically (1–22) then X, Y, others."""
    name = chrom.replace("chr", "").replace("Chr", "")
    if name.isdigit():
        return (int(name), "")
    elif name.upper() == "X":
        return (23, "")
    elif name.upper() == "Y":
        return (24, "")
    else:
        # Anything nonstandard goes last
        return (25, name)


def compress_all_results(main_output_dir):
    """
    Create a single .tar.gz archive containing all chromosome outputs and global summary.
    The archive will be named after the main output directory.
    """
    main_output_dir = Path(main_output_dir)
    global_summary = main_output_dir / "global_summary.txt"

    # Use the name of the main directory for the tar file
    tar_name = f"{main_output_dir.name}.tar.gz"
    final_archive = main_output_dir / tar_name

    print(f"\nCreating single archive: {final_archive}")

    with tarfile.open(final_archive, "w:gz") as tar:
        # Add global summary
        if global_summary.exists():
            tar.add(global_summary, arcname="global_summary.txt")

        # Add per-chromosome outputs
        for subdir in sorted(main_output_dir.glob("*_out")):
            chrom = subdir.name.replace("_out", "")
            stats_file = subdir / "component_statistics.txt"
            graphs_dir = subdir / "graphs"

            if stats_file.exists():
                tar.add(stats_file, arcname=f"{chrom}/component_statistics.txt")
            if graphs_dir.exists():
                tar.add(graphs_dir, arcname=f"{chrom}/graphs")

    print(f"All results successfully packaged into: {final_archive}")



def main(main_output_dir):
    main_output_dir = Path(main_output_dir)
    global_file = main_output_dir / "global_summary.txt"

    all_summaries = []

    for subdir in sorted(main_output_dir.glob("*_out")):
        chrom = subdir.name.replace("_out", "")
        stats_file = subdir / "component_statistics.txt"
        if not stats_file.exists():
            print(f"Warning: Missing {stats_file}, skipping.")
            continue

        print(f"Processing {chrom} ...")
        summary = summarize_component_stats(stats_file)
        summary["chromosome"] = chrom
        all_summaries.append(summary)

    if not all_summaries:
        print("No summaries found. Exiting.")
        return

    import pandas as pd
    df_all = pd.DataFrame(all_summaries)

    # --- Sort chromosomes properly ---
    df_all = df_all.sort_values(by="chromosome", key=lambda x: x.map(chrom_sort_key))

    # --- Reorder columns: put chromosome first ---
    desired_col_order = ["chromosome"] + [c for c in df_all.columns if c != "chromosome"]
    df_all = df_all[desired_col_order]

    # --- Write global summary ---
    df_all.to_csv(global_file, sep="\t", index=False)
    print(f"\nGlobal summary written to: {global_file}")




if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Summarize component_statistics.txt files for all chromosomes into one global summary."
    )
    parser.add_argument("main_output_dir", help="Main directory containing *_out folders")
    args = parser.parse_args()

    main(args.main_output_dir)
    compress_all_results(args.main_output_dir)

