#!/usr/bin/env python3
"""
Run postprocessing steps on the main program output: annotate (if --vcfs),
report_condensed_timing_chains (always), evaluate (if --tree), then
summarize_and_compress and phase-block plotting (always). Only runs
annotate/evaluate when their required inputs are provided.
"""

import argparse
import subprocess
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(
        description="Run postprocessing on main output: annotate (optional), report_condensed_timing_chains, evaluate (optional), then summarize_and_compress and plot_max_phase_block_by_chr (always)."
    )
    parser.add_argument(
        "--outdir",
        required=True,
        type=Path,
        help="Main output directory (used by all steps that run).",
    )
    parser.add_argument(
        "--vcfs",
        type=Path,
        default=None,
        help="Directory containing VCF/VCF.GZ files; if provided, run annotate_source_vcf.",
    )
    parser.add_argument(
        "--cn-bed-hp1",
        type=Path,
        default=None,
        help="Copy-number BED segments file for haplotype 1.",
    )
    parser.add_argument(
        "--cn-bed-hp2",
        type=Path,
        default=None,
        help="Copy-number BED segments file for haplotype 2.",
    )
    parser.add_argument(
        "--tree",
        type=Path,
        default=None,
        help="Tree file (parent<TAB>child); if provided, run evaluate_graphs.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands that would be run without executing them.",
    )
    args = parser.parse_args()

    postprocess_dir = Path(__file__).resolve().parent
    evo_dir = postprocess_dir.parent
    annotate_script = postprocess_dir / "annotate_source_vcf.py"
    timing_chains_script = postprocess_dir / "report_condensed_timing_chains.py"
    evaluate_script = postprocess_dir / "evaluate_graphs.py"
    summarize_script = evo_dir / "summarize_and_compress.py"
    plot_script = postprocess_dir / "plot_max_phase_block_by_chr.py"

    outdir = args.outdir.resolve()

    steps = []
    if (args.cn_bed_hp1 is None) != (args.cn_bed_hp2 is None):
        parser.error("Provide both --cn-bed-hp1 and --cn-bed-hp2 together.")

    if args.vcfs is not None or args.cn_bed_hp1 is not None:
        annotate_cmd = [sys.executable, str(annotate_script), "--main_out", str(outdir)]
        if args.vcfs is not None:
            annotate_cmd.extend(["--vcfs", str(args.vcfs.resolve())])
        if args.cn_bed_hp1 is not None and args.cn_bed_hp2 is not None:
            annotate_cmd.extend(["--cn-bed-hp1", str(args.cn_bed_hp1.resolve())])
            annotate_cmd.extend(["--cn-bed-hp2", str(args.cn_bed_hp2.resolve())])
        steps.append(
            (
                "annotate_source_vcf",
                annotate_cmd,
            )
        )
    steps.append(
        ("report_condensed_timing_chains", [sys.executable, str(timing_chains_script), "--outdir", str(outdir)]),
    )
    if args.tree is not None:
        steps.append(
            (
                "evaluate_graphs",
                [sys.executable, str(evaluate_script), "--outdir", str(outdir), "--tree", str(args.tree.resolve())],
            )
        )
    steps.append(
        ("summarize_and_compress", [sys.executable, str(summarize_script), str(outdir)])
    )
    steps.append(
        (
            "plot_max_phase_block_by_chr",
            [
                sys.executable,
                str(plot_script),
                str(outdir / "global_summary.txt"),
            ],
        )
    )

    for i, (name, argv) in enumerate(steps, start=1):
        print(f"Running ({i}/{len(steps)}): {name}")
        if args.dry_run:
            print("  ", " ".join(argv))
            continue
        result = subprocess.run(argv)
        if result.returncode != 0:
            print(f"Error: {name} exited with code {result.returncode}", file=sys.stderr)
            sys.exit(result.returncode)

    if args.dry_run:
        print("(dry-run; no commands executed)")
    else:
        print("Done.")


# tree /data/KolmogorovLab/donmeza2/cancer_phasing_data/mouse/tree.tsv
# example command: python postprocess/run_postprocessing.py --outdir /data/KolmogorovLab/donmeza2/cancer_phasing_data/mouse/output --vcfs /data/KolmogorovLab/donmeza2/cancer_phasing_data/mouse/vcfs --tree /data/KolmogorovLab/donmeza2/cancer_phasing_data/mouse/tree.tsv
if __name__ == "__main__":
    main()
