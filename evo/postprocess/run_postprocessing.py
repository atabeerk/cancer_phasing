#!/usr/bin/env python3
"""
Run postprocessing steps on main output:
  1) annotate_source_vcf (optional; if --vcfs and/or --cn-bed)
  2) report_condensed_timing_chains (always)
  3) evaluate_graphs (optional; if --tree)
  4) plot_edge_eval_vcfpair_heatmap (optional; if --tree)
  5) summarize_and_compress (always)
  6) plot_timing_cn_interactive (always; CN overlays optional)
  7) plot_max_phase_block_by_chr (always)
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


_SUMMARY_PREFIXES = (
    "Total ",
    "Rank distribution",
    "Wrote ",
)


def _extract_summary(text: str) -> str:
    """
    Return only summary lines plus indented detail lines that follow a summary
    header (e.g. per-rank counts under "Rank distribution ...").
    """
    out: list[str] = []
    in_block = False
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped:
            in_block = False
            continue
        if any(stripped.startswith(p) for p in _SUMMARY_PREFIXES):
            out.append(stripped)
            in_block = True
        elif in_block and line.startswith("  "):
            out.append(stripped)
        else:
            in_block = False
    return "\n".join(out)


def _tail_lines(text: str, n: int = 25) -> str:
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if not lines:
        return ""
    if len(lines) <= n:
        return "\n".join(lines)
    return "\n".join(lines[-n:])


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Run postprocessing on main output: annotate (optional), "
            "report_condensed_timing_chains, evaluate (optional), then "
            "summarize_and_compress and plotting."
        )
    )
    parser.add_argument("--outdir", required=True, type=Path, help="Main output directory.")
    parser.add_argument("--vcfs", type=Path, default=None, help="VCF directory for source annotation.")
    parser.add_argument(
        "--cn-bed",
        type=Path,
        default=None,
        help="Wakhan integer_profile.bed containing HP1 and HP2 copy number.",
    )
    parser.add_argument("--tree", type=Path, default=None, help="Parent-child TSV for evaluation.")
    parser.add_argument("--dry-run", action="store_true", help="Print commands and exit.")
    args = parser.parse_args()

    has_vcfs = args.vcfs is not None
    has_tree = args.tree is not None
    has_cn = args.cn_bed is not None

    postprocess_dir = Path(__file__).resolve().parent
    evo_dir = postprocess_dir.parent
    annotate_script = postprocess_dir / "annotate_source_vcf.py"
    timing_script = postprocess_dir / "report_condensed_timing_chains.py"
    evaluate_script = postprocess_dir / "evaluate_graphs.py"
    eval_heatmap_script = postprocess_dir / "plot_edge_eval_vcfpair_heatmap.py"
    summarize_script = evo_dir / "summarize_and_compress.py"
    max_span_plot_script = postprocess_dir / "plot_max_phase_block_by_chr.py"
    interactive_plot_script = postprocess_dir / "plot_timing_cn_interactive.py"

    outdir = args.outdir.resolve()
    steps: list[tuple[str, list[str]]] = []

    if has_vcfs or has_cn:
        cmd = [sys.executable, str(annotate_script), "--main_out", str(outdir)]
        if has_vcfs:
            cmd.extend(["--vcfs", str(args.vcfs.resolve())])
        if has_cn:
            cmd.extend(["--cn-bed", str(args.cn_bed.resolve())])
        steps.append(("annotate_source_vcf", cmd))

    steps.append(("report_condensed_timing_chains", [sys.executable, str(timing_script), "--outdir", str(outdir)]))

    if has_tree:
        steps.append(
            (
                "evaluate_graphs",
                [sys.executable, str(evaluate_script), "--outdir", str(outdir), "--tree", str(args.tree.resolve())],
            )
        )
    if has_tree:
        steps.append(
            (
                "plot_edge_eval_vcfpair_heatmap",
                [
                    sys.executable,
                    str(eval_heatmap_script),
                    "--outdir",
                    str(outdir),
                    "--tree",
                    str(args.tree.resolve()),
                ],
            )
        )

    steps.append(("summarize_and_compress", [sys.executable, str(summarize_script), str(outdir)]))

    interactive_cmd = [sys.executable, str(interactive_plot_script), "--outdir", str(outdir)]
    if has_cn:
        interactive_cmd.extend(["--cn-bed", str(args.cn_bed.resolve())])
    steps.append(("plot_timing_cn_interactive", interactive_cmd))

    steps.append(
        (
            "plot_max_phase_block_by_chr",
            [sys.executable, str(max_span_plot_script), str(outdir / "global_summary.txt")],
        )
    )

    succeeded = 0
    for i, (name, argv) in enumerate(steps, start=1):
        print(f"Running ({i}/{len(steps)}): {name}")
        if args.dry_run:
            print(" ", " ".join(argv))
            continue

        result = subprocess.run(argv, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Error: {name} exited with code {result.returncode}", file=sys.stderr)
            err_tail = _tail_lines(result.stderr or "")
            out_tail = _tail_lines(result.stdout or "")
            if err_tail:
                print("--- stderr (tail) ---", file=sys.stderr)
                print(err_tail, file=sys.stderr)
            elif out_tail:
                print("--- stdout (tail) ---", file=sys.stderr)
                print(out_tail, file=sys.stderr)
            sys.exit(result.returncode)

        summary = _extract_summary(result.stdout or "")
        if summary:
            print(summary)
        succeeded += 1
        print(f"Completed: {name}")

    if args.dry_run:
        print("(dry-run; no commands executed)")
    else:
        print(f"Done. Completed {succeeded}/{len(steps)} steps successfully.")


if __name__ == "__main__":
    main()
