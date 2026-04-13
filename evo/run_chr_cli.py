#!/usr/bin/env python3

import argparse
import os
import shutil
import subprocess


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Run evo pipeline per chromosome, aggregate genome-level stats, "
            "and write run/chromosome logs."
        )
    )
    parser.add_argument("--vcf", dest="input_vcf_gz", required=True, help="Input VCF.gz")
    parser.add_argument("--bam", dest="haplotagged_bam", required=True, help="Input haplotagged BAM")
    parser.add_argument("--output-dir", dest="output_dir", required=True, help="Output directory")
    parser.add_argument("--vcf-sample-name", dest="vcf_sample_name", default="", help="Optional VCF sample name")
    parser.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="Number of chromosomes to process in parallel (default: 1)",
    )
    # Optional pass-through arguments for postprocessing.
    parser.add_argument("--vcfs", type=str, default=None, help="Postprocess: VCF directory for source annotation.")
    parser.add_argument("--cn-bed-hp1", type=str, default=None, help="Postprocess: HP1 copy-number BED.")
    parser.add_argument("--cn-bed-hp2", type=str, default=None, help="Postprocess: HP2 copy-number BED.")
    parser.add_argument("--tree", type=str, default=None, help="Postprocess: parent-child TSV for evaluation.")

    args = parser.parse_args()
    if args.jobs <= 0:
        parser.error("--jobs must be >= 1")
    if (args.cn_bed_hp1 is None) != (args.cn_bed_hp2 is None):
        parser.error("Provide both --cn-bed-hp1 and --cn-bed-hp2 together.")
    return args


def preflight_requirements(main_program: str) -> None:
    """Fail fast on missing executables before starting any chromosome work."""
    required_bins = ("samtools", "bcftools")
    missing = [b for b in required_bins if shutil.which(b) is None]
    if missing:
        raise RuntimeError(
            "Missing required executables on PATH: " + ", ".join(missing)
            + ". Activate the conda environment from evo/environment.yml."
        )

    if not os.path.exists(main_program):
        raise RuntimeError(f"Required program not found: {main_program}")
    probe = subprocess.run(
        [main_program, "--help"],
        capture_output=True,
        text=True,
    )
    if probe.returncode != 0:
        detail = (probe.stderr or probe.stdout or "").strip()
        raise RuntimeError(
            f"Failed to execute {main_program} --help (exit {probe.returncode}). "
            f"This can indicate missing runtime libraries (e.g., htslib).\n{detail}"
        )
