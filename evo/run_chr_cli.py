#!/usr/bin/env python3

import argparse
import os
import shutil
import subprocess


def normalize_chrom_name(chrom):
    chrom = str(chrom).strip()
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    return chrom.upper()


def parse_exclude_chroms(values):
    excluded = set()
    for value in values or []:
        for chrom in value.split(","):
            chrom = chrom.strip()
            if chrom:
                excluded.add(normalize_chrom_name(chrom))
    return excluded


def chrom_is_excluded(chrom, excluded_chroms):
    return normalize_chrom_name(chrom) in excluded_chroms


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Run evo pipeline per chromosome, aggregate genome-level stats, "
            "and write run/chromosome logs."
        )
    )
    parser.add_argument(
        "--somatic-vcf",
        dest="input_vcf_gz",
        required=True,
        help="Input somatic VCF.gz",
    )
    parser.add_argument(
        "--bam",
        dest="haplotagged_bams",
        action="append",
        required=True,
        help="Input tumor BAM. Repeat --bam to provide multiple BAMs.",
    )
    parser.add_argument("--output-dir", dest="output_dir", required=True, help="Output directory")
    parser.add_argument(
        "--somatic-sample-name",
        dest="vcf_sample_name",
        default="",
        help="Somatic VCF sample name; required when the VCF has multiple samples.",
    )
    parser.add_argument(
        "--germline-vcf",
        dest="germline_vcf_gz",
        default=None,
        help="Optional phased germline VCF.gz used to infer read haplotypes.",
    )
    parser.add_argument(
        "--germline-vcf-sample-name",
        dest="germline_vcf_sample_name",
        default="",
        help="Required sample name when --germline-vcf is provided.",
    )
    parser.add_argument(
        "--exclude-regions-bed",
        default=None,
        help=(
            "Optional BED file whose intervals are excluded from the somatic "
            "mutation set."
        ),
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="Number of chromosomes to process in parallel (default: 1)",
    )
    parser.add_argument(
        "--exclude-chrom",
        action="append",
        default=[],
        help=(
            "Chromosome to exclude from processing. Can be repeated or comma-separated; "
            "chr prefixes are optional, e.g. --exclude-chrom chr7 --exclude-chrom 9."
        ),
    )
    parser.add_argument(
        "--divergent-same-hp",
        action="store_true",
        help=(
            "Only retain divergent relationships whose two mutations are both "
            "HP1 or both HP2."
        ),
    )
    # Optional pass-through arguments for postprocessing.
    parser.add_argument("--vcfs", type=str, default=None, help="Postprocess: VCF directory for source annotation.")
    parser.add_argument(
        "--cn-bed",
        type=str,
        default=None,
        help="Postprocess: Wakhan integer_profile.bed containing HP1 and HP2 copy number.",
    )
    parser.add_argument("--tree", type=str, default=None, help="Postprocess: parent-child TSV for evaluation.")

    args = parser.parse_args()
    if args.jobs <= 0:
        parser.error("--jobs must be >= 1")
    if args.germline_vcf_gz and not args.germline_vcf_sample_name:
        parser.error("--germline-vcf-sample-name is required with --germline-vcf.")
    if args.germline_vcf_sample_name and not args.germline_vcf_gz:
        parser.error("--germline-vcf-sample-name requires --germline-vcf.")
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
