#!/usr/bin/env python3

import argparse
import os
import re
import shutil
import subprocess
import sys


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
        help="Input somatic VCF.gz; sites are candidates and sample FORMAT fields are unused.",
    )
    parser.add_argument(
        "--bam",
        dest="bams",
        action="append",
        required=True,
        help="Input tumor BAM; existing HP tags are ignored. Repeat --bam for multiple BAMs.",
    )
    parser.add_argument("--output-dir", dest="output_dir", required=True, help="Output directory")
    parser.add_argument(
        "--germline-vcf",
        dest="germline_vcf_gz",
        default=None,
        help="Optional phased germline VCF.gz for read haplotagging; otherwise reads are unassigned.",
    )
    parser.add_argument(
        "--germline-vcf-sample-name",
        dest="germline_vcf_sample_name",
        default="",
        help="Required sample name when --germline-vcf is provided.",
    )
    parser.add_argument(
        "--exclude-regions-bed",
        action="append",
        default=[],
        help=(
            "BED file whose intervals are excluded from the somatic mutation set. "
            "Repeat this option to exclude the union of multiple BED files."
        ),
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="Number of chromosomes to process in parallel (default: 1)",
    )
    parser.add_argument(
        "--top-components",
        type=int,
        default=50,
        help="Number of components in each genome-wide top-N graph (default: 50)",
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
            "HP1 or both HP2. Requires --germline-vcf."
        ),
    )
    parser.add_argument(
        "--window-size", type=int, choices=range(2, 6), default=4,
        help="Maximum mutations per shared-tree window (2–5; default: 4).",
    )
    parser.add_argument(
        "--window-stride", type=int, default=3,
        help="Number of mutations to advance between windows (default: 3).",
    )
    # Optional copy-number labels for final outputs.
    parser.add_argument(
        "--cn-bed",
        type=str,
        default=None,
        help="Annotate outputs with HP1/HP2 copy numbers from Wakhan integer_profile.bed; VAFs are unchanged.",
    )

    args = parser.parse_args()
    if args.window_stride < 1:
        parser.error("--window-stride must be >= 1")
    if args.jobs <= 0:
        parser.error("--jobs must be >= 1")
    if args.top_components <= 0:
        parser.error("--top-components must be >= 1")
    if args.divergent_same_hp and not args.germline_vcf_gz:
        parser.error("--divergent-same-hp requires --germline-vcf.")
    if args.germline_vcf_gz and not args.germline_vcf_sample_name:
        parser.error("--germline-vcf-sample-name is required with --germline-vcf.")
    if args.germline_vcf_sample_name and not args.germline_vcf_gz:
        parser.error("--germline-vcf-sample-name requires --germline-vcf.")
    return args


def parse_and_validate_args(main_program: str = "./main") -> argparse.Namespace:
    """Return normalized arguments after all startup validation succeeds."""
    args = parse_args()
    args.main_program = main_program
    args.output_dir = os.path.abspath(args.output_dir)
    if args.cn_bed:
        args.cn_bed = os.path.abspath(args.cn_bed)
    args.excluded_chroms = parse_exclude_chroms(args.exclude_chrom)
    try:
        preflight_requirements(main_program, args)
    except RuntimeError as exc:
        sys.exit(f"Error: {exc}")
    return args


def get_chromosomes_from_bam(bam_path):
    """Return (chrom, length) pairs for standard chromosomes only (chr or no-chr)."""
    cmd = ["samtools", "idxstats", bam_path]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)

    chromosomes = []
    for line in result.stdout.strip().split("\n"):
        fields = line.split("\t")
        if len(fields) < 2:
            continue

        chrom = fields[0]
        length = int(fields[1])

        # Normalize: drop optional "chr" prefix for checking
        c = chrom[3:] if chrom.startswith("chr") else chrom

        # Keep only 1-22, X, Y, M, MT
        if re.fullmatch(r"(?:[1-9]|1[0-9]|2[0-2]|X|Y|M|MT)", c) is None:
            print(f"Skipping non-standard contig {chrom}")
            continue

        chromosomes.append((chrom, length))

    return chromosomes


def get_vcf_samples(vcf_path):
    """Return sample names from a VCF header in declared order."""
    result = subprocess.run(
        ["bcftools", "query", "-l", vcf_path],
        capture_output=True,
        text=True,
        check=True,
    )
    return [line.strip() for line in result.stdout.splitlines() if line.strip()]


def validate_vcf_chromosomes(vcf_path, chromosomes, excluded_chroms):
    """Match populated, processed VCF contigs to BAM names without requiring every contig."""
    result = subprocess.run(
        ["bcftools", "index", "--stats", vcf_path],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"Could not read the VCF index for {vcf_path}: {result.stderr.strip()}"
        )
    bam_names = {chrom for chrom, _ in chromosomes}
    for line in result.stdout.splitlines():
        chrom = line.split("\t", 1)[0]
        normalized = normalize_chrom_name(chrom)
        if chrom_is_excluded(chrom, excluded_chroms):
            continue
        if re.fullmatch(r"(?:[1-9]|1[0-9]|2[0-2]|X|Y|M|MT)", normalized) is None:
            continue
        if chrom not in bam_names:
            raise RuntimeError(
                f"VCF chromosome '{chrom}' in {vcf_path} does not match any BAM sequence name. "
                "VCF and BAM chromosome names must match exactly, including any 'chr' prefix."
            )


def validate_exclusion_bed(path):
    """Validate BED intervals once, before launching chromosome workers."""
    with open(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.split()
            try:
                if any(re.fullmatch(r"[+-]?[0-9]+", value) is None for value in fields[1:3]):
                    raise ValueError
                start, end = int(fields[1]), int(fields[2])
                if start < 0 or end <= start or end > 2**63 - 1:
                    raise ValueError
            except (IndexError, ValueError):
                raise RuntimeError(f"invalid BED interval at {path}:{line_number}") from None


def preflight_requirements(main_program: str, args: argparse.Namespace) -> None:
    """Validate dependencies, input files, samples, and chromosome selection."""
    if os.path.lexists(args.output_dir):
        if not os.path.isdir(args.output_dir):
            raise RuntimeError(f"output path is not a directory: {args.output_dir}")
        with os.scandir(args.output_dir) as entries:
            if next(entries, None) is not None:
                raise RuntimeError(f"output directory must be empty: {args.output_dir}")

    required_bins = ("samtools", "bcftools")
    missing = [b for b in required_bins if shutil.which(b) is None]
    if missing:
        raise RuntimeError(
            "Missing required executables on PATH: " + ", ".join(missing)
            + ". Activate the conda environment from the repository-root environment.yml."
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

    required_files = [args.input_vcf_gz, *args.bams]
    if args.germline_vcf_gz:
        required_files.append(args.germline_vcf_gz)
    if args.cn_bed:
        required_files.append(args.cn_bed)
    if args.exclude_regions_bed:
        required_files.extend(args.exclude_regions_bed)
    for path in required_files:
        if not os.path.exists(path):
            raise RuntimeError(f"required file not found: {path}")

    for path in args.exclude_regions_bed or []:
        validate_exclusion_bed(path)

    if args.germline_vcf_gz:
        germline_samples = get_vcf_samples(args.germline_vcf_gz)
        if args.germline_vcf_sample_name not in germline_samples:
            raise RuntimeError(
                "germline VCF sample "
                f"'{args.germline_vcf_sample_name}' was not found."
            )

    print("Extracting chromosome list from BAM header using samtools...")
    args.detected_chromosomes = get_chromosomes_from_bam(args.bams[0])
    for bam_path in args.bams[1:]:
        if get_chromosomes_from_bam(bam_path) != args.detected_chromosomes:
            raise RuntimeError(
                "BAM sequence names and lengths do not match the first "
                f"input BAM: {bam_path}"
            )
    args.chromosomes = [
        (chrom, length)
        for chrom, length in args.detected_chromosomes
        if not chrom_is_excluded(chrom, args.excluded_chroms)
    ]
    args.excluded_detected_chromosomes = [
        (chrom, length)
        for chrom, length in args.detected_chromosomes
        if chrom_is_excluded(chrom, args.excluded_chroms)
    ]
    print(f"Detected {len(args.detected_chromosomes)} chromosomes.")
    if args.excluded_chroms:
        requested = ", ".join(sorted(args.excluded_chroms))
        matched = ", ".join(chrom for chrom, _ in args.excluded_detected_chromosomes) or "none"
        print(f"Excluded chromosomes requested: {requested}")
        print(f"Excluded chromosomes found in BAM: {matched}")
    if not args.chromosomes:
        raise RuntimeError("no chromosomes left after applying --exclude-chrom.")
    for vcf_path in (args.input_vcf_gz, args.germline_vcf_gz):
        if vcf_path:
            validate_vcf_chromosomes(vcf_path, args.chromosomes, args.excluded_chroms)
