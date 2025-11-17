#!/usr/bin/env python3
import os
import sys
import subprocess

def run_cmd(cmd):
    """Run a shell command and stream its output."""
    print(f"\n[Running] {' '.join(cmd)}")
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    for line in process.stdout:
        print(line, end="")
    process.wait()
    if process.returncode != 0:
        raise RuntimeError(f"Command failed ({process.returncode}): {' '.join(cmd)}")

def get_chromosomes_from_bam(bam_path):
    """Return (chrom, length) pairs for standard chromosomes only."""
    import subprocess

    # List of standard chromosomes
    standard_chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]

    # Get BAM idxstats
    cmd = ["samtools", "idxstats", bam_path]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)

    chromosomes = []
    for line in result.stdout.strip().split("\n"):
        fields = line.split("\t")
        if len(fields) < 2:
            continue
        chrom, length = fields[0], int(fields[1])
        if chrom not in standard_chroms:
            print(f"Skipping non-standard contig {chrom}")
            continue
        chromosomes.append((chrom, length))

    return chromosomes


def main():
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <input_vcf.gz> <haplotagged_bam> <output_dir>")
        sys.exit(1)

    vcf_path = sys.argv[1]
    bam_path = sys.argv[2]
    outdir = os.path.abspath(sys.argv[3])

    main_program = "./snv_occur"

    # Validate inputs
    for f in [vcf_path, bam_path, main_program]:
        if not os.path.exists(f):
            sys.exit(f"Error: required file not found: {f}")

    os.makedirs(outdir, exist_ok=True)

    # Derive genome name
    genome_name = os.path.basename(bam_path).replace(".bam", "")
    print(f"Genome name: {genome_name}")

    # Get chromosomes
    print("Extracting chromosome list from BAM header using samtools...")
    chromosomes = get_chromosomes_from_bam(bam_path)
    print(f"Detected {len(chromosomes)} chromosomes.")

    for chrom, length in chromosomes:
        region = f"{chrom}:1-{length}"

        pre_dir = os.path.join(outdir, f"{chrom}_pre")
        out_dir = os.path.join(outdir, f"{chrom}_out")
        os.makedirs(pre_dir, exist_ok=True)
        os.makedirs(out_dir, exist_ok=True)

        print(f"\n=== Processing {chrom} ({region}) ===")
        print(f"Preprocessing directory: {pre_dir}")
        print(f"Main program output directory: {out_dir}")

        # --- Step 1: Split BAM for this chromosome with filters ---
        chrom_bam = os.path.join(pre_dir, f"{chrom}.bam")
        # Exclude unmapped, secondary, supplementary reads; MAPQ < 20
        run_cmd([
            "samtools", "view", "-b", "-F", "2308", "-q", "20",
            bam_path, region, "-o", chrom_bam
        ])
        run_cmd(["samtools", "index", chrom_bam])

        # --- Step 2: Split VCF for this chromosome with PASS filter ---
        chrom_vcf = os.path.join(pre_dir, f"{chrom}.vcf.gz")
        run_cmd([
            "bcftools", "view", "-Oz", "-r", chrom, "-f", "PASS",
            vcf_path, "-o", chrom_vcf
        ])
        run_cmd(["bcftools", "index", chrom_vcf])

        # --- Step 3: Run main program ---
        run_cmd([main_program, chrom_vcf, chrom_bam, out_dir])

        print(f"=== Done {chrom} ===")

    print(f"\nAll chromosomes processed for {genome_name}. Results saved in {outdir}")

if __name__ == "__main__":
    main()
