#!/usr/bin/env python3
import os
import sys
import subprocess
import glob

def run_cmd(cmd):
    """Run shell command and stream its output."""
    print(f"\n[Running] {' '.join(cmd)}")
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    for line in process.stdout:
        print(line, end="")
    process.wait()
    if process.returncode != 0:
        raise RuntimeError(f"Command failed ({process.returncode}): {' '.join(cmd)}")

def get_chromosomes_from_bam(bam_path):
    """Return (chrom, length) pairs extracted via `samtools idxstats`."""
    cmd = ["samtools", "idxstats", bam_path]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    chromosomes = []
    for line in result.stdout.strip().split("\n"):
        fields = line.split("\t")
        if len(fields) < 2 or fields[0] == "*" or not fields[0]:
            continue
        chrom, length = fields[0], int(fields[1])
        chromosomes.append((chrom, length))
    return chromosomes

def main():
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <input_vcf.gz> <haplotagged_bam> <output_dir>")
        sys.exit(1)

    vcf_path = sys.argv[1]
    bam_path = sys.argv[2]
    outdir = sys.argv[3]

    preprocess_script_dir = "/Users/donmeza2/Documents/cancer_phasing/data/cancer_phasing_data"
    preprocess_script = "run_phasing_pipeline.sh"
    main_program = "./snv_occur"

    # --- Validate inputs ---
    for f in [vcf_path, bam_path, main_program]:
        if not os.path.exists(f):
            sys.exit(f"Error: required file not found: {f}")

    os.makedirs(outdir, exist_ok=True)

    # --- Derive genome name ---
    genome_name = os.path.basename(bam_path).replace(".bam", "")
    print(f"Genome name: {genome_name}")

    # --- Get chromosomes from BAM via samtools ---
    print("Extracting chromosome list from BAM header using samtools...")
    chromosomes = get_chromosomes_from_bam(bam_path)
    print(f"Detected {len(chromosomes)} chromosomes.")

    # --- Process each chromosome ---
    for chrom, length in chromosomes:
        region = f"{chrom}:1-{length}"

        pre_dir = os.path.join(outdir, f"{genome_name}_{chrom}_pre")
        result_dir = os.path.join(outdir, f"{genome_name}_{chrom}_results")
        os.makedirs(pre_dir, exist_ok=True)
        os.makedirs(result_dir, exist_ok=True)

        print(f"\n=== Processing {chrom} ({region}) ===")
        print(f"Preprocessing output → {pre_dir}")
        print(f"Main program output → {result_dir}")

        # Step 1: Run preprocessing
        cwd = os.getcwd()
        os.chdir(preprocess_script_dir)
        run_cmd(["bash", preprocess_script, vcf_path, bam_path, region, pre_dir])
        os.chdir(cwd)

        # Step 2: Locate merged BAM + VCF
        merged_bam_list = glob.glob(os.path.join(pre_dir, "*.merged.bam"))
        merged_vcf_list = glob.glob(os.path.join(pre_dir, "*.merged.vcf.gz"))

        if not merged_bam_list:
            raise FileNotFoundError(f"No merged BAM found in {pre_dir}")
        if not merged_vcf_list:
            raise FileNotFoundError(f"No merged VCF found in {pre_dir}")

        merged_bam = merged_bam_list[0]
        merged_vcf = merged_vcf_list[0]

        print(f"Found merged BAM: {os.path.basename(merged_bam)}")
        print(f"Found merged VCF: {os.path.basename(merged_vcf)}")

        # Step 3: Run main program
        run_cmd([main_program, merged_bam, merged_vcf, result_dir])

        print(f"=== Done {chrom} ===")

    print(f"\nAll chromosomes processed for {genome_name}. Results saved in {outdir}")

if __name__ == "__main__":
    main()
