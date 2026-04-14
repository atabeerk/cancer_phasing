# Phasing Tumor Clones By Timing Somatic Mutations In Bulk Long-Read Sequencing

This repo contains an actively developed tool for phasing somatic SNVs into clonal haplotypes based on their co-occurrence patterns in the reads.
Required inputs are:

1. A `.vcf.gz` file containing somatic only variants.
2. A `.bam` file where the tumor reads aligned against a reference (e.g., grch38)

The program will produce the following outputs under the specified directory:

1. Condensed and uncondensed mutation graphs for each region (10 Mb by default) in `.json` format under `chr*_out/graphs` folders. We recommend [Cytoscape](https://cytoscape.org) for visualizating these files.
2. `global_summary.txt` will contain genome-level summary of connected components in the mutation graph.
3. `logs/` directory contains the run logs and additional summary of detected mutations and their relationships.
4. Additional interactive plots in the `.html` format based on provided optional parameters (in progress).

## Prerequisites

- `conda` or `mamba`
- A C++17-capable compiler

## Installation

Clone the repository and enter it:

```bash
git clone https://github.com/atabeerk/cancer_phasing.git
cd cancer_phasing/evo
```

Create the environment with either tool:

```bash
conda env create -n cancer-phasing -f environment.yml
# or
mamba env create -n cancer-phasing -f environment.yml
```

Activate the environment:

```bash
conda activate cancer-phasing
```

## Build

```bash
cd evo
make build
```

Dependency checks will be run by the program and fail early at startup before processing starts:

- `run_genome_by_chr.py` checks required executables (`samtools`, `bcftools`) and tries `./main --help`

## Quickstart Usage

* ### Run chromosome-parallel wrapper (Recommended)

This wrapper does the required preprocessing, runs the main program with each chromosome in parallel (set `--jobs 1` for serial) and applies available postprocessing. 

```bash
python run_genome_by_chr.py \
  --vcf /path/to/input.vcf.gz \
  --bam /path/to/input.haplotagged.bam \
  --output-dir /path/to/output_dir \
  --jobs 4
```

* ### Alternatively, run `evo/main` directly (not recommended)

```bash
./main \
  --vcf /path/to/input.vcf.gz \
  --bam /path/to/input.haplotagged.bam \
  --output-dir /path/to/output_dir \
```

## Troubleshooting

- `Error: required executable not found on PATH: samtools` (from `./main`) or missing `samtools/bcftools` (from `run_genome_by_chr.py`):
  - Ensure the environment is active: `conda activate cancer-phasing`
- `Failed to execute ./main --help ... missing runtime libraries (e.g., htslib)` (from `run_genome_by_chr.py` preflight):
  - Reinstall environment from `environment.yml`
  - Rebuild `main` after activating the environment: `cd evo && make build`
- Compile-time errors about `htslib/vcf.h` or `-lhts`:
  - Activate the conda/mamba environment first so `CONDA_PREFIX` is set
  - Or set `HTSLIB_PREFIX` manually before running `make build`

## Make Targets

From `evo/`:

- `make help`
- `make build`
- `make test-smoke`
- `make clean`
