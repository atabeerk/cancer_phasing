# EvoPhase: Phasing tumor mutations with bulk long-read sequencing

EvoPhase infers local evolutionary relationships between somatic SNVs from their
co-occurrence on long-read tumor sequencing. It fits candidate mutation trees with
expectation maximization (EM), detects evolutionary relationship between each mutation pair and builds chromosome-level phase graphs from these
pairwise relationships.


## Inputs

Required:

- A coordinate-sorted, bgzip-compressed and indexed somatic `.vcf.gz` file.
  Sample fields (including GT) are ignored; read counts and VAFs come from the BAMs.
- One or more coordinate-sorted and indexed tumor `.bam` files aligned to the
  same reference, with matching sequence names and lengths. VCF chromosome names
  must match the BAMs exactly (for example, `chr1` and `1` do not match).

Optional:

- A sorted, bgzip-compressed and indexed phased germline `.vcf.gz` for haplotagging somatic mutations. `--germline-vcf-sample-name` is required with this input.
- One or more BED files defining genomic regions (and somatic mutations in these regions) to exclude.
- A Wakhan `integer_profile.bed` for haplotype copy-number annotations.


## Installation and build

Install `conda` or `mamba`, then clone the repository:

```bash
git clone https://github.com/atabeerk/cancer_phasing.git
cd cancer_phasing
```

Create and activate the environment from the repository root:

```bash
conda env create -n cancer-phasing -f environment.yml
conda activate cancer-phasing
```

Alternatively, use `mamba env create -n cancer-phasing -f environment.yml`.
The environment provides Python, a C++17 compiler, `make`, `samtools`,
`bcftools`, HTSlib, zlib, NumPy and pysam.

Build:

```bash
cd evo
make build
```

Keep the environment active when building and running the program.
The remaining commands below are run from `evo/`.

### Small end-to-end test

From `evo/`, with the environment active, run:

```bash
make test
```

This runs the bundled BAM, somatic VCF and phased germline VCF in `toy_data/`
through the full pipeline and compares results with `expected_output.tsv` in
that same folder. The example has five nearby mutations, 50 reads covering all
five, and cooccurring, ordered and divergent relationships. The command also
checks read counts, VAFs and germline haplotagging, and prints the output directory.
See `toy_data/README.md` for the input design and expected relationships.

From `evo/`, use `make help`, `make build`, `make test`, `make test-smoke` or `make clean`.

## Quickstart

```bash
python run_genome_by_chr.py \
  --somatic-vcf /path/to/somatic.vcf.gz \
  --bam /path/to/tumor.bam \
  --output-dir /path/to/new_output_dir \
  --jobs 4
```

A new (empty) output directory is required for each run.

To haplotag reads from phased germline markers:

```bash
python run_genome_by_chr.py \
  --somatic-vcf /path/to/somatic.vcf.gz \
  --bam /path/to/tumor.bam \
  --germline-vcf /path/to/phased_germline.vcf.gz \
  --germline-vcf-sample-name NORMAL \
  --output-dir /path/to/new_output_dir \
  --jobs 4
```

To retain divergent relationships only when both mutations are on the same germline haplotype assignments (i.e., both HP1 or both HP2), also add `--divergent-same-hp`. This option requires
`--germline-vcf`.

### Options

Run `python run_genome_by_chr.py --help` for the complete CLI.

| Option | Meaning | Default |
| --- | --- | --- |
| `--germline-vcf` | Phased germline VCF for read haplotagging. | None |
| `--germline-vcf-sample-name` | Germline sample; required with `--germline-vcf`. | None |
| `--jobs` | Number of chromosomes processed in parallel. | `1` |
| `--window-size` | Maximum mutations per fitting window; supported range `2`–`5`. | `4` |
| `--window-stride` | Number of mutations to advance between regular windows; positive integer. | `3` |
| `--divergent-same-hp` | Require matching HP1/HP2 assignments for divergent proposals; requires a germline VCF. | Off |
| `--exclude-regions-bed` | Exclude somatic mutations overlapping BED intervals; repeat for multiple files. | None |
| `--exclude-chrom` | Exclude chromosomes; repeat or use comma-separated names, with optional `chr` prefixes. | None |
| `--top-components` | Number of components in each genome-wide top-N graph. | `50` |
| `--cn-bed` | Annotate final outputs with haplotype copy numbers. | None |

For copy-number annotation, supply a tab-separated Wakhan `integer_profile.bed`
with these nine columns:

```text
chr  start  end  hp1_coverage  hp1_copynumber_state  hp1_confidence  hp2_coverage  hp2_copynumber_state  hp2_confidence
```

Copy numbers annotate final graph nodes and component spans. They do not adjust
VAFs or alter inferred relationships or accepted edges.


## Outputs

Paths below are relative to `--output-dir`. `<chrom>` is the chromosome name
from the BAM header.

| Output | Contents |
| --- | --- |
| `<chrom>_out/graphs/<chrom>.json` | Mutation graph in Cytoscape JSON format. |
| `<chrom>_out/graphs/<chrom>_condensed.json` | Graph with cooccurring mutations grouped into nodes. |
| `<chrom>_out/graphs/<chrom>_simplified.json` | Simplified graph for visualization. |
| `<chrom>_out/component_statistics.txt` | Component sizes, spans and mutation statistics; includes CN annotations when requested. |
| `<chrom>_out/inference/` | Pair results, fitted windows and summaries. `conflicting_pairs.tsv` is produced only when there are pairs covered by multiple windows and contains rows only for conflicts. |
| `global_summary.txt` | Per-chromosome component summaries, including chromosomes with no components. |
| `topN_graphs/` | Genome-wide graphs of the largest selected components. |
| `timing_chains.tsv` | Complete root-to-leaf timing paths with at least three condensed nodes. |
| `haplotagged_snvs.tsv` and `<chrom>_out/haplotagged_snvs.tsv` | Retained SNVs, read counts, haplotype assignments and haplotype-specific ALT fractions. **Produced only when `--germline-vcf` is provided.** Assignments can remain UNKNOWN. |
| `logs/` | Run and chromosome logs, mutation/edge summaries, and Python runtime measurements. |
| `<output_directory_name>.tar.gz` | Final graphs, inference results, summaries, logs and optional haplotagged/CN tables for the selected chromosomes, prepared for export.  |

Phase graphs produced in JSON format can be visualized with [Cytoscape](https://cytoscape.org).
