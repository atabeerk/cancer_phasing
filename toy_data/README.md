# Toy end-to-end test

This folder contains all input files and `expected_output.tsv`. From `evo/`,
with the normal EvoPhase environment active, run:

```bash
make test
```

The command runs the complete pipeline with germline haplotagging,
`--divergent-same-hp`, window size 4 and stride 3. It compares all ten accepted
pairs with the expected table, including relationship, direction, four allele
counts, coverage and VAFs. It also checks the final graph, haplotagged SNV counts
and absence of rejected edges. Output is retained in a printed temporary directory.
To choose a new or empty output directory, run from the repository root:

```bash
python tests/run_toy_test.py --output-dir /path/to/new_toy_output
```

The synthetic reference is `chr1` (1,000 bases of A). Fifty primary reads each
align at position 1 with CIGAR `250M`, MAPQ 60 and base quality 40. There are no
BAM haplotype tags. All reads cover the five A>G somatic sites at positions
100, 120, 140, 160 and 180, called A, B, C, D and E below.

| ALT pattern ABCDE | Reads |
| --- | ---: |
| 00000 | 25 |
| 11000 | 5 |
| 11100 | 5 |
| 11110 | 10 |
| 11001 | 5 |

A and B cooccur. They precede C, D and E; C precedes D; E diverges from C and D.
This gives one cooccurring, seven ordered (`timing`) and two divergent pairs.
The VAFs are 50%, 50%, 30%, 20% and 10%. Expected counts and relationships were
calculated from this table, independently of the program output.

Two phased A>G germline markers at positions 20 and 40 have sample `NORMAL`
and GT `1|0`. All 25 reads carrying any somatic ALT carry both germline ALTs
(HP1); the other 25 carry both reference alleles (HP2). All five somatic sites
should therefore be assigned HP1. The somatic VCF is sites-only.

The bundled files are ready to run. To regenerate the inputs with pysam, run
`python tests/generate_toy_fixture.py` from the repository root. This does not
regenerate or overwrite the expected output. This small, error-free example
checks the pipeline wiring and basic inference; it is not an accuracy benchmark.
