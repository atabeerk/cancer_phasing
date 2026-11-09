# Production shared-tree inference

The chromosome worker invokes `infer_evo_relationships.py`. It calls the model
through Python functions and passes the resulting proposals to the standard
`evo/graph_ops` builder.

| Module | Responsibility |
| --- | --- |
| `infer_evo_relationships.py` | Load retained mutations and annotations, run fitting, apply support/HP gates, and emit relationship files. |
| `read_observations.py` | Read prepared chromosome pileups, count eligible pairs, and project observations into windows. |
| `shared_tree_model.py` | Select connected windows, enumerate local trees, and compare fitted BIC scores. |
| `joint_clone_model.py` | Represent local genotypes, calculate likelihoods, and fit each candidate tree's mixture using EM. |
| `fit_window_trees.py` | Coordinate chromosome fitting and write pair results, fitted windows, and conflict diagnostics. |
| `relationship_calls.py` | Convert pair rows into calls and implement the required support gate and BIC-margin ordering. |

Run the pipeline through `evo/run_genome_by_chr.py`. These six modules are
production code. The model libraries have no
standalone fitting or benchmark CLI, truth input, or truth-scoring code.

Windows default to at most four mutations, advancing three mutations at a time.
Use `--window-size` (2–5) and `--window-stride` (positive integer) to override.
Connected supplemental windows cover eligible pairs missed by the regular slide.
Each pair is assigned one covering window deterministically before fitting.

Each prepared BAM alignment record is an observation; primary and supplementary
records remain separate. Reads covering only one window mutation also contribute
to the likelihood. HP annotations are excluded from fitting. The optional
`--divergent-same-hp` gate retains divergent proposals only when both mutations
have matching HP1 or HP2 assignments.

Candidate eligibility is `AA + AR + RA >= 3 OR AA >= 2`. Before graph construction,
proposals must have at least one observation of every defining allele pattern:

| Relationship | Required patterns |
| --- | --- |
| Cooccurring | AA |
| Left before right | AA and AR |
| Right before left | AA and RA |
| Divergent | AR and RA |

The fixed-tree fit uses at most 50 EM iterations, a 1% base-error floor, and a
1% outlier component. In `groups.tsv`, `best_converged` indicates whether the
selected tree met the existing EM stopping tolerance (1) or reached the iteration
limit (0). `em_converged` and `em_iteration_limited` count completed candidate
fits; `em_pruned` counts fits stopped early by a safe BIC bound and is excluded
from those two counts. Reaching the iteration limit does not discard a fit.
Pair `prob_*` fields are relative BIC weights. Exported
`RELIABILITY` is half the gap between the two lowest relationship BICs and is
used for graph ordering; it is not a calibrated correctness probability.

`pair_results.tsv`, `groups.tsv`, and `summary.json` preserve the production
inference results and diagnostics. `conflicting_pairs.tsv` reports each covering
window's relationship and BIC scores for conflicting pairs. Conflict reporting
is disabled when no eligible pair appears together in multiple windows.

Tests named `test_*.py` exercise the production code and are not runtime inputs.
