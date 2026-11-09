"""Fit chromosome mutation windows and save pair calls and conflict diagnostics."""
import csv
import json
import time
from collections import Counter, defaultdict
from itertools import combinations

from read_observations import (
    candidate_key, candidates_from_observations, project_reads,
    read_observations, passes_inference_gate,
)
from joint_clone_model import inferred_relation
from shared_tree_model import (
    RELATIONS, assign_pairs_to_windows, build_windows,
    fit_shared_tree_relations, probabilities_from_bics,
)


def annotate_window_relationships(pair_rows, pair_opinions):
    """Flag disagreements without using the opinions to change pair selection."""
    overlap = Counter()
    multiplicity = Counter()
    conflicts = []
    for key, row in pair_rows.items():
        opinions = pair_opinions[key]
        assigned = [opinion for opinion in opinions if opinion["is_assigned"]]
        if len(assigned) != 1:
            raise RuntimeError(f"pair must have exactly one assigned-window opinion: {key}")
        if assigned[0]["relationship"] != row["best_pair_relationship"]:
            raise RuntimeError(f"assigned-window relationship mismatch: {key}")
        counts = Counter(opinion["relationship"] for opinion in opinions)
        count = len(opinions)
        conflict = len(counts) > 1
        row.update({
            "covering_group_count": count,
            "covering_group_ids": ",".join(str(opinion["group_id"]) for opinion in opinions),
            "distinct_group_relationships": len(counts),
            "covering_relationship_counts": json.dumps(counts, sort_keys=True),
            "group_relationship_conflict": int(conflict),
        })
        multiplicity[count] += 1
        overlap["pairs"] += 1
        overlap["single_coverage_pairs" if count == 1 else "multiple_coverage_pairs"] += 1
        if conflict:
            conflicts.append(key)
            overlap["conflicting_pairs"] += 1
        elif count > 1:
            overlap["unanimous_multiple_pairs"] += 1
    summary = {field: overlap[field] for field in (
        "pairs", "single_coverage_pairs", "multiple_coverage_pairs",
        "unanimous_multiple_pairs", "conflicting_pairs",
    )}
    summary["coverage_multiplicity"] = {str(count): pairs for count, pairs in sorted(multiplicity.items())}
    return summary, sorted(conflicts)


def write_window_conflicts(output_dir, pair_rows, pair_opinions, conflict_keys, overlap):
    """For conflicting pairs, report each window's relationship and BIC scores."""
    if not overlap["multiple_coverage_pairs"]:
        return None
    filename = "conflicting_pairs.tsv"
    fields = [
        "chrom", "left_pos", "right_pos", "assigned_group_id", "group_id",
        "is_assigned", "relationship", "group_best_bic",
    ] + [f"bic_{relation}" for relation in RELATIONS]
    with (output_dir / filename).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for key in conflict_keys:
            row = pair_rows[key]
            common = {field: row[field] for field in (
                "chrom", "left_pos", "right_pos", "assigned_group_id",
            )}
            for opinion in pair_opinions[key]:
                writer.writerow({**common, **opinion})
    return filename


def fit_window_trees(mutations, args):
    """Fit one chromosome, save diagnostics, and return its pair-result rows."""
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise ValueError(f"output directory is not empty: {args.output_dir}")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    observations, reads_by_position, observation_stats = read_observations(
        args.pileup_dir, mutations, args.evo_pileup_bam_count,
    )
    candidate_pool = candidates_from_observations(
        observations,
        min(args.minimum_informative_reads, args.minimum_alt_alt_reads),
    )
    candidates = [
        candidate for candidate in candidate_pool
        if passes_inference_gate(
            candidate,
            args.minimum_informative_reads,
            args.minimum_alt_alt_reads,
        )
    ]
    windows, component_sizes = build_windows(
        candidates, args.window_size, args.window_stride,
    )
    assignments = assign_pairs_to_windows(candidates, windows)
    candidate_by_pair = {
        ((candidate.chrom, candidate.left), (candidate.chrom, candidate.right)):
        candidate for candidate in candidates
    }
    assignments_by_window = defaultdict(list)
    for pair, window_index in assignments.items():
        assignments_by_window[window_index].append(pair)
    # Preserve saved column names for existing result readers.
    pair_fields = [
        "chrom", "left_pos", "right_pos", "alt_alt", "alt_ref", "ref_alt",
        "ref_ref", "informative_reads", "spanning_reads",
        "best_pair_relationship", "assigned_group_id", "group_size",
        "group_positions", "group_best_bic", "covering_group_count",
        "covering_group_ids", "distinct_group_relationships",
        "covering_relationship_counts", "group_relationship_conflict",
    ] + [f"bic_{relation}" for relation in RELATIONS] + [
        f"prob_{relation}" for relation in RELATIONS
    ]
    window_fields = [
        "group_id", "chrom", "positions", "mutation_count", "physical_reads",
        "unique_reads", "eligible_pairs", "assigned_pairs", "best_bic",
        "best_groups", "best_parents", "trees", "bound_pruned", "em_pruned",
        "em_fits", "em_iterations", "em_converged", "em_iteration_limited",
        "best_converged", "elapsed_seconds",
    ]
    pair_rows = {}
    pair_opinions = defaultdict(list)
    total_fit_seconds = 0.0
    with (args.output_dir / "groups.tsv").open("w", newline="") as window_handle:
        window_writer = csv.DictWriter(
            window_handle, fieldnames=window_fields, delimiter="\t"
        )
        window_writer.writeheader()
        for window_index, window in enumerate(windows):
            window_started = time.perf_counter()
            local_index = {position: index for index, position in enumerate(window)}
            eligible_pairs = [
                pair for pair in combinations(window, 2)
                if pair in candidate_by_pair
            ]
            pair_indices = [
                (local_index[left], local_index[right])
                for left, right in eligible_pairs
            ]
            reads = project_reads(window, observations, reads_by_position)
            best, relation_bics, stats = fit_shared_tree_relations(
                reads, len(window), pair_indices, args.em_iterations,
            )
            elapsed = time.perf_counter() - window_started
            total_fit_seconds += elapsed
            assigned = set(assignments_by_window.get(window_index, ()))
            window_positions = ",".join(str(position[1]) for position in window)
            best_groups = json.dumps(best.state.groups)
            best_parents = json.dumps(best.state.parents)
            for pair, indices, bics in zip(
                eligible_pairs, pair_indices, relation_bics
            ):
                candidate = candidate_by_pair[pair]
                best_relation = inferred_relation(best.state, *indices)
                bic_fields = {f"bic_{relation}": bics[relation] for relation in RELATIONS}
                opinion = {
                    "group_id": window_index,
                    "is_assigned": int(pair in assigned),
                    "relationship": best_relation,
                    "group_best_bic": best.bic,
                }
                opinion.update(bic_fields)
                pair_opinions[candidate_key(candidate)].append(opinion)
                if pair not in assigned:
                    continue
                probabilities = probabilities_from_bics(bics)
                row = {
                    "chrom": candidate.chrom,
                    "left_pos": candidate.left,
                    "right_pos": candidate.right,
                    "alt_alt": candidate.alt_alt,
                    "alt_ref": candidate.alt_ref,
                    "ref_alt": candidate.ref_alt,
                    "ref_ref": candidate.ref_ref,
                    "informative_reads": candidate.informative_reads,
                    "spanning_reads": candidate.total,
                    "best_pair_relationship": best_relation,
                    "assigned_group_id": window_index,
                    "group_size": len(window),
                    "group_positions": window_positions,
                    "group_best_bic": best.bic,
                }
                row.update(bic_fields)
                row.update({
                    f"prob_{relation}": probabilities[relation]
                    for relation in RELATIONS
                })
                pair_rows[candidate_key(candidate)] = row
            window_writer.writerow({
                "group_id": window_index,
                "chrom": window[0][0],
                "positions": window_positions,
                "mutation_count": len(window),
                "physical_reads": stats["reads"],
                "unique_reads": stats["unique_reads"],
                "eligible_pairs": len(eligible_pairs),
                "assigned_pairs": len(assigned),
                "best_bic": best.bic,
                "best_groups": best_groups,
                "best_parents": best_parents,
                "trees": stats["trees"],
                "bound_pruned": stats["bound_pruned"],
                "em_pruned": stats["em_pruned"],
                "em_fits": stats["em_fits"],
                "em_iterations": stats["em_iterations"],
                "em_converged": stats["em_converged"],
                "em_iteration_limited": stats["em_iteration_limited"],
                "best_converged": int(best.converged),
                "elapsed_seconds": elapsed,
            })
            if (window_index + 1) % 100 == 0 or window_index + 1 == len(windows):
                print(
                    f"windows={window_index + 1}/{len(windows)} "
                    f"pairs={len(pair_rows)}/{len(candidates)} "
                    f"fit_seconds={total_fit_seconds:.1f}",
                    flush=True,
                )
    if len(pair_rows) != len(candidates):
        raise RuntimeError(
            f"wrote {len(pair_rows)} pair calls for {len(candidates)} candidates"
        )
    overlap_summary, conflict_keys = annotate_window_relationships(
        pair_rows, pair_opinions,
    )
    with (args.output_dir / "pair_results.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=pair_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(pair_rows[key] for key in sorted(pair_rows))

    conflict_log = write_window_conflicts(
        args.output_dir, pair_rows, pair_opinions, conflict_keys, overlap_summary,
    )

    summary = {
        "chromosome": args.chrom,
        "configuration": {
            "method": "exhaustive_shared_tree",
            "maximum_mutations": args.window_size,
            "group_selection": "stride",
            "group_stride": args.window_stride,
            "large_component_windows": f"genomic_order_stride_{args.window_stride}",
            "pair_assignment": "largest_group_then_endpoint_centrality",
            "minimum_informative_reads": args.minimum_informative_reads,
            "minimum_alt_alt_reads": args.minimum_alt_alt_reads,
            "connectivity": (
                "ALT_ALT + ALT_REF + REF_ALT >= minimum_informative_reads "
                "OR ALT_ALT >= minimum_alt_alt_reads"
            ),
            "read_weighting": "one_equal_weight_per_physical_read",
            "read_scope": "all reads overlapping at least one group mutation",
            "supplementary_alignment_handling": "inherited_from_prepared_pileup",
            "em_iterations": args.em_iterations,
            "outlier_rate": 0.01,
            "error_floor": 0.01,
        },
        "mutations_in_vcf": len(mutations),
        "candidate_pool_pairs": len(candidate_pool),
        "connectivity_gate_rejected_pairs": len(candidate_pool) - len(candidates),
        "eligible_pairs": len(candidates),
        "connected_mutations": sum(component_sizes),
        "components": len(component_sizes),
        "component_size_counts": dict(Counter(component_sizes)),
        "groups": len(windows),
        "group_size_counts": dict(Counter(map(len, windows))),
        "assigned_pairs": len(pair_rows),
        "group_relationship_overlap": overlap_summary,
        "conflicting_pair_log": conflict_log,
        "input": observation_stats,
        "fit_seconds": total_fit_seconds,
        "runtime_seconds": time.perf_counter() - started,
    }
    (args.output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))
    return pair_rows
