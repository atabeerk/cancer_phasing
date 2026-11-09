"""Sliding and supplemental mutation windows with exhaustive shared-tree fitting."""

import math
from collections import defaultdict
from functools import lru_cache
from itertools import combinations, product

import numpy as np

from joint_clone_model import (
    GenotypeLikelihoodCache, canonical_tree, fit_tree,
    is_ancestor, node_for_mutation, optimistic_bic_lower_bound,
)


RELATIONS = (
    "cooccurring",
    "left_before_right",
    "right_before_left",
    "divergent",
)


def build_windows(candidates, window_size, window_stride):
    """Cover every eligible pair with connected mutation windows.

    Use genomic windows at the requested stride within large components, then
    add supplemental windows where regular windows miss eligible pairs.
    Return windows and component sizes; singleton mutations produce no window.
    """
    if not 2 <= window_size <= 5:
        raise ValueError("window_size must be between two and five")
    if window_stride < 1:
        raise ValueError("window stride must be positive")
    adjacency = defaultdict(set)
    edges = []
    for candidate in candidates:
        left = candidate.chrom, candidate.left
        right = candidate.chrom, candidate.right
        adjacency[left].add(right)
        adjacency[right].add(left)
        edges.append((left, right))
    components = []
    unseen = set(adjacency)
    while unseen:
        start = min(unseen)
        unseen.remove(start)
        stack, component = [start], []
        while stack:
            vertex = stack.pop()
            component.append(vertex)
            for neighbor in adjacency[vertex]:
                if neighbor in unseen:
                    unseen.remove(neighbor)
                    stack.append(neighbor)
        components.append(tuple(sorted(component)))

    def is_connected(window):
        remaining = set(window)
        stack = [min(remaining)]
        remaining.remove(stack[0])
        while stack:
            for neighbor in adjacency[stack.pop()] & remaining:
                remaining.remove(neighbor)
                stack.append(neighbor)
        return not remaining

    windows = set()
    covered = set()

    def add(window):
        window = tuple(sorted(window))
        windows.add(window)
        for left, right in combinations(window, 2):
            if right in adjacency[left]:
                covered.add((left, right))

    for component in components:
        if len(component) <= window_size:
            add(component)
        else:
            last_offset = len(component) - window_size
            offsets = list(range(0, last_offset + 1, window_stride))
            if offsets[-1] != last_offset:
                offsets.append(last_offset)
            for offset in offsets:
                window = component[offset:offset + window_size]
                if is_connected(window):
                    add(window)
    for edge in sorted(edges):
        if edge in covered:
            continue
        selected = set(edge)
        while len(selected) < window_size:
            frontier = set().union(*(adjacency[node] for node in selected)) - selected
            if not frontier:
                break
            neighbor = min(frontier, key=lambda node: (
                min(abs(node[1] - other[1]) for other in selected), node,
            ))
            selected.add(neighbor)
        add(selected)
    return sorted(windows), sorted(map(len, components))


def assign_pairs_to_windows(candidates, windows):
    """Choose one deterministic covering window for every eligible pair.

    Prefer a larger window, then the window where both endpoints are most
    central. This prevents overlapping windows from casting duplicate votes.
    """
    candidate_keys = {
        ((candidate.chrom, candidate.left), (candidate.chrom, candidate.right))
        for candidate in candidates
    }
    covering = defaultdict(list)
    for window_index, window in enumerate(windows):
        for left, right in combinations(window, 2):
            if (left, right) in candidate_keys:
                covering[(left, right)].append(window_index)
    assignments = {}
    for pair in sorted(candidate_keys):
        choices = covering.get(pair, ())
        if not choices:
            raise ValueError(f"eligible pair is not covered by a window: {pair}")

        def rank(window_index):
            window = windows[window_index]
            center = (len(window) - 1) / 2.0
            endpoint_indices = (window.index(pair[0]), window.index(pair[1]))
            distances = tuple(abs(index - center) for index in endpoint_indices)
            return (
                -len(window),
                sum(distances),
                max(distances),
                window[-1][1] - window[0][1],
                window,
            )

        assignments[pair] = min(choices, key=rank)
    return assignments


def _partitions(items):
    if not items:
        yield ()
        return
    first, rest = items[0], items[1:]
    for partition in _partitions(rest):
        yield ((first,),) + partition
        for index in range(len(partition)):
            yield partition[:index] + ((first,) + partition[index],) + partition[index + 1:]


@lru_cache(maxsize=None)
def _rooted_parents(node_count):
    parents = []
    choices = [tuple(parent for parent in range(node_count + 1) if parent != node)
               for node in range(1, node_count + 1)]
    for assignments in product(*choices):
        values = (-1,) + assignments
        valid = True
        for node in range(1, node_count + 1):
            visited = set()
            while node != 0:
                if node in visited:
                    valid = False
                    break
                visited.add(node)
                node = values[node]
            if not valid:
                break
        if valid:
            parents.append(values)
    return tuple(parents)


@lru_cache(maxsize=4)
def exhaustive_trees(mutation_count):
    """Enumerate partitions into acquisition nodes and every rooted topology."""
    if not 2 <= mutation_count <= 5:
        raise ValueError("shared trees support two through five mutations")
    trees = []
    for partition in _partitions(tuple(range(mutation_count))):
        for parents in _rooted_parents(len(partition)):
            trees.append(canonical_tree(((),) + partition, parents))
    # Simple trees first give safe BIC pruning an early incumbent.
    return tuple(sorted(trees, key=lambda tree: (
        len(tree.groups), tree.groups, tree.parents,
    )))


@lru_cache(maxsize=4)
def tree_relation_signatures(mutation_count):
    """Cache every pair relationship for every enumerated tree."""
    pairs = tuple(combinations(range(mutation_count), 2))
    signatures = []
    for tree in exhaustive_trees(mutation_count):
        mutation_nodes = node_for_mutation(tree)
        relations = []
        for left, right in pairs:
            left_node, right_node = mutation_nodes[left], mutation_nodes[right]
            if left_node == right_node:
                relation = "cooccurring"
            elif is_ancestor(tree, left_node, right_node):
                relation = "left_before_right"
            elif is_ancestor(tree, right_node, left_node):
                relation = "right_before_left"
            else:
                relation = "divergent"
            relations.append(relation)
        signatures.append(tuple(relations))
    return pairs, tuple(signatures)


@lru_cache(maxsize=4)
def tree_genotype_indices(mutation_count):
    """Return dense-genotype column indices for every enumerated tree."""
    powers = 1 << np.arange(mutation_count, dtype=np.int64)
    return tuple(
        tuple(
            int(index)
            for index in tree.genotype_matrix(mutation_count).astype(
                np.int64
            ) @ powers
        )
        for tree in exhaustive_trees(mutation_count)
    )


def dense_genotype_likelihoods(cache, mutation_count):
    """Compute every possible genotype likelihood once for a local window."""
    masks = np.arange(1 << mutation_count, dtype=np.uint32)
    bits = np.arange(mutation_count, dtype=np.uint32)
    genotypes = ((masks[:, np.newaxis] >> bits) & 1).astype(np.bool_)
    return cache.log_likelihood_matrix(genotypes)


def ranked_tree_likelihoods(cache, mutation_count, signatures):
    """Order trees by their certified optimistic BIC lower bounds."""
    trees = exhaustive_trees(mutation_count)
    genotype_indices = tree_genotype_indices(mutation_count)
    dense_likelihoods = dense_genotype_likelihoods(cache, mutation_count)
    ranked = []
    for original_index, (tree, signature, indices) in enumerate(
        zip(trees, signatures, genotype_indices)
    ):
        clone_likelihoods = dense_likelihoods[:, indices]
        lower_bound = optimistic_bic_lower_bound(
            tree,
            mutation_count,
            0.01,
            cache,
            clone_log_likelihood=clone_likelihoods,
        )
        ranked.append(
            (lower_bound, original_index, tree, signature, indices)
        )
    ranked.sort(key=lambda item: (item[0], item[1]))
    return ranked, dense_likelihoods


def probabilities_from_bics(bics):
    minimum = min(bics.values())
    weights = {
        relation: math.exp(-0.5 * min(1400.0, bic - minimum))
        for relation, bic in bics.items()
    }
    denominator = sum(weights.values())
    return {relation: weight / denominator for relation, weight in weights.items()}


def fit_shared_tree_relations(
    reads, mutation_count, pair_indices, max_iterations=50,
):
    """Fit every topology once and retain relation-specific BIC minima.

    The returned best tree supplies the joint labels. For each requested pair,
    the four BIC values are the best shared-tree fits constrained to that
    relationship. A tree is pruned only when it cannot improve any applicable
    pair/relationship bucket.
    """
    if not reads:
        raise ValueError("at least one read is required")
    all_pairs, signatures = tree_relation_signatures(mutation_count)
    pair_to_slot = {pair: slot for slot, pair in enumerate(all_pairs)}
    pair_indices = tuple(tuple(sorted(pair)) for pair in pair_indices)
    if len(set(pair_indices)) != len(pair_indices):
        raise ValueError("pair indices must be unique")
    slots = []
    for pair in pair_indices:
        if pair not in pair_to_slot:
            raise ValueError(f"invalid mutation pair: {pair}")
        slots.append(pair_to_slot[pair])

    cache = GenotypeLikelihoodCache(reads, mutation_count)
    tree_entries, dense_likelihoods = ranked_tree_likelihoods(
        cache, mutation_count, signatures
    )
    best = None
    best_original_index = None
    relation_bics = [dict.fromkeys(RELATIONS, math.inf) for _ in pair_indices]
    stats = dict(
        trees=0, bound_pruned=0, em_pruned=0, em_fits=0, em_iterations=0,
        reads=len(reads), unique_reads=cache.unique_read_count,
        em_converged=0, em_iteration_limited=0,
    )
    for lower_bound, original_index, tree, signature, genotype_indices in tree_entries:
        stats["trees"] += 1
        relations = tuple(signature[slot] for slot in slots)
        thresholds = [
            relation_bics[index][relation]
            for index, relation in enumerate(relations)
        ]
        cutoff = max(thresholds) if thresholds and all(
            math.isfinite(value) for value in thresholds
        ) else None
        if cutoff is not None and lower_bound > cutoff:
            stats["bound_pruned"] += 1
            continue
        clone_likelihoods = np.ascontiguousarray(
            dense_likelihoods[:, genotype_indices]
        )
        fit = fit_tree(
            tree, reads, mutation_count, 0.01, max_iterations,
            likelihood_cache=cache,
            bic_cutoff=(
                math.nextafter(cutoff, math.inf)
                if cutoff is not None else None
            ),
            clone_log_likelihood=clone_likelihoods,
        )
        stats["em_fits"] += 1
        stats["em_iterations"] += fit.em_iterations
        if fit.certified_pruned:
            stats["em_pruned"] += 1
            continue
        stats["em_converged" if fit.converged else "em_iteration_limited"] += 1
        if (
            best is None
            or fit.bic < best.bic
            or (
                fit.bic == best.bic
                and original_index < best_original_index
            )
        ):
            best = fit
            best_original_index = original_index
        for index, relation in enumerate(relations):
            relation_bics[index][relation] = min(
                relation_bics[index][relation], fit.bic
            )
    if best is None:
        raise RuntimeError("no shared tree was fitted")
    for bics in relation_bics:
        if not all(math.isfinite(value) for value in bics.values()):
            raise RuntimeError("not every pair relationship received a fit")
    return best, relation_bics, stats
