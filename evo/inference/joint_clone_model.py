"""Read likelihoods, fixed-tree EM fitting, and mutation-tree utilities."""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np


@dataclass(frozen=True)
class ReadData:
    mutation_indices: np.ndarray
    alleles: np.ndarray
    errors: np.ndarray


@dataclass(frozen=True)
class TreeState:
    """Rooted infinite-sites tree of distinct local genotypes.

    Node zero is the all-reference local state. Each other node owns the
    mutations acquired on its incoming branch. Every node is a possible
    sampled local clone state. Parents precede their children, as produced
    by canonical_tree().
    """

    groups: Tuple[Tuple[int, ...], ...]
    parents: Tuple[int, ...]

    def validate(self, mutation_count: int) -> None:
        if len(self.parents) != len(self.groups):
            raise ValueError("each node must have one parent entry")
        if not self.groups or self.groups[0] or self.parents[0] != -1:
            raise ValueError("node zero must be an empty root with parent -1")
        seen: List[int] = []
        for node in range(1, len(self.groups)):
            if not self.groups[node]:
                raise ValueError("non-root nodes cannot be empty")
            parent = self.parents[node]
            if not 0 <= parent < node:
                raise ValueError(f"parent {parent} must precede node {node}")
            seen.extend(self.groups[node])
        if sorted(seen) != list(range(mutation_count)):
            raise ValueError("each mutation must occur exactly once")

    def genotype_matrix(self, mutation_count: int) -> np.ndarray:
        self.validate(mutation_count)
        genotypes = np.zeros((len(self.groups), mutation_count), dtype=np.bool_)
        for node in range(1, len(self.groups)):
            genotypes[node] = genotypes[self.parents[node]]
            genotypes[node, list(self.groups[node])] = True
        return genotypes


@dataclass
class Fit:
    state: TreeState
    log_likelihood: float
    bic: float
    mixture: np.ndarray
    em_iterations: int
    certified_pruned: bool = False
    converged: Optional[bool] = None  # None for fits stopped by safe pruning.


class GenotypeLikelihoodCache:
    """Compress identical reads and build the shared genotype likelihood table."""

    def __init__(self, reads: Sequence[ReadData], mutation_count: int) -> None:
        self.reads = tuple(reads)
        self.mutation_count = mutation_count
        unique_reads: List[ReadData] = []
        read_counts: List[int] = []
        signature_to_index: Dict[Tuple[bytes, bytes, bytes], int] = {}
        for read in self.reads:
            signature = (
                read.mutation_indices.tobytes(),
                read.alleles.tobytes(),
                read.errors.tobytes(),
            )
            index = signature_to_index.get(signature)
            if index is None:
                signature_to_index[signature] = len(unique_reads)
                unique_reads.append(read)
                read_counts.append(1)
            else:
                read_counts[index] += 1
        self.unique_reads = tuple(unique_reads)
        self.read_weights = np.asarray(read_counts, dtype=np.float64)
        self.outlier_log_likelihood = np.asarray(
            [
                -len(read.mutation_indices) * math.log(2.0)
                for read in self.unique_reads
            ],
            dtype=np.float64,
        )

    @property
    def unique_read_count(self) -> int:
        return len(self.unique_reads)

    def log_likelihood_matrix(self, genotypes: np.ndarray) -> np.ndarray:
        genotypes = np.asarray(genotypes, dtype=np.bool_)
        if genotypes.ndim != 2 or genotypes.shape[1] != self.mutation_count:
            raise ValueError(
                "genotypes must be a matrix with one column per mutation"
            )

        values = np.empty(
            (self.unique_read_count, len(genotypes)), dtype=np.float64
        )
        for read_index, read in enumerate(self.unique_reads):
            expected = genotypes[:, read.mutation_indices]
            matches = expected == read.alleles[np.newaxis, :]
            probabilities = np.where(matches, 1.0 - read.errors, read.errors)
            values[read_index] = np.log(probabilities).sum(axis=1)
        return values


def canonical_tree(groups: Sequence[Sequence[int]], parents: Sequence[int]) -> TreeState:
    children: Dict[int, List[int]] = {node: [] for node in range(len(groups))}
    for node in range(1, len(groups)):
        children[parents[node]].append(node)
    for values in children.values():
        values.sort(key=lambda node: tuple(sorted(groups[node])))
    order: List[int] = []

    def visit(node: int) -> None:
        order.append(node)
        for child in children[node]:
            visit(child)

    visit(0)
    if len(order) != len(groups):
        raise ValueError("tree is cyclic or disconnected")
    remap = {old: new for new, old in enumerate(order)}
    return TreeState(
        tuple(tuple(sorted(groups[old])) for old in order),
        tuple(-1 if old == 0 else remap[parents[old]] for old in order),
    )


def _prepare_clone_log_likelihood(
    state: TreeState,
    mutation_count: int,
    likelihood_cache: GenotypeLikelihoodCache,
    clone_log_likelihood: Optional[np.ndarray],
) -> np.ndarray:
    if likelihood_cache.mutation_count != mutation_count:
        raise ValueError("likelihood cache does not match the mutation count")
    if clone_log_likelihood is None:
        return likelihood_cache.log_likelihood_matrix(
            state.genotype_matrix(mutation_count)
        )
    clone_log_likelihood = np.asarray(clone_log_likelihood, dtype=np.float64)
    expected_shape = (likelihood_cache.unique_read_count, len(state.groups))
    if clone_log_likelihood.shape != expected_shape:
        raise ValueError(
            "clone likelihood matrix has shape "
            f"{clone_log_likelihood.shape}; expected {expected_shape}"
        )
    return clone_log_likelihood


def optimistic_bic_lower_bound(
    state: TreeState,
    mutation_count: int,
    outlier_rate: float,
    likelihood_cache: GenotypeLikelihoodCache,
    clone_log_likelihood: Optional[np.ndarray] = None,
) -> float:
    """Return an optimistic BIC obtained by letting each read pick its best clone."""
    clone_log_likelihood = _prepare_clone_log_likelihood(
        state, mutation_count, likelihood_cache, clone_log_likelihood,
    )
    best_clone_log_likelihood = np.max(clone_log_likelihood, axis=1)
    log_regular_weight = math.log1p(-outlier_rate) if outlier_rate < 1.0 else -np.inf
    log_outlier_weight = math.log(outlier_rate) if outlier_rate > 0.0 else -np.inf
    read_log_likelihood = np.logaddexp(
        log_regular_weight + best_clone_log_likelihood,
        log_outlier_weight + likelihood_cache.outlier_log_likelihood,
    )
    upper_log_likelihood = float(
        (read_log_likelihood * likelihood_cache.read_weights).sum()
    )
    parameter_count = len(state.groups) - 1
    return (
        -2.0 * upper_log_likelihood
        + parameter_count * math.log(max(2, len(likelihood_cache.reads)))
    )


def fit_tree(
    state: TreeState,
    reads: Sequence[ReadData],
    mutation_count: int,
    outlier_rate: float,
    max_iterations: int,
    tolerance: float = 1e-8,
    likelihood_cache: Optional[GenotypeLikelihoodCache] = None,
    bic_cutoff: Optional[float] = None,
    clone_log_likelihood: Optional[np.ndarray] = None,
) -> Fit:
    if likelihood_cache is None:
        likelihood_cache = GenotypeLikelihoodCache(reads, mutation_count)
    elif len(likelihood_cache.reads) != len(reads):
        raise ValueError("likelihood cache does not match the fit inputs")
    clone_log_likelihood = _prepare_clone_log_likelihood(
        state, mutation_count, likelihood_cache, clone_log_likelihood,
    )
    outlier_log_likelihood = likelihood_cache.outlier_log_likelihood
    read_weights = likelihood_cache.read_weights
    parameter_count = len(state.groups) - 1
    bic_penalty = parameter_count * math.log(max(2, len(reads)))

    mixture = np.full(len(state.groups), 1.0 / len(state.groups), dtype=np.float64)
    log_regular_weight = math.log1p(-outlier_rate) if outlier_rate < 1.0 else -np.inf
    log_outlier_weight = math.log(outlier_rate) if outlier_rate > 0.0 else -np.inf
    previous_log_likelihood = -np.inf
    iteration = 0
    converged = False
    for iteration in range(1, max_iterations + 1):
        clone_component_log_likelihood = (
            clone_log_likelihood + np.log(np.maximum(mixture, 1e-300))
        )
        clone_mixture_log_likelihood = np.logaddexp.reduce(
            clone_component_log_likelihood, axis=1
        )
        read_log_likelihood = np.logaddexp(
            log_regular_weight + clone_mixture_log_likelihood,
            log_outlier_weight + outlier_log_likelihood,
        )
        log_likelihood = float((read_log_likelihood * read_weights).sum())
        responsibilities = np.exp(
            log_regular_weight + clone_component_log_likelihood
            - read_log_likelihood[:, np.newaxis]
        )
        raw_counts = (
            responsibilities * read_weights[:, np.newaxis]
        ).sum(axis=0)
        counts = raw_counts + 1e-6
        updated_mixture = counts / counts.sum()
        if bic_cutoff is not None:
            gradient = raw_counts / np.maximum(mixture, 1e-300)
            duality_gap = max(
                0.0,
                float(np.max(gradient) - np.dot(gradient, mixture)),
            )
            upper_log_likelihood = log_likelihood + duality_gap
            bic_lower_bound = -2.0 * upper_log_likelihood + bic_penalty
            if bic_lower_bound >= bic_cutoff:
                bic = -2.0 * log_likelihood + bic_penalty
                return Fit(
                    state,
                    log_likelihood,
                    bic,
                    mixture,
                    iteration,
                    certified_pruned=True,
                )
        if (
            abs(log_likelihood - previous_log_likelihood) < tolerance
            and np.max(abs(updated_mixture - mixture)) < tolerance
        ):
            mixture = updated_mixture
            converged = True
            break
        mixture = updated_mixture
        previous_log_likelihood = log_likelihood

    # Score the final mixture, which was updated after the last E-step.
    clone_component_log_likelihood = (
        clone_log_likelihood + np.log(np.maximum(mixture, 1e-300))
    )
    clone_mixture_log_likelihood = np.logaddexp.reduce(
        clone_component_log_likelihood, axis=1
    )
    read_log_likelihood = np.logaddexp(
        log_regular_weight + clone_mixture_log_likelihood,
        log_outlier_weight + outlier_log_likelihood,
    )
    log_likelihood = float((read_log_likelihood * read_weights).sum())
    bic = -2.0 * log_likelihood + bic_penalty
    return Fit(state, log_likelihood, bic, mixture, iteration, converged=converged)


def node_for_mutation(state: TreeState) -> Dict[int, int]:
    return {
        mutation: node
        for node, group in enumerate(state.groups)
        for mutation in group
    }


def is_ancestor(state: TreeState, ancestor: int, node: int) -> bool:
    while node != 0:
        node = state.parents[node]
        if node == ancestor:
            return True
    return ancestor == 0


def inferred_relation(state: TreeState, left: int, right: int) -> str:
    mutation_nodes = node_for_mutation(state)
    left_node, right_node = mutation_nodes[left], mutation_nodes[right]
    if left_node == right_node:
        return "cooccurring"
    if is_ancestor(state, left_node, right_node):
        return "left_before_right"
    if is_ancestor(state, right_node, left_node):
        return "right_before_left"
    return "divergent"
