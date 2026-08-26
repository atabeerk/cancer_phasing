#!/usr/bin/env python3

import argparse
import os
import sys
import time
from collections import Counter

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
if CURRENT_DIR not in sys.path:
    sys.path.insert(0, CURRENT_DIR)

from builder import GraphBuilder
from fast_builder import FastGraphBuilder
from parsing import load_edges_from_base


def accepted_key(edge):
    return (
        edge.chrom,
        edge.u,
        edge.v,
        edge.relation,
        edge.reliability,
        edge.alt_alt,
        edge.alt_ref,
        edge.ref_alt,
        edge.ref_ref,
    )


def cluster_sets(builder):
    return sorted(tuple(sorted(nodes)) for nodes in builder.cluster_nodes.values())


def run_builder(builder_cls, edges):
    builder = builder_cls()
    results = []
    t0 = time.monotonic()
    for edge in edges:
        results.append(builder.try_add_edge(edge))
    elapsed = time.monotonic() - t0
    return builder, results, elapsed


def compare_builders(base_path, limit):
    t0 = time.monotonic()
    edges = load_edges_from_base(base_path)
    load_elapsed = time.monotonic() - t0
    if limit is not None:
        edges = edges[:limit]

    legacy, legacy_results, legacy_elapsed = run_builder(GraphBuilder, edges)
    fast, fast_results, fast_elapsed = run_builder(FastGraphBuilder, edges)

    checks = {
        "results": legacy_results == fast_results,
        "accepted_edges": [accepted_key(e) for e in legacy.edges] == [accepted_key(e) for e in fast.edges],
        "rejection_reasons": [r["reason"] for r in legacy.inconsistencies]
        == [r["reason"] for r in fast.inconsistencies],
        "cluster_sets": cluster_sets(legacy) == cluster_sets(fast),
        "nodes": legacy.nodes == fast.nodes,
    }

    print(f"base={base_path}")
    print(f"limit={limit if limit is not None else 'all'}")
    print(f"edges_compared={len(edges)}")
    print(f"load_sort_seconds={load_elapsed:.3f}")
    print(f"legacy_seconds={legacy_elapsed:.3f}")
    print(f"fast_seconds={fast_elapsed:.3f}")
    if fast_elapsed > 0:
        print(f"speedup={legacy_elapsed / fast_elapsed:.2f}x")
    print(f"legacy_accepted={len(legacy.edges)}")
    print(f"fast_accepted={len(fast.edges)}")
    print(f"legacy_rejected={len(legacy.inconsistencies)}")
    print(f"fast_rejected={len(fast.inconsistencies)}")
    print(f"legacy_reason_counts={dict(Counter(r['reason'] for r in legacy.inconsistencies))}")
    print(f"fast_reason_counts={dict(Counter(r['reason'] for r in fast.inconsistencies))}")
    for name, ok in checks.items():
        print(f"check_{name}={'OK' if ok else 'FAIL'}")

    if not all(checks.values()):
        return 1
    return 0


def main():
    parser = argparse.ArgumentParser(description="Compare legacy and fast graph builders on one chunk base.")
    parser.add_argument("base_path", help="Chunk base path without relation suffix or .txt")
    parser.add_argument("--limit", type=int, default=None, help="Compare only the first N sorted edges.")
    args = parser.parse_args()
    raise SystemExit(compare_builders(args.base_path, args.limit))


if __name__ == "__main__":
    main()
