#!/usr/bin/env python3
"""
Scan condensed graphs under a main output directory, find timing-edge chains
(A→B→C or longer, 3–20 nodes), write detailed chain data to a TSV and print
a summary (chain count and length distribution) to stdout.
"""

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path


MAX_CHAIN_NODES = 20


def find_condensed_graph_files(main_out_dir: Path) -> list[Path]:
    """Recursively find condensed graphs: chunk_*_condensed.json"""
    out: list[Path] = []
    for p in main_out_dir.rglob("chunk_*_condensed.json"):
        out.append(p)
    return sorted(out)


def load_graph(path: Path) -> tuple[list, list]:
    """Load nodes and edges from a graph JSON (elements.nodes / elements.edges)."""
    with path.open("r", encoding="utf-8") as f:
        obj = json.load(f)
    if isinstance(obj, dict) and "elements" in obj and isinstance(obj["elements"], dict):
        nodes = obj["elements"].get("nodes", []) or []
        edges = obj["elements"].get("edges", []) or []
    else:
        nodes = obj.get("nodes", []) or []
        edges = obj.get("edges", []) or []
    if not isinstance(nodes, list) or not isinstance(edges, list):
        raise ValueError(f"{path}: unexpected JSON schema for nodes/edges")
    return nodes, edges


def format_members(members) -> str:
    """Stringify node.data['members'] for TSV."""
    if members is None:
        return ""
    if isinstance(members, list):
        return ",".join(str(m) for m in members)
    return str(members)


def find_chains(successors: dict[str, list[str]], node_data: dict[str, dict]) -> set[tuple[str, ...]]:
    """
    Enumerate all simple paths of 3 to MAX_CHAIN_NODES nodes (2 to MAX_CHAIN_NODES-1 edges)
    in the timing digraph. Returns set of paths, each path a tuple of node ids in order.
    """
    chains: set[tuple[str, ...]] = set()

    def dfs(path: list[str], visited: set[str]) -> None:
        if len(path) >= 3:
            chains.add(tuple(path))
        if len(path) >= MAX_CHAIN_NODES:
            return
        u = path[-1]
        for v in successors.get(u, []):
            if v not in visited and v in node_data:
                dfs(path + [v], visited | {v})

    for start in node_data:
        dfs([start], {start})

    return chains


def process_file(
    path: Path,
    node_data: dict[str, dict],
    successors: dict[str, list[str]],
) -> list[tuple[str, ...]]:
    """Load one condensed graph, build timing digraph, return list of chains (ordered)."""
    nodes, edges = load_graph(path)
    node_data.clear()
    successors.clear()
    for n in nodes:
        data = n.get("data", {}) or {}
        nid = data.get("id")
        if nid is None:
            continue
        node_data[str(nid)] = data
    successors.clear()
    for e in edges:
        data = e.get("data", {}) or {}
        if data.get("relation") != "timing":
            continue
        src = data.get("source")
        tgt = data.get("target")
        if src is None or tgt is None:
            continue
        src, tgt = str(src), str(tgt)
        if src in node_data and tgt in node_data:
            successors[src].append(tgt)
    chains = find_chains(successors, node_data)
    return sorted(chains)  # deterministic order for TSV


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Report timing-edge chains (3–20 nodes) in condensed graphs; write details to <outdir>/timing_chains.tsv and print summary to stdout.",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        type=Path,
        help="Main output directory containing condensed graph files.",
    )
    args = parser.parse_args()
    outdir = args.outdir.resolve()
    outpath = outdir / "timing_chains.tsv"

    condensed_files = find_condensed_graph_files(outdir)
    if not condensed_files:
        print("No condensed graphs found.", file=sys.stderr)
        with outpath.open("w", encoding="utf-8") as f:
            f.write("file\tchain_index\tnode_rank\tnode_id\tmembers\tsource_vcf\n")
        return

    node_data: dict[str, dict] = {}
    successors: dict[str, list[str]] = defaultdict(list)
    length_counts: dict[int, int] = defaultdict(int)

    with outpath.open("w", encoding="utf-8") as tsv:
        tsv.write("file\tchain_index\tnode_rank\tnode_id\tmembers\tsource_vcf\n")
        for path in condensed_files:
            try:
                chains = process_file(path, node_data, successors)
            except Exception as e:
                print(f"Error processing {path}: {e}", file=sys.stderr)
                continue
            try:
                file_rel = path.relative_to(outdir)
            except ValueError:
                file_rel = path
            file_str = str(file_rel)
            for i, chain in enumerate(chains):
                length_counts[len(chain)] += 1
                # Re-load node_data for this file (process_file cleared it)
                _, _ = load_graph(path)
                for rank, nid in enumerate(chain):
                    data = node_data.get(str(nid), {})
                    members = format_members(data.get("members"))
                    source_vcf = data.get("source_vcf") or ""
                    # Escape tabs in string fields
                    members_esc = members.replace("\t", " ")
                    source_vcf_esc = (source_vcf or "").replace("\t", " ")
                    tsv.write(f"{file_str}\t{i}\t{rank}\t{nid}\t{members_esc}\t{source_vcf_esc}\n")

    total = sum(length_counts.values())
    print(f"Total chains detected: {total}")
    if total > 0:
        print("Chains by length (nodes):")
        for k in sorted(length_counts):
            print(f"  {k}: {length_counts[k]}")
    else:
        print("Chains by length (nodes): (none)")


if __name__ == "__main__":
    main()
