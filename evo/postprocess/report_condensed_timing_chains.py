#!/usr/bin/env python3
"""
Find timing chains in condensed graphs, assign DAG ranks, write:
  - timing_chains.tsv
  - node_ranks.tsv
and write rank fields back to both condensed and uncondensed graph JSON files.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict, deque
from pathlib import Path
from typing import Optional


MAX_CHAIN_NODES = 20


def find_condensed_graph_files(main_out_dir: Path) -> list[Path]:
    return sorted(main_out_dir.rglob("chunk_*_condensed.json"))


def load_graph(path: Path) -> tuple[dict, list, list]:
    with path.open("r", encoding="utf-8") as f:
        obj = json.load(f)
    if isinstance(obj, dict) and isinstance(obj.get("elements"), dict):
        nodes = obj["elements"].get("nodes", []) or []
        edges = obj["elements"].get("edges", []) or []
    else:
        nodes = obj.get("nodes", []) or []
        edges = obj.get("edges", []) or []
    if not isinstance(nodes, list) or not isinstance(edges, list):
        raise ValueError(f"{path}: unexpected nodes/edges schema")
    return obj, nodes, edges


def format_members(members: object) -> str:
    if members is None:
        return ""
    if isinstance(members, list):
        return ",".join(str(m) for m in members)
    return str(members)


def find_chains(successors: dict[str, list[str]], node_data: dict[str, dict]) -> set[tuple[str, ...]]:
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


def compute_ranks(node_data: dict[str, dict], successors: dict[str, list[str]]) -> dict[str, Optional[int]]:
    """
    rank = 1 + max(rank of predecessors)
    roots (in timing subgraph) get rank 1
    nodes with no timing edges get None
    """
    timing_nodes: set[str] = set()
    for src, tgts in successors.items():
        for tgt in tgts:
            if src in node_data and tgt in node_data:
                timing_nodes.add(src)
                timing_nodes.add(tgt)

    in_degree: dict[str, int] = {n: 0 for n in timing_nodes}
    for src in timing_nodes:
        for tgt in successors.get(src, []):
            if tgt in timing_nodes:
                in_degree[tgt] += 1

    ranks: dict[str, Optional[int]] = {}
    q: deque[str] = deque()
    for n in timing_nodes:
        if in_degree[n] == 0:
            ranks[n] = 1
            q.append(n)

    while q:
        u = q.popleft()
        for v in successors.get(u, []):
            if v not in timing_nodes:
                continue
            cand = (ranks[u] or 0) + 1
            if ranks.get(v) is None or cand > (ranks[v] or 0):
                ranks[v] = cand
            in_degree[v] -= 1
            if in_degree[v] == 0:
                q.append(v)

    for nid in node_data:
        if nid not in ranks:
            ranks[nid] = None
    return ranks


def write_ranks_to_uncondensed_graph(condensed_path: Path, node_data: dict[str, dict], ranks: dict[str, Optional[int]]) -> None:
    name = condensed_path.name
    if not name.endswith("_condensed.json"):
        return
    uncondensed = condensed_path.with_name(name[: -len("_condensed.json")] + ".json")
    if not uncondensed.exists():
        return

    pos_to_rank: dict[str, Optional[int]] = {}
    for nid, data in node_data.items():
        r = ranks.get(nid)
        members = data.get("members")
        if isinstance(members, list) and members:
            for m in members:
                pos_to_rank[str(m)] = r
        else:
            pos_to_rank[str(nid)] = r

    with uncondensed.open("r", encoding="utf-8") as f:
        obj = json.load(f)
    if isinstance(obj, dict) and isinstance(obj.get("elements"), dict):
        nodes = obj["elements"].get("nodes", []) or []
    else:
        nodes = obj.get("nodes", []) or []

    for n in nodes:
        data = n.get("data", {})
        nid = str(data.get("id", ""))
        data["rank"] = pos_to_rank.get(nid)
        n["data"] = data

    with uncondensed.open("w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2)


def process_file(path: Path, node_data: dict[str, dict], successors: dict[str, list[str]]) -> tuple[list[tuple[str, ...]], dict]:
    obj, nodes, edges = load_graph(path)
    node_data.clear()
    successors.clear()

    for n in nodes:
        data = n.get("data", {}) or {}
        nid = data.get("id")
        if nid is not None:
            node_data[str(nid)] = data

    for e in edges:
        data = e.get("data", {}) or {}
        if data.get("relation") != "timing":
            continue
        src = data.get("source")
        tgt = data.get("target")
        if src is None or tgt is None:
            continue
        s = str(src)
        t = str(tgt)
        if s in node_data and t in node_data:
            successors[s].append(t)

    chains = sorted(find_chains(successors, node_data))
    return chains, obj


def main() -> None:
    ap = argparse.ArgumentParser(
        description=(
            "Report timing chains in condensed graphs, write timing_chains.tsv, "
            "assign ranks and write node_ranks.tsv, and annotate rank fields in graph JSONs."
        )
    )
    ap.add_argument("--outdir", required=True, type=Path, help="Main output directory.")
    args = ap.parse_args()
    outdir = args.outdir.resolve()

    chains_out = outdir / "timing_chains.tsv"
    ranks_out = outdir / "node_ranks.tsv"
    condensed_files = find_condensed_graph_files(outdir)

    if not condensed_files:
        print("No condensed graphs found.", file=sys.stderr)
        with chains_out.open("w", encoding="utf-8") as f:
            f.write("file\tchain_index\tnode_rank\tnode_id\tmembers\tsource_vcf\n")
        with ranks_out.open("w", encoding="utf-8") as f:
            f.write("file\tnode_id\trank\tchrom\tstart\tend\tmember_count\tmembers\n")
        return

    node_data: dict[str, dict] = {}
    successors: dict[str, list[str]] = defaultdict(list)
    length_counts: dict[int, int] = defaultdict(int)
    rank_distribution: dict[str, int] = defaultdict(int)
    total_ranked = 0
    total_unranked = 0
    all_rank_rows: list[tuple[str, str, Optional[int], str, int, int, int, str]] = []

    with chains_out.open("w", encoding="utf-8") as tsv:
        tsv.write("file\tchain_index\tnode_rank\tnode_id\tmembers\tsource_vcf\n")
        for path in condensed_files:
            try:
                chains, obj = process_file(path, node_data, successors)
            except Exception as e:
                print(f"Error processing {path}: {e}", file=sys.stderr)
                continue

            ranks = compute_ranks(node_data, successors)
            for nid, r in ranks.items():
                if nid in node_data:
                    node_data[nid]["rank"] = r

            with path.open("w", encoding="utf-8") as f:
                json.dump(obj, f, indent=2)
            write_ranks_to_uncondensed_graph(path, node_data, ranks)

            try:
                file_rel = str(path.relative_to(outdir))
            except ValueError:
                file_rel = str(path)

            for nid, data in node_data.items():
                r = ranks.get(nid)
                chrom = str(data.get("chrom", ""))
                members = data.get("members")
                positions: list[int] = []
                if isinstance(members, list):
                    for m in members:
                        try:
                            positions.append(int(m))
                        except (TypeError, ValueError):
                            continue
                else:
                    try:
                        positions.append(int(nid))
                    except (TypeError, ValueError):
                        pass
                start = min(positions) if positions else 0
                end = max(positions) if positions else 0
                member_count = len(positions) if positions else 1
                members_str = ",".join(str(p) for p in sorted(positions))
                all_rank_rows.append((file_rel, str(nid), r, chrom, start, end, member_count, members_str))

                if r is None:
                    rank_distribution["N/A"] += member_count
                    total_unranked += member_count
                else:
                    rank_distribution[str(r)] += member_count
                    total_ranked += member_count

            for i, chain in enumerate(chains):
                length_counts[len(chain)] += 1
                for rank_in_chain, nid in enumerate(chain):
                    data = node_data.get(str(nid), {})
                    members = format_members(data.get("members"))
                    source_vcf = str(data.get("source_vcf") or "")
                    tsv.write(
                        f"{file_rel}\t{i}\t{rank_in_chain}\t{nid}\t"
                        f"{members.replace(chr(9), ' ')}\t{source_vcf.replace(chr(9), ' ')}\n"
                    )

    with ranks_out.open("w", encoding="utf-8") as f:
        f.write("file\tnode_id\trank\tchrom\tstart\tend\tmember_count\tmembers\n")
        for file_rel, nid, r, chrom, start, end, n, members_str in all_rank_rows:
            r_str = "N/A" if r is None else str(r)
            f.write(f"{file_rel}\t{nid}\t{r_str}\t{chrom}\t{start}\t{end}\t{n}\t{members_str}\n")

    total_chains = sum(length_counts.values())
    print(f"Total chains detected: {total_chains}")
    if total_chains > 0:
        print("Chains by length (nodes):")
        for k in sorted(length_counts):
            print(f"  {k}: {length_counts[k]}")
    else:
        print("Chains by length (nodes): (none)")

    print("Rank distribution (uncondensed mutations):")
    na = rank_distribution.pop("N/A", 0)
    print(f"  N/A: {na}")
    for k in sorted(rank_distribution.keys(), key=int):
        print(f"  {k}: {rank_distribution[k]}")
    print(f"Total ranked: {total_ranked}")
    print(f"Total unranked (N/A): {total_unranked}")


if __name__ == "__main__":
    main()

