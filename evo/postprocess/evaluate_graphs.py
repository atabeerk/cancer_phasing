#!/usr/bin/env python3
"""
evaluate_edges.py

Evaluate timing/cooccurrence edges in uncondensed chunk graphs against a ground-truth tree.

Directory layout expected:
  OUTDIR/
    12_out/graphs/*.json
    12_pre/...
    chr3_out/graphs/**/*.json
    chr3_pre/...
    ...

We ignore any graph JSON whose filename ends with "_condensed.json".

Node entries are expected to include:
  node["data"]["source_vcf"] = "N1.fixed.merged.withcontigs.sorted"
We map SNV -> tree node label by taking the prefix before the first dot (e.g., "N1").

Tree input is parent<TAB>child (no header).

Rules:
- Timing edge A->B is CORRECT iff node(A) is a STRICT ancestor of node(B).
  Timing between SNPs from the same node is an error.
- Cooccurrence edge is CORRECT iff node(A) == node(B). (Undirected)

UNKNOWN edge:
- If either endpoint cannot be mapped to a tree node (missing node, missing source_vcf, label not in tree)
  then the edge is UNKNOWN and excluded from the "inconsistent edge detail" files and node-pair counts.

Outputs (only these):
- edge_eval_summary.tsv (per chromosome + ALL, includes % consistent over known edges)
- edge_eval_inconsistent_timing_edges.tsv (edge-level details)
- edge_eval_inconsistent_cooccur_edges.tsv (edge-level details)
- edge_eval_timing_nodepair_counts.tsv (genome-wide node-pair counts, ordered)
- edge_eval_cooccur_nodepair_counts.tsv (genome-wide node-pair counts, unordered)
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple
from collections import defaultdict, Counter


# ---------------- Sorting utilities ----------------

def chrom_sort_key(chrom: str) -> Tuple[int, int, str]:
    """
    Natural chromosome order:
      chr1..chr22 (or 1..22), then chrX, chrY, chrMT/M/mito, then anything else.
    Returns a tuple used for sorting.
    """
    c = chrom.strip()
    c_lower = c.lower()
    if c_lower.startswith("chr"):
        c_core = c[3:]
    else:
        c_core = c

    c_core_lower = c_core.lower()

    # numeric chromosomes
    if c_core_lower.isdigit():
        return (0, int(c_core_lower), "")

    # special chromosomes
    special_rank = {
        "x": (1, 23, "X"),
        "y": (1, 24, "Y"),
        "mt": (1, 25, "MT"),
        "m": (1, 25, "MT"),
        "mito": (1, 25, "MT"),
        "mitochondria": (1, 25, "MT"),
    }
    if c_core_lower in special_rank:
        grp, num, label = special_rank[c_core_lower]
        return (grp, num, label)

    # fallback: lexicographic after known
    return (2, 10**9, c)


def nlabel_sort_key(n: str) -> Tuple[int, int, str]:
    """
    Sort N labels like N1, N2, ... N10 naturally.
    """
    s = n.strip()
    if len(s) >= 2 and s[0].upper() == "N" and s[1:].isdigit():
        return (0, int(s[1:]), "")
    return (1, 10**9, s)


# ---------------- CLI ----------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--outdir", required=True, help="Main output directory containing *_out folders.")
    p.add_argument("--tree", required=True, help="Tree file (parent<TAB>child).")
    p.add_argument("--out-prefix", default="edge_eval", help="Prefix for output files.")
    return p.parse_args()


# ---------------- Tree utilities ----------------

def read_tree_parent_child_tsv(path: Path) -> Tuple[Dict[str, str], Set[str]]:
    parent_of: Dict[str, str] = {}
    nodes: Set[str] = set()

    with path.open("r", encoding="utf-8") as f:
        for line_no, line in enumerate(f, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) != 2:
                raise ValueError(f"Tree parse error at line {line_no}: expected 2 tab-separated columns, got: {line!r}")
            parent, child = parts[0].strip(), parts[1].strip()
            if not parent or not child:
                raise ValueError(f"Tree parse error at line {line_no}: empty parent/child in: {line!r}")
            if child in parent_of and parent_of[child] != parent:
                raise ValueError(f"Tree conflict: child {child} has multiple parents: {parent_of[child]} and {parent}")
            parent_of[child] = parent
            nodes.add(parent)
            nodes.add(child)

    return parent_of, nodes


def build_depth_cache(parent_of: Dict[str, str]) -> Dict[str, int]:
    depth: Dict[str, int] = {}

    def get_depth(n: str, stack: Set[str]) -> int:
        if n in depth:
            return depth[n]
        if n in stack:
            raise ValueError(f"Cycle detected in tree at node {n}")
        stack.add(n)
        if n not in parent_of:
            d = 0
        else:
            d = get_depth(parent_of[n], stack) + 1
        stack.remove(n)
        depth[n] = d
        return d

    all_nodes = set(parent_of.keys()) | set(parent_of.values())
    for n in all_nodes:
        get_depth(n, set())
    return depth


def is_ancestor_strict(anc: str, desc: str, parent_of: Dict[str, str], depth: Dict[str, int]) -> bool:
    if anc == desc:
        return False
    if anc not in depth or desc not in depth:
        return False
    if depth[anc] >= depth[desc]:
        return False

    cur = desc
    while cur in parent_of and depth[cur] > depth[anc]:
        cur = parent_of[cur]
    return cur == anc


# ---------------- Graph parsing utilities ----------------

def load_graph_json(path: Path) -> Tuple[List[dict], List[dict]]:
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


def node_label_from_source_vcf(source_vcf: Optional[str]) -> Optional[str]:
    if not source_vcf or not isinstance(source_vcf, str):
        return None
    return source_vcf.split(".", 1)[0].strip() or None


def chromosome_from_out_folder(out_folder_name: str) -> str:
    return out_folder_name[:-4] if out_folder_name.endswith("_out") else out_folder_name


def should_skip_json(path: Path) -> bool:
    return path.name.endswith("_condensed.json")


# ---------------- Data structures ----------------

class ChrAgg:
    def __init__(self) -> None:
        self.timing = Counter()   # total/correct/incorrect/unknown
        self.cooccur = Counter()  # total/correct/incorrect/unknown


# ---------------- Evaluation ----------------

def eval_graph_file(
    json_path: Path,
    fallback_chr_label: str,
    parent_of: Dict[str, str],
    tree_nodes: Set[str],
    depth: Dict[str, int],
    agg_by_chr: Dict[str, ChrAgg],
    agg_all: ChrAgg,
    # edge-level inconsistent details
    inconsistent_timing_rows: List[Tuple],
    inconsistent_cooccur_rows: List[Tuple],
    # genome-wide node-pair counts (known only)
    timing_pair_counts: Dict[Tuple[str, str], Counter],
    cooccur_pair_counts: Dict[Tuple[str, str], Counter],
) -> None:
    nodes, edges = load_graph_json(json_path)

    # Map graph node id -> (chrom, position, tree_label)
    id_to_info: Dict[str, Tuple[str, int, Optional[str]]] = {}

    for n in nodes:
        data = n.get("data", {}) if isinstance(n, dict) else {}
        nid = data.get("id")
        if nid is None:
            continue

        # chrom: prefer node chrom field, else folder-derived
        chrom = str(data.get("chrom", fallback_chr_label))
        # position: prefer explicit position, else parse id if numeric
        pos_val = data.get("position", None)
        if isinstance(pos_val, int):
            pos = pos_val
        else:
            # try to parse
            try:
                pos = int(pos_val) if pos_val is not None else int(nid)
            except Exception:
                # if position truly unavailable, use -1
                pos = -1

        label = node_label_from_source_vcf(data.get("source_vcf"))
        id_to_info[str(nid)] = (chrom, pos, label)

    for e in edges:
        data = e.get("data", {}) if isinstance(e, dict) else {}
        relation = data.get("relation")
        if relation not in ("timing", "cooccurring"):
            continue

        src_id = data.get("source")
        dst_id = data.get("target")
        if src_id is None or dst_id is None:
            continue

        src_info = id_to_info.get(str(src_id))
        dst_info = id_to_info.get(str(dst_id))

        # If either endpoint node is missing in node list => unknown
        if src_info is None or dst_info is None:
            if relation == "timing":
                agg_by_chr[fallback_chr_label].timing["total"] += 1
                agg_by_chr[fallback_chr_label].timing["unknown"] += 1
                agg_all.timing["total"] += 1
                agg_all.timing["unknown"] += 1
            else:
                agg_by_chr[fallback_chr_label].cooccur["total"] += 1
                agg_by_chr[fallback_chr_label].cooccur["unknown"] += 1
                agg_all.cooccur["total"] += 1
                agg_all.cooccur["unknown"] += 1
            continue

        src_chr, src_pos, src_label = src_info
        dst_chr, dst_pos, dst_label = dst_info

        # choose chromosome bucket:
        # if both endpoints have a chromosome and they match, use that; else fallback to folder label
        chr_bucket = src_chr if (src_chr == dst_chr and src_chr is not None) else fallback_chr_label

        known = (
            src_label is not None and dst_label is not None
            and src_label in tree_nodes and dst_label in tree_nodes
        )

        if relation == "timing":
            agg_by_chr[chr_bucket].timing["total"] += 1
            agg_all.timing["total"] += 1

            if not known:
                agg_by_chr[chr_bucket].timing["unknown"] += 1
                agg_all.timing["unknown"] += 1
                continue

            # known labels
            if src_label == dst_label:
                # incorrect timing within same node
                agg_by_chr[chr_bucket].timing["incorrect"] += 1
                agg_all.timing["incorrect"] += 1
                timing_pair_counts[(src_label, dst_label)]["incorrect"] += 1

                inconsistent_timing_rows.append((
                    chr_bucket, src_chr, src_pos, src_label,
                    dst_chr, dst_pos, dst_label,
                    "same_node",
                ))
            else:
                if is_ancestor_strict(src_label, dst_label, parent_of, depth):
                    agg_by_chr[chr_bucket].timing["correct"] += 1
                    agg_all.timing["correct"] += 1
                    timing_pair_counts[(src_label, dst_label)]["correct"] += 1
                else:
                    agg_by_chr[chr_bucket].timing["incorrect"] += 1
                    agg_all.timing["incorrect"] += 1
                    timing_pair_counts[(src_label, dst_label)]["incorrect"] += 1

                    if is_ancestor_strict(dst_label, src_label, parent_of, depth):
                        reason = "reverse_direction"
                    else:
                        reason = "no_ancestry_relation"

                    inconsistent_timing_rows.append((
                        chr_bucket, src_chr, src_pos, src_label,
                        dst_chr, dst_pos, dst_label,
                        reason,
                    ))

        else:  # cooccurring
            agg_by_chr[chr_bucket].cooccur["total"] += 1
            agg_all.cooccur["total"] += 1

            if not known:
                agg_by_chr[chr_bucket].cooccur["unknown"] += 1
                agg_all.cooccur["unknown"] += 1
                continue

            # unordered pair key for counts
            a, b = sorted((src_label, dst_label), key=nlabel_sort_key)

            if src_label == dst_label:
                agg_by_chr[chr_bucket].cooccur["correct"] += 1
                agg_all.cooccur["correct"] += 1
                cooccur_pair_counts[(a, b)]["correct"] += 1
            else:
                agg_by_chr[chr_bucket].cooccur["incorrect"] += 1
                agg_all.cooccur["incorrect"] += 1
                cooccur_pair_counts[(a, b)]["incorrect"] += 1

                inconsistent_cooccur_rows.append((
                    chr_bucket, src_chr, src_pos, src_label,
                    dst_chr, dst_pos, dst_label,
                ))


# ---------------- Output writing ----------------

def pct(correct: int, total: int, unknown: int) -> str:
    known = total - unknown
    if known <= 0:
        return "NA"
    return f"{(100.0 * correct / known):.2f}"


def write_summary_tsv(path: Path, agg_by_chr: Dict[str, ChrAgg], agg_all: ChrAgg) -> None:
    hdr = [
        "chrom",
        "timing_total", "timing_correct", "timing_incorrect", "timing_unknown", "timing_pct_consistent",
        "cooccur_total", "cooccur_correct", "cooccur_incorrect", "cooccur_unknown", "cooccur_pct_consistent",
    ]
    with path.open("w", encoding="utf-8") as f:
        f.write("\t".join(hdr) + "\n")

        for chrom in sorted(agg_by_chr.keys(), key=chrom_sort_key):
            a = agg_by_chr[chrom]
            f.write("\t".join([
                chrom,
                str(a.timing["total"]), str(a.timing["correct"]), str(a.timing["incorrect"]), str(a.timing["unknown"]),
                pct(a.timing["correct"], a.timing["total"], a.timing["unknown"]),
                str(a.cooccur["total"]), str(a.cooccur["correct"]), str(a.cooccur["incorrect"]), str(a.cooccur["unknown"]),
                pct(a.cooccur["correct"], a.cooccur["total"], a.cooccur["unknown"]),
            ]) + "\n")

        A = agg_all
        f.write("\t".join([
            "ALL",
            str(A.timing["total"]), str(A.timing["correct"]), str(A.timing["incorrect"]), str(A.timing["unknown"]),
            pct(A.timing["correct"], A.timing["total"], A.timing["unknown"]),
            str(A.cooccur["total"]), str(A.cooccur["correct"]), str(A.cooccur["incorrect"]), str(A.cooccur["unknown"]),
            pct(A.cooccur["correct"], A.cooccur["total"], A.cooccur["unknown"]),
        ]) + "\n")


def write_inconsistent_timing_edges(path: Path, rows: List[Tuple]) -> None:
    hdr = [
        "chrom_bucket",
        "src_chrom", "src_pos", "src_tree_node",
        "dst_chrom", "dst_pos", "dst_tree_node",
        "reason",
    ]
    # sort by chromosome then coordinates
    rows_sorted = sorted(
        rows,
        key=lambda r: (chrom_sort_key(str(r[0])), int(r[2]), int(r[5]), nlabel_sort_key(str(r[3])), nlabel_sort_key(str(r[6]))),
    )
    with path.open("w", encoding="utf-8") as f:
        f.write("\t".join(hdr) + "\n")
        for r in rows_sorted:
            f.write("\t".join(map(str, r)) + "\n")


def write_inconsistent_cooccur_edges(path: Path, rows: List[Tuple]) -> None:
    hdr = [
        "chrom_bucket",
        "src_chrom", "src_pos", "src_tree_node",
        "dst_chrom", "dst_pos", "dst_tree_node",
    ]
    rows_sorted = sorted(
        rows,
        key=lambda r: (chrom_sort_key(str(r[0])), int(r[2]), int(r[5]), nlabel_sort_key(str(r[3])), nlabel_sort_key(str(r[6]))),
    )
    with path.open("w", encoding="utf-8") as f:
        f.write("\t".join(hdr) + "\n")
        for r in rows_sorted:
            f.write("\t".join(map(str, r)) + "\n")


def write_timing_nodepair_counts(path: Path, counts: Dict[Tuple[str, str], Counter]) -> None:
    hdr = ["src_tree_node", "dst_tree_node", "consistent", "inconsistent"]
    rows = []
    for (src, dst), c in counts.items():
        rows.append((src, dst, int(c["correct"]), int(c["incorrect"])))
    rows_sorted = sorted(rows, key=lambda r: (nlabel_sort_key(r[0]), nlabel_sort_key(r[1])))
    with path.open("w", encoding="utf-8") as f:
        f.write("\t".join(hdr) + "\n")
        for r in rows_sorted:
            f.write("\t".join(map(str, r)) + "\n")


def write_cooccur_nodepair_counts(path: Path, counts: Dict[Tuple[str, str], Counter]) -> None:
    hdr = ["tree_node_a", "tree_node_b", "consistent", "inconsistent"]
    rows = []
    for (a, b), c in counts.items():
        rows.append((a, b, int(c["correct"]), int(c["incorrect"])))
    rows_sorted = sorted(rows, key=lambda r: (nlabel_sort_key(r[0]), nlabel_sort_key(r[1])))
    with path.open("w", encoding="utf-8") as f:
        f.write("\t".join(hdr) + "\n")
        for r in rows_sorted:
            f.write("\t".join(map(str, r)) + "\n")


# ---------------- Main ----------------

def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    tree_path = Path(args.tree)

    parent_of, tree_nodes = read_tree_parent_child_tsv(tree_path)
    depth = build_depth_cache(parent_of)

    agg_by_chr: Dict[str, ChrAgg] = defaultdict(ChrAgg)
    agg_all = ChrAgg()

    inconsistent_timing_rows: List[Tuple] = []
    inconsistent_cooccur_rows: List[Tuple] = []

    # genome-wide nodepair counts (known edges only)
    timing_pair_counts: Dict[Tuple[str, str], Counter] = defaultdict(Counter)   # (src, dst)
    cooccur_pair_counts: Dict[Tuple[str, str], Counter] = defaultdict(Counter) # (a, b) unordered

    out_folders = [p for p in outdir.iterdir() if p.is_dir() and p.name.endswith("_out")]
    if not out_folders:
        raise SystemExit(f"No *_out folders found under {outdir}")

    for of in sorted(out_folders, key=lambda p: chrom_sort_key(chromosome_from_out_folder(p.name))):
        chr_label = chromosome_from_out_folder(of.name)
        graphs_dir = of / "graphs"
        if not graphs_dir.is_dir():
            continue

        for jp in sorted(graphs_dir.rglob("*.json")):
            if should_skip_json(jp):
                continue
            try:
                eval_graph_file(
                    json_path=jp,
                    fallback_chr_label=chr_label,
                    parent_of=parent_of,
                    tree_nodes=tree_nodes,
                    depth=depth,
                    agg_by_chr=agg_by_chr,
                    agg_all=agg_all,
                    inconsistent_timing_rows=inconsistent_timing_rows,
                    inconsistent_cooccur_rows=inconsistent_cooccur_rows,
                    timing_pair_counts=timing_pair_counts,
                    cooccur_pair_counts=cooccur_pair_counts,
                )
            except json.JSONDecodeError as e:
                raise SystemExit(f"JSON parse error in {jp}: {e}") from e
            except Exception as e:
                raise SystemExit(f"Error processing {jp}: {e}") from e

    prefix = args.out_prefix

    summary_path = outdir / f"{prefix}_summary.tsv"
    incons_timing_path = outdir / f"{prefix}_inconsistent_timing_edges.tsv"
    incons_cooccur_path = outdir / f"{prefix}_inconsistent_cooccur_edges.tsv"
    timing_nodepair_path = outdir / f"{prefix}_timing_nodepair_counts.tsv"
    cooccur_nodepair_path = outdir / f"{prefix}_cooccur_nodepair_counts.tsv"

    write_summary_tsv(summary_path, agg_by_chr, agg_all)
    write_inconsistent_timing_edges(incons_timing_path, inconsistent_timing_rows)
    write_inconsistent_cooccur_edges(incons_cooccur_path, inconsistent_cooccur_rows)
    write_timing_nodepair_counts(timing_nodepair_path, timing_pair_counts)
    write_cooccur_nodepair_counts(cooccur_nodepair_path, cooccur_pair_counts)

    print("Wrote:")
    print(f"  {summary_path}")
    print(f"  {incons_timing_path}")
    print(f"  {incons_cooccur_path}")
    print(f"  {timing_nodepair_path}")
    print(f"  {cooccur_nodepair_path}")


if __name__ == "__main__":
    main()