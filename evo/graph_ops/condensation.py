from typing import Tuple, Dict, Set, List
import networkx as nx
from pathlib import Path
import json
import csv


Edge = Tuple[str, str, str]  # (chrom, src, dst)


def condense_graph(nodes: Set[str],
                   directed_edges: Set[Edge],
                   directed_loss_edges: Set[Edge],
                   cooccurring_edges: Set[Edge],
                   cooccurring_loss_edges: Set[Edge],
                   divergent_edges: Set[Edge]) -> Tuple[List[Dict], Dict[str, Set[Tuple[str, str, str]]]]:
    """Collapse cooccurrence clusters into CNs and aggregate edges between them.

    Returns (cn_list, cn_edges_dict) where cn_edges_dict maps edge_type -> set of (chrom, src_cn, dst_cn).
    """
    G_co = nx.Graph()
    G_co.add_nodes_from(nodes)
    G_co.add_edges_from([(s, t) for _, s, t in cooccurring_edges])
    G_co.add_edges_from([(s, t) for _, s, t in cooccurring_loss_edges])

    cn_list = []
    node_to_cn: Dict[str, str] = {}
    for cn_id, cluster in enumerate(nx.connected_components(G_co)):
        cluster = list(cluster)
        cn_name = f"CN{cn_id}"
        span = None
        if all(n.isdigit() for n in cluster):
            positions = [int(n) for n in cluster]
            span = max(positions) - min(positions) if positions else 0

        cn_list.append({
            "id": cn_name,
            "members": cluster,
            "n_nodes": len(cluster),
            "span": span,
            "mixed_edges": False
        })
        for n in cluster:
            node_to_cn[n] = cn_name

    # Flag CNs with internal directed/divergent edges
    for edge_set in [directed_edges, directed_loss_edges, divergent_edges]:
        for _, s, t in edge_set:
            if s in node_to_cn and t in node_to_cn and node_to_cn[s] == node_to_cn[t]:
                cn_name = node_to_cn[s]
                for c in cn_list:
                    if c["id"] == cn_name:
                        c["mixed_edges"] = True

    # Build CN-level edges
    cn_edges: Dict[str, Set[Tuple[str, str, str]]] = {"directed": set(), "directed_loss": set(), "divergent": set()}

    def add_cn_edge(edge_type: str, s: str, t: str, chrom: str):
        src, dst = node_to_cn.get(s), node_to_cn.get(t)
        if not src or not dst or src == dst:
            return
        cn_edges[edge_type].add((chrom, src, dst))

    for chrom, s, t in directed_edges:
        add_cn_edge("directed", s, t, chrom)
    for chrom, s, t in directed_loss_edges:
        add_cn_edge("directed_loss", s, t, chrom)
    for chrom, s, t in divergent_edges:
        add_cn_edge("divergent", s, t, chrom)

    return cn_list, cn_edges


def write_condensed_json_and_tsv(cn_nodes: Dict[str, Set[str]],
                                 cn_edges: List[Tuple[str, str, str]],
                                 cn_flags: Dict[str, bool],
                                 out_base: Path) -> None:
    """Write condensed JSON and TSV stats based on CN-level nodes and edges.

    cn_nodes: mapping CN_id -> set(members)
    cn_edges: list of tuples (src_cn, dst_cn, edge_type)
    cn_flags: mapping CN_id -> mixed_edges boolean
    out_base: Path-like (base path, e.g. Path(out_dir) / base)
    """
    elements = {"nodes": [], "edges": []}

    for cn_id, members in cn_nodes.items():
        elements["nodes"].append({
            "data": {"id": cn_id, "label": cn_id,
                     "size": len(members), "mixed_edges": cn_flags.get(cn_id, False)}
        })

    for src, tgt, edge_type in cn_edges:
        elements["edges"].append({
            "data": {"source": src, "target": tgt, "edge_type": edge_type,
                     "directed": edge_type in ("directed", "directed_loss")}
        })

    json_path = f"{out_base}_condensed.json"
    tsv_path = f"{out_base}_condensed.stats.tsv"

    with open(json_path, "w") as jf:
        json.dump({"elements": elements}, jf, indent=2)

    with open(tsv_path, "w", newline="") as tf:
        writer = csv.writer(tf, delimiter="\t")
        writer.writerow(["CN_id", "num_nodes", "span_bp", "nodes_list", "mixed_edges"])

        for cn_id, members in cn_nodes.items():
            positions = sorted(int(n) for n in members if n.isdigit())
            span = max(positions) - min(positions) if positions else 0
            writer.writerow([cn_id, len(members), span, ",".join(sorted(members)), cn_flags.get(cn_id, False)])

    print(f"Condensed JSON/TSV written: {json_path}, {tsv_path}")
