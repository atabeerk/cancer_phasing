import os
import json
import csv
import argparse
from pathlib import Path
from collections import defaultdict
import networkx as nx


def read_edges_from_file(file_path, reverse=False):
    """Read edges, optionally reversed, with TOTAL and ORIGINAL metadata."""
    edges = set()
    nodes = set()
    edge_data = {}  # key: (chrom, pos1, pos2), value: {"original": str, "total": int}

    with open(file_path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            chrom, pos1, pos2 = parts[0], parts[1], parts[2]
            if reverse:
                pos1, pos2 = pos2, pos1
            total = None
            original = None
            for token in parts[3:]:
                if token.startswith("TOTAL="):
                    total = int(token.split("=")[1])
                elif token.startswith("ORIGINAL="):
                    original = token.split("=")[1]
            edges.add((chrom, pos1, pos2))
            nodes.update([pos1, pos2])
            edge_data[(chrom, pos1, pos2)] = {"total": total, "original": original}
    return nodes, edges, edge_data


def build_graphs(chunk_dir):
    """Group files by chunk and type; build nodes, edges, and metadata."""
    chunk_dir = Path(chunk_dir)
    graphs = {}
    grouped = defaultdict(dict)

    for f in chunk_dir.glob("chunk_*_*.txt"):
        name = f.stem
        if "_snp1_before_snp2_loss" in name:
            grouped[name.replace("_snp1_before_snp2_loss", "")]["forward_loss"] = f
        elif "_snp2_before_snp1_loss" in name:
            grouped[name.replace("_snp2_before_snp1_loss", "")]["reverse_loss"] = f
        elif "_snp1_before_snp2" in name:
            grouped[name.replace("_snp1_before_snp2", "")]["forward"] = f
        elif "_snp2_before_snp1" in name:
            grouped[name.replace("_snp2_before_snp1", "")]["reverse"] = f
        elif "_cooccurring" in name:
            grouped[name.replace("_cooccurring", "")]["cooccurring"] = f
        elif "_divergent" in name:
            grouped[name.replace("_divergent", "")]["divergent"] = f

    for base, files in grouped.items():
        nodes, directed_edges, directed_loss_edges = set(), set(), set()
        cooccurring_edges, divergent_edges = set(), set()
        edge_meta = {}

        # Directed edges
        for key, rev in [("forward", False), ("reverse", True)]:
            if key in files:
                n, e, meta = read_edges_from_file(files[key], reverse=rev)
                nodes |= n
                directed_edges |= e
                edge_meta.update(meta)
        # Directed loss
        for key, rev in [("forward_loss", False), ("reverse_loss", True)]:
            if key in files:
                n, e, meta = read_edges_from_file(files[key], reverse=rev)
                nodes |= n
                directed_loss_edges |= e
                edge_meta.update(meta)
        # Cooccurring
        if "cooccurring" in files:
            n, e, meta = read_edges_from_file(files["cooccurring"])
            nodes |= n
            cooccurring_edges |= e
            edge_meta.update(meta)
        # Divergent (only if both nodes exist)
        if "divergent" in files:
            n, e, meta = read_edges_from_file(files["divergent"])
            for chrom, s, t in e:
                if s in nodes and t in nodes:
                    divergent_edges.add((chrom, s, t))
                    edge_meta[(chrom, s, t)] = meta.get((chrom, s, t), {})

        graphs[base] = {
            "nodes": nodes,
            "directed": directed_edges,
            "directed_loss": directed_loss_edges,
            "cooccurring": cooccurring_edges,
            "divergent": divergent_edges,
            "edge_meta": edge_meta
        }
    return graphs


def compute_component_statistics(nodes, directed_edges, directed_loss_edges,
                                 cooccurring_edges, divergent_edges):
    G_all = nx.Graph()
    G_all.add_nodes_from(nodes)
    G_all.add_edges_from([(s, t) for _, s, t in directed_edges])
    G_all.add_edges_from([(s, t) for _, s, t in directed_loss_edges])
    G_all.add_edges_from([(s, t) for _, s, t in cooccurring_edges])
    G_all.add_edges_from([(s, t) for _, s, t in divergent_edges])

    component_stats, node_to_component = [], {}
    for comp_id, comp_nodes in enumerate(nx.connected_components(G_all)):
        comp_nodes = list(comp_nodes)
        positions = sorted(int(n) for n in comp_nodes if n.isdigit())
        span = max(positions) - min(positions) if positions else 0
        G_co = nx.Graph()
        for chrom, s, t in cooccurring_edges:
            if s in comp_nodes and t in comp_nodes:
                G_co.add_edge(s, t)
        for n in comp_nodes:
            if n not in G_co:
                G_co.add_node(n)
        haplotypes = nx.number_connected_components(G_co)
        for n in comp_nodes:
            node_to_component[n] = comp_id
        component_stats.append({
            "component_id": comp_id,
            "nodes": comp_nodes,
            "span_bp": span,
            "haplotypes": haplotypes
        })
    return component_stats, node_to_component


def write_cytoscape_json(nodes, directed_edges, directed_loss_edges,
                         cooccurring_edges, divergent_edges,
                         component_stats, node_to_component,
                         edge_meta, out_file):
    elements = {"nodes": [], "edges": []}
    for n in sorted(nodes):
        elements["nodes"].append({"data": {"id": n, "label": n,
                                           "component_id": node_to_component.get(n, -1)}})
    def add_edges(edges, edge_type, directed_flag):
        for chrom, s, t in edges:
            meta = edge_meta.get((chrom, s, t), {})
            elements["edges"].append({
                "data": {
                    "source": s,
                    "target": t,
                    "chrom": chrom,
                    "directed": directed_flag,
                    "edge_type": edge_type,
                    "total": meta.get("total"),
                    "original": meta.get("original")
                }
            })
    add_edges(directed_edges, "directed", True)
    add_edges(directed_loss_edges, "directed_loss", True)
    add_edges(cooccurring_edges, "cooccurring", False)
    add_edges(divergent_edges, "divergent", False)
    with open(out_file, "w") as out:
        json.dump({"elements": elements, "component_stats": component_stats}, out, indent=2)


def write_component_statistics(out_dir, base, component_stats):
    stats_path = Path(out_dir) / "component_statistics.txt"
    write_header = not stats_path.exists()
    with open(stats_path, "a") as f:
        if write_header:
            f.write("base\tcomponent_id\tnum_nodes\tspan_bp\thaplotypes\tnodes\n")
        for comp in component_stats:
            nodes_str = ",".join(sorted(comp["nodes"]))
            f.write(f"{base}\t{comp['component_id']}\t{len(comp['nodes'])}\t"
                    f"{comp['span_bp']}\t{comp['haplotypes']}\t{nodes_str}\n")


def condense_graph(nodes, directed_edges, directed_loss_edges,
                   cooccurring_edges, divergent_edges):
    """Collapse cooccurrence clusters into CNs and aggregate edges between them."""
    G_co = nx.Graph()
    G_co.add_nodes_from(nodes)
    G_co.add_edges_from([(s, t) for _, s, t in cooccurring_edges])

    cn_list, node_to_cn = [], {}
    for cn_id, cluster in enumerate(nx.connected_components(G_co)):
        cluster = list(cluster)
        cn_name = f"CN{cn_id}"
        cn_list.append({
            "id": cn_name,
            "members": cluster,
            "n_nodes": len(cluster),
            "span": max(map(int, cluster)) - min(map(int, cluster)) if all(n.isdigit() for n in cluster) else None,
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
    cn_edges = {"directed": set(), "directed_loss": set(), "divergent": set()}
    def add_cn_edge(edge_type, s, t, chrom):
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


def write_condensed_json_and_tsv(cn_nodes, cn_edges, cn_flags, out_base):
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
    with open(f"{out_base}_condensed.json", "w") as jf:
        json.dump({"elements": elements}, jf, indent=2)

    with open(f"{out_base}_condensed.stats.tsv", "w", newline="") as tf:
        writer = csv.writer(tf, delimiter="\t")
        writer.writerow(["CN_id", "num_nodes", "span_bp", "nodes_list", "mixed_edges"])
        for cn_id, members in cn_nodes.items():
            positions = sorted(int(n) for n in members if n.isdigit())
            span = max(positions) - min(positions) if positions else 0
            writer.writerow([cn_id, len(members), span, ",".join(sorted(members)), cn_flags.get(cn_id, False)])
    print(f"Condensed JSON/TSV written: {out_base}_condensed.*")


def process_graph(base, g, out_dir):
    nodes, directed, directed_loss = g["nodes"], g["directed"], g["directed_loss"]
    cooccurring, divergent, edge_meta = g["cooccurring"], g["divergent"], g["edge_meta"]

    # --- Step 1: Compute component stats and write uncondensed JSON ---
    component_stats, node_to_component = compute_component_statistics(
        nodes, directed, directed_loss, cooccurring, divergent
    )
    out_file = Path(out_dir) / f"{base}.json"
    write_cytoscape_json(
        nodes, directed, directed_loss, cooccurring, divergent,
        component_stats, node_to_component, edge_meta, out_file
    )
    write_component_statistics(out_dir, base, component_stats)
    print(f"{out_file}: {len(nodes)} nodes, {len(component_stats)} components")

    # --- Step 2: Condense CN graph ---
    cn_list, cn_edges_dict = condense_graph(nodes, directed, directed_loss, cooccurring, divergent)
    cn_flags = {c["id"]: c["mixed_edges"] for c in cn_list}
    cn_nodes = {c["id"]: set(c["members"]) for c in cn_list}

    # Combine edges into simple list (chrom, src, dst, edge_type)
    cn_edges = []
    for edge_type, edges in cn_edges_dict.items():
        for chrom, src, dst in edges:
            cn_edges.append((src, dst, edge_type))

    write_condensed_json_and_tsv(cn_nodes, cn_edges, cn_flags, Path(out_dir) / base)


def main(chunk_dir, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    graphs = build_graphs(chunk_dir)
    for base, g in graphs.items():
        process_graph(base, g, out_dir)
    print(f"Component statistics written to {Path(out_dir) / 'component_statistics.txt'}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate Cytoscape JSON DAGs including TOTAL/ORIGINAL and condensed CN graphs."
    )
    parser.add_argument("chunk_dir", help="Directory containing chunk *_before_* and related files")
    parser.add_argument("--outdir", default="cytoscape_chunks", help="Output directory for JSON files")
    args = parser.parse_args()
    main(args.chunk_dir, args.outdir)
