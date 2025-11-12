import os
import json
import argparse
from pathlib import Path
from collections import defaultdict

import networkx as nx


def read_edges_from_file(file_path, reverse=False):
    """
    Reads formatted lines like:
      chr1 100 200 ALT_ALT=10 ALT_REF=5 REF_ALT=7 REF_REF=9 TOTAL=31 ORIGINAL=12/6/8/10
    Returns:
      nodes, edges
    Each edge now includes extra metadata: original_counts and total.
    """
    edges = set()
    nodes = set()
    edge_data = {}  # key: (chrom, pos1, pos2), value: {"original": str, "total": int}

    with open(file_path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue

            parts = line.strip().split()
            if len(parts) < 4:
                continue

            chrom, pos1, pos2 = parts[0], parts[1], parts[2]
            if reverse:
                pos1, pos2 = pos2, pos1

            # Extract fields
            total = None
            original = None
            for token in parts[3:]:
                if token.startswith("TOTAL="):
                    total = int(token.split("=")[1])
                elif token.startswith("ORIGINAL="):
                    original = token.split("=")[1]

            edges.add((chrom, pos1, pos2))
            nodes.update([pos1, pos2])
            edge_data[(chrom, pos1, pos2)] = {
                "original": original,
                "total": total
            }

    return nodes, edges, edge_data


def build_graphs(chunk_dir):
    chunk_dir = Path(chunk_dir)
    graphs = {}

    grouped = defaultdict(dict)
    for f in chunk_dir.glob("chunk_*_*.txt"):
        name = f.stem
        if "_snp1_before_snp2_loss" in name:
            base = name.replace("_snp1_before_snp2_loss", "")
            grouped[base]["forward_loss"] = f
        elif "_snp2_before_snp1_loss" in name:
            base = name.replace("_snp2_before_snp1_loss", "")
            grouped[base]["reverse_loss"] = f
        elif "_snp1_before_snp2" in name:
            base = name.replace("_snp1_before_snp2", "")
            grouped[base]["forward"] = f
        elif "_snp2_before_snp1" in name:
            base = name.replace("_snp2_before_snp1", "")
            grouped[base]["reverse"] = f
        elif "_cooccurring" in name:
            base = name.replace("_cooccurring", "")
            grouped[base]["cooccurring"] = f
        elif "_divergent" in name:
            base = name.replace("_divergent", "")
            grouped[base]["divergent"] = f

    for base, files in grouped.items():
        nodes = set()
        directed_edges = set()
        directed_loss_edges = set()
        cooccurring_edges = set()
        divergent_edges = set()

        edge_meta = {}  # combined edge metadata for all types

        # Directed strong edges
        for key, reverse in [("forward", False), ("reverse", True)]:
            if key in files:
                n, e, meta = read_edges_from_file(files[key], reverse=reverse)
                nodes |= n
                directed_edges |= {(chrom, s, t) for chrom, s, t in e}
                edge_meta.update(meta)

        # Directed loss edges
        for key, reverse in [("forward_loss", False), ("reverse_loss", True)]:
            if key in files:
                n, e, meta = read_edges_from_file(files[key], reverse=reverse)
                nodes |= n
                directed_loss_edges |= {(chrom, s, t) for chrom, s, t in e}
                edge_meta.update(meta)

        # Cooccurring edges
        if "cooccurring" in files:
            n, e, meta = read_edges_from_file(files["cooccurring"])
            nodes |= n
            cooccurring_edges |= e
            edge_meta.update(meta)

        # Divergent edges
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

    component_stats = []
    node_to_component = {}

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
                         cooccurring_edges, divergent_edges, component_stats,
                         node_to_component, edge_meta, out_file):
    """Write combined graph with TOTAL and ORIGINAL counts per edge."""
    elements = {"nodes": [], "edges": []}

    # Nodes
    for n in sorted(nodes):
        elements["nodes"].append({
            "data": {
                "id": n,
                "label": n,
                "component_id": node_to_component.get(n, -1)
            }
        })

    # Helper to add edges with metadata
    def add_edges(edge_list, edge_type, directed_flag):
        for chrom, s, t in edge_list:
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

    # Add all edge types
    add_edges(directed_edges, "directed", True)
    add_edges(directed_loss_edges, "directed_loss", True)
    add_edges(cooccurring_edges, "cooccurring", False)
    add_edges(divergent_edges, "divergent", False)

    with open(out_file, "w") as out:
        json.dump({
            "elements": elements,
            "component_stats": component_stats
        }, out, indent=2)


def write_component_statistics(out_dir, base, component_stats):
    parent_dir = Path(out_dir).parent
    stats_path = parent_dir / "component_statistics.txt"

    write_header = not stats_path.exists()
    with open(stats_path, "a") as stats_out:
        if write_header:
            stats_out.write(
                "base\tcomponent_id\tnum_nodes\tspan_bp\thaplotypes\tnodes\n"
            )
        for comp in component_stats:
            node_list = ",".join(sorted(comp["nodes"]))
            stats_out.write(
                f"{base}\t{comp['component_id']}\t{len(comp['nodes'])}\t"
                f"{comp['span_bp']}\t{comp['haplotypes']}\t{node_list}\n"
            )


def main(chunk_dir, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    graphs = build_graphs(chunk_dir)

    for base, g in graphs.items():
        nodes = g["nodes"]
        directed = g["directed"]
        directed_loss = g["directed_loss"]
        cooccurring = g["cooccurring"]
        divergent = g["divergent"]
        edge_meta = g["edge_meta"]

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

    print(f"Component statistics written to {Path(out_dir) / 'component_statistics.txt'}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate Cytoscape JSON DAGs including TOTAL and ORIGINAL edge data."
    )
    parser.add_argument("chunk_dir", help="Directory containing chunk *_before_* and related files")
    parser.add_argument("--outdir", default="cytoscape_chunks", help="Output directory for JSON files")
    args = parser.parse_args()

    main(args.chunk_dir, args.outdir)
