import os
import json
from pathlib import Path
from collections import defaultdict


def read_edges_from_file(file_path, reverse=False):
    """Reads (chrom, pos1, pos2) → returns node set and edge set (tuples of str)."""
    edges = set()
    nodes = set()

    with open(file_path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue

            chrom, pos1, pos2 = parts[0], parts[1], parts[2]
            if reverse:
                pos1, pos2 = pos2, pos1
            edges.add((chrom, pos1, pos2))
            nodes.update([pos1, pos2])

    return nodes, edges


def build_graphs(chunk_dir):
    """
    For each chunk:
      - Read *_before_* and *_loss_* files to create directed edges
      - Read *_cooccurring.txt to create undirected edges (include all)
      - Read *_divergent.txt to create divergent edges (red)
    """
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

        # Directed strong edges
        for key, reverse in [("forward", False), ("reverse", True)]:
            if key in files:
                n, e = read_edges_from_file(files[key], reverse=reverse)
                nodes |= n
                directed_edges |= {(chrom, s, t) for chrom, s, t in e}

        # Directed loss edges
        for key, reverse in [("forward_loss", False), ("reverse_loss", True)]:
            if key in files:
                n, e = read_edges_from_file(files[key], reverse=reverse)
                nodes |= n
                directed_loss_edges |= {(chrom, s, t) for chrom, s, t in e}

        # Cooccurring edges (include all)
        if "cooccurring" in files:
            n, e = read_edges_from_file(files["cooccurring"], reverse=False)
            nodes |= n
            cooccurring_edges |= e

        # Divergent edges (only if both exist in directed nodes)
        if "divergent" in files:
            n, e = read_edges_from_file(files["divergent"], reverse=False)
            for chrom, s, t in e:
                if s in nodes and t in nodes:
                    divergent_edges.add((chrom, s, t))

        graphs[base] = {
            "nodes": nodes,
            "directed": directed_edges,
            "directed_loss": directed_loss_edges,
            "cooccurring": cooccurring_edges,
            "divergent": divergent_edges
        }

    return graphs


import networkx as nx


import networkx as nx

def compute_component_statistics(nodes, directed_edges, directed_loss_edges,
                                 cooccurring_edges, divergent_edges):
    """Compute per-component statistics: region span and estimated haplotypes."""
    # Build full undirected graph (for component detection)
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

        # --- Compute genomic span ---
        positions = sorted(int(n) for n in comp_nodes if n.isdigit())
        span = max(positions) - min(positions) if positions else 0

        # --- Build cooccurrence-based subgraph ---
        G_co = nx.Graph()
        for chrom, s, t in cooccurring_edges:
            if s in comp_nodes and t in comp_nodes:
                G_co.add_edge(s, t)

        # --- Add singleton nodes (no cooccurrence edges) ---
        for n in comp_nodes:
            if n not in G_co:
                G_co.add_node(n)

        # --- Haplotype count = connected components in cooccurrence subgraph ---
        haplotypes = nx.number_connected_components(G_co)

        # --- Record component mapping and stats ---
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
                         node_to_component, out_file):
    """Write combined graph with component metadata as Cytoscape-compatible JSON."""
    elements = {
        "nodes": [],
        "edges": []
    }

    # Include node component metadata
    for n in sorted(nodes):
        elements["nodes"].append({
            "data": {
                "id": n,
                "label": n,
                "component_id": node_to_component.get(n, -1)
            }
        })

    # Directed edges
    for chrom, s, t in directed_edges:
        elements["edges"].append({
            "data": {
                "source": s, "target": t,
                "chrom": chrom,
                "directed": True,
                "edge_type": "directed"
            }
        })

    # Directed loss
    for chrom, s, t in directed_loss_edges:
        elements["edges"].append({
            "data": {
                "source": s, "target": t,
                "chrom": chrom,
                "directed": True,
                "edge_type": "directed_loss"
            }
        })

    # Cooccurring
    for chrom, s, t in cooccurring_edges:
        elements["edges"].append({
            "data": {
                "source": s, "target": t,
                "chrom": chrom,
                "directed": False,
                "edge_type": "cooccurring"
            }
        })

    # Divergent
    for chrom, s, t in divergent_edges:
        elements["edges"].append({
            "data": {
                "source": s, "target": t,
                "chrom": chrom,
                "directed": False,
                "edge_type": "divergent"
            }
        })

    # Write everything (including stats)
    with open(out_file, "w") as out:
        json.dump({
            "elements": elements,
            "component_stats": component_stats
        }, out, indent=2)


def write_component_statistics(out_dir, base, component_stats):
    """
    Append per-component statistics to 'component_statistics.txt' located
    alongside the output directory.

    Columns:
        base, component_id, num_nodes, span_bp, haplotypes, nodes
    """
    # Get the parent directory of the output directory
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

        # Compute per-component stats and node-component mapping
        component_stats, node_to_component = compute_component_statistics(
            nodes, directed, directed_loss, cooccurring, divergent
        )

        # Write Cytoscape JSON
        out_file = Path(out_dir) / f"{base}.json"
        write_cytoscape_json(
            nodes,
            directed,
            directed_loss,
            cooccurring,
            divergent,
            component_stats,
            node_to_component,
            out_file
        )

        # New modularized call
        write_component_statistics(out_dir, base, component_stats)

        print(f"{out_file}: {len(nodes)} nodes, {len(component_stats)} components")

    print(f"Component statistics written to {Path(out_dir) / 'component_statistics.txt'}")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Generate Cytoscape JSON DAGs including *loss* and *divergent* edges."
    )
    parser.add_argument("chunk_dir", help="Directory containing chunk *_before_* and related files")
    parser.add_argument("--outdir", default="cytoscape_chunks", help="Output directory for JSON files")
    args = parser.parse_args()

    main(args.chunk_dir, args.outdir)
