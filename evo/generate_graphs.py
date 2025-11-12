import os
import json
from pathlib import Path
from collections import defaultdict
import networkx as nx
import csv

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


def compute_component_statistics(nodes, directed_edges, directed_loss_edges,
                                 cooccurring_edges, divergent_edges):
    """Compute per-component statistics: region span and estimated haplotypes."""
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

        # Build cooccurrence-only subgraph within this component
        G_co = nx.Graph()
        for chrom, s, t in cooccurring_edges:
            if s in comp_nodes and t in comp_nodes:
                G_co.add_edge(s, t)

        haplotypes = (
            nx.number_connected_components(G_co)
            if len(G_co.nodes) > 0
            else len(comp_nodes)
        )

        for n in comp_nodes:
            node_to_component[n] = comp_id

        component_stats.append({
            "component_id": comp_id,
            "nodes": comp_nodes,
            "span_bp": span,
            "haplotypes": haplotypes
        })

    return component_stats, node_to_component


def condense_graph(nodes, directed_edges, directed_loss_edges,
                   cooccurring_edges, divergent_edges):
    """Collapse cooccurrence clusters into CNs and aggregate edges between them."""
    G_co = nx.Graph()
    G_co.add_edges_from([(s, t) for _, s, t in cooccurring_edges])
    G_co.add_nodes_from(nodes)

    cn_list = []
    node_to_cn = {}

    for cn_id, cluster in enumerate(nx.connected_components(G_co)):
        cluster = list(cluster)
        cn_name = f"CN{cn_id}"
        cn_list.append({
            "id": cn_name,
            "members": cluster,
            "n_nodes": len(cluster),
            "span": (
                max(map(int, cluster)) - min(map(int, cluster))
                if all(n.isdigit() for n in cluster) else None
            ),
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
    cn_edges = {
        "directed": set(),
        "directed_loss": set(),
        "divergent": set()
    }

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
    """
    Writes:
      - Condensed Cytoscape JSON: <out_base>_condensed.json
      - Condensed stats TSV: <out_base>_condensed.stats.tsv
    cn_nodes: dict of CN_id -> set of original nodes
    cn_edges: list of tuples (source_CN, target_CN, edge_type)
    cn_flags: dict CN_id -> bool (True if CN has mixed edges)
    out_base: Path or str, base filename
    """
    # --- JSON for Cytoscape ---
    elements = {"nodes": [], "edges": []}
    
    # Add CN nodes
    for cn_id, members in cn_nodes.items():
        elements["nodes"].append({
            "data": {
                "id": cn_id,
                "label": cn_id,
                "size": len(members),
                "mixed_edges": cn_flags.get(cn_id, False)
            }
        })
    
    # Add edges
    for src, tgt, edge_type in cn_edges:
        elements["edges"].append({
            "data": {
                "source": src,
                "target": tgt,
                "edge_type": edge_type,
                "directed": edge_type in ("directed", "directed_loss")
            }
        })
    
    json_file = f"{out_base}_condensed.json"
    with open(json_file, "w") as jf:
        json.dump({"elements": elements}, jf, indent=2)
    
    # --- TSV stats ---
    tsv_file = f"{out_base}_condensed.stats.tsv"
    with open(tsv_file, "w", newline="") as tf:
        writer = csv.writer(tf, delimiter="\t")
        writer.writerow(["CN_id", "num_nodes", "span_bp", "nodes_list", "mixed_edges"])
        for cn_id, members in cn_nodes.items():
            positions = sorted(int(n) for n in members if n.isdigit())
            span = max(positions) - min(positions) if positions else 0
            writer.writerow([
                cn_id,
                len(members),
                span,
                ",".join(sorted(members)),
                cn_flags.get(cn_id, False)
            ])
    
    print(f"Condensed graph JSON: {json_file}, stats TSV: {tsv_file}")


def main(chunk_dir, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    graphs = build_graphs(chunk_dir)

    for base, g in graphs.items():
        nodes = g["nodes"]
        directed = g["directed"]
        directed_loss = g["directed_loss"]
        cooccurring = g["cooccurring"]
        divergent = g["divergent"]

        # --- Step 1: write uncondensed graph JSON ---
        out_file = Path(out_dir) / f"{base}.json"
        # Compute per-component stats and node-component mapping
        component_stats, node_to_component = compute_component_statistics(
            nodes, directed, directed_loss, cooccurring, divergent
        )

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
        print(f"{out_file}: {len(nodes)} nodes, {len(component_stats)} components (uncondensed)")

        # --- Step 2: compute CNs (condensed graph) ---
        # Use cooccurrence edges to form clusters (CNs)
        G_co = nx.Graph()
        G_co.add_nodes_from(nodes)
        for _, s, t in cooccurring:
            G_co.add_edge(s, t)
        # Each connected component in G_co becomes a CN
        cn_nodes = {}  # CN_id -> set of member nodes
        cn_flags = {}  # CN_id -> bool (True if any mixed edges inside)
        for cn_id, comp_nodes in enumerate(nx.connected_components(G_co)):
            comp_nodes = set(comp_nodes)
            cn_nodes[f"CN{cn_id}"] = comp_nodes

            # Check for mixed edges: any directed, directed_loss, or divergent inside CN
            mixed = False
            for s in comp_nodes:
                for t in comp_nodes:
                    if ((s, t) in [(x[1], x[2]) for x in directed] or
                        (s, t) in [(x[1], x[2]) for x in directed_loss] or
                        (s, t) in [(x[1], x[2]) for x in divergent]):
                        mixed = True
                        break
                if mixed:
                    break
            cn_flags[f"CN{cn_id}"] = mixed

        # --- Step 3: compute edges between CNs ---
        cn_edges = set()
        node_to_cn = {}
        for cn_id, members in cn_nodes.items():
            for n in members:
                node_to_cn[n] = cn_id

        # All inter-CN edges (directed, directed_loss, divergent)
        for edge_type, edge_set in [("directed", directed), ("directed_loss", directed_loss), ("divergent", divergent)]:
            for _, s, t in edge_set:
                if s in node_to_cn and t in node_to_cn:
                    cn_s, cn_t = node_to_cn[s], node_to_cn[t]
                    if cn_s != cn_t:
                        cn_edges.add((cn_s, cn_t, edge_type))

        # --- Step 4: write condensed graph JSON + TSV ---
        out_base = Path(out_dir) / base
        write_condensed_json_and_tsv(cn_nodes, cn_edges, cn_flags, out_base)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Generate Cytoscape JSON DAGs and condensed CN-level graphs."
    )
    parser.add_argument("chunk_dir", help="Directory containing chunk *_before_* and related files")
    parser.add_argument("--outdir", default="cytoscape_chunks", help="Output directory for JSON files")
    args = parser.parse_args()

    main(args.chunk_dir, args.outdir)
