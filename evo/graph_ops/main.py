import os
from pathlib import Path
import argparse

from graph_building import build_graphs
from graph_analysis import compute_component_statistics
from io_utils import parse_pileup_depth, write_component_statistics
from cytoscape_export import write_cytoscape_json
from condensation import condense_graph, write_condensed_json_and_tsv


def process_graph(base: str, g, out_dir: str):
    """Process a single GraphData object: analyze, export, condense."""
    nodes = g.nodes
    directed = g.directed
    directed_loss = g.directed_loss
    cooccurring = g.cooccurring
    cooccurring_loss = g.cooccurring_loss
    divergent = g.divergent
    edge_meta = g.edge_meta

    # Step 1: Compute component stats & write uncondensed JSON
    component_stats, node_to_component = compute_component_statistics(g)

    # Determine pileup path (best-effort -- keep old naming convention)
    try:
        chrom, start = base.split("_")[1:3]
    except Exception:
        chrom, start = ("unknown_chr", "0")

    pileup_file = Path(out_dir).parent / "pileup_files" / f"mpileup_{chrom}_{start}.out"
    pileup_depths = parse_pileup_depth(pileup_file)

    out_file = Path(out_dir) / f"{base}.json"
    write_cytoscape_json(
        nodes, directed, directed_loss, cooccurring, cooccurring_loss,
        divergent, component_stats, node_to_component, edge_meta, out_file,
        pileup_depths
    )
    write_component_statistics(Path(out_dir), base, component_stats)
    print(f"{out_file}: {len(nodes)} nodes, {len(component_stats)} components")

    # Step 2: Condense CN graph
    cn_list, cn_edges_dict = condense_graph(nodes, directed, directed_loss, cooccurring, cooccurring_loss, divergent)
    cn_flags = {c["id"]: c["mixed_edges"] for c in cn_list}
    cn_nodes = {c["id"]: set(c["members"]) for c in cn_list}

    # Flatten CN edges into required format (src, dst, edge_type)
    cn_edges = []
    for edge_type, edges in cn_edges_dict.items():
        for chrom, src, dst in edges:
            cn_edges.append((src, dst, edge_type))

    write_condensed_json_and_tsv(cn_nodes, cn_edges, cn_flags, Path(out_dir) / base)


def main(chunk_dir: str, out_dir: str):
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
