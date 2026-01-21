import os
import sys
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
if CURRENT_DIR not in sys.path:
    sys.path.insert(0, CURRENT_DIR)

from parsing import find_chunk_bases, load_edges_from_base
from builder import GraphBuilder
from export_graph import export_cytoscape_json, export_condensed_cytoscape_json, write_inconsistency_log, write_accepted_edges
from component_stats import compute_component_statistics_rows, append_component_statistics_tsv



def main():
    parser = argparse.ArgumentParser(
        description=(
            "Build consistent SNV graphs for all chunk files in a directory.\n"
            "Each chunk's per-relation files ( *_cooccurring.txt etc. ) must "
            "already exist in the same directory."
        )
    )
    parser.add_argument(
        "chunk_dir",
        help="Directory containing chunk_*.txt and their per-relation outputs.",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Directory where Cytoscape JSON and inconsistency logs will be written.",
    )
    args = parser.parse_args()

    chunk_dir = os.path.abspath(args.chunk_dir)
    outdir = os.path.abspath(args.outdir)
    os.makedirs(outdir, exist_ok=True)

    # Per-chromosome component statistics next to the graphs/ folder.
    chrom_out_dir = os.path.dirname(outdir) if os.path.basename(outdir) == "graphs" else outdir
    component_stats_path = os.path.join(chrom_out_dir, "component_statistics.txt")
    
    # Start fresh each run to avoid duplicated rows when re-running.
    if os.path.exists(component_stats_path):
        os.remove(component_stats_path)

    bases = find_chunk_bases(chunk_dir)
    if not bases:
        print(f"[graph_ops] No chunk base files found in {chunk_dir}")
        return

    print(f"[graph_ops] Found {len(bases)} chunk(s) in {chunk_dir}")

    for base in bases:
        basename = os.path.basename(base)
        print(f"[graph_ops] Processing chunk base: {basename}")

        edges = load_edges_from_base(base)

        builder = GraphBuilder()
        for e in edges:
            builder.try_add_edge(e)

        # Skip chunks with no accepted edges; still could emit empty graphs if desired.
        json_path = os.path.join(outdir, f"{basename}.json")
        condensed_json_path = os.path.join(outdir, f"{basename}_condensed.json")
        log_path = os.path.join(outdir, f"{basename}_inconsistencies.tsv")
        accepted_path = os.path.join(outdir, f"{basename}_accepted_edges.txt")

        export_cytoscape_json(builder, json_path, name=basename)
        export_condensed_cytoscape_json(
            builder,
            condensed_json_path,
            name=basename + "_condensed",
        )
        write_inconsistency_log(builder, log_path)
        write_accepted_edges(builder, accepted_path)

        rows = compute_component_statistics_rows(builder, chunk_base=basename)
        append_component_statistics_tsv(rows, component_stats_path)

        print(
            f"[graph_ops]   -> wrote graph JSON: {json_path}\n"
            f"condensed JSON: {condensed_json_path}\n"
            f"inconsistencies: {log_path}\n"
            f"accepted edges: {accepted_path}\n"
            f"(nodes: {len(builder.nodes)}, edges: {len(builder.edges)})"
        )


if __name__ == "__main__":
    main()
