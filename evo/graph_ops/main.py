import os
import sys
import argparse

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
if CURRENT_DIR not in sys.path:
    sys.path.insert(0, CURRENT_DIR)

from parsing import find_chunk_bases, load_edges_from_base
from builder import GraphBuilder
from export_graph import export_cytoscape_json, write_inconsistency_log


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
        log_path = os.path.join(outdir, f"{basename}_inconsistencies.tsv")

        export_cytoscape_json(builder, json_path, name=basename)
        write_inconsistency_log(builder, log_path)

        print(
            f"[graph_ops]   -> wrote graph JSON: {json_path}, "
            f"inconsistencies: {log_path} "
            f"(nodes: {len(builder.nodes)}, edges: {len(builder.edges)})"
        )


if __name__ == "__main__":
    main()
