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
    parser.add_argument(
        "--total-snvs",
        type=int,
        default=0,
        help="Total number of SNVs in the filtered VCF (for connectivity reporting).",
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

    all_relation_positions = set()
    all_graph_nodes = set()
    all_timing_nodes = set()
    total_accepted_edges = 0
    total_timing_edges = 0

    for base in bases:
        basename = os.path.basename(base)
        print(f"[graph_ops] Processing chunk base: {basename}")

        edges = load_edges_from_base(base)

        for e in edges:
            all_relation_positions.add(e.u)
            all_relation_positions.add(e.v)

        builder = GraphBuilder()
        for e in edges:
            builder.try_add_edge(e)

        all_graph_nodes.update(builder.nodes)
        total_accepted_edges += len(builder.edges)
        for e in builder.edges:
            if e.relation == "timing":
                all_timing_nodes.add(e.u)
                all_timing_nodes.add(e.v)
                total_timing_edges += 1

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

    # ---- Mutation connectivity statistics ----
    total_snvs = args.total_snvs
    n_in_relations = len(all_relation_positions)
    n_in_graph = len(all_graph_nodes)
    n_with_timing = len(all_timing_nodes)
    n_orphaned = n_in_relations - n_in_graph
    n_singleton = total_snvs - n_in_relations if total_snvs > 0 else 0

    def _pct(num, denom):
        return f"{100.0 * num / denom:.1f}%" if denom > 0 else "N/A"

    lines = ["\n=== Mutation connectivity statistics ==="]
    if total_snvs > 0:
        lines.append(f"Total mutations (filtered VCF):        {total_snvs}")
    lines.append(
        f"  In pairwise relations:               {n_in_relations}"
        + (f"  ({_pct(n_in_relations, total_snvs)} of total)" if total_snvs > 0 else "")
    )
    lines.append(
        f"    In graph (accepted edges):          {n_in_graph}"
        + (f"  ({_pct(n_in_graph, n_in_relations)} of connected)" if n_in_relations > 0 else "")
    )
    lines.append(
        f"      - with timing edge:              {n_with_timing}"
        + (f"  ({_pct(n_with_timing, n_in_graph)} of in-graph)" if n_in_graph > 0 else "")
    )
    lines.append(
        f"    Orphaned (all edges rejected):      {n_orphaned}"
        + (f"  ({_pct(n_orphaned, n_in_relations)} of connected)" if n_in_relations > 0 else "")
    )
    if total_snvs > 0:
        lines.append(
            f"  Singleton (no pairwise relations):    {n_singleton}"
            f"  ({_pct(n_singleton, total_snvs)} of total)"
        )
    lines.append(f"Total accepted edges:                  {total_accepted_edges}")
    lines.append(f"  - timing edges:                      {total_timing_edges}")

    summary_text = "\n".join(lines)
    print(summary_text)

    stats_path = os.path.join(chrom_out_dir, "mutation_stats.txt")
    with open(stats_path, "w") as f:
        f.write(summary_text + "\n")
    print(f"[graph_ops] Wrote mutation stats: {stats_path}")

    tsv_path = os.path.join(chrom_out_dir, "mutation_stats.tsv")
    with open(tsv_path, "w") as f:
        f.write("metric\tcount\n")
        f.write(f"total_snvs\t{total_snvs}\n")
        f.write(f"in_relations\t{n_in_relations}\n")
        f.write(f"in_graph\t{n_in_graph}\n")
        f.write(f"with_timing\t{n_with_timing}\n")
        f.write(f"orphaned\t{n_orphaned}\n")
        f.write(f"singleton\t{n_singleton}\n")
        f.write(f"accepted_edges\t{total_accepted_edges}\n")
        f.write(f"timing_edges\t{total_timing_edges}\n")


if __name__ == "__main__":
    main()
