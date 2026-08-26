import os
import sys
import argparse
import csv
import time
from bisect import bisect_right

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
if CURRENT_DIR not in sys.path:
    sys.path.insert(0, CURRENT_DIR)

from parsing import find_chunk_bases, group_chunk_bases_by_chromosome, load_edges_from_base
from builder import GraphBuilder
from fast_builder import FastGraphBuilder
from export_graph import export_cytoscape_json, export_condensed_cytoscape_json, write_inconsistency_log, write_accepted_edges
from component_stats import compute_component_statistics_rows, append_component_statistics_tsv


RELATIONS = ("cooccurring", "timing", "divergent")


def _partition_component_rows(rows, chunk_records):
    """
    Start from the original consecutive batch ranges. If a connected component
    crosses a batch end, extend that output region to the component's maximum
    coordinate and start the next region at the following coordinate.
    """
    batch_ranges = sorted((int(start), int(end)) for _base, start, end in chunk_records)
    if not batch_ranges:
        return []

    batch_ends = [end for _start, end in batch_ranges]
    ordered_rows = sorted(
        rows,
        key=lambda row: (
            int(row["min_pos"]),
            int(row["max_pos"]),
            int(row["component_id"]),
        ),
    )

    partitions = []
    row_index = 0
    batch_index = 0
    region_start = batch_ranges[0][0]

    while batch_index < len(batch_ranges):
        region_end = batch_ranges[batch_index][1]
        partition_rows = []

        # region_end may grow while rows are consumed. That closure ensures a
        # whole component is kept and no later component begins inside a region
        # assigned to an earlier file.
        while (
            row_index < len(ordered_rows)
            and int(ordered_rows[row_index]["min_pos"]) <= region_end
        ):
            row = ordered_rows[row_index]
            if int(row["min_pos"]) < region_start:
                raise AssertionError("Connected component assigned to two regions")
            partition_rows.append(row)
            region_end = max(region_end, int(row["max_pos"]))
            row_index += 1

        partitions.append((region_start, region_end, partition_rows))
        region_start = region_end + 1
        batch_index = bisect_right(batch_ends, region_end)

    if row_index != len(ordered_rows):
        raise ValueError("Graph components extend beyond the final batch range")

    return partitions


def _component_partition_nodes(rows):
    nodes = set()
    for row in rows:
        nodes.update(int(token) for token in str(row["nodes"]).split(",") if token)
    return nodes


def _new_edge_haplotype_counters():
    return {
        relation: {
            "total": 0,
            "same_haplotype": 0,
            "different_haplotype": 0,
            "at_least_one_unknown": 0,
            "hp1_hp1": 0,
            "hp2_hp2": 0,
        }
        for relation in RELATIONS
    }


def _normalize_hp(hp):
    if hp is None:
        return "UNKNOWN"
    s = str(hp).strip().upper()
    if s in {"HP1", "1", "H1"}:
        return "HP1"
    if s in {"HP2", "2", "H2"}:
        return "HP2"
    # Per requirement, treat MIXED as UNKNOWN.
    return "UNKNOWN"


def _update_edge_haplotype_counters(counters, edge):
    rel = edge.relation
    if rel not in counters:
        return
    rec = counters[rel]
    rec["total"] += 1
    hu = _normalize_hp(edge.hap_u)
    hv = _normalize_hp(edge.hap_v)
    if hu == "UNKNOWN" or hv == "UNKNOWN":
        rec["at_least_one_unknown"] += 1
        return
    if hu == hv:
        rec["same_haplotype"] += 1
        if hu == "HP1":
            rec["hp1_hp1"] += 1
        elif hu == "HP2":
            rec["hp2_hp2"] += 1
    else:
        rec["different_haplotype"] += 1


def _load_snv_assignment_counts(chrom_out_dir):
    """
    Load per-SNV assignment counts from haplotagged_snvs.tsv.
    Returns counts over all filtered SNVs in this chromosome run.
    """
    counts = {"HP1": 0, "HP2": 0, "MIXED": 0, "UNKNOWN": 0}
    path = os.path.join(chrom_out_dir, "haplotagged_snvs.tsv")
    if not os.path.exists(path):
        return counts

    try:
        with open(path, "r", newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                label = str(row.get("HP_ASSIGNMENT", "")).strip().upper()
                if label in counts:
                    counts[label] += 1
                else:
                    counts["UNKNOWN"] += 1
    except Exception as exc:
        print(f"[graph_ops] Warning: failed to parse {path}: {exc}")
    return counts


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Build chromosome-level consistent SNV graphs from chunk files.\n"
            "Each chunk's per-relation files ( *_cooccurring.txt etc. ) must "
            "already exist in the same directory. All chromosome edges are globally sorted."
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
    parser.add_argument(
        "--builder",
        choices=("fast", "legacy"),
        default="fast",
        help="Graph consistency builder implementation to use.",
    )
    parser.add_argument(
        "--progress-interval",
        type=int,
        default=50000,
        help="Print graph-build progress every N attempted edges; use 0 to disable.",
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

    chromosome_groups = group_chunk_bases_by_chromosome(bases)

    builder_cls = FastGraphBuilder if args.builder == "fast" else GraphBuilder
    print(
        f"[graph_ops] Found {len(bases)} chunk(s) across "
        f"{len(chromosome_groups)} chromosome(s) in {chunk_dir}", flush=True
    )
    print(f"[graph_ops] Builder: {args.builder}", flush=True)

    all_relation_positions = set()
    all_graph_nodes = set()
    all_timing_nodes = set()
    total_accepted_edges = 0
    total_timing_edges = 0
    edge_haplotype_counters = _new_edge_haplotype_counters()

    for chrom, chunk_records in chromosome_groups.items():
        chrom_t0 = time.monotonic()
        core_start = min(record[1] for record in chunk_records)
        core_end = max(record[2] for record in chunk_records)
        basename = f"chunk_{chrom}_{core_start}_{core_end}"
        print(
            f"[graph_ops] Processing chromosome: {chrom} "
            f"chunks={len(chunk_records)} output_base={basename}", flush=True
        )

        edges = []
        for base, _start, _end in chunk_records:
            edges.extend(load_edges_from_base(base))
        edges.sort(key=lambda e: e.reliability, reverse=True)
        print(
            f"[graph_ops]   loaded_globally_sorted_edges={len(edges)} "
            f"elapsed={time.monotonic() - chrom_t0:.3f}s",
            flush=True,
        )

        for e in edges:
            if e.chrom != chrom:
                raise ValueError(
                    f"Edge chromosome {e.chrom!r} does not match chunk chromosome {chrom!r}: "
                    f"{e.source_file}"
                )
            all_relation_positions.add((e.chrom, e.u))
            all_relation_positions.add((e.chrom, e.v))

        builder = builder_cls()
        build_t0 = time.monotonic()
        accepted = 0
        for idx, e in enumerate(edges, start=1):
            if builder.try_add_edge(e):
                accepted += 1
            if args.progress_interval > 0 and idx % args.progress_interval == 0:
                print(
                    f"[graph_ops]   build_progress chromosome={chrom} "
                    f"attempted={idx}/{len(edges)} accepted={accepted} "
                    f"rejected={idx - accepted} elapsed={time.monotonic() - build_t0:.3f}s",
                    flush=True,
                )
        print(
            f"[graph_ops]   build_done chromosome={chrom} attempted={len(edges)} "
            f"accepted={accepted} rejected={len(edges) - accepted} "
            f"elapsed={time.monotonic() - build_t0:.3f}s",
            flush=True,
        )

        all_graph_nodes.update((chrom, node) for node in builder.nodes)
        total_accepted_edges += len(builder.edges)
        for e in builder.edges:
            _update_edge_haplotype_counters(edge_haplotype_counters, e)
            if e.relation == "timing":
                all_timing_nodes.add((e.chrom, e.u))
                all_timing_nodes.add((e.chrom, e.v))
                total_timing_edges += 1

        log_path = os.path.join(outdir, f"{basename}_inconsistencies.tsv")
        accepted_path = os.path.join(outdir, f"{basename}_accepted_edges.txt")

        export_t0 = time.monotonic()
        write_inconsistency_log(builder, log_path)
        write_accepted_edges(builder, accepted_path)

        component_rows = compute_component_statistics_rows(
            builder,
            chunk_base=basename,
            chromosome=chrom,
        )
        partitions = _partition_component_rows(component_rows, chunk_records)
        partition_specs = []
        for partition_start, partition_end, partition_rows in partitions:
            partition_nodes = _component_partition_nodes(partition_rows)
            partition_base = f"chunk_{chrom}_{partition_start}_{partition_end}"
            partition_specs.append(
                (partition_base, partition_nodes, partition_rows)
            )

        partition_paths = []

        for partition_base, partition_nodes, partition_rows in partition_specs:
            json_path = os.path.join(outdir, f"{partition_base}.json")
            condensed_json_path = os.path.join(
                outdir,
                f"{partition_base}_condensed.json",
            )

            export_cytoscape_json(
                builder,
                json_path,
                name=partition_base,
                node_subset=partition_nodes,
            )
            export_condensed_cytoscape_json(
                builder,
                condensed_json_path,
                name=partition_base + "_condensed",
                node_subset=partition_nodes,
            )

            stats_rows = []
            for component_id, row in enumerate(partition_rows):
                output_row = dict(row)
                output_row["chunk_base"] = partition_base
                output_row["component_id"] = component_id
                stats_rows.append(output_row)
            append_component_statistics_tsv(stats_rows, component_stats_path)
            partition_paths.append((json_path, condensed_json_path))

        print(
            f"[graph_ops]   export_done chromosome={chrom} "
            f"json_partitions={len(partition_paths)} "
            f"elapsed={time.monotonic() - export_t0:.3f}s",
            flush=True,
        )

        print(
            f"[graph_ops]   -> wrote {len(partition_paths)} paired graph JSON partition(s)\n"
            f"inconsistencies: {log_path}\n"
            f"accepted edges: {accepted_path}\n"
            f"(nodes: {len(builder.nodes)}, edges: {len(builder.edges)})"
            f"\nchromosome_elapsed={time.monotonic() - chrom_t0:.3f}s",
            flush=True,
        )

    # ---- Mutation connectivity statistics ----
    total_snvs = args.total_snvs
    n_in_relations = len(all_relation_positions)
    n_in_graph = len(all_graph_nodes)
    n_with_timing = len(all_timing_nodes)
    n_orphaned = n_in_relations - n_in_graph
    n_singleton = total_snvs - n_in_relations if total_snvs > 0 else 0
    snv_hp_counts = _load_snv_assignment_counts(chrom_out_dir)

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
        f"    In graph (has at least one accepted edge): {n_in_graph}"
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
            f"  Singleton (has no pairwise relations): {n_singleton}"
            f"  ({_pct(n_singleton, total_snvs)} of total)"
        )
    lines.append(
        f"  HP1 assigned:                         {snv_hp_counts['HP1']}"
        + (f"  ({_pct(snv_hp_counts['HP1'], total_snvs)} of total)" if total_snvs > 0 else "")
    )
    lines.append(
        f"  HP2 assigned:                         {snv_hp_counts['HP2']}"
        + (f"  ({_pct(snv_hp_counts['HP2'], total_snvs)} of total)" if total_snvs > 0 else "")
    )
    lines.append(
        f"  MIXED assigned:                       {snv_hp_counts['MIXED']}"
        + (f"  ({_pct(snv_hp_counts['MIXED'], total_snvs)} of total)" if total_snvs > 0 else "")
    )
    lines.append(
        f"  UNKNOWN assigned:                     {snv_hp_counts['UNKNOWN']}"
        + (f"  ({_pct(snv_hp_counts['UNKNOWN'], total_snvs)} of total)" if total_snvs > 0 else "")
    )
    lines.append(f"Total accepted edges:                  {total_accepted_edges}")
    lines.append(f"  - timing edges:                      {total_timing_edges}")

    summary_text = "\n".join(lines)
    print(summary_text)
    # Machine-readable lines for chromosome logs (parsed by run_genome_by_chr.py).
    print(
        "[graph_ops_mutation_stats] "
        f"total_snvs={total_snvs} "
        f"in_relations={n_in_relations} "
        f"in_graph={n_in_graph} "
        f"with_timing={n_with_timing} "
        f"orphaned={n_orphaned} "
        f"singleton={n_singleton} "
        f"accepted_edges={total_accepted_edges} "
        f"timing_edges={total_timing_edges} "
        f"hp1={snv_hp_counts['HP1']} "
        f"hp2={snv_hp_counts['HP2']} "
        f"mixed={snv_hp_counts['MIXED']} "
        f"unknown={snv_hp_counts['UNKNOWN']}"
    )
    for relation in RELATIONS:
        rec = edge_haplotype_counters[relation]
        print(
            "[graph_ops_edge_hap] "
            f"relation={relation} "
            f"total={rec['total']} "
            f"same_haplotype={rec['same_haplotype']} "
            f"different_haplotype={rec['different_haplotype']} "
            f"at_least_one_unknown={rec['at_least_one_unknown']} "
            f"hp1_hp1={rec['hp1_hp1']} "
            f"hp2_hp2={rec['hp2_hp2']}"
        )


if __name__ == "__main__":
    main()
