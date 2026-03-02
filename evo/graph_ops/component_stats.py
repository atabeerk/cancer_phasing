# graph_ops/component_stats.py

from __future__ import annotations

from typing import Dict, List, Any, Set, Optional
import csv

from builder import GraphBuilder


def _edge_support(e) -> int:
    """Total read support as used elsewhere in exports."""
    return int(e.alt_alt) + int(e.alt_ref) + int(e.ref_alt) + int(e.ref_ref)


def compute_component_statistics_rows(
    builder: GraphBuilder,
    chunk_base: str,
    chromosome: Optional[str] = None,
) -> List[Dict[str, Any]]:
    """
    Compute one row per connected component on the UNCONDENSED SNV graph.

    Connected components are computed on the undirected graph induced by ALL accepted edges.
    Timing direction is ignored for connectivity (we use builder.adj which is undirected).

    Haplotype statistics reuse the existing cooccurrence union-find inside GraphBuilder:
      - haplotypes = number of distinct cooccurrence cluster reps inside the CC
      - multi_node_haplotypes = number of those reps whose cooccurrence cluster size > 1
    """
    nodes = sorted(builder.nodes)
    if not nodes:
        return []

    # Infer chromosome if not given
    if chromosome is None:
        if builder.edges:
            chromosome = builder.edges[0].chrom
        else:
            chromosome = "NA"

    # --- Connected components on full (undirected) graph ---
    visited: Set[int] = set()
    components: List[List[int]] = []

    for start in nodes:
        if start in visited:
            continue
        stack = [start]
        visited.add(start)
        comp_nodes: List[int] = []

        while stack:
            u = stack.pop()
            comp_nodes.append(u)
            # builder.adj is undirected adjacency over all accepted edges
            for v in builder.adj.get(u, {}).keys():
                if v not in visited:
                    visited.add(v)
                    stack.append(v)

        components.append(comp_nodes)

    node_to_comp: Dict[int, int] = {}
    for cid, comp_nodes in enumerate(components):
        for n in comp_nodes:
            node_to_comp[n] = cid

    # --- Precompute cooccurrence cluster rep + sizes (reuse builder's DSU) ---
    node_to_cluster: Dict[int, int] = {n: builder.get_cluster_rep(n) for n in nodes}
    cluster_size: Dict[int, int] = builder.get_cluster_sizes()
    
    # --- Accumulate edge-based stats per component in one pass ---
    num_comps = len(components)
    timing_edges = [0] * num_comps
    coocc_edges = [0] * num_comps
    div_edges = [0] * num_comps
    edge_count = [0] * num_comps
    support_sum = [0] * num_comps

    for e in builder.edges:
        cu = node_to_comp.get(e.u)
        cv = node_to_comp.get(e.v)
        if cu is None or cv is None:
            continue
        # since this is an edge, cu and cv should be the same CC; guard anyway
        if cu != cv:
            continue

        cid = cu
        edge_count[cid] += 1
        support_sum[cid] += _edge_support(e)

        if e.relation == "timing":
            timing_edges[cid] += 1
        elif e.relation == "cooccurring":
            coocc_edges[cid] += 1
        elif e.relation == "divergent":
            div_edges[cid] += 1

    # --- Build final rows ---
    rows: List[Dict[str, Any]] = []
    for cid, comp_nodes in enumerate(components):
        comp_nodes_sorted = sorted(comp_nodes)
        min_pos = comp_nodes_sorted[0]
        max_pos = comp_nodes_sorted[-1]
        span_bp = max_pos - min_pos

        reps_in_cc = {node_to_cluster[n] for n in comp_nodes_sorted}
        haplotypes = len(reps_in_cc)
        multi_node_haplotypes = sum(1 for r in reps_in_cc if cluster_size.get(r, 1) > 1)

        avg_support = (
            (support_sum[cid] / edge_count[cid]) if edge_count[cid] > 0 else float("nan")
        )

        rows.append(
            {
                "chromosome": chromosome,
                "chunk_base": chunk_base,
                "component_id": cid,

                "min_pos": min_pos,
                "max_pos": max_pos,
                "span_bp": span_bp,
                "num_nodes": len(comp_nodes_sorted),

                "timing_edges": timing_edges[cid],
                "cooccurrence_edges": coocc_edges[cid],
                "divergent_edges": div_edges[cid],
                "avg_read_support": avg_support,

                "haplotypes": haplotypes,
                "multi_node_haplotypes": multi_node_haplotypes,

                # Copy-number annotation fields are populated during postprocessing
                # when a CN BED file is provided.
                "cn_copy_number_state": "NA",
                "cn_segment_start": "NA",
                "cn_segment_end": "NA",
                "cn_overlap_bp": "NA",
                "cn_segment_coverage": "NA",
                "cn_segment_confidence": "NA",
                "cn_haplotype_source": "NA",

                "nodes": ",".join(str(x) for x in comp_nodes_sorted),
            }
        )

    return rows


def append_component_statistics_tsv(rows: List[Dict[str, Any]], out_path: str) -> None:
    """Append rows to a TSV, writing a header if the file does not yet exist."""
    if not rows:
        return

    fieldnames = list(rows[0].keys())
    file_exists = False
    try:
        with open(out_path, "r"):
            file_exists = True
    except FileNotFoundError:
        file_exists = False

    with open(out_path, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        if not file_exists:
            w.writeheader()
        for r in rows:
            w.writerow(r)
