# graph_ops/export_graph.py

import json
import csv
from collections import Counter
from typing import Dict, Iterable, List, Optional, Set, Tuple

from builder import GraphBuilder


def _median(values: List[float]) -> Optional[float]:
    """Return the median of a non-empty list; otherwise None."""
    if not values:
        return None
    vals = sorted(values)
    n = len(vals)
    mid = n // 2
    if n % 2 == 1:
        return float(vals[mid])
    return float((vals[mid - 1] + vals[mid]) / 2.0)


def _lower_median(values: List[float]) -> Optional[float]:
    """Return median using lower middle value for even-length inputs."""
    if not values:
        return None
    vals = sorted(values)
    return float(vals[(len(vals) - 1) // 2])


def _component_layout_positions(
    node_ids: List[int],
    edges: Iterable[Tuple[int, int]],
    node_x_sort: Dict[int, float],
    node_y_sort: Dict[int, float],
    node_x_gap: float = 80.0,
    node_y_gap: float = 70.0,
    component_gap: float = 220.0,
) -> Dict[int, Dict[str, float]]:
    """
    Compute Cytoscape positions with these rules:
      1) Genomic coordinate and VAF are used only for sorting.
      2) Actual x/y positions are synthetic, rank-based, and evenly spaced.
      3) Connected components are ordered left-to-right by component minimum
         genomic coordinate key and separated by a fixed gap.
    """
    adj: Dict[int, set] = {n: set() for n in node_ids}
    for u, v in edges:
        if u in adj and v in adj:
            adj[u].add(v)
            adj[v].add(u)

    # Build connected components.
    components: List[List[int]] = []
    seen = set()
    for n in sorted(node_ids):
        if n in seen:
            continue
        stack = [n]
        seen.add(n)
        comp = []
        while stack:
            cur = stack.pop()
            comp.append(cur)
            for nxt in adj.get(cur, ()):
                if nxt not in seen:
                    seen.add(nxt)
                    stack.append(nxt)
        components.append(comp)

    # Sort components by their minimum coordinate, then by min node id for stability.
    components.sort(
        key=lambda comp: (
            min(node_x_sort[m] for m in comp),
            min(comp),
        )
    )

    positions: Dict[int, Dict[str, float]] = {}
    comp_start_x = 0.0
    for comp in components:
        x_order = sorted(comp, key=lambda n: (node_x_sort[n], node_y_sort[n], n))
        # Bin VAF values into 0.05-width buckets so each bucket shares a y-level.
        # Example bins: [0.00,0.05), [0.05,0.10), ...
        y_bin = {n: int(node_y_sort[n] / 0.05) for n in comp}
        present_bins = sorted(set(y_bin.values()))
        bin_rank = {b: i for i, b in enumerate(present_bins)}
        y_mid = (len(present_bins) - 1) / 2.0 if present_bins else 0.0

        for i, n in enumerate(x_order):
            positions[n] = {
                "x": float(comp_start_x + i * node_x_gap),
                # Cytoscape screen coordinates increase downward in y, so negate
                # the bin offset to place higher VAF bins visually higher.
                "y": float((y_mid - bin_rank[y_bin[n]]) * node_y_gap),
            }

        comp_width = max(0.0, (len(comp) - 1) * node_x_gap)
        comp_start_x += comp_width + component_gap

    return positions


def export_cytoscape_json(
    builder: GraphBuilder,
    out_path: str,
    name: Optional[str] = None,
    node_subset: Optional[Set[int]] = None,
) -> None:
    """
    Export the graph to a Cytoscape-compatible JSON file.

    Node schema:
        {
          "data": {
            "id": "<pos>",
            "chrom": "<chr>",   # chromosome label for this SNV/node
            "position": <pos>,
            "vaf": <float>,     # VAF for this SNV/node
            "haplotype": "<HP1|HP2|UNKNOWN|MIXED>"
          }
        }

    Edge schema:
        {
          "data": {
            "id": "<u>-<v>-<index>",
            "source": "<u>",
            "target": "<v>",
            "relation": "timing" | "cooccurring" | "divergent",
            "reliability": <float>,
            "directed": <bool>,         # True only for timing
            "read_counts": "x/y/z/t",   # ALT_ALT/ALT_REF/REF_ALT/REF_REF
            "label": "x/y/z/t | rel=R"  # convenience string for visualization
          }
        }
    """
    if name is None:
        name = "chunk_graph"

    selected_nodes = (
        set(builder.nodes)
        if node_subset is None
        else set(builder.nodes) & set(node_subset)
    )
    selected_edges = [
        e for e in builder.edges
        if e.u in selected_nodes and e.v in selected_nodes
    ]
    sorted_nodes = sorted(selected_nodes)
    node_x_sort = {pos: float(pos) for pos in sorted_nodes}
    node_y_sort = {
        pos: float(builder.node_vaf[pos]) if builder.node_vaf.get(pos) is not None else 0.0
        for pos in sorted_nodes
    }
    uncondensed_positions = _component_layout_positions(
        node_ids=sorted_nodes,
        edges=((e.u, e.v) for e in selected_edges),
        node_x_sort=node_x_sort,
        node_y_sort=node_y_sort,
    )

    # Nodes: positions + VAF
    nodes_json = []
    for pos in sorted_nodes:
        # Keep position as string for visualization tools (e.g., Cytoscape)
        # to avoid scientific-notation formatting of large genomic coordinates.
        data = {"id": str(pos), "position": str(pos)}
        chrom = builder.node_chrom.get(pos)
        if chrom is not None:
            data["chrom"] = chrom
        vaf = builder.node_vaf.get(pos)
        if vaf is not None:
            data["vaf"] = vaf
        data["haplotype"] = builder.node_haplotype.get(pos, "UNKNOWN")
        data["hp_reads"] = builder.node_hp_reads.get(pos, "UNKNOWN")
        nodes_json.append({
            "data": data,
            "position": uncondensed_positions[pos],
        })

    # Edges with directed flag and read counts / label
    edges_json = []
    for idx, e in enumerate(selected_edges):
        edge_id = f"{e.u}-{e.v}-{idx}"

        # Build "x/y/z/t" read count string
        read_counts_str = f"{e.alt_alt}/{e.alt_ref}/{e.ref_alt}/{e.ref_ref}"

        # Build label string combining reads + reliability
        label_str = f"{read_counts_str} | rel={e.reliability:.3f}"

        edges_json.append({
            "data": {
                "id": edge_id,
                "source": str(e.u),
                "target": str(e.v),
                "relation": e.relation,           # "cooccurring" / "timing" / "divergent"
                "reliability": e.reliability,
                "directed": (e.relation == "timing"),
                "read_counts": read_counts_str,
                "label": label_str,
            }
        })

    graph = {
        "data": {"name": name},
        "elements": {"nodes": nodes_json, "edges": edges_json},
    }

    with open(out_path, "w") as f:
        json.dump(graph, f, indent=2)


def export_condensed_cytoscape_json(
    builder: GraphBuilder,
    out_path: str,
    name: Optional[str] = None,
    node_subset: Optional[Set[int]] = None,
) -> None:
    """
    Export a condensed graph where:
      - Each node is a cooccurrence cluster (union–find component).
      - Cooccurring edges are internal and removed.
      - Timing/divergent edges between SNVs induce edges between clusters.

    Node schema (per cluster):
        {
          "data": {
            "id": "<cluster_rep>",
            "chrom": "<chr>",                        # chromosome for the cluster
            "cluster_id": <cluster_rep>,            # numeric representative
            "position": <one member SNV position>,  # representative position
            "members": [pos1, pos2, ...],           # all SNVs in this cluster
            "cluster_size": <int>,
            "vaf": <float>,                         # mean VAF over members (if available)
            "haplotype": "<majority haplotype>",
            "mixed_haplotypes": <bool>,
            "haplotype_counts": {"HP1": n1, ...}
          }
        }

    Edge schema (between clusters):
        {
          "data": {
            "id": "<u>-<v>-<index>",
            "source": "<cluster_rep_u>",
            "target": "<cluster_rep_v>",
            "relation": "timing" | "divergent",
            "reliability": <float>,        # max reliability of supporting SNV–SNV edges
            "directed": <bool>,            # True only for timing
            "read_counts": "x/y/z/t",      # summed ALT_ALT/ALT_REF/REF_ALT/REF_REF
            "label": "x/y/z/t | rel=R | n=N",
            "support_edge_count": <int>,   # how many SNV–SNV edges collapsed into this
            "support_edge_details": [      # per supporting uncondensed edge
              {
                "source": "<uncondensed_u>",
                "target": "<uncondensed_v>",
                "read_counts": "x/y/z/t",
                "reliability": <float>
              },
              ...
            ],
          }
        }
    """
    if name is None:
        name = "chunk_graph_condensed"

    selected_nodes = (
        set(builder.nodes)
        if node_subset is None
        else set(builder.nodes) & set(node_subset)
    )
    selected_edges = [
        e for e in builder.edges
        if e.u in selected_nodes and e.v in selected_nodes
    ]

    # ---- Build cluster mapping ----
    # builder.cluster_nodes: rep -> set(nodes)  (union–find components)
    # builder.nodes: SNVs that actually appear in accepted edges. :contentReference[oaicite:4]{index=4} :contentReference[oaicite:5]{index=5}
    clusters_used = {}
    for rep, members in builder.cluster_nodes.items():
        used = members & selected_nodes
        if used:
            clusters_used[rep] = used

    # Map each SNV to its cluster representative
    node_to_cluster = {}
    for rep, members in clusters_used.items():
        for n in members:
            node_to_cluster[n] = rep

    # ---- Condensed nodes ----
    nodes_json = []
    for rep in sorted(clusters_used.keys()):
        members = sorted(clusters_used[rep])

        # Aggregate VAFs over members (if available); condensed values are medians.
        vafs = [
            builder.node_vaf.get(n)
            for n in members
            if n in builder.node_vaf and builder.node_vaf.get(n) is not None
        ]
        median_vaf = _median([float(v) for v in vafs]) if vafs else None
        median_position = _lower_median([float(m) for m in members]) if members else None

        data = {
            "id": str(rep),
            # Store large integer-like fields as strings to prevent scientific
            # notation in downstream graph viewers.
            "members": [str(m) for m in members],
            "cluster_size": len(members),
        }

        # All members should be on the same chromosome; store it for convenience.
        chrom = builder.node_chrom.get(members[0]) if members else None
        if chrom is not None:
            data["chrom"] = chrom

        # Store median coordinate in the canonical "position" field.
        # Keep as string to avoid scientific notation in downstream viewers.
        if median_position is not None:
            data["position"] = str(int(median_position))
        else:
            data["position"] = str(members[0])
        data["position_min"] = str(members[0])
        if median_vaf is not None:
            data["vaf"] = median_vaf

        # Condensed node haplotype summary:
        # choose majority haplotype; explicitly flag if mixed within cluster.
        member_haps = [builder.node_haplotype.get(n, "UNKNOWN") for n in members]
        hap_counts = Counter(member_haps)
        top_hap, top_count = sorted(hap_counts.items(), key=lambda kv: (-kv[1], kv[0]))[0]
        informative_haps = {h for h in hap_counts if h not in {"UNKNOWN", "MIXED"}}
        hp1_count = hap_counts.get("HP1", 0)
        hp2_count = hap_counts.get("HP2", 0)
        unk_count = len(members) - hp1_count - hp2_count
        data["haplotype"] = top_hap
        data["mixed_haplotypes"] = ("MIXED" in hap_counts) or (len(informative_haps) > 1)
        data["hp_counts(hp1/hp2/unk)"] = f"{hp1_count}/{hp2_count}/{unk_count}"

        nodes_json.append({"data": data})

    # ---- Aggregate edges between clusters ----
    # We skip cooccurring edges entirely (they are internal to clusters now).
    # For timing edges, direction is preserved (cluster_u -> cluster_v).
    # For divergent edges, treat as undirected (canonical ordering).
    agg_edges = {}

    for e in selected_edges:  # accepted edges only :contentReference[oaicite:6]{index=6}
        if e.relation == "cooccurring":
            continue  # internal to clusters in the condensed view

        cu = node_to_cluster.get(e.u)
        cv = node_to_cluster.get(e.v)
        if cu is None or cv is None:
            # Edge involves a node that didn't make it into the condensed graph;
            # should be rare, but we guard anyway.
            continue
        if cu == cv:
            # Defensive: timing/divergent edges inside a cluster should already
            # be forbidden by GraphBuilder invariants. :contentReference[oaicite:7]{index=7}
            continue

        if e.relation == "timing":
            key = (cu, cv, "timing")
            directed = True
            rep_u, rep_v = cu, cv
        else:  # "divergent"
            # Undirected; canonical ordering of cluster reps
            rep_u, rep_v = (cu, cv) if cu <= cv else (cv, cu)
            key = (rep_u, rep_v, "divergent")
            directed = False

        read_counts_str = f"{e.alt_alt}/{e.alt_ref}/{e.ref_alt}/{e.ref_ref}"
        source_id = str(e.u)
        target_id = str(e.v)
        support_detail = {
            "source": source_id,
            "target": target_id,
            "read_counts": read_counts_str,
            "reliability": e.reliability,
        }

        if key not in agg_edges:
            agg_edges[key] = {
                "u": rep_u,
                "v": rep_v,
                "relation": e.relation,
                "directed": directed,
                "alt_alt": e.alt_alt,
                "alt_ref": e.alt_ref,
                "ref_alt": e.ref_alt,
                "ref_ref": e.ref_ref,
                "reliability": e.reliability,   # start with this edge's reliability
                "support_edge_count": 1,
                "support_edge_details": [support_detail],
            }
        else:
            rec = agg_edges[key]
            rec["alt_alt"] += e.alt_alt
            rec["alt_ref"] += e.alt_ref
            rec["ref_alt"] += e.ref_alt
            rec["ref_ref"] += e.ref_ref
            rec["support_edge_count"] += 1
            rec["support_edge_details"].append(support_detail)
            if e.reliability > rec["reliability"]:
                rec["reliability"] = e.reliability

    # ---- Build Cytoscape edges JSON ----
    edges_json = []
    for idx, rec in enumerate(agg_edges.values()):
        edge_id = f"{rec['u']}-{rec['v']}-{idx}"
        read_counts_str = f"{rec['alt_alt']}/{rec['alt_ref']}/{rec['ref_alt']}/{rec['ref_ref']}"
        label_str = (
            f"{read_counts_str} | rel={rec['reliability']:.3f} "
            f"| n={rec['support_edge_count']}"
        )

        edges_json.append({
            "data": {
                "id": edge_id,
                "source": str(rec["u"]),
                "target": str(rec["v"]),
                "relation": rec["relation"],
                "reliability": rec["reliability"],
                "directed": rec["directed"],
                "read_counts": read_counts_str,
                "label": label_str,
                "support_edge_count": rec["support_edge_count"],
                "support_edge_details": rec["support_edge_details"],
            }
        })

    # Compute condensed-layout positions using median coordinate and median VAF.
    condensed_node_ids = sorted(clusters_used.keys())
    condensed_x_sort: Dict[int, float] = {}
    condensed_y_sort: Dict[int, float] = {}
    node_index = {
        int(node_obj["data"]["id"]): node_obj
        for node_obj in nodes_json
    }
    for rep in condensed_node_ids:
        d = node_index[rep]["data"]
        condensed_x_sort[rep] = float(d["position"])
        condensed_y_sort[rep] = float(d["vaf"]) if d.get("vaf") is not None else 0.0

    condensed_positions = _component_layout_positions(
        node_ids=condensed_node_ids,
        edges=((int(rec["u"]), int(rec["v"])) for rec in agg_edges.values()),
        node_x_sort=condensed_x_sort,
        node_y_sort=condensed_y_sort,
    )
    for node_obj in nodes_json:
        rep = int(node_obj["data"]["id"])
        node_obj["position"] = condensed_positions[rep]

    graph = {
        "data": {"name": name},
        "elements": {"nodes": nodes_json, "edges": edges_json},
    }

    with open(out_path, "w") as f:
        json.dump(graph, f, indent=2)


def _in_region(u: int, v: int, region_start: Optional[int], region_end: Optional[int]) -> bool:
    if region_start is None and region_end is None:
        return True
    if region_start is None or region_end is None:
        raise ValueError("region_start and region_end must be provided together")
    return region_start <= int(u) <= region_end and region_start <= int(v) <= region_end


def write_inconsistency_log(
    builder: GraphBuilder,
    out_path: str,
    region_start: Optional[int] = None,
    region_end: Optional[int] = None,
) -> int:
    """Write rejected edges, optionally restricted to one inclusive region."""
    fieldnames = [
        "u",
        "v",
        "relation",
        "reliability",
        "read_counts",
        "reason",
        "conflict_u",
        "conflict_v",
        "conflict_relation_label",
        "conflict_reliability",
    ]

    written = 0
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for rec in builder.inconsistencies:
            if not _in_region(rec["u"], rec["v"], region_start, region_end):
                continue
            writer.writerow(rec)
            written += 1
    return written


def write_accepted_edges(
    builder: GraphBuilder,
    out_path: str,
    region_start: Optional[int] = None,
    region_end: Optional[int] = None,
) -> int:
    """Write accepted edges, optionally restricted to one inclusive region."""
    written = 0
    with open(out_path, "w") as f:
        for e in builder.edges:
            if not _in_region(e.u, e.v, region_start, region_end):
                continue
            total = e.alt_alt + e.alt_ref + e.ref_alt + e.ref_ref
            line = (
                f"{e.chrom} {e.u} {e.v} "
                f"RELATION={e.relation} "
                f"VAF1={e.vaf_u} "
                f"VAF2={e.vaf_v} "
                f"ALT_ALT={e.alt_alt} "
                f"ALT_REF={e.alt_ref} "
                f"REF_ALT={e.ref_alt} "
                f"REF_REF={e.ref_ref} "
                f"TOTAL={total} "
                f"HP_READS1={(e.hp_reads_u if e.hp_reads_u is not None else 'UNKNOWN')} "
                f"HP_READS2={(e.hp_reads_v if e.hp_reads_v is not None else 'UNKNOWN')} "
                f"RELIABILITY={e.reliability} "
                f"BEST_SCORE={e.best_score} "
                f"MARGIN={e.margin}"
            )
            f.write(line + "\n")
            written += 1
    return written
