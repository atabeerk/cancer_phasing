# graph_ops/export_graph.py

import json
import csv
from collections import Counter
from typing import Optional

from builder import GraphBuilder


def export_cytoscape_json(
    builder: GraphBuilder,
    out_path: str,
    name: Optional[str] = None,
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
            "loss": <bool>,
            "reliability": <float>,
            "directed": <bool>,         # True only for timing
            "read_counts": "x/y/z/t",   # ALT_ALT/ALT_REF/REF_ALT/REF_REF
            "label": "x/y/z/t | rel=R"  # convenience string for visualization
          }
        }
    """
    if name is None:
        name = "chunk_graph"

    # Nodes: positions + VAF
    nodes_json = []
    for pos in sorted(builder.nodes):
        data = {"id": str(pos), "position": pos}
        chrom = builder.node_chrom.get(pos)
        if chrom is not None:
            data["chrom"] = chrom
        vaf = builder.node_vaf.get(pos)
        if vaf is not None:
            data["vaf"] = vaf
        data["haplotype"] = builder.node_haplotype.get(pos, "UNKNOWN")
        nodes_json.append({"data": data})

    # Edges with directed flag and read counts / label
    edges_json = []
    for idx, e in enumerate(builder.edges):
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
                "loss": e.loss,
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
            "loss": <bool>,
            "reliability": <float>,        # max reliability of supporting SNV–SNV edges
            "directed": <bool>,            # True only for timing
            "read_counts": "x/y/z/t",      # summed ALT_ALT/ALT_REF/REF_ALT/REF_REF
            "label": "x/y/z/t | rel=R | n=N",
            "support_edge_count": <int>,   # how many SNV–SNV edges collapsed into this
          }
        }
    """
    if name is None:
        name = "chunk_graph_condensed"

    # ---- Build cluster mapping ----
    # builder.cluster_nodes: rep -> set(nodes)  (union–find components)
    # builder.nodes: SNVs that actually appear in accepted edges. :contentReference[oaicite:4]{index=4} :contentReference[oaicite:5]{index=5}
    clusters_used = {}
    for rep, members in builder.cluster_nodes.items():
        used = members & builder.nodes
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

        # Aggregate VAFs over members (if available)
        vafs = [
            builder.node_vaf.get(n)
            for n in members
            if n in builder.node_vaf and builder.node_vaf.get(n) is not None
        ]
        mean_vaf = sum(vafs) / len(vafs) if vafs else None

        data = {
            "id": str(rep),
            "cluster_id": rep,
            "members": members,
            "cluster_size": len(members),
        }

        # All members should be on the same chromosome; store it for convenience.
        chrom = builder.node_chrom.get(members[0]) if members else None
        if chrom is not None:
            data["chrom"] = chrom

        # Use the smallest position as a representative "position" for layout
        data["position"] = members[0]
        if mean_vaf is not None:
            data["vaf"] = mean_vaf

        # Condensed node haplotype summary:
        # choose majority haplotype; explicitly flag if mixed within cluster.
        member_haps = [builder.node_haplotype.get(n, "UNKNOWN") for n in members]
        hap_counts = Counter(member_haps)
        top_hap, top_count = sorted(hap_counts.items(), key=lambda kv: (-kv[1], kv[0]))[0]
        informative_haps = {h for h in hap_counts if h not in {"UNKNOWN", "MIXED"}}
        data["haplotype"] = top_hap
        data["mixed_haplotypes"] = ("MIXED" in hap_counts) or (len(informative_haps) > 1)
        data["haplotype_counts"] = dict(sorted(hap_counts.items()))
        data["haplotype_majority_fraction"] = (top_count / len(members)) if members else 0.0

        nodes_json.append({"data": data})

    # ---- Aggregate edges between clusters ----
    # We skip cooccurring edges entirely (they are internal to clusters now).
    # For timing edges, direction is preserved (cluster_u -> cluster_v).
    # For divergent edges, treat as undirected (canonical ordering).
    agg_edges = {}

    for e in builder.edges:  # accepted edges only :contentReference[oaicite:6]{index=6}
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
            key = (cu, cv, "timing", e.loss)
            directed = True
            rep_u, rep_v = cu, cv
        else:  # "divergent"
            # Undirected; canonical ordering of cluster reps
            rep_u, rep_v = (cu, cv) if cu <= cv else (cv, cu)
            key = (rep_u, rep_v, "divergent", e.loss)
            directed = False

        if key not in agg_edges:
            agg_edges[key] = {
                "u": rep_u,
                "v": rep_v,
                "relation": e.relation,
                "loss": e.loss,
                "directed": directed,
                "alt_alt": e.alt_alt,
                "alt_ref": e.alt_ref,
                "ref_alt": e.ref_alt,
                "ref_ref": e.ref_ref,
                "reliability": e.reliability,   # start with this edge's reliability
                "support_edge_count": 1,
            }
        else:
            rec = agg_edges[key]
            rec["alt_alt"] += e.alt_alt
            rec["alt_ref"] += e.alt_ref
            rec["ref_alt"] += e.ref_alt
            rec["ref_ref"] += e.ref_ref
            rec["support_edge_count"] += 1
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
                "loss": rec["loss"],
                "reliability": rec["reliability"],
                "directed": rec["directed"],
                "read_counts": read_counts_str,
                "label": label_str,
                "support_edge_count": rec["support_edge_count"],
            }
        })

    graph = {
        "data": {"name": name},
        "elements": {"nodes": nodes_json, "edges": edges_json},
    }

    with open(out_path, "w") as f:
        json.dump(graph, f, indent=2)


def write_inconsistency_log(builder: GraphBuilder, out_path: str) -> None:
    """
    Write a TSV log of rejected edges and reasons.

    Columns:
      u, v                  : endpoints of the skipped edge
      relation, loss        : relation type + loss flag of the skipped edge
      reliability           : reliability of the skipped edge
      read_counts           : "ALT_ALT/ALT_REF/REF_ALT/REF_REF" for the skipped edge
      reason                : reason for skipping
      conflict_u, conflict_v: endpoints of the conflicting existing edge (if any)
      conflict_relation_label: relation+loss for the conflicting edge (e.g. "timing_loss")
      conflict_reliability  : reliability of the conflicting edge

    Notable reason values include:
      - pair_type_conflict
      - pair_direction_conflict
      - cluster_merge_noncooccurring
      - cluster_internal_timing
      - cluster_internal_divergent
      - timing_cycle
      - cluster_pair_relation_conflict
      - cluster_pair_timing_direction_conflict
    """
    fieldnames = [
        "u",
        "v",
        "relation",
        "loss",
        "reliability",
        "read_counts",
        "reason",
        "conflict_u",
        "conflict_v",
        "conflict_relation_label",
        "conflict_reliability",
    ]

    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for rec in builder.inconsistencies:
            writer.writerow(rec)


def write_accepted_edges(builder: GraphBuilder, out_path: str) -> None:
    """
    Write all accepted edges (those actually added to the graph)
    in a single file, with fields similar to the C++ text output plus:

      RELATION=<cooccurring|timing|divergent>
      LOSS=<true|false>

    Format (space-separated):

      chr pos1 pos2 RELATION=... LOSS=... VAF1=... VAF2=...
      ALT_ALT=... ALT_REF=... REF_ALT=... REF_REF=... TOTAL=...
      RELIABILITY=... BEST_SCORE=... MARGIN=...

    where pos1=u, pos2=v (the oriented endpoints used in the graph).
    """
    with open(out_path, "w") as f:
        for e in builder.edges:
            total = e.alt_alt + e.alt_ref + e.ref_alt + e.ref_ref
            loss_str = "true" if e.loss else "false"

            line = (
                f"{e.chrom} {e.u} {e.v} "
                f"RELATION={e.relation} "
                f"LOSS={loss_str} "
                f"VAF1={e.vaf_u} "
                f"VAF2={e.vaf_v} "
                f"ALT_ALT={e.alt_alt} "
                f"ALT_REF={e.alt_ref} "
                f"REF_ALT={e.ref_alt} "
                f"REF_REF={e.ref_ref} "
                f"TOTAL={total} "
                f"RELIABILITY={e.reliability} "
                f"BEST_SCORE={e.best_score} "
                f"MARGIN={e.margin}"
            )
            f.write(line + "\n")