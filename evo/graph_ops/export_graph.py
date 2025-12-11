# graph_ops/export_graph.py

import json
import csv
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
        { "data": { "id": "<pos>", "position": <pos> } }

    Edge schema:
        {
          "data": {
            "id": "<u>-<v>-<index>",
            "source": "<u>",
            "target": "<v>",
            "relation": "timing" | "cooccurring" | "divergent",
            "loss": <bool>,
            "reliability": <float>,
            "directed": <bool>   # True only for timing relations
          }
        }
    """
    if name is None:
        name = "chunk_graph"

    # Nodes: positions only
    nodes_json = [
        {"data": {"id": str(pos), "position": pos}}
        for pos in sorted(builder.nodes)
    ]

    # Edges with directed flag
    edges_json = []
    for idx, e in enumerate(builder.edges):
        edge_id = f"{e.u}-{e.v}-{idx}"
        edges_json.append({
            "data": {
                "id": edge_id,
                "source": str(e.u),
                "target": str(e.v),
                "relation": e.relation,         # "cooccurring" / "timing" / "divergent"
                "loss": e.loss,
                "reliability": e.reliability,
                "directed": (e.relation == "timing"),
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
    For timing cycles, 'details' states that there exists a path from target to source.
    """
    fieldnames = [
        "u",
        "v",
        "relation",
        "loss",
        "reliability",
        "reason",
        "conflict_u",
        "conflict_v",
        "conflict_relation",
        "conflict_loss",
        "conflict_reliability",
        "details",
    ]
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for rec in builder.inconsistencies:
            writer.writerow(rec)
