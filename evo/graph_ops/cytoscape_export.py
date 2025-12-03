from typing import Dict
from pathlib import Path
import json


def write_cytoscape_json(nodes, directed_edges, directed_loss_edges,
                         cooccurring_edges, cooccurring_loss_edges, 
                         divergent_edges, component_stats, node_to_component,
                         edge_meta, out_file: Path, pileup_depths: Dict[str, int]):
    """Write a Cytoscape-format JSON file containing nodes, edges, and component stats."""
    elements = {"nodes": [], "edges": []}

    for n in sorted(nodes, key=lambda x: (not x.isdigit(), x)):
        depth = 0
        if pileup_depths:
            depth = pileup_depths.get(n, 0)

        elements["nodes"].append({
            "data": {
                "id": n,
                "label": n,
                "component_id": node_to_component.get(n, -1),
                "depth": depth
            }
        })

    def add_edges(edges, edge_type, directed_flag):
        for chrom, s, t in edges:
            meta = edge_meta.get((chrom, s, t), {})
            elements["edges"].append({
                "data": {
                    "source": s,
                    "target": t,
                    "chrom": chrom,
                    "directed": directed_flag,
                    "edge_type": edge_type,
                    "total": meta.get("total"),
                    "original": meta.get("original")
                }
            })

    add_edges(directed_edges, "directed", True)
    add_edges(directed_loss_edges, "directed_loss", True)
    add_edges(cooccurring_edges, "cooccurring", False)
    add_edges(cooccurring_loss_edges, "cooccurring_loss", False)
    add_edges(divergent_edges, "divergent", False)

    out_file = Path(out_file)
    with open(out_file, "w") as out:
        json.dump({"elements": elements, "component_stats": component_stats}, out, indent=2)
