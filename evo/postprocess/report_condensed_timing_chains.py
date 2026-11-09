"""Export complete root-to-leaf timing paths from final condensed graphs."""

import csv
import json
from graphlib import TopologicalSorter
from pathlib import Path


def find_chains(successors, node_data):
    """Yield maximal timing paths of at least three nodes, without a length cap."""
    adjacency = {node: sorted(set(successors.get(node, ()))) for node in node_data}
    if any(target not in node_data for targets in adjacency.values() for target in targets):
        raise ValueError("Timing edge references a missing node")
    # Validate the entire graph, including components without a root.
    TopologicalSorter(adjacency).prepare()
    targets = {target for neighbors in adjacency.values() for target in neighbors}
    for root in sorted(set(node_data) - targets):
        stack = [(root,)]
        while stack:
            path = stack.pop()
            neighbors = adjacency[path[-1]]
            if neighbors:
                stack.extend(path + (target,) for target in reversed(neighbors))
            elif len(path) >= 3:
                yield path


def write_timing_chains(outdir, chromosomes):
    """Write chain order without calculating ranks or modifying graph JSONs."""
    outdir = Path(outdir)
    output = outdir / "timing_chains.tsv"
    count = 0
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["file", "chain_index", "chain_position", "node_id", "members"])
        for chrom, _length in chromosomes:
            path = outdir / f"{chrom}_out" / "graphs" / f"{chrom}_condensed.json"
            if not path.exists():
                continue
            with path.open(encoding="utf-8") as graph_file:
                graph = json.load(graph_file)["elements"]
            nodes = {str(node["data"]["id"]): node["data"] for node in graph["nodes"]}
            successors = {node: [] for node in nodes}
            for edge in graph["edges"]:
                data = edge["data"]
                if data["relation"] == "timing":
                    successors[str(data["source"])].append(str(data["target"]))
            for chain_index, chain in enumerate(find_chains(successors, nodes)):
                count += 1
                for position, node in enumerate(chain):
                    members = ",".join(str(member) for member in nodes[node]["members"])
                    writer.writerow([str(path.relative_to(outdir)), chain_index, position, node, members])
    return output, count
