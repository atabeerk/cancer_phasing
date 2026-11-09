#!/usr/bin/env python3

import csv
import os

RELATIONS = ("cooccurring", "timing", "divergent")


def merged_interval_coverage_bp(intervals):
    """Return covered bp after merging inclusive [start, end] intervals."""
    if not intervals:
        return 0

    intervals = sorted(intervals)
    merged_start, merged_end = intervals[0]
    covered = 0

    for start, end in intervals[1:]:
        if start <= merged_end + 1:
            merged_end = max(merged_end, end)
        else:
            covered += merged_end - merged_start + 1
            merged_start, merged_end = start, end

    covered += merged_end - merged_start + 1
    return covered


def component_span_stats_from_stats(component_stats_path):
    """
    Compute two connected-component span metrics from component_statistics.txt:
      1) summed_component_bp: sum of all component spans (counts overlaps multiple times)
      2) unique_component_bp: covered bp after merging overlapping component intervals
         (counts overlaps once)
    """
    if not os.path.exists(component_stats_path):
        return 0, 0

    intervals = []
    summed_component_bp = 0
    with open(component_stats_path, "r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            try:
                start = int(row["min_pos"])
                end = int(row["max_pos"])
            except (KeyError, TypeError, ValueError) as error:
                raise ValueError(f"Invalid component coordinates in {component_stats_path}") from error
            if end < start:
                raise ValueError(f"Reversed component interval in {component_stats_path}: {start}-{end}")
            intervals.append((start, end))
            # Inclusive genomic coordinates: [start, end]
            summed_component_bp += (end - start + 1)

    unique_component_bp = merged_interval_coverage_bp(intervals)
    return summed_component_bp, unique_component_bp


def new_edge_haplotype_summary():
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


MUTATION_STAT_FIELDS = (
    "total_snvs", "in_relations", "in_graph", "with_timing", "orphaned", "singleton",
    "accepted_edges", "timing_edges", "hp1", "hp2", "mixed", "unknown",
)


def _parse_kv_tokens(line):
    return dict(token.split("=", 1) for token in line.split() if "=" in token)


def read_chromosome_log_stats(chrom_log_path):
    """Read current mutation, edge and haplotag status records in one pass."""
    mutations = {}
    edges = new_edge_haplotype_summary()
    haplotag_status = {}
    with open(chrom_log_path, encoding="utf-8") as handle:
        for line in handle:
            prefix = line.split(maxsplit=1)[0] if line.strip() else ""
            if prefix == "[graph_ops_mutation_stats]":
                fields = _parse_kv_tokens(line)
                mutations = {key: int(fields[key]) for key in MUTATION_STAT_FIELDS}
            elif prefix == "[graph_ops_edge_hap]":
                fields = _parse_kv_tokens(line)
                relation = fields["relation"]
                edges[relation] = {key: int(fields[key]) for key in edges[relation]}
            elif prefix == "[haplotag_status]":
                haplotag_status = _parse_kv_tokens(line)
    return mutations, edges, haplotag_status


def add_edge_haplotype_summary(total_summary, chrom_summary):
    for relation in RELATIONS:
        for key in total_summary[relation]:
            total_summary[relation][key] += chrom_summary[relation][key]


def parse_haplotag_detected_value(haplotag_status):
    value = haplotag_status.get("haplotagged")
    return int(value) if value in ("0", "1") else None
