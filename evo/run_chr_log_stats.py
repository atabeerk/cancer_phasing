#!/usr/bin/env python3

import csv
import os
import re

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
            except (KeyError, TypeError, ValueError):
                continue
            if end < start:
                continue
            intervals.append((start, end))
            # Inclusive genomic coordinates: [start, end]
            summed_component_bp += (end - start + 1)

    unique_component_bp = merged_interval_coverage_bp(intervals)
    return summed_component_bp, unique_component_bp


def new_edge_haplotype_summary():
    return {
        relation: {
            "total": 0,
            "without_loss": 0,
            "with_loss": 0,
            "same_haplotype": 0,
            "different_haplotype": 0,
            "at_least_one_unknown": 0,
            "hp1_hp1": 0,
            "hp2_hp2": 0,
        }
        for relation in RELATIONS
    }


def _parse_kv_tokens(line):
    out = {}
    for tok in line.split():
        if "=" not in tok:
            continue
        k, v = tok.split("=", 1)
        out[k] = v
    return out


def extract_mutation_stats_from_chrom_log(chrom_log_path):
    """Parse [graph_ops_mutation_stats] from a chromosome log."""
    stats = {}
    if not os.path.exists(chrom_log_path):
        return stats
    with open(chrom_log_path, "r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s.startswith("[graph_ops_mutation_stats]"):
                continue
            kv = _parse_kv_tokens(s)
            for key in (
                "total_snvs",
                "in_relations",
                "in_graph",
                "with_timing",
                "orphaned",
                "singleton",
                "accepted_edges",
                "timing_edges",
                "hp1",
                "hp2",
                "mixed",
                "unknown",
            ):
                try:
                    stats[key] = int(kv.get(key, "0"))
                except ValueError:
                    stats[key] = 0
    return stats


def extract_edge_haplotype_summary_from_chrom_log(chrom_log_path):
    """Parse [graph_ops_edge_hap] entries from a chromosome log."""
    summary = new_edge_haplotype_summary()
    if not os.path.exists(chrom_log_path):
        return summary
    with open(chrom_log_path, "r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s.startswith("[graph_ops_edge_hap]"):
                continue
            kv = _parse_kv_tokens(s)
            relation = kv.get("relation", "")
            if relation not in summary:
                continue
            for key in summary[relation]:
                try:
                    summary[relation][key] = int(kv.get(key, "0"))
                except ValueError:
                    summary[relation][key] = 0
    return summary


def add_edge_haplotype_summary(total_summary, chrom_summary):
    for relation in RELATIONS:
        for key in total_summary[relation]:
            total_summary[relation][key] += chrom_summary[relation][key]


def extract_max_span_scan_from_chrom_log(chrom_log_path):
    """
    Parse a '[max_span_scan]' line written by ./main from a chromosome log.
    Returns a dict of parsed key/value fields plus the raw line.
    """
    if not os.path.exists(chrom_log_path):
        return {}
    last_line = ""
    with open(chrom_log_path, "r") as f:
        for line in f:
            if line.startswith("[max_span_scan]"):
                last_line = line.strip()
    if not last_line:
        return {}

    out = {"raw": last_line}
    for tok in last_line.split():
        if "=" not in tok:
            continue
        k, v = tok.split("=", 1)
        out[k] = v
    return out


def extract_haplotag_status_from_chrom_log(chrom_log_path):
    """
    Parse BAM haplotag detection status from a chromosome log.
    Returns a dict with keys:
      - haplotagged: "0" or "1" when available
      - source: text source when available
      - raw_status: raw "Using BAM haplotagged status = ..." line when available
      - has_hp_tag: parsed from [haplotag_detect] line when available
    """
    if not os.path.exists(chrom_log_path):
        return {}

    status_line = ""
    detect_line = ""
    with open(chrom_log_path, "r") as f:
        for line in f:
            s = line.strip()
            if s.startswith("Using BAM haplotagged status ="):
                status_line = s
            if s.startswith("[haplotag_detect]"):
                detect_line = s

    out = {}
    if status_line:
        out["raw_status"] = status_line
        m = re.search(r"Using BAM haplotagged status\s*=\s*([01])\s*\(([^)]+)\)", status_line)
        if m:
            out["haplotagged"] = m.group(1)
            out["source"] = m.group(2).strip()
        else:
            m2 = re.search(r"Using BAM haplotagged status\s*=\s*([01])", status_line)
            if m2:
                out["haplotagged"] = m2.group(1)

    if detect_line:
        for tok in detect_line.split():
            if "=" not in tok:
                continue
            k, v = tok.split("=", 1)
            out[k] = v

    if "haplotagged" not in out and "has_hp_tag" in out:
        out["haplotagged"] = out["has_hp_tag"]
        out.setdefault("source", "auto-detected")

    return out


def parse_haplotag_detected_value(haplotag_status):
    """Return 0/1 when parsable, else None."""
    v = haplotag_status.get("haplotagged")
    if v in ("0", "1"):
        return int(v)
    return None
