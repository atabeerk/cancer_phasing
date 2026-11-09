#!/usr/bin/env python3
"""Replace Evo's pair classifier with shared-tree inference, retaining its pipeline."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from contextlib import ExitStack
import csv
import json
from pathlib import Path

from read_observations import read_mutations
from fit_window_trees import fit_window_trees
from relationship_calls import (
    calls_from_rows, greedy_order, passes_relationship_support,
)

SUFFIXES = {
    "cooccurring": "_cooccurring.txt",
    "left_before_right": "_snp1_before_snp2.txt",
    "right_before_left": "_snp2_before_snp1.txt",
    "divergent": "_divergent.txt",
}

FIT_SETTINGS = {
    "minimum_informative_reads": 3,
    "minimum_alt_alt_reads": 2,
    "em_iterations": 50,
}


def annotation_fields(annotation):
    hp1 = int(annotation["HP1_SUPPORTING"])
    hp2 = int(annotation["HP2_SUPPORTING"])
    alt = hp1 + hp2 + int(annotation["UNKNOWN_SUPPORTING"])
    ref = sum(int(annotation[name]) for name in (
        "HP1_UNSUPPORTING", "HP2_UNSUPPORTING", "UNKNOWN_UNSUPPORTING",
    ))
    label = annotation["HP_ASSIGNMENT"]
    hvaf = annotation.get(f"{label}_HVAF", "NA") if label in ("HP1", "HP2") else "NA"
    return {
        "vaf": alt / (alt + ref) if alt + ref else 0.0,
        "haplotype_vaf": hvaf,
        "haplotype": label,
        "hp_reads": f"{hp1}/{hp2}",
    }


def format_call(call, annotations):
    left = annotations[(call.chrom, call.left)]
    right = annotations[(call.chrom, call.right)]
    bic_margin = call.bic_margin
    fields = {
        "VAF1": left["vaf"], "VAF2": right["vaf"],
        "HVAF1": left["haplotype_vaf"], "HVAF2": right["haplotype_vaf"],
        "HAP1": left["haplotype"], "HAP2": right["haplotype"],
        "HP_READS1": left["hp_reads"], "HP_READS2": right["hp_reads"],
        "ALT_ALT": call.alt_alt, "ALT_REF": call.alt_ref,
        "REF_ALT": call.ref_alt, "REF_REF": call.ref_ref,
        "TOTAL": call.spanning_reads,
        # Evo sorts by RELIABILITY; use the agreed BIC margin as its priority.
        "RELIABILITY": bic_margin,
        "BEST_SCORE": -call.bics[call.best_relation], "MARGIN": bic_margin,
    }
    return f"{call.chrom} {call.left} {call.right} " + " ".join(
        f"{name}={value:.17g}" if isinstance(value, float) else f"{name}={value}"
        for name, value in fields.items()
    ) + "\n"


def emit_relationships(calls, chrom, annotations, chunk_dir,
                       divergent_same_hp=False):
    formatted_annotations = {
        key: annotation_fields(annotation)
        for key, annotation in annotations.items() if key[0] == chrom
    }
    positions = [position for _, position in formatted_annotations]
    base = chunk_dir / f"chunk_{chrom}_{min(positions)}_{max(positions)}"
    counts = Counter(inferred_pairs=len(calls))
    with ExitStack() as stack:
        canonical = stack.enter_context(Path(str(base) + ".txt").open("w"))
        outputs = {
            relation: stack.enter_context(Path(str(base) + suffix).open("w"))
            for relation, suffix in SUFFIXES.items()
        }
        for call in greedy_order(calls.values()):
            formatted = format_call(call, formatted_annotations)
            canonical.write(formatted)
            if divergent_same_hp and call.best_relation == "divergent":
                hp_left = formatted_annotations[(call.chrom, call.left)]["haplotype"]
                hp_right = formatted_annotations[(call.chrom, call.right)]["haplotype"]
                if hp_left != hp_right or hp_left not in ("HP1", "HP2"):
                    continue
            if not passes_relationship_support(call):
                continue
            outputs[call.best_relation].write(formatted)
            counts[call.best_relation] += 1
            counts["graph_candidate_pairs"] += 1
    return dict(counts)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--evo-pileup-bam-count", type=int, required=True)
    parser.add_argument("--divergent-same-hp", action="store_true")
    parser.add_argument("--window-size", type=int, choices=range(2, 6), default=4)
    parser.add_argument("--window-stride", type=int, default=3)
    args = parser.parse_args()
    if args.window_stride < 1:
        parser.error("--window-stride must be >= 1")
    if args.evo_pileup_bam_count < 1:
        parser.error("--evo-pileup-bam-count must be positive")
    return args


def load_retained_mutations(vcf_path, annotation_path):
    with annotation_path.open() as handle:
        annotations = {
            (row["CHR"], int(row["POSITION"])): row
            for row in csv.DictReader(handle, delimiter="\t")
        }
    # Use exactly the SNVs retained by exclusion and somatic support filtering.
    mutations = read_mutations(vcf_path, wanted=list(annotations))
    return mutations, annotations


def fit_chromosome(args, chrom, mutations):
    settings = argparse.Namespace(
        pileup_dir=args.output_dir / "pileup_files",
        output_dir=args.output_dir / "inference" / chrom,
        chrom=chrom,
        window_size=args.window_size,
        window_stride=args.window_stride,
        evo_pileup_bam_count=args.evo_pileup_bam_count,
        **FIT_SETTINGS,
    )
    pair_rows = fit_window_trees(mutations, settings)
    return calls_from_rows(pair_rows.values())


def main():
    args = parse_args()
    mutations, annotations = load_retained_mutations(
        args.vcf, args.output_dir / "snv_counts.tsv"
    )
    by_chrom = defaultdict(list)
    for key in sorted(mutations):
        by_chrom[key[0]].append(key)
    inference = args.output_dir / "inference"
    inference.mkdir()
    summary = {
        "method": f"shared_tree_{args.window_size}ms{args.window_stride}",
        "maximum_mutations": args.window_size, "group_stride": args.window_stride,
        "group_selection": "stride", **FIT_SETTINGS,
        "relationship_support_gate": True, "greedy_priority": "bic-margin",
        "haplotypes_used_in_tree_fitting": False,
        "divergent_same_hp": args.divergent_same_hp,
        "chromosomes": {},
    }
    for chrom, keys in by_chrom.items():
        calls = fit_chromosome(args, chrom, {key: mutations[key] for key in keys})
        counts = emit_relationships(
            calls, chrom, annotations,
            args.output_dir / "chunk_files", args.divergent_same_hp,
        )
        summary["chromosomes"][chrom] = counts
        print(f"[shared_tree_relationships] chromosome={chrom} " + " ".join(
            f"{key}={value}" for key, value in counts.items()
        ), flush=True)
    (inference / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
