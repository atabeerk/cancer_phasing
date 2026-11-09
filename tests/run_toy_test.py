"""Run the bundled five-mutation example through the complete production pipeline."""

import argparse
import csv
import difflib
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile


REPO = Path(__file__).resolve().parents[1]
FIXTURE = REPO / "toy_data"
FIELDS = ["CHR", "POSITION1", "POSITION2", "RELATION", "ALT_ALT", "ALT_REF",
          "REF_ALT", "REF_REF", "TOTAL", "VAF1", "VAF2"]


def check_equal(label, actual, expected):
    if actual != expected:
        raise ValueError(f"{label}: expected {expected!r}, got {actual!r}")


def check_output(output):
    with (FIXTURE / "expected_output.tsv").open() as handle:
        expected = list(csv.DictReader(handle, delimiter="\t"))
    actual = []
    graph_dir = output / "chr1_out" / "graphs"
    for line in (graph_dir / "chr1_accepted_edges.txt").read_text().splitlines():
        chrom, first, second, *values = line.split()
        row = dict(item.split("=", 1) for item in values)
        row.update(CHR=chrom, POSITION1=first, POSITION2=second)
        actual.append({field: row[field] for field in FIELDS})

    # Ignore line order, floating-point presentation, BIC scores and log timings.
    def normalized(rows):
        result = io.StringIO()
        writer = csv.writer(result, delimiter="\t", lineterminator="\n")
        writer.writerow(FIELDS)
        for row in sorted(rows, key=lambda r: (int(r["POSITION1"]), int(r["POSITION2"]))):
            writer.writerow(f"{float(row[field]):.6f}" if field.startswith("VAF")
                            else row[field] for field in FIELDS)
        return result.getvalue().splitlines(keepends=True)

    difference = "".join(difflib.unified_diff(normalized(expected), normalized(actual),
                                            fromfile="expected_output.tsv", tofile="actual output"))
    if difference:
        raise ValueError("Accepted edges differ:\n" + difference)
    check_equal("Rejected edges", (graph_dir / "chr1_rejected_edges.txt").read_text().strip(), "")

    # Verify that the final graph contains exactly the expected edges and counts.
    graph = json.loads((graph_dir / "chr1.json").read_text())["elements"]
    expected_edges = sorted((r["POSITION1"], r["POSITION2"], r["RELATION"],
                             "/".join(r[k] for k in ("ALT_ALT", "ALT_REF", "REF_ALT", "REF_REF")))
                            for r in expected)
    actual_edges = sorted((e["data"]["source"], e["data"]["target"], e["data"]["relation"],
                           e["data"]["read_counts"]) for e in graph["edges"])
    check_equal("Graph edges", actual_edges, expected_edges)
    vafs = {r[f"POSITION{i}"]: float(r[f"VAF{i}"]) for r in expected for i in (1, 2)}
    check_equal("Graph nodes", sorted(n["data"]["id"] for n in graph["nodes"]), sorted(vafs))
    for node in graph["nodes"]:
        data = node["data"]
        check_equal(f"Node {data['id']} VAF", round(float(data["vaf"]), 6), vafs[data["id"]])
        check_equal(f"Node {data['id']} haplotype", data["haplotype"], "HP1")

    # All ALT reads are HP1; the 25 all-reference reads are HP2.
    with (output / "haplotagged_snvs.tsv").open() as handle:
        snvs = list(csv.DictReader(handle, delimiter="\t"))
    check_equal("Haplotagged positions", sorted(r["POSITION"] for r in snvs), sorted(vafs))
    for row in snvs:
        position = row["POSITION"]
        supporting = round(50 * vafs[position])
        values = {"CHR": "chr1", "TOTAL_COVERAGE": "50", "HP_ASSIGNMENT": "HP1",
                  "HP1_SUPPORTING": str(supporting), "HP1_UNSUPPORTING": str(25 - supporting),
                  "HP2_SUPPORTING": "0", "HP2_UNSUPPORTING": "25",
                  "UNKNOWN_SUPPORTING": "0", "UNKNOWN_UNSUPPORTING": "0",
                  "HP1_HVAF": f"{supporting / 25:.6f}", "HP2_HVAF": "0.000000"}
        for field, value in values.items():
            check_equal(f"SNV {position} {field}", row[field], value)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path,
                        help="New or empty output directory; defaults to a retained temporary directory.")
    args = parser.parse_args()
    output = args.output_dir.resolve() if args.output_dir else Path(tempfile.mkdtemp(prefix="evophase-toy-"))
    print(f"Toy test output: {output}", flush=True)
    subprocess.run([
        sys.executable, "run_genome_by_chr.py",
        "--somatic-vcf", str(FIXTURE / "somatic.vcf.gz"),
        "--bam", str(FIXTURE / "reads.bam"),
        "--germline-vcf", str(FIXTURE / "germline.vcf.gz"),
        "--germline-vcf-sample-name", "NORMAL", "--divergent-same-hp",
        "--window-size", "4", "--window-stride", "3", "--jobs", "1",
        "--output-dir", str(output),
    ], cwd=REPO / "evo", check=True)
    check_output(output)
    print("Toy test passed: 5 SNVs, 50 reads, 1 cooccurring / 7 ordered / 2 divergent pairs.")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError) as error:
        print(f"Toy test failed: {error}", file=sys.stderr)
        sys.exit(1)
