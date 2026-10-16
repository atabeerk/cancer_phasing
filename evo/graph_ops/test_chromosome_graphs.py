import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


GRAPH_OPS_DIR = Path(__file__).resolve().parent
if str(GRAPH_OPS_DIR) not in sys.path:
    sys.path.insert(0, str(GRAPH_OPS_DIR))

from parsing import group_chunk_bases_by_chromosome, parse_chunk_base
from main import merge_chunk_outputs_by_chromosome


def relation_line(chrom, pos1, pos2, reliability):
    return (
        f"{chrom} {pos1} {pos2} "
        "VAF1=0.2 VAF2=0.2 HAP1=HP1 HAP2=HP1 "
        "HP_READS1=3/0 HP_READS2=3/0 "
        "ALT_ALT=3 ALT_REF=0 REF_ALT=0 REF_REF=3 TOTAL=6 "
        f"RELIABILITY={reliability} BEST_SCORE=1.0 MARGIN=1.0\n"
    )


class ChunkGroupingTests(unittest.TestCase):
    def test_parse_chunk_base_allows_underscores_in_chromosome(self):
        self.assertEqual(
            parse_chunk_base("/tmp/chunk_chrUn_gl000220_10_10000009"),
            ("chrUn_gl000220", 10, 10000009),
        )

    def test_grouping_orders_cores_by_coordinate(self):
        groups = group_chunk_bases_by_chromosome([
            "/tmp/chunk_chr2_10000001_20000000",
            "/tmp/chunk_chr1_1_10000000",
            "/tmp/chunk_chr2_1_10000000",
        ])
        self.assertEqual(list(groups), ["chr1", "chr2"])
        self.assertEqual(
            [record[1] for record in groups["chr2"]],
            [1, 10000001],
        )


    def test_merge_chunk_outputs_preserves_suffixes_and_genomic_order(self):
        with tempfile.TemporaryDirectory() as tmp:
            chunks = Path(tmp)
            left = chunks / "chunk_chrT_1_10000000"
            right = chunks / "chunk_chrT_10000001_20000000"
            left.with_suffix(".txt").write_text("raw-left\n", encoding="utf-8")
            right.with_suffix(".txt").write_text("raw-right\n", encoding="utf-8")
            Path(str(left) + "_cooccurring.txt").write_text(
                "co-left\n", encoding="utf-8"
            )
            Path(str(right) + "_cooccurring.txt").write_text(
                "co-right\n", encoding="utf-8"
            )
            Path(str(right) + "_errors.txt").write_text(
                "error-right\n", encoding="utf-8"
            )
            unrelated = chunks / "notes.txt"
            unrelated.write_text("keep\n", encoding="utf-8")

            merged = merge_chunk_outputs_by_chromosome(
                str(chunks), [str(right), str(left)]
            )

            merged_base = chunks / "chrT"
            self.assertEqual(merged, [str(merged_base)])
            self.assertEqual(
                merged_base.with_suffix(".txt").read_text(encoding="utf-8"),
                "raw-left\nraw-right\n",
            )
            self.assertEqual(
                Path(str(merged_base) + "_cooccurring.txt").read_text(
                    encoding="utf-8"
                ),
                "co-left\nco-right\n",
            )
            self.assertEqual(
                Path(str(merged_base) + "_errors.txt").read_text(encoding="utf-8"),
                "error-right\n",
            )
            self.assertFalse(left.with_suffix(".txt").exists())
            self.assertFalse(right.with_suffix(".txt").exists())
            self.assertFalse(Path(str(left) + "_cooccurring.txt").exists())
            self.assertFalse(Path(str(right) + "_cooccurring.txt").exists())
            self.assertEqual(unrelated.read_text(encoding="utf-8"), "keep\n")


    def test_merge_single_chunk_uses_chromosome_base(self):
        with tempfile.TemporaryDirectory() as tmp:
            chunks = Path(tmp)
            source = chunks / "chunk_chrT_1_10000000"
            source.with_suffix(".txt").write_text("raw\n", encoding="utf-8")
            Path(str(source) + "_timing.txt").write_text(
                "timing\n", encoding="utf-8"
            )

            merged = merge_chunk_outputs_by_chromosome(
                str(chunks), [str(source)]
            )

            merged_base = chunks / "chrT"
            self.assertEqual(merged, [str(merged_base)])
            self.assertEqual(
                merged_base.with_suffix(".txt").read_text(encoding="utf-8"),
                "raw\n",
            )
            self.assertEqual(
                Path(str(merged_base) + "_timing.txt").read_text(
                    encoding="utf-8"
                ),
                "timing\n",
            )
            self.assertFalse(source.with_suffix(".txt").exists())
            self.assertFalse(Path(str(source) + "_timing.txt").exists())



class ChromosomeGraphIntegrationTests(unittest.TestCase):
    def test_edges_from_adjacent_cores_form_one_chromosome_graph(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            chunks = root / "chunks"
            graphs = root / "graphs"
            chunks.mkdir()

            left = chunks / "chunk_chrT_1_10000000"
            right = chunks / "chunk_chrT_10000001_20000000"
            left.with_suffix(".txt").write_text("", encoding="utf-8")
            right.with_suffix(".txt").write_text("", encoding="utf-8")
            Path(str(left) + "_cooccurring.txt").write_text(
                relation_line("chrT", 9999999, 10000001, 0.9),
                encoding="utf-8",
            )
            Path(str(right) + "_cooccurring.txt").write_text(
                relation_line("chrT", 10000001, 10000002, 0.8),
                encoding="utf-8",
            )

            result = subprocess.run(
                [
                    sys.executable,
                    str(GRAPH_OPS_DIR / "main.py"),
                    str(chunks),
                    "--outdir",
                    str(graphs),
                    "--total-snvs",
                    "3",
                    "--progress-interval",
                    "0",
                ],
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

            graph_paths = sorted(
                path for path in graphs.glob("*.json")
                if not path.name.endswith("_condensed.json")
            )
            self.assertEqual(
                [path.name for path in graph_paths],
                ["chrT.json"],
            )
            graph = json.loads(graph_paths[0].read_text(encoding="utf-8"))
            node_ids = {
                int(node["data"]["id"])
                for node in graph["elements"]["nodes"]
            }
            self.assertEqual(node_ids, {9999999, 10000001, 10000002})
            self.assertGreaterEqual(min(node_ids), 1)
            self.assertLessEqual(max(node_ids), 20000000)

            stats_path = root / "component_statistics.txt"
            with stats_path.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 1)
            self.assertEqual(int(rows[0]["num_nodes"]), 3)
            self.assertEqual(rows[0]["chunk_base"], "chrT")

            condensed_paths = sorted(graphs.glob("*_condensed.json"))
            self.assertEqual(
                [path.name for path in condensed_paths],
                ["chrT_condensed.json"],
            )
            chunk_paths = sorted(path.name for path in chunks.glob("*.txt"))
            self.assertEqual(
                chunk_paths,
                [
                    "chrT.txt",
                    "chrT_cooccurring.txt",
                ],
            )

    def test_region_edge_files_match_graph(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            chunks = root / "chunks"
            graphs = root / "graphs"
            chunks.mkdir()

            base = chunks / "chunk_chrT_1_10000000"
            base.with_suffix(".txt").write_text("", encoding="utf-8")
            Path(str(base) + "_cooccurring.txt").write_text(
                relation_line("chrT", 1, 2, 0.9),
                encoding="utf-8",
            )
            Path(str(base) + "_divergent.txt").write_text(
                relation_line("chrT", 1, 2, 0.8),
                encoding="utf-8",
            )

            result = subprocess.run(
                [
                    sys.executable,
                    str(GRAPH_OPS_DIR / "main.py"),
                    str(chunks),
                    "--outdir",
                    str(graphs),
                    "--progress-interval",
                    "0",
                ],
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

            base_name = "chrT"
            graph = json.loads(
                (graphs / f"{base_name}.json").read_text(encoding="utf-8")
            )
            graph_edges = graph["elements"]["edges"]
            self.assertEqual(len(graph_edges), 1)

            accepted_lines = (graphs / f"{base_name}_accepted_edges.txt").read_text(
                encoding="utf-8"
            ).splitlines()
            self.assertEqual(len(accepted_lines), 1)
            accepted_fields = accepted_lines[0].split()
            self.assertEqual(
                {int(accepted_fields[1]), int(accepted_fields[2])},
                {1, 2},
            )

            with (graphs / f"{base_name}_inconsistencies.tsv").open(
                newline="", encoding="utf-8"
            ) as handle:
                inconsistencies = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(inconsistencies), 1)
            self.assertEqual(
                {int(inconsistencies[0]["u"]), int(inconsistencies[0]["v"])},
                {1, 2},
            )

    def test_disconnected_components_are_written_to_one_chromosome_graph(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            chunks = root / "chunks"
            graphs = root / "graphs"
            chunks.mkdir()

            left = chunks / "chunk_chrT_1_10000000"
            right = chunks / "chunk_chrT_10000001_20000000"
            left.with_suffix(".txt").write_text("", encoding="utf-8")
            right.with_suffix(".txt").write_text("", encoding="utf-8")
            Path(str(left) + "_cooccurring.txt").write_text(
                relation_line("chrT", 1, 2, 0.9),
                encoding="utf-8",
            )
            Path(str(right) + "_cooccurring.txt").write_text(
                relation_line("chrT", 15000000, 15000001, 0.8),
                encoding="utf-8",
            )

            result = subprocess.run(
                [
                    sys.executable,
                    str(GRAPH_OPS_DIR / "main.py"),
                    str(chunks),
                    "--outdir",
                    str(graphs),
                    "--progress-interval",
                    "0",
                ],
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

            graph_paths = sorted(
                path for path in graphs.glob("*.json")
                if not path.name.endswith("_condensed.json")
            )
            self.assertEqual(
                [path.name for path in graph_paths],
                ["chrT.json"],
            )
            graph = json.loads(graph_paths[0].read_text(encoding="utf-8"))
            self.assertEqual(
                {
                    int(node["data"]["id"])
                    for node in graph["elements"]["nodes"]
                },
                {1, 2, 15000000, 15000001},
            )

    def test_component_wider_than_target_is_not_split(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            chunks = root / "chunks"
            graphs = root / "graphs"
            chunks.mkdir()

            left = chunks / "chunk_chrT_1_10000000"
            right = chunks / "chunk_chrT_10000001_20000000"
            left.with_suffix(".txt").write_text("", encoding="utf-8")
            right.with_suffix(".txt").write_text("", encoding="utf-8")
            Path(str(left) + "_cooccurring.txt").write_text(
                relation_line("chrT", 1, 12000001, 0.9),
                encoding="utf-8",
            )

            result = subprocess.run(
                [
                    sys.executable,
                    str(GRAPH_OPS_DIR / "main.py"),
                    str(chunks),
                    "--outdir",
                    str(graphs),
                    "--progress-interval",
                    "0",
                ],
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

            graph_paths = sorted(
                path for path in graphs.glob("*.json")
                if not path.name.endswith("_condensed.json")
            )
            self.assertEqual(
                [path.name for path in graph_paths],
                ["chrT.json"],
            )
            graph = json.loads(graph_paths[0].read_text(encoding="utf-8"))
            self.assertEqual(
                {
                    int(node["data"]["id"])
                    for node in graph["elements"]["nodes"]
                },
                {1, 12000001},
            )

    def test_empty_chromosome_graph_still_writes_json_pair(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            chunks = root / "chunks"
            graphs = root / "graphs"
            chunks.mkdir()

            base = chunks / "chunk_chrT_1_10000000"
            base.with_suffix(".txt").write_text("", encoding="utf-8")

            result = subprocess.run(
                [
                    sys.executable,
                    str(GRAPH_OPS_DIR / "main.py"),
                    str(chunks),
                    "--outdir",
                    str(graphs),
                    "--progress-interval",
                    "0",
                ],
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

            graph_path = graphs / "chrT.json"
            condensed_path = graphs / "chrT_condensed.json"
            self.assertTrue(graph_path.exists())
            self.assertTrue(condensed_path.exists())
            graph = json.loads(graph_path.read_text(encoding="utf-8"))
            condensed = json.loads(condensed_path.read_text(encoding="utf-8"))
            self.assertEqual(graph["elements"], {"nodes": [], "edges": []})
            self.assertEqual(condensed["elements"], {"nodes": [], "edges": []})



if __name__ == "__main__":
    unittest.main()
