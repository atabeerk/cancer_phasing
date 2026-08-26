import unittest

from builder import GraphBuilder
from fast_builder import FastGraphBuilder
from models import Edge


def edge(u, v, relation, reliability=1.0):
    return Edge(
        chrom="chrT",
        u=u,
        v=v,
        relation=relation,
        reliability=reliability,
        alt_alt=1,
        alt_ref=1,
        ref_alt=1,
        ref_ref=1,
        vaf_u=0.2,
        vaf_v=0.3,
    )


def accepted_key(e):
    return (e.u, e.v, e.relation, e.reliability)


def cluster_sets(builder):
    return sorted(tuple(sorted(nodes)) for nodes in builder.cluster_nodes.values())


class FastBuilderEquivalenceTests(unittest.TestCase):
    def assert_equivalent(self, edges):
        legacy = GraphBuilder()
        fast = FastGraphBuilder()
        legacy_results = [legacy.try_add_edge(e) for e in edges]
        fast_results = [fast.try_add_edge(e) for e in edges]

        self.assertEqual(legacy_results, fast_results)
        self.assertEqual([accepted_key(e) for e in legacy.edges], [accepted_key(e) for e in fast.edges])
        self.assertEqual(
            [r["reason"] for r in legacy.inconsistencies],
            [r["reason"] for r in fast.inconsistencies],
        )
        self.assertEqual(cluster_sets(legacy), cluster_sets(fast))
        self.assertEqual(legacy.nodes, fast.nodes)
        return legacy, fast, legacy_results

    def assert_rejected_with_reason(self, edges, expected_results, expected_reason):
        legacy, fast, results = self.assert_equivalent(edges)
        self.assertEqual(results, expected_results)
        self.assertEqual([r["reason"] for r in legacy.inconsistencies], [expected_reason])
        self.assertEqual([r["reason"] for r in fast.inconsistencies], [expected_reason])

    def test_cooccurring_rejects_internal_timing_and_divergent(self):
        self.assert_equivalent([
            edge(1, 2, "cooccurring", 3.0),
            edge(1, 2, "timing", 2.0),
            edge(1, 2, "divergent", 1.0),
        ])

    def test_direct_timing_rejects_divergent_and_cooccurring(self):
        self.assert_equivalent([
            edge(1, 2, "timing", 3.0),
            edge(1, 2, "divergent", 2.0),
            edge(1, 2, "cooccurring", 1.0),
        ])

    def test_transitive_timing_rejects_divergent_and_cooccurring(self):
        self.assert_equivalent([
            edge(1, 2, "timing", 5.0),
            edge(2, 3, "timing", 4.0),
            edge(1, 3, "divergent", 3.0),
            edge(1, 3, "cooccurring", 2.0),
        ])

    def test_long_transitive_timing_rejects_cooccurring_cycle(self):
        self.assert_equivalent([
            edge(1, 2, "timing", 5.0),
            edge(2, 3, "timing", 4.0),
            edge(3, 4, "timing", 3.0),
            edge(1, 4, "cooccurring", 2.0),
        ])

    def test_timing_cycle(self):
        self.assert_equivalent([
            edge(1, 2, "timing", 3.0),
            edge(2, 3, "timing", 2.0),
            edge(3, 1, "timing", 1.0),
        ])

    def test_timing_implies_divergent_conflict(self):
        self.assert_equivalent([
            edge(1, 4, "divergent", 5.0),
            edge(2, 3, "timing", 4.0),
            edge(1, 2, "timing", 3.0),
            edge(3, 4, "timing", 2.0),
        ])

    def test_divergent_clusters_cannot_share_direct_timing_descendant(self):
        self.assert_rejected_with_reason(
            [
                edge(1, 2, "divergent", 3.0),
                edge(1, 3, "timing", 2.0),
                edge(2, 3, "timing", 1.0),
            ],
            [True, True, False],
            "cluster_pair_divergent_shared_timing_descendant_conflict",
        )

    def test_divergent_edge_rejected_after_timing_convergence(self):
        self.assert_rejected_with_reason(
            [
                edge(1, 3, "timing", 3.0),
                edge(2, 3, "timing", 2.0),
                edge(1, 2, "divergent", 1.0),
            ],
            [True, True, False],
            "cluster_pair_divergent_shared_timing_descendant_conflict",
        )

    def test_divergent_clusters_cannot_share_transitive_timing_descendant(self):
        self.assert_rejected_with_reason(
            [
                edge(1, 2, "divergent", 5.0),
                edge(1, 3, "timing", 4.0),
                edge(3, 5, "timing", 3.0),
                edge(2, 4, "timing", 2.0),
                edge(4, 5, "timing", 1.0),
            ],
            [True, True, True, True, False],
            "cluster_pair_divergent_shared_timing_descendant_conflict",
        )

    def test_rule_applies_to_cooccurring_clusters(self):
        self.assert_rejected_with_reason(
            [
                edge(1, 10, "cooccurring", 5.0),
                edge(2, 20, "cooccurring", 4.0),
                edge(1, 2, "divergent", 3.0),
                edge(10, 3, "timing", 2.0),
                edge(20, 3, "timing", 1.0),
            ],
            [True, True, True, True, False],
            "cluster_pair_divergent_shared_timing_descendant_conflict",
        )

    def test_cooccurring_merge_cannot_create_shared_descendant(self):
        self.assert_rejected_with_reason(
            [
                edge(1, 2, "divergent", 4.0),
                edge(1, 3, "timing", 3.0),
                edge(2, 4, "timing", 2.0),
                edge(3, 4, "cooccurring", 1.0),
            ],
            [True, True, True, False],
            "cluster_merge_divergent_shared_timing_descendant_conflict",
        )

    def test_divergent_clusters_may_have_separate_descendants(self):
        legacy, fast, results = self.assert_equivalent([
            edge(1, 2, "divergent", 3.0),
            edge(1, 3, "timing", 2.0),
            edge(2, 4, "timing", 1.0),
        ])
        self.assertEqual(results, [True, True, True])
        self.assertEqual(legacy.inconsistencies, [])
        self.assertEqual(fast.inconsistencies, [])

    def test_divergent_clusters_may_share_timing_ancestor(self):
        legacy, fast, results = self.assert_equivalent([
            edge(3, 1, "timing", 3.0),
            edge(3, 2, "timing", 2.0),
            edge(1, 2, "divergent", 1.0),
        ])
        self.assertEqual(results, [True, True, True])
        self.assertEqual(legacy.inconsistencies, [])
        self.assertEqual(fast.inconsistencies, [])

    def test_cluster_boundary_direction_conflict(self):
        self.assert_equivalent([
            edge(1, 2, "cooccurring", 6.0),
            edge(3, 4, "cooccurring", 5.0),
            edge(1, 3, "timing", 4.0),
            edge(4, 2, "timing", 3.0),
        ])

    def test_cluster_boundary_relation_conflict(self):
        self.assert_equivalent([
            edge(1, 2, "cooccurring", 6.0),
            edge(3, 4, "cooccurring", 5.0),
            edge(1, 3, "divergent", 4.0),
            edge(2, 4, "timing", 3.0),
        ])


if __name__ == "__main__":
    unittest.main()
