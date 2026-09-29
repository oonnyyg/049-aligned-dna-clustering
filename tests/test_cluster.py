import unittest
from fractions import Fraction

from alignment_cluster.cluster import upgma
from alignment_cluster.distance import build_distance_matrix
from alignment_cluster.fasta import parse_fasta


def cluster(text):
    records = parse_fasta(text)
    ids, pairs = build_distance_matrix(records)
    return upgma(ids, pairs)


def collect(node):
    if node.is_leaf:
        return node.identifier
    return [collect(child) for child in node.children]


class TestUpgma(unittest.TestCase):
    def test_all_leaves_and_zero_distance_kept(self):
        root, merges = cluster(">a\nACGT\n>b\nACGT\n>c\nTCGA\n")
        self.assertEqual(len(merges), 2)
        leaves = sorted(root.members)
        self.assertEqual(leaves, ["a", "b", "c"])
        self.assertEqual(merges[0].distance, Fraction(0))
        self.assertEqual(merges[0].height, Fraction(0))
        # identical pair merges first, then the cluster joins c
        self.assertEqual(merges[0].left_members, ("a",))
        self.assertEqual(merges[0].right_members, ("b",))
        self.assertEqual(merges[1].left_members, ("a", "b"))
        self.assertEqual(merges[1].right_members, ("c",))

    def test_tie_breaking_by_min_identifier_pair(self):
        # a/b and c/d differ on 1 of 2 sites; cross pairs differ on 2 of 2.
        text = ">a\nAA\n>b\nAC\n>c\nGG\n>d\nGT\n"
        root, merges = cluster(text)
        first_pairs = [
            (m.left_members[0], m.right_members[0]) for m in merges[:2]
        ]
        self.assertEqual(first_pairs, [("a", "b"), ("c", "d")])
        # final cluster distance is the mean of four cross pairs (all 1.0)
        self.assertEqual(merges[2].distance, Fraction(1))
        self.assertEqual(merges[2].height, Fraction(1, 2))
        # children ordered by their minimum identifier
        self.assertEqual(root.children[0].min_id, "a")
        self.assertEqual(root.children[1].min_id, "c")

    def test_unequal_cluster_average_uses_all_leaf_pairs(self):
        # Synthetic pairs: a,b,c merge internally first; d stays alone.
        # {a,b,c} vs d must average three distances (size-weighted):
        # (1/2 + 1/4 + 1) / 3 = 7/12, not a two-cluster equal-weight mean.
        from alignment_cluster.distance import PairComparison

        def pair(a, b, value):
            value = Fraction(value)
            return PairComparison(a, b, value.denominator,
                                  value.numerator, value)

        pairs = {
            ("a", "b"): pair("a", "b", 0),
            ("a", "c"): pair("a", "c", 0),
            ("b", "c"): pair("b", "c", 0),
            ("a", "d"): pair("a", "d", "1/2"),
            ("b", "d"): pair("b", "d", "1/4"),
            ("c", "d"): pair("c", "d", 1),
        }
        _, merges = upgma(["a", "b", "c", "d"], pairs)
        final = merges[-1]
        self.assertEqual(final.left_members, ("a", "b", "c"))
        self.assertEqual(final.right_members, ("d",))
        self.assertEqual(final.distance, Fraction(7, 12))
        self.assertEqual(final.height, Fraction(7, 24))

    def test_input_reordering_does_not_change_tree(self):
        text1 = ">a\nACGT\n>b\nACGA\n>c\nTCGT\n>d\nTCGA\n"
        text2 = ">d\nTCGA\n>b\nACGA\n>a\nACGT\n>c\nTCGT\n"
        root1, merges1 = cluster(text1)
        root2, merges2 = cluster(text2)
        self.assertEqual(collect(root1), collect(root2))
        self.assertEqual(
            [(m.distance, m.height) for m in merges1],
            [(m.distance, m.height) for m in merges2],
        )

    def test_branch_lengths_are_height_differences(self):
        root, merges = cluster(">a\nAC\n>b\nAG\n>c\nTC\n")

        def check(node):
            for child in node.children:
                self.assertEqual(node.height - child.height, child.branch_length)
                self.assertGreaterEqual(child.branch_length, 0)
                check(child)

        check(root)
        self.assertEqual(root.branch_length, 0)


if __name__ == "__main__":
    unittest.main()
