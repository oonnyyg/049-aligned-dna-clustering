import unittest
from fractions import Fraction

from alignment_cluster.distance import (
    DistanceError,
    build_distance_matrix,
)
from alignment_cluster.fasta import parse_fasta


class TestDistances(unittest.TestCase):
    def _matrix(self, text):
        return build_distance_matrix(parse_fasta(text))

    def test_pairwise_masking_is_pair_local(self):
        # column layout:
        # 0: a=A b=N c=A  -> (a,b) skipped, (a,c) comp equal, (b,c) skipped
        # 1: a=C b=C c=N  -> (a,b) comp equal, (a,c) skipped, (b,c) skipped
        # 2: a=G b=- c=T  -> (a,b) skipped, (a,c) comp diff,  (b,c) skipped
        # 3: a=T b=A c=T  -> all pairs comparable; (a,b) diff, others equal
        text = ">a\nACGT\n>b\nNC-A\n>c\nANTT\n"
        ids, pairs = self._matrix(text)
        self.assertEqual(ids, ["a", "b", "c"])
        ab, ac, bc = pairs[("a", "b")], pairs[("a", "c")], pairs[("b", "c")]
        self.assertEqual((ab.comparable, ab.differences), (2, 1))
        self.assertEqual(ab.distance, Fraction(1, 2))
        self.assertEqual((ac.comparable, ac.differences), (3, 1))
        self.assertEqual(ac.distance, Fraction(1, 3))
        # b,c share only column 3: A vs T -> one differing comparable site
        self.assertEqual((bc.comparable, bc.differences), (1, 1))
        self.assertEqual(bc.distance, Fraction(1, 1))

    def test_zero_distance_pairs_kept(self):
        ids, pairs = self._matrix(">a\nACGT\n>b\nACGT\n>c\nACGN\n")
        self.assertEqual(pairs[("a", "b")].distance, 0)
        self.assertEqual(pairs[("a", "c")].distance, 0)
        self.assertEqual(ids, ["a", "b", "c"])

    def test_no_comparable_sites_fails_and_names_pair(self):
        text = ">a\nNNNN\n>b\n----\n>c\nACGT\n"
        with self.assertRaisesRegex(DistanceError, r"a/b"):
            self._matrix(text)

    def test_all_incomparable_pairs_reported(self):
        text = ">a\nNN\n>b\n--\n>c\nNN\n"
        with self.assertRaises(DistanceError) as ctx:
            self._matrix(text)
        message = str(ctx.exception)
        self.assertIn("a/b", message)
        self.assertIn("a/c", message)
        self.assertIn("b/c", message)


if __name__ == "__main__":
    unittest.main()
