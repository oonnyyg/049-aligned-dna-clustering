"""Automated tests for the aligned DNA clustering pipeline."""

from __future__ import annotations

import csv
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path

from alignment_cluster.clustering import upgma
from alignment_cluster.distances import DistanceError, build_comparisons
from alignment_cluster.fasta import FastaError, parse_fasta
from alignment_cluster.__main__ import main

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def analyze(text: str):
    records = parse_fasta(dedent(text))
    identifiers, comparisons, distances = build_comparisons(records)
    merges, tree = upgma(identifiers, distances)
    return records, identifiers, comparisons, distances, merges, tree


def dedent(text: str) -> str:
    return "\n".join(line.strip() for line in text.strip().splitlines())


class FastaParsingTests(unittest.TestCase):
    def test_multiline_whitespace_uppercase_and_description(self):
        records = parse_fasta(dedent("""
            >s1 hello world
             ac g t
            N - C
            >s2
            ACGTN-C
        """))
        self.assertEqual(records[0].identifier, "s1")
        self.assertEqual(records[0].description, "hello world")
        self.assertEqual(records[0].sequence, "ACGTN-C")
        self.assertEqual(records[1].description, "")

    def test_sequence_before_header_rejected(self):
        with self.assertRaises(FastaError):
            parse_fasta("ACGT\n>s1\nACGT\n>s2\nACGT\n")

    def test_duplicate_identifier_rejected(self):
        with self.assertRaises(FastaError):
            parse_fasta(">s1\nACGT\n>s1 dup\nACGT\n")

    def test_length_mismatch_rejected(self):
        with self.assertRaises(FastaError):
            parse_fasta(">s1\nACGT\n>s2\nACG\n")

    def test_invalid_character_rejected(self):
        with self.assertRaises(FastaError):
            parse_fasta(">s1\nACGX\n>s2\nACGT\n")

    def test_too_few_or_empty_sequences_rejected(self):
        with self.assertRaises(FastaError):
            parse_fasta(">s1\nACGT\n")
        with self.assertRaises(FastaError):
            parse_fasta(">s1\n\n>s2\nACGT\n")


class DistanceTests(unittest.TestCase):
    def test_n_and_gaps_excluded_per_pair(self):
        # Sites: N in s2 excludes site 0 only for pairs containing s2;
        # gap in s3 excludes site 3 only for pairs containing s3.
        _, _, comparisons, _, _, _ = analyze("""
            >s1
            ACGT
            >s2
            NCGT
            >s3
            ACG-
        """)
        by_pair = {(c.identifier_a, c.identifier_b): c for c in comparisons}
        self.assertEqual(by_pair[("s1", "s2")].comparable_sites, 3)
        self.assertEqual(by_pair[("s1", "s3")].comparable_sites, 3)
        self.assertEqual(by_pair[("s2", "s3")].comparable_sites, 2)
        for comparison in comparisons:
            self.assertEqual(comparison.differing_sites, 0)
            self.assertEqual(comparison.distance, 0.0)

    def test_distance_is_differing_over_comparable(self):
        _, _, comparisons, _, _, _ = analyze("""
            >s1
            ACGT
            >s2
            CCGN
        """)
        comparison = comparisons[0]
        self.assertEqual(comparison.comparable_sites, 3)
        self.assertEqual(comparison.differing_sites, 1)
        self.assertAlmostEqual(comparison.distance, 1 / 3)

    def test_no_comparable_sites_fails_with_identifiers(self):
        with self.assertRaises(DistanceError) as context:
            analyze("""
                >s1
                NNNN
                >s2
                ----
            """)
        message = str(context.exception)
        self.assertIn("s1", message)
        self.assertIn("s2", message)


class ClusteringTests(unittest.TestCase):
    def test_zero_distance_keeps_every_leaf(self):
        _, identifiers, _, _, merges, tree = analyze("""
            >s1
            ACGT
            >s2
            ACGT
            >s3
            ACGT
        """)
        leaves = [identifier for identifier in identifiers]
        self.assertEqual(sorted(tree.members), leaves)
        self.assertEqual(len(merges), 2)
        self.assertTrue(all(merge["inter_cluster_distance"] == 0.0 for merge in merges))

    def test_ties_break_by_min_identifier_pair(self):
        # Every pair differs at every site: all distances tie at 1.0.
        _, _, _, _, merges, tree = analyze("""
            >s3
            AAAA
            >s1
            CCCC
            >s2
            GGGG
        """)
        first = merges[0]
        self.assertEqual(first["child_ids"], ["s1", "s2"])
        self.assertAlmostEqual(first["inter_cluster_distance"], 1.0)
        self.assertEqual(merges[1]["child_ids"], ["merge_1", "s3"])
        self.assertEqual(tree.children[0].id, "merge_1")

    def test_unequal_cluster_sizes_use_mean_over_all_pairs(self):
        # s1/s2 identical; s3 and s4 each differ from them at one distinct
        # site of four, and from each other at two sites.
        _, _, _, distances, merges, _ = analyze("""
            >s1
            ACGT
            >s2
            ACGT
            >s3
            CCGT
            >s4
            ATGT
        """)
        # First merge: s1 + s2 at distance zero.
        self.assertEqual(merges[0]["child_ids"], ["s1", "s2"])
        # Distance {s1,s2} -> s3 = mean(0.25, 0.25) = 0.25 (two pairs,
        # not an equal-weight average of cluster "representatives").
        expected = (
            distances[("s1", "s4")] + distances[("s2", "s4")]
        ) / 2
        self.assertAlmostEqual(expected, 0.25)
        second = merges[1]
        self.assertEqual(second["child_ids"], ["merge_1", "s3"])
        self.assertAlmostEqual(second["inter_cluster_distance"], 0.25)

    def test_heights_and_branch_lengths(self):
        _, _, _, _, merges, tree = analyze("""
            >s1
            ACGT
            >s2
            CCGT
        """)
        self.assertAlmostEqual(tree.height, merges[0]["inter_cluster_distance"] / 2)
        self.assertAlmostEqual(tree.children[0].branch_length, tree.height)
        self.assertEqual(tree.children[0].height, 0.0)


class InputOrderIndependenceTests(unittest.TestCase):
    def test_reordered_input_gives_identical_tree(self):
        text_a = dedent("""
            >s3
            ACGT
            >s1
            CCGT
            >s2
            AGGT
        """)
        order_b = [">s2", "AGGT", ">s3", "ACGT", ">s1", "CCGT"]
        text_b = "\n".join(order_b)

        result_a = analyze(text_a)
        result_b = analyze(text_b)
        self.assertEqual(
            json.dumps(result_a[4], sort_keys=True),
            json.dumps(result_b[4], sort_keys=True),
        )
        self.assertEqual(
            json.dumps(result_a[5].to_dict(), sort_keys=True),
            json.dumps(result_b[5].to_dict(), sort_keys=True),
        )


class CommandLineTests(unittest.TestCase):
    def test_end_to_end_reports(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            fasta = tmp_path / "input.fasta"
            out = tmp_path / "reports"
            fasta.write_text(
                ">a desc\nACGT\n>b\nNCGT\n>c\nACG-\n", encoding="utf-8"
            )
            code = main([str(fasta), "-o", str(out)])
            self.assertEqual(code, 0)

            with (out / "distance_matrix.csv").open(encoding="utf-8") as handle:
                rows = list(csv.reader(handle))
            self.assertEqual(rows[0], ["", "a", "b", "c"])
            self.assertEqual(rows[1][0], "a")
            self.assertEqual(rows[1][1], "0")

            with (out / "pairwise_comparisons.csv").open(encoding="utf-8") as handle:
                pair_rows = list(csv.DictReader(handle))
            self.assertEqual(len(pair_rows), 3)
            self.assertEqual(pair_rows[0]["identifier_a"], "a")
            self.assertEqual(pair_rows[0]["identifier_b"], "b")

            payload = json.loads((out / "tree.json").read_text(encoding="utf-8"))
            self.assertEqual(payload["leaf_order"], ["a", "b", "c"])
            self.assertEqual(len(payload["merges"]), 2)

    def test_invalid_input_writes_nothing_and_returns_nonzero(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            fasta = tmp_path / "bad.fasta"
            out = tmp_path / "reports"
            out.mkdir()
            sentinel = out / "keep.txt"
            sentinel.write_text("untouched", encoding="utf-8")
            fasta.write_text(">a\nNNN\n>b\n---\n", encoding="utf-8")
            code = main([str(fasta), "-o", str(out)])
            self.assertNotEqual(code, 0)
            self.assertEqual(set(p.name for p in out.iterdir()), {"keep.txt"})

    def test_example_file_runs(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            code = main(
                [str(ROOT / "examples" / "example.fasta"), "-o", str(out)]
            )
            self.assertEqual(code, 0)
            self.assertTrue((out / "tree.json").exists())


if __name__ == "__main__":
    unittest.main()
