import csv
import json
import subprocess
import sys
from pathlib import Path

import unittest


REPO = Path(__file__).resolve().parents[1]


VALID_FASTA = (
    ">alpha sample one\nACGT-N\n"
    ">beta sample two\nAC-TNN\n"
    ">gamma sample three\nTCGT-A\n"
)


class TestCli(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(__file__).resolve().parent / "_tmp"
        self.tmp.mkdir(exist_ok=True)
        self.fasta = self.tmp / "in.fasta"
        self.out = self.tmp / "reports"
        self.fasta.write_text(VALID_FASTA, encoding="utf-8")

    def tearDown(self):
        import shutil

        shutil.rmtree(self.tmp, ignore_errors=True)

    def _run(self):
        return subprocess.run(
            [sys.executable, "-m", "alignment_cluster",
             str(self.fasta), "-o", str(self.out)],
            cwd=REPO,
            capture_output=True,
            text=True,
        )

    def test_success_outputs(self):
        result = self._run()
        self.assertEqual(result.returncode, 0, result.stderr)
        names = sorted(p.name for p in self.out.iterdir())
        self.assertEqual(
            names,
            ["cluster_tree.json", "distance_matrix.csv",
             "pairwise_comparisons.csv"],
        )
        with open(self.out / "distance_matrix.csv", newline="",
                  encoding="utf-8") as handle:
            rows = list(csv.reader(handle))
        self.assertEqual(rows[0], ["", "alpha", "beta", "gamma"])
        self.assertEqual(rows[1][0], "alpha")
        self.assertEqual(rows[1][1], "0")
        self.assertEqual(rows[2][2], "0")
        self.assertEqual(rows[3][3], "0")
        self.assertEqual(rows[1][2], rows[2][1])

        tree = json.loads((self.out / "cluster_tree.json").read_text("utf-8"))
        self.assertEqual(len(tree["merges"]), 2)
        self.assertEqual(tree["identifiers"], ["alpha", "beta", "gamma"])
        self.assertEqual(tree["tree"]["members"],
                         ["alpha", "beta", "gamma"])

    def test_invalid_input_nonzero_and_no_output_dir(self):
        self.fasta.write_text(">only\nACGT\n", encoding="utf-8")
        result = self._run()
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(self.out.exists())

    def test_incomparable_pair_nonzero_and_existing_reports_untouched(self):
        stale = self.out
        stale.mkdir(parents=True)
        marker = stale / "distance_matrix.csv"
        marker.write_text("STALE", encoding="utf-8")
        self.fasta.write_text(">a\nNN\n>b\n--\n", encoding="utf-8")
        result = self._run()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("a/b", result.stderr)
        self.assertEqual(marker.read_text(encoding="utf-8"), "STALE")


if __name__ == "__main__":
    unittest.main()
