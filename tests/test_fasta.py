import unittest

from alignment_cluster.fasta import FastaError, parse_fasta


class TestParseFasta(unittest.TestCase):
    def test_multiline_blank_lines_uppercase_and_description(self):
        text = ">s1 first sample\n\nacgt\n  a C \n>s2 other\nACGTAC\n"
        records = parse_fasta(text)
        self.assertEqual(len(records), 2)
        self.assertEqual(records[0].identifier, "s1")
        self.assertEqual(records[0].description, "first sample")
        self.assertEqual(records[0].sequence, "ACGTAC")
        self.assertEqual(records[1].identifier, "s2")
        self.assertEqual(records[1].description, "other")
        self.assertEqual(records[1].sequence, "ACGTAC")

    def test_bom_accepted(self):
        records = parse_fasta("\ufeff>a\nACGT\n>b\nACGT\n")
        self.assertEqual([r.identifier for r in records], ["a", "b"])

    def test_allowed_gap_and_n(self):
        records = parse_fasta(">a\nAC-GTN\n>b\nNNC-GT\n")
        self.assertEqual(records[0].sequence, "AC-GTN")

    def test_sequence_before_header_rejected(self):
        with self.assertRaisesRegex(FastaError, "before any header"):
            parse_fasta("ACGT\n>a\nACGT\n>b\nACGT\n")

    def test_duplicate_identifier_rejected(self):
        with self.assertRaisesRegex(FastaError, "duplicate identifier"):
            parse_fasta(">x\nACGT\n>x again\nACGT\n")

    def test_empty_header_rejected(self):
        with self.assertRaisesRegex(FastaError, "no identifier"):
            parse_fasta(">   \nACGT\n>b\nACGT\n")

    def test_invalid_character_rejected(self):
        with self.assertRaisesRegex(FastaError, "invalid character"):
            parse_fasta(">a\nACGX\n>b\nACGT\n")

    def test_length_mismatch_rejected(self):
        with self.assertRaisesRegex(FastaError, "lengths differ"):
            parse_fasta(">a\nACG\n>b\nACGT\n")

    def test_too_few_and_empty_rejected(self):
        with self.assertRaisesRegex(FastaError, "at least two"):
            parse_fasta(">a\nACGT\n")
        with self.assertRaisesRegex(FastaError, "empty"):
            parse_fasta(">a\n\n>b\nACGT\n")


if __name__ == "__main__":
    unittest.main()
