"""Aligned DNA sequence local-similarity analysis tools.

Only the Python 3.12 standard library is used.  The package converts an
already aligned FASTA file into a reviewable distance matrix, pairwise
comparison table and UPGMA hierarchical clustering tree.
"""

from .fasta import FastaError, SequenceRecord, parse_fasta

__all__ = ["FastaError", "SequenceRecord", "parse_fasta"]
