"""Aligned DNA sequence pairwise-distance and UPGMA clustering toolkit.

Only the Python 3.12 standard library is used.
"""

from .fasta import FastaError, SequenceRecord, parse_fasta
from .distance import PairComparison, build_distance_matrix
from .cluster import ClusterNode, MergeRecord, upgma

__all__ = [
    "FastaError",
    "SequenceRecord",
    "parse_fasta",
    "PairComparison",
    "build_distance_matrix",
    "ClusterNode",
    "MergeRecord",
    "upgma",
]
