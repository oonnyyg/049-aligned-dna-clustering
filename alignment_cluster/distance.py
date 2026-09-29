"""Pairwise comparable-site distances for aligned DNA sequences."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from .fasta import SequenceRecord


@dataclass(frozen=True)
class PairComparison:
    """One unordered sequence pair.

    ``comparable`` counts columns where both sequences have A/C/G/T.
    ``differences`` counts comparable columns where the bases differ.
    ``distance`` = differences / comparable (exact Fraction).
    """

    id_a: str
    id_b: str
    comparable: int
    differences: int
    distance: Fraction


class DistanceError(ValueError):
    """Raised when a sequence pair has no comparable sites."""


def compare_pair(a: SequenceRecord, b: SequenceRecord) -> tuple[int, int]:
    """Return (comparable_sites, differing_sites) for one aligned pair.

    Columns where either sequence has N or a gap are excluded for this
    pair only; they remain available to every other pair.
    """

    comparable = 0
    differences = 0
    for base_a, base_b in zip(a.sequence, b.sequence, strict=True):
        if base_a in "ACGT" and base_b in "ACGT":
            comparable += 1
            if base_a != base_b:
                differences += 1
    return comparable, differences


def build_distance_matrix(
    records: list[SequenceRecord],
) -> tuple[list[str], dict[tuple[str, str], PairComparison]]:
    """Build every unordered pair comparison and a symmetric distance map.

    Returns the lexicographically sorted identifiers and a dictionary keyed
    by sorted identifier pairs. Raises ``DistanceError`` naming every pair
    that has zero comparable sites; the whole analysis fails in that case.
    """

    identifiers = sorted(r.identifier for r in records)
    by_id = {r.identifier: r for r in records}
    comparisons: dict[tuple[str, str], PairComparison] = {}
    incomparable: list[str] = []

    for i, id_a in enumerate(identifiers):
        for id_b in identifiers[i + 1 :]:
            comparable, differences = compare_pair(by_id[id_a], by_id[id_b])
            if comparable == 0:
                incomparable.append(f"{id_a}/{id_b}")
                continue
            comparisons[(id_a, id_b)] = PairComparison(
                id_a=id_a,
                id_b=id_b,
                comparable=comparable,
                differences=differences,
                distance=Fraction(differences, comparable),
            )

    if incomparable:
        raise DistanceError(
            "no comparable sites for pair(s): " + ", ".join(incomparable)
        )

    return identifiers, comparisons
