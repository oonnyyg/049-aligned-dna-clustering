"""Pairwise comparison and distance matrix construction.

For a pair of aligned sequences only sites where *both* sequences contain
one of A, C, G or T are compared.  N and gap characters are excluded for
that pair alone; such sites may still be used by other sequence pairs.
"""

from __future__ import annotations

from dataclasses import dataclass

from .fasta import SequenceRecord

COMPARABLE = frozenset("ACGT")


class DistanceError(ValueError):
    """Raised when a sequence pair has no comparable sites."""


@dataclass(frozen=True)
class PairComparison:
    """Result for one unordered sequence pair."""

    identifier_a: str
    identifier_b: str
    comparable_sites: int
    differing_sites: int
    distance: float


def compare_pair(first: SequenceRecord, second: SequenceRecord) -> PairComparison:
    """Compare two equal-length aligned sequences.

    Raises:
        DistanceError: When the pair shares no ACGT site.
    """

    comparable = 0
    differing = 0
    for base_a, base_b in zip(first.sequence, second.sequence, strict=True):
        if base_a in COMPARABLE and base_b in COMPARABLE:
            comparable += 1
            if base_a != base_b:
                differing += 1
    if comparable == 0:
        raise DistanceError(
            "no comparable sites (both bases in A/C/G/T) between "
            f"{first.identifier!r} and {second.identifier!r}"
        )
    return PairComparison(
        identifier_a=first.identifier,
        identifier_b=second.identifier,
        comparable_sites=comparable,
        differing_sites=differing,
        distance=differing / comparable,
    )


def build_comparisons(
    records: list[SequenceRecord],
) -> tuple[list[str], list[PairComparison], dict[tuple[str, str], float]]:
    """Compare every unordered pair.

    Returns:
        A tuple of identifiers in lexicographic order, pair comparisons
        ordered by lexicographic identifier pair, and a symmetric distance
        mapping keyed by ordered identifier pairs.

    Raises:
        DistanceError: When any distinct pair has no comparable sites.
    """

    ordered = sorted(records, key=lambda record: record.identifier)
    by_id = {record.identifier: record for record in ordered}
    identifiers = [record.identifier for record in ordered]

    comparisons: list[PairComparison] = []
    distances: dict[tuple[str, str], float] = {}
    for identifier in identifiers:
        distances[(identifier, identifier)] = 0.0

    for i, id_a in enumerate(identifiers):
        for id_b in identifiers[i + 1 :]:
            comparison = compare_pair(by_id[id_a], by_id[id_b])
            comparisons.append(comparison)
            distances[(id_a, id_b)] = comparison.distance
            distances[(id_b, id_a)] = comparison.distance

    return identifiers, comparisons, distances
