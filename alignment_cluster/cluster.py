"""UPGMA-style average-linkage hierarchical clustering.

Cluster distance is the arithmetic mean of the original pairwise distances
over every cross-cluster leaf pair, so clusters of different sizes are
weighted by size automatically. Exact ``Fraction`` arithmetic is used so
merge decisions are never driven by rounded numbers.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction

from .distance import PairComparison


@dataclass
class ClusterNode:
    """A node of the rooted clustering tree."""

    min_id: str
    height: Fraction
    members: tuple[str, ...]
    identifier: str | None = None
    children: list["ClusterNode"] = field(default_factory=list)
    branch_length: Fraction = Fraction(0)

    @property
    def is_leaf(self) -> bool:
        return self.identifier is not None


@dataclass(frozen=True)
class MergeRecord:
    """One merge event, in the order it happened."""

    step: int
    left_members: tuple[str, ...]
    right_members: tuple[str, ...]
    distance: Fraction
    height: Fraction


def _cluster_distance(
    members_a: tuple[str, ...],
    members_b: tuple[str, ...],
    pair_map: dict[tuple[str, str], PairComparison],
) -> Fraction:
    """Mean of all original leaf-pair distances across two clusters."""

    total = Fraction(0)
    for id_a in members_a:
        for id_b in members_b:
            key = (id_a, id_b) if id_a < id_b else (id_b, id_a)
            total += pair_map[key].distance
    return total / (len(members_a) * len(members_b))


def upgma(
    identifiers: list[str],
    pair_map: dict[tuple[str, str], PairComparison],
) -> tuple[ClusterNode, list[MergeRecord]]:
    """Run average-linkage UPGMA and return (root, merge records)."""

    clusters: list[ClusterNode] = [
        ClusterNode(
            min_id=identifier,
            height=Fraction(0),
            members=(identifier,),
            identifier=identifier,
        )
        for identifier in identifiers
    ]
    merges: list[MergeRecord] = []

    for step in range(1, len(identifiers)):
        best: tuple[Fraction, str, str, int, int] | None = None
        for i, cluster_a in enumerate(clusters):
            for j in range(i + 1, len(clusters)):
                cluster_b = clusters[j]
                distance = _cluster_distance(
                    cluster_a.members, cluster_b.members, pair_map
                )
                candidate = (distance, cluster_a.min_id, cluster_b.min_id, i, j)
                if best is None or candidate < best:
                    best = candidate

        assert best is not None
        distance, _, _, i, j = best
        child_b = clusters.pop(j)
        child_a = clusters.pop(i)
        # Children are ordered by their smallest identifier.
        if child_a.min_id > child_b.min_id:
            child_a, child_b = child_b, child_a

        height = distance / 2
        child_a.branch_length = height - child_a.height
        child_b.branch_length = height - child_b.height
        merged = ClusterNode(
            min_id=child_a.min_id,
            height=height,
            members=tuple(sorted((*child_a.members, *child_b.members))),
            children=[child_a, child_b],
        )
        clusters.append(merged)
        merges.append(
            MergeRecord(
                step=step,
                left_members=child_a.members,
                right_members=child_b.members,
                distance=distance,
                height=height,
            )
        )

    root = clusters[0]
    root.branch_length = Fraction(0)
    return root, merges
