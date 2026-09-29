"""Serialize the distance matrix, pairwise table and clustering tree."""

from __future__ import annotations

import csv
import io
import json
from fractions import Fraction

from .cluster import ClusterNode, MergeRecord
from .distance import PairComparison


def _fraction_dict(value: Fraction) -> dict[str, object]:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "value": float(value),
    }


def render_distance_matrix_csv(
    identifiers: list[str],
    pair_map: dict[tuple[str, str], PairComparison],
) -> str:
    """Render the symmetric, zero-diagonal distance matrix as CSV text."""

    buffer = io.StringIO(newline="")
    writer = csv.writer(buffer, lineterminator="\n")
    writer.writerow([""] + identifiers)
    for id_a in identifiers:
        row = [id_a]
        for id_b in identifiers:
            if id_a == id_b:
                row.append("0")
            else:
                key = (id_a, id_b) if id_a < id_b else (id_b, id_a)
                row.append(repr(float(pair_map[key].distance)))
        writer.writerow(row)
    return buffer.getvalue()


def render_pairwise_csv(
    pair_map: dict[tuple[str, str], PairComparison],
) -> str:
    """Render one row per unordered pair, ordered by identifier pair."""

    buffer = io.StringIO(newline="")
    writer = csv.writer(buffer, lineterminator="\n")
    writer.writerow(
        ["id_a", "id_b", "comparable_sites", "differences", "distance"]
    )
    for key in sorted(pair_map):
        comparison = pair_map[key]
        writer.writerow(
            [
                comparison.id_a,
                comparison.id_b,
                comparison.comparable,
                comparison.differences,
                repr(float(comparison.distance)),
            ]
        )
    return buffer.getvalue()


def _node_to_dict(node: ClusterNode) -> dict[str, object]:
    if node.is_leaf:
        return {
            "type": "leaf",
            "identifier": node.identifier,
            "height": _fraction_dict(node.height),
            "branch_length": _fraction_dict(node.branch_length),
        }
    return {
        "type": "internal",
        "min_identifier": node.min_id,
        "height": _fraction_dict(node.height),
        "branch_length": _fraction_dict(node.branch_length),
        "members": list(node.members),
        "children": [_node_to_dict(child) for child in node.children],
    }


def _merge_to_dict(merge: MergeRecord) -> dict[str, object]:
    return {
        "step": merge.step,
        "left": {"members": list(merge.left_members),
                 "min_identifier": merge.left_members[0]},
        "right": {"members": list(merge.right_members),
                  "min_identifier": merge.right_members[0]},
        "cluster_distance": _fraction_dict(merge.distance),
        "node_height": _fraction_dict(merge.height),
    }


def render_tree_json(
    identifiers: list[str],
    pair_map: dict[tuple[str, str], PairComparison],
    root: ClusterNode,
    merges: list[MergeRecord],
) -> str:
    """Render merge history plus the final tree as JSON text."""

    payload = {
        "method": "UPGMA (unweighted pair group method with arithmetic mean)",
        "cluster_distance": (
            "arithmetic mean of original pairwise distances over all "
            "cross-cluster leaf pairs"
        ),
        "tie_breaking": (
            "lexicographic ordered pair of each cluster's minimum identifier"
        ),
        "identifiers": identifiers,
        "pairwise_distances": {
            f"{id_a}|{id_b}": _fraction_dict(comparison.distance)
            for (id_a, id_b), comparison in sorted(pair_map.items())
        },
        "merges": [_merge_to_dict(merge) for merge in merges],
        "tree": _node_to_dict(root),
    }
    return json.dumps(payload, indent=2, ensure_ascii=False) + "\n"
