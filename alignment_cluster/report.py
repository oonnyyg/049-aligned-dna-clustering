"""Write reviewable CSV and JSON reports."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from .clustering import TreeNode
from .distances import PairComparison
from .fasta import SequenceRecord

MATRIX_FILENAME = "distance_matrix.csv"
PAIRS_FILENAME = "pairwise_comparisons.csv"
TREE_FILENAME = "tree.json"


def _format_number(value: float) -> str:
    """Format a distance without losing exact values needed for review."""

    if value == 0.0:
        return "0"
    return repr(value)


def write_distance_matrix(
    path: Path,
    identifiers: list[str],
    distances: dict[tuple[str, str], float],
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([""] + identifiers)
        for row_id in identifiers:
            writer.writerow(
                [row_id]
                + [_format_number(distances[(row_id, col_id)]) for col_id in identifiers]
            )


def write_pairwise_comparisons(
    path: Path, comparisons: list[PairComparison]
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "identifier_a",
                "identifier_b",
                "comparable_sites",
                "differing_sites",
                "distance",
            ]
        )
        for comparison in comparisons:
            writer.writerow(
                [
                    comparison.identifier_a,
                    comparison.identifier_b,
                    comparison.comparable_sites,
                    comparison.differing_sites,
                    _format_number(comparison.distance),
                ]
            )


def build_tree_payload(
    records: list[SequenceRecord],
    identifiers: list[str],
    alignment_length: int,
    merge_records: list[dict],
    root: TreeNode,
) -> dict:
    descriptions = {
        record.identifier: record.description for record in records
    }
    return {
        "schema_version": 1,
        "method": "UPGMA",
        "distance_definition": (
            "differing sites divided by sites where both sequences are in "
            "A/C/G/T; N and gap sites are excluded per pair"
        ),
        "height_definition": "half of the merging inter-cluster distance",
        "tie_breaking": (
            "lexicographic order of the two clusters' minimum identifiers"
        ),
        "alignment_length": alignment_length,
        "leaf_order": identifiers,
        "leaves": [
            {
                "id": identifier,
                "description": descriptions.get(identifier, ""),
                "height": 0.0,
            }
            for identifier in identifiers
        ],
        "merges": merge_records,
        "tree": root.to_dict(),
    }


def write_tree_json(path: Path, payload: dict) -> None:
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
