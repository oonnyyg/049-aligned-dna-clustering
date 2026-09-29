"""Command line entry point: ``python -m alignment_cluster``."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .clustering import upgma
from .distances import DistanceError, build_comparisons
from .fasta import FastaError, parse_fasta
from .report import (
    MATRIX_FILENAME,
    PAIRS_FILENAME,
    TREE_FILENAME,
    build_tree_payload,
    write_distance_matrix,
    write_pairwise_comparisons,
    write_tree_json,
)


def create_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m alignment_cluster",
        description=(
            "Build a pairwise distance matrix and UPGMA tree from an "
            "already aligned DNA FASTA file."
        ),
    )
    parser.add_argument("fasta", type=Path, help="input FASTA file")
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        required=True,
        help=(
            "directory for distance_matrix.csv, pairwise_comparisons.csv "
            "and tree.json"
        ),
    )
    return parser


def run(fasta_path: Path, output_dir: Path) -> None:
    """Validate, compute and then write all reports.

    Reports are written only after every computation succeeds, so invalid
    input never creates or overwrites report files.
    """

    try:
        text = fasta_path.read_text(encoding="utf-8")
    except OSError as exc:
        raise FastaError(f"cannot read {fasta_path}: {exc}") from exc

    records = parse_fasta(text)
    identifiers, comparisons, distances = build_comparisons(records)
    merge_records, root = upgma(identifiers, distances)
    payload = build_tree_payload(
        records=records,
        identifiers=identifiers,
        alignment_length=len(records[0].sequence),
        merge_records=merge_records,
        root=root,
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    write_distance_matrix(output_dir / MATRIX_FILENAME, identifiers, distances)
    write_pairwise_comparisons(output_dir / PAIRS_FILENAME, comparisons)
    write_tree_json(output_dir / TREE_FILENAME, payload)


def main(argv: list[str] | None = None) -> int:
    parser = create_parser()
    args = parser.parse_args(argv)
    try:
        run(args.fasta, args.output_dir)
    except (FastaError, DistanceError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
