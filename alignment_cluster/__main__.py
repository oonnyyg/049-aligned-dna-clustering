"""Command-line entry point: ``python -m alignment_cluster``."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .cluster import upgma
from .distance import DistanceError, build_distance_matrix
from .fasta import FastaError, parse_fasta
from .report import (
    render_distance_matrix_csv,
    render_pairwise_csv,
    render_tree_json,
)


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="python -m alignment_cluster",
        description=(
            "Compute pairwise comparable-site distances and an UPGMA tree "
            "for an aligned FASTA file (standard library only)."
        ),
    )
    parser.add_argument("fasta", type=Path, help="input aligned FASTA file")
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        required=True,
        help="directory for distance_matrix.csv, pairwise_comparisons.csv "
        "and cluster_tree.json",
    )
    return parser.parse_args(argv)


def run(fasta_path: Path, output_dir: Path) -> None:
    """Validate everything before writing any report file."""

    try:
        text = fasta_path.read_text(encoding="utf-8")
    except OSError as exc:
        raise FastaError(f"cannot read {fasta_path}: {exc}") from exc

    records = parse_fasta(text)
    identifiers, pair_map = build_distance_matrix(records)
    root, merges = upgma(identifiers, pair_map)

    outputs = {
        "distance_matrix.csv": render_distance_matrix_csv(identifiers, pair_map),
        "pairwise_comparisons.csv": render_pairwise_csv(pair_map),
        "cluster_tree.json": render_tree_json(
            identifiers, pair_map, root, merges
        ),
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    for name, content in outputs.items():
        (output_dir / name).write_text(content, encoding="utf-8", newline="")


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    try:
        run(args.fasta, args.output_dir)
    except (FastaError, DistanceError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    print(
        "wrote distance_matrix.csv, pairwise_comparisons.csv and "
        f"cluster_tree.json in {args.output_dir}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
