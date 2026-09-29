"""Strict FASTA parser for already-aligned DNA sequences."""

from __future__ import annotations

from dataclasses import dataclass


class FastaError(ValueError):
    """Raised when the FASTA input violates the required format."""


@dataclass(frozen=True)
class SequenceRecord:
    """One FASTA record.

    ``identifier`` is the first non-whitespace token after ``>``;
    ``description`` is the remainder of the header line;
    ``sequence`` is the aligned sequence in upper case with whitespace removed.
    """

    identifier: str
    description: str
    sequence: str


_ALLOWED = frozenset("ACGTN-")


def parse_fasta(text: str) -> list[SequenceRecord]:
    """Parse FASTA ``text`` into records with strict validation.

    Rules:
    * every header line starts with ``>``;
    * the first non-whitespace token after ``>`` is a unique identifier;
    * sequence lines may span multiple lines and blank lines are allowed;
    * in-sequence whitespace is ignored and letters are upper-cased;
    * only ``A`` ``C`` ``G`` ``T`` ``N`` and ``-`` are accepted;
    * at least two non-empty, equal-length sequences are required;
    * sequence data before any header is rejected.
    """

    records: list[SequenceRecord] = []
    current_id: str | None = None
    current_desc = ""
    current_parts: list[str] = []

    def finish() -> None:
        if current_id is not None:
            records.append(
                SequenceRecord(current_id, current_desc, "".join(current_parts))
            )

    if text.startswith("\ufeff"):
        text = text[1:]
    for line_no, raw_line in enumerate(text.splitlines(), start=1):
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            finish()
            rest = line[1:]
            stripped = rest.strip()
            if not stripped:
                raise FastaError(f"line {line_no}: header has no identifier")
            parts = stripped.split(None, 1)
            identifier = parts[0]
            description = parts[1] if len(parts) > 1 else ""
            if any(r.identifier == identifier for r in records):
                raise FastaError(f"line {line_no}: duplicate identifier {identifier!r}")
            current_id = identifier
            current_desc = description
            current_parts = []
            continue
        if current_id is None:
            raise FastaError(f"line {line_no}: sequence data appears before any header")
        segment = "".join(raw_line.split()).upper()
        bad = {ch for ch in segment if ch not in _ALLOWED}
        if bad:
            raise FastaError(
                f"line {line_no}: invalid character(s) {sorted(bad)!r}; "
                "only A, C, G, T, N and '-' are allowed"
            )
        current_parts.append(segment)

    finish()

    if len(records) < 2:
        raise FastaError("at least two sequences are required")

    expected = len(records[0].sequence)
    if expected == 0:
        raise FastaError(f"sequence {records[0].identifier!r} is empty")
    for record in records[1:]:
        if len(record.sequence) == 0:
            raise FastaError(f"sequence {record.identifier!r} is empty")
        if len(record.sequence) != expected:
            raise FastaError(
                "sequence lengths differ: "
                f"{records[0].identifier!r} has {expected}, "
                f"{record.identifier!r} has {len(record.sequence)}"
            )

    return records
