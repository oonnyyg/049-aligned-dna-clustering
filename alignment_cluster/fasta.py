"""Strict FASTA parsing for already aligned DNA sequences.

Accepted sequence characters are A, C, G, T, N and the gap character
``-``.  Sequences may span multiple lines and blank lines are allowed.
Whitespace inside sequence lines is ignored and letters are upper-cased.
"""

from __future__ import annotations

from dataclasses import dataclass

ALLOWED_CHARS = frozenset("ACGTN-")


class FastaError(ValueError):
    """Raised when the FASTA input does not satisfy the requirements."""


@dataclass(frozen=True)
class SequenceRecord:
    """One FASTA record.

    Attributes:
        identifier: First non-whitespace token after the ``>`` marker.
        description: Remaining header text (may be empty).
        sequence: Aligned sequence with whitespace removed, upper-cased.
    """

    identifier: str
    description: str
    sequence: str


def parse_fasta(text: str) -> list[SequenceRecord]:
    """Parse FASTA text into ordered records.

    Args:
        text: Complete FASTA file content.

    Returns:
        Records in file order.

    Raises:
        FastaError: On any structural or content violation.
    """

    records: list[SequenceRecord] = []
    seen: set[str] = set()
    current_id: str | None = None
    current_desc: str = ""
    current_parts: list[str] = []

    def finish() -> None:
        if current_id is None:
            return
        sequence = "".join(current_parts)
        if not sequence:
            raise FastaError(f"sequence for identifier {current_id!r} is empty")
        records.append(
            SequenceRecord(
                identifier=current_id,
                description=current_desc,
                sequence=sequence,
            )
        )

    for line_number, raw_line in enumerate(text.splitlines(), start=1):
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            finish()
            header = line[1:].strip()
            if not header:
                raise FastaError(f"line {line_number}: header has no identifier")
            parts = header.split(None, 1)
            identifier = parts[0]
            description = parts[1].strip() if len(parts) > 1 else ""
            if identifier in seen:
                raise FastaError(
                    f"line {line_number}: duplicate identifier {identifier!r}"
                )
            seen.add(identifier)
            current_id = identifier
            current_desc = description
            current_parts = []
            continue
        if current_id is None:
            raise FastaError(
                f"line {line_number}: sequence data appears before any header"
            )
        fragment = "".join(raw_line.split()).upper()
        invalid = {ch for ch in fragment if ch not in ALLOWED_CHARS}
        if invalid:
            raise FastaError(
                "line "
                f"{line_number}: invalid character(s) "
                f"{''.join(sorted(invalid))!r} for identifier {current_id!r}; "
                "only A, C, G, T, N and '-' are accepted"
            )
        current_parts.append(fragment)

    finish()

    if len(records) < 2:
        raise FastaError("at least two non-empty aligned sequences are required")

    expected_length = len(records[0].sequence)
    for record in records[1:]:
        if len(record.sequence) != expected_length:
            raise FastaError(
                "alignment length mismatch: "
                f"{records[0].identifier!r} has {expected_length} site(s) but "
                f"{record.identifier!r} has {len(record.sequence)}"
            )

    return records
