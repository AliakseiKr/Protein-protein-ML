"""Sequence region entities."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


@dataclass(frozen=True)
class SequenceRegion:
    """A bounded region inside a protein sequence.

    Coordinates use biological one-based inclusive indexing.
    """

    protein_id: str
    region_type: str
    start: int
    end: int
    sequence: str
    description: str | None = None
    source: str | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.protein_id.strip():
            raise ValueError("protein_id must not be empty.")

        if not self.region_type.strip():
            raise ValueError("region_type must not be empty.")

        if self.start < 1:
            raise ValueError("Region start must be at least 1.")

        if self.end < self.start:
            raise ValueError(
                "Region end must be greater than or equal to start."
            )

        expected_length = self.end - self.start + 1

        if len(self.sequence) != expected_length:
            raise ValueError(
                "Region sequence length does not match its coordinates: "
                f"expected {expected_length}, received {len(self.sequence)}."
            )

    @property
    def length(self) -> int:
        """Return the region length."""
        return self.end - self.start + 1

    def relative_position(self, protein_length: int) -> float:
        """Return the region midpoint as a fraction of protein length."""
        if protein_length <= 0:
            raise ValueError("protein_length must be greater than zero.")

        midpoint = (self.start + self.end) / 2

        return midpoint / protein_length

    def position_category(self, protein_length: int) -> str:
        """Classify the region as N-terminal, middle, or C-terminal."""
        relative_position = self.relative_position(protein_length)

        if relative_position < 1 / 3:
            return "n_terminal"

        if relative_position < 2 / 3:
            return "middle"

        return "c_terminal"

    def to_dict(self) -> dict[str, Any]:
        """Convert the region to a serializable dictionary."""
        return {
            "protein_id": self.protein_id,
            "region_type": self.region_type,
            "start": self.start,
            "end": self.end,
            "sequence": self.sequence,
            "length": self.length,
            "description": self.description,
            "source": self.source,
            "metadata": self.metadata,
        }
