"""Protein motif entities."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from protein_graph.entities.region import SequenceRegion


@dataclass(frozen=True)
class MotifRecord(SequenceRegion):
    """A concrete motif occurrence inside a protein sequence."""

    motif_id: str | None = None
    motif_class: str | None = None
    pattern: str | None = None

    matched_sequence: str | None = None

    recognized_by: tuple[str, ...] = field(
        default_factory=tuple
    )
    interaction_role: str | None = None

    evidence_source: str | None = None
    evidence_level: str | None = None
    confidence: float | None = None

    def __post_init__(self) -> None:
        super().__post_init__()

        if self.region_type != "motif":
            raise ValueError(
                "MotifRecord region_type must be 'motif'."
            )

        if not self.motif_id:
            raise ValueError(
                "MotifRecord requires motif_id."
            )

        matched_sequence = (
            self.matched_sequence or self.sequence
        )

        if matched_sequence != self.sequence:
            raise ValueError(
                "matched_sequence must equal sequence."
            )

        if (
            self.confidence is not None
            and not 0.0 <= self.confidence <= 1.0
        ):
            raise ValueError(
                "confidence must be between 0.0 and 1.0."
            )

        object.__setattr__(
            self,
            "matched_sequence",
            matched_sequence,
        )

        object.__setattr__(
            self,
            "recognized_by",
            tuple(self.recognized_by),
        )

    def to_dict(self) -> dict[str, Any]:
        """Convert the motif to a serializable dictionary."""
        data = super().to_dict()

        data.update(
            {
                "motif_id": self.motif_id,
                "motif_class": self.motif_class,
                "pattern": self.pattern,
                "matched_sequence": (
                    self.matched_sequence
                ),
                "recognized_by": list(
                    self.recognized_by
                ),
                "interaction_role": (
                    self.interaction_role
                ),
                "evidence_source": (
                    self.evidence_source
                ),
                "evidence_level": (
                    self.evidence_level
                ),
                "confidence": self.confidence,
            }
        )

        return data
