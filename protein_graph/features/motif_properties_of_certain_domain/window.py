"""0-based, half-open Python-index window helper.

Distinct from protein_graph.entities.region.SequenceRegion (and its
subclasses DomainRecord / MotifRecord), which store 1-based inclusive
biological coordinates plus identity (protein_id, region_type, sequence).
WindowInterval carries no sequence and no identity - it exists only to
express local Python-slice windows and interval overlap for the motif
context/combined feature extractors in this package.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class WindowInterval:
    """A 0-based, half-open Python-index window.

    Callers holding biological coordinates (1-based inclusive) should
    convert explicitly using `from_biological`.
    """

    start: int
    end: int

    def __post_init__(self) -> None:
        if self.start < 0:
            raise ValueError("start must be >= 0")

        if self.end <= self.start:
            raise ValueError(
                "end must be greater than start"
            )

    def overlaps(
        self,
        other: WindowInterval,
    ) -> bool:
        """Return whether this window overlaps another window."""
        return (
            self.start < other.end
            and other.start < self.end
        )

    @classmethod
    def from_biological(
        cls,
        start: int,
        end: int,
    ) -> WindowInterval:
        """Convert 1-based inclusive biological coordinates
        (entities.region.SequenceRegion convention) to a 0-based
        half-open Python window.
        """
        if start < 1:
            raise ValueError(
                "biological start must be >= 1"
            )

        if end < start:
            raise ValueError(
                "biological end must be >= start"
            )

        return cls(start - 1, end)
