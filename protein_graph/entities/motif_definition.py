"""Motif definition entity."""

from __future__ import annotations

import re
from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class MotifDefinition:
    """A motif class definition obtained from a knowledge source.

    Parameters
    ----------
    motif_id:
        Stable identifier of the motif definition.

    motif_class:
        Human-readable motif class or functional name.

    pattern:
        Regular expression used to find motif occurrences in a
        protein sequence.

    recognized_by:
        Domain classes or other interaction partners known to
        recognize the motif.

    interaction_role:
        Optional biological role associated with the motif.

    source:
        Name of the knowledge source from which the definition
        originated.

    evidence_level:
        Optional source-specific evidence classification.
    """

    motif_id: str
    motif_class: str
    pattern: str

    recognized_by: tuple[str, ...] = ()

    interaction_role: str | None = None
    source: str = "unknown"
    evidence_level: str | None = None

    def __post_init__(self) -> None:
        """Validate and normalize the motif definition."""
        motif_id = self.motif_id.strip()
        motif_class = self.motif_class.strip()
        pattern = self.pattern.strip()
        source = self.source.strip()

        if not motif_id:
            raise ValueError("motif_id cannot be empty.")

        if not motif_class:
            raise ValueError("motif_class cannot be empty.")

        if not pattern:
            raise ValueError("pattern cannot be empty.")

        if not source:
            raise ValueError("source cannot be empty.")

        try:
            re.compile(pattern)
        except re.error as error:
            raise ValueError(
                f"Invalid motif regex {pattern!r}."
            ) from error

        recognized_by = self._normalize_recognized_by(
            self.recognized_by
        )

        interaction_role = self._normalize_optional_string(
            self.interaction_role
        )
        evidence_level = self._normalize_optional_string(
            self.evidence_level
        )

        object.__setattr__(self, "motif_id", motif_id)
        object.__setattr__(self, "motif_class", motif_class)
        object.__setattr__(self, "pattern", pattern)
        object.__setattr__(self, "source", source)
        object.__setattr__(
            self,
            "recognized_by",
            recognized_by,
        )
        object.__setattr__(
            self,
            "interaction_role",
            interaction_role,
        )
        object.__setattr__(
            self,
            "evidence_level",
            evidence_level,
        )

    @staticmethod
    def _normalize_recognized_by(
        values: tuple[str, ...],
    ) -> tuple[str, ...]:
        """Normalize recognizing-domain names."""
        if isinstance(values, str):
            raise TypeError(
                "recognized_by must be a tuple of strings, "
                "not a single string."
            )

        normalized: list[str] = []

        for value in values:
            if not isinstance(value, str):
                raise TypeError(
                    "recognized_by values must be strings."
                )

            cleaned = value.strip()

            if cleaned and cleaned not in normalized:
                normalized.append(cleaned)

        return tuple(normalized)

    @staticmethod
    def _normalize_optional_string(
        value: str | None,
    ) -> str | None:
        """Strip an optional string and convert blanks to None."""
        if value is None:
            return None

        if not isinstance(value, str):
            raise TypeError(
                "Optional motif definition fields must be "
                "strings or None."
            )

        cleaned = value.strip()
        return cleaned or None
