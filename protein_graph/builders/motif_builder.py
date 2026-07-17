"""Build motif entities from protein sequences."""

from __future__ import annotations

from collections.abc import Iterable
import re
from typing import Any

from protein_graph.entities.motif import MotifRecord
from protein_graph.entities.motif_definition import (
    MotifDefinition,
)


MatchResult = tuple[int, int, str]


class MotifBuilder:
    """Find motif occurrences in protein sequences."""

    def __init__(
        self,
        definitions: Iterable[MotifDefinition],
        *,
        allow_overlapping: bool = True,
    ) -> None:
        self.definitions = tuple(definitions)
        self.allow_overlapping = allow_overlapping

        motif_ids = [
            definition.motif_id
            for definition in self.definitions
        ]

        if len(motif_ids) != len(set(motif_ids)):
            raise ValueError(
                "Motif definitions contain duplicate motif_id values."
            )

        self._compiled_patterns = {
            definition.motif_id: self._compile_pattern(
                definition
            )
            for definition in self.definitions
        }

    def build(
        self,
        protein: object,
    ) -> list[MotifRecord]:
        """Find all known motif occurrences in one protein."""
        protein_id = self._required_attribute(
            protein,
            "protein_id",
        )

        sequence = self._required_attribute(
            protein,
            "sequence",
        ).upper()

        motifs: list[MotifRecord] = []

        for definition in self.definitions:
            pattern = self._compiled_patterns[
                definition.motif_id
            ]

            for start, end, matched_sequence in (
                self._find_matches(
                    pattern,
                    sequence,
                )
            ):
                motifs.append(
                    MotifRecord(
                        protein_id=protein_id,
                        region_type="motif",
                        start=start + 1,
                        end=end,
                        sequence=matched_sequence,
                        description=(
                            definition.motif_class
                        ),
                        source=definition.source,
                        metadata={
                            "pattern": (
                                definition.pattern
                            ),
                        },
                        motif_id=definition.motif_id,
                        motif_class=(
                            definition.motif_class
                        ),
                        pattern=definition.pattern,
                        matched_sequence=(
                            matched_sequence
                        ),
                        recognized_by=(
                            definition.recognized_by
                        ),
                        interaction_role=(
                            definition.interaction_role
                        ),
                        evidence_source=(
                            definition.source
                        ),
                        evidence_level=(
                            definition.evidence_level
                        ),
                    )
                )

        return motifs

    def build_many(
        self,
        protein: dict[str, Any],
    ) -> list[MotifRecord]:
        """Build motifs from a canonical protein dictionary."""
        return self.build(protein)

    def _find_matches(
        self,
        pattern: re.Pattern[str],
        sequence: str,
    ) -> Iterable[MatchResult]:
        if not self.allow_overlapping:
            for match in pattern.finditer(sequence):
                yield (
                    match.start(),
                    match.end(),
                    match.group(0),
                )

            return

        lookahead = re.compile(
            f"(?=({pattern.pattern}))",
            flags=pattern.flags,
        )

        for match in lookahead.finditer(sequence):
            matched_sequence = match.group(1)

            if not matched_sequence:
                continue

            yield (
                match.start(1),
                match.end(1),
                matched_sequence,
            )

    @staticmethod
    def _compile_pattern(
        definition: MotifDefinition,
    ) -> re.Pattern[str]:
        try:
            return re.compile(
                definition.pattern,
                flags=re.IGNORECASE,
            )
        except re.error as error:
            raise ValueError(
                f"Invalid regex for motif "
                f"{definition.motif_id!r}: "
                f"{definition.pattern!r}"
            ) from error

    @staticmethod
    def _required_attribute(
        entity: object,
        name: str,
    ) -> str:
        if isinstance(entity, dict):
            value = entity.get(name)
        else:
            value = getattr(entity, name, None)

        if not isinstance(value, str) or not value.strip():
            raise ValueError(
                f"Protein must provide a non-empty "
                f"{name!r}."
            )

        return value.strip()

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}("
            f"definitions={len(self.definitions)}, "
            f"allow_overlapping={self.allow_overlapping}"
            f")"
        )
