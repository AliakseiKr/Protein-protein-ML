"""Builder for motif definition entities."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from typing import Any

from protein_graph.entities.motif_definition import (
    MotifDefinition,
)


class MotifDefinitionBuilder:
    """Build MotifDefinition entities from canonical records.

    The builder is independent of any particular external database.
    Data sources such as ELM or PROSITE must first normalize their
    records into the canonical field names expected here.
    """

    def build(
        self,
        record: Mapping[str, Any],
    ) -> MotifDefinition:
        """Build one motif definition.

        Parameters
        ----------
        record:
            Canonical motif-definition record.

        Returns
        -------
        MotifDefinition
            Validated immutable motif definition.
        """
        motif_id = self._required_string(
            record,
            "motif_id",
        )
        motif_class = self._required_string(
            record,
            "motif_class",
        )
        pattern = self._required_string(
            record,
            "pattern",
        )

        return MotifDefinition(
            motif_id=motif_id,
            motif_class=motif_class,
            pattern=pattern,
            recognized_by=self._string_tuple(
                record.get("recognized_by")
            ),
            interaction_role=self._optional_string(
                record.get("interaction_role"),
                field_name="interaction_role",
            ),
            source=(
                self._optional_string(
                    record.get("source"),
                    field_name="source",
                )
                or "unknown"
            ),
            evidence_level=self._optional_string(
                record.get("evidence_level"),
                field_name="evidence_level",
            ),
        )

    def build_many(
        self,
        records: Iterable[Mapping[str, Any]],
    ) -> list[MotifDefinition]:
        """Build multiple unique motif definitions."""
        definitions: list[MotifDefinition] = []
        seen_ids: set[str] = set()

        for record in records:
            definition = self.build(record)

            if definition.motif_id in seen_ids:
                raise ValueError(
                    "Duplicate motif definition ID: "
                    f"{definition.motif_id!r}."
                )

            seen_ids.add(definition.motif_id)
            definitions.append(definition)

        return definitions

    @staticmethod
    def _required_string(
        record: Mapping[str, Any],
        field_name: str,
    ) -> str:
        """Read a required non-empty string."""
        value = record.get(field_name)

        if not isinstance(value, str) or not value.strip():
            raise ValueError(
                "Motif definition record requires a "
                f"non-empty {field_name!r}."
            )

        return value.strip()

    @staticmethod
    def _optional_string(
        value: Any,
        *,
        field_name: str,
    ) -> str | None:
        """Normalize an optional string value."""
        if value is None:
            return None

        if not isinstance(value, str):
            raise TypeError(
                f"{field_name} must be a string or None."
            )

        cleaned = value.strip()
        return cleaned or None

    @staticmethod
    def _string_tuple(
        value: Any,
    ) -> tuple[str, ...]:
        """Normalize a collection of recognizing-domain names."""
        if value is None:
            return ()

        if isinstance(value, str):
            raw_values = value.split(",")
        elif isinstance(value, (list, tuple, set)):
            raw_values = value
        else:
            raise TypeError(
                "recognized_by must be a string or an "
                "iterable of strings."
            )

        normalized: list[str] = []

        for item in raw_values:
            if not isinstance(item, str):
                raise TypeError(
                    "recognized_by values must be strings."
                )

            cleaned = item.strip()

            if cleaned and cleaned not in normalized:
                normalized.append(cleaned)

        return tuple(normalized)

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}()"
