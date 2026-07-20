"""Protein-level ELM motif profile extraction."""

from __future__ import annotations

import re
from collections import Counter
from typing import Any

import pandas as pd

from protein_graph.features.base import FeatureExtractor


class ELMMotifExtractor(FeatureExtractor):
    """Extract an aggregated ELM motif profile from a protein sequence.

    The extractor scans the complete protein sequence using regular
    expressions from an ELM motif table.

    It reports:

    - total number of motif occurrences;
    - number of unique motif classes;
    - IDs of detected ELM motifs;
    - classes of detected ELM motifs;
    - number of occurrences for each motif class.
    """

    name = "elm_motifs"

    REQUIRED_COLUMNS = {
        "elm_id",
        "regex",
    }

    def __init__(self, motifs: pd.DataFrame) -> None:
        """Initialize the extractor from an ELM motif table."""
        missing_columns = (
            self.REQUIRED_COLUMNS - set(motifs.columns)
        )

        if missing_columns:
            missing = ", ".join(
                sorted(missing_columns)
            )

            raise ValueError(
                "ELM motif table is missing required columns: "
                f"{missing}."
            )

        self.motifs = motifs.copy()

    def extract(
        self,
        protein: dict[str, Any],
    ) -> dict[str, Any]:
        """Extract the aggregated ELM motif profile."""
        sequence = str(
            protein.get("sequence", "")
        ).strip().upper()

        if not sequence:
            raise ValueError(
                "Protein sequence is empty."
            )

        total_hits = 0
        matched_motif_ids: set[str] = set()
        matched_motif_classes: set[str] = set()
        motif_class_counts: Counter[str] = Counter()

        for _, row in self.motifs.iterrows():
            elm_id = self._normalize_value(
                row["elm_id"]
            )
            pattern = self._normalize_value(
                row["regex"]
            )

            if not elm_id or not pattern:
                continue

            motif_class = self._get_motif_class(
                row=row,
                elm_id=elm_id,
            )

            try:
                matches = list(
                    re.finditer(
                        pattern,
                        sequence,
                    )
                )
            except re.error:
                continue

            if not matches:
                continue

            number_of_hits = len(matches)

            total_hits += number_of_hits
            matched_motif_ids.add(elm_id)
            matched_motif_classes.add(
                motif_class
            )
            motif_class_counts[
                motif_class
            ] += number_of_hits

        return {
            "elm_total_motif_hits": total_hits,
            "elm_unique_motif_classes": len(
                matched_motif_classes
            ),
            "elm_matched_motif_ids": sorted(
                matched_motif_ids
            ),
            "elm_matched_motif_classes": sorted(
                matched_motif_classes
            ),
            "elm_motif_class_counts": dict(
                sorted(
                    motif_class_counts.items()
                )
            ),
        }

    @staticmethod
    def _normalize_value(
        value: Any,
    ) -> str:
        """Convert a table value to a normalized string."""
        if pd.isna(value):
            return ""

        return str(value).strip()

    def _get_motif_class(
        self,
        row: pd.Series,
        elm_id: str,
    ) -> str:
        """Return an explicit motif class or derive it from ELM ID."""
        if "motif_class" in row.index:
            motif_class = self._normalize_value(
                row["motif_class"]
            )

            if motif_class:
                return motif_class

        return self._derive_motif_class(
            elm_id
        )

    @staticmethod
    def _derive_motif_class(
        elm_id: str,
    ) -> str:
        """Derive a broader ELM class from its identifier.

        Examples:

        LIG_SH3_1          -> LIG_SH3
        MOD_CK2_1          -> MOD_CK2
        DOC_MAPK_gen_1     -> DOC_MAPK
        TRG_NLS_MonoExtN_4 -> TRG_NLS

        This fallback is used only when the input table does not contain
        a populated ``motif_class`` column.
        """
        parts = elm_id.split("_")

        if len(parts) < 2:
            return elm_id

        return "_".join(
            parts[:2]
        )
