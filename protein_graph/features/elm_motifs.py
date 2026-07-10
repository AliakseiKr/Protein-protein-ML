"""ELM motif profile extractor."""

from __future__ import annotations

import re
from typing import Any

import pandas as pd

from protein_graph.features.base import FeatureExtractor


class ELMMotifExtractor(FeatureExtractor):
    """Detect ELM motif profile in protein sequences."""

    name = "elm_motifs"

    def __init__(self, motifs: pd.DataFrame) -> None:
        self.motifs = motifs

    def extract(self, protein: dict[str, Any]) -> dict[str, Any]:
        sequence = protein.get("sequence", "").upper()

        if not sequence:
            raise ValueError("Protein sequence is empty.")

        matched_motif_ids: list[str] = []
        total_hits = 0

        for _, row in self.motifs.iterrows():
            elm_id = str(row["elm_id"])
            pattern = str(row["regex"])

            try:
                matches = re.findall(pattern, sequence)
            except re.error:
                continue

            if matches:
                matched_motif_ids.append(elm_id)
                total_hits += len(matches)

        return {
            "elm_total_motif_hits": total_hits,
            "elm_unique_motif_classes": len(set(matched_motif_ids)),
            "elm_matched_motif_ids": sorted(set(matched_motif_ids)),
        }
