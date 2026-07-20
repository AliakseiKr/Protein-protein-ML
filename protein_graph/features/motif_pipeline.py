"""Pipeline for extracting motif-level features."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from typing import Any, Protocol

import pandas as pd

from protein_graph.entities.motif import MotifRecord
from protein_graph.features.base import FeatureMap
from protein_graph.features.elm_motifs import (
    MotifAnnotationFeatureExtractor,
    MotifFeatureContext,
    MotifGeometryFeatureExtractor,
)


class MotifFeatureExtractor(Protocol):
    """Interface implemented by motif-level feature extractors."""

    def extract(
        self,
        motif: MotifRecord,
        context: MotifFeatureContext,
    ) -> FeatureMap:
        """Extract features from one motif occurrence."""
        ...


class MotifFeaturePipeline:
    """Combine multiple motif-level feature extractors."""

    def __init__(
        self,
        extractors: Iterable[
            MotifFeatureExtractor
        ] | None = None,
    ) -> None:
        if extractors is None:
            extractors = (
                MotifGeometryFeatureExtractor(),
                MotifAnnotationFeatureExtractor(),
            )

        self.extractors = tuple(extractors)

    def build_from_motifs(
        self,
        motifs: Iterable[MotifRecord],
        protein_length: int,
    ) -> pd.DataFrame:
        """Build one feature-table row per motif."""
        context = MotifFeatureContext(
            protein_length=protein_length
        )

        rows = [
            self._build_row(
                motif=motif,
                context=context,
            )
            for motif in motifs
        ]

        return pd.DataFrame(rows)

    def _build_row(
        self,
        motif: MotifRecord,
        context: MotifFeatureContext,
    ) -> dict[str, Any]:
        row = self._motif_metadata(motif)

        for extractor in self.extractors:
            features = extractor.extract(
                motif=motif,
                context=context,
            )

            self._merge_features(
                row=row,
                features=features,
                extractor=extractor,
            )

        return row

    @staticmethod
    def _merge_features(
        row: dict[str, Any],
        features: Mapping[str, Any],
        extractor: MotifFeatureExtractor,
    ) -> None:
        duplicate_keys = (
            set(row) & set(features)
        )

        if duplicate_keys:
            names = ", ".join(
                sorted(duplicate_keys)
            )

            raise ValueError(
                "Duplicate motif feature keys produced by "
                f"{type(extractor).__name__}: {names}."
            )

        row.update(features)

    @staticmethod
    def _motif_metadata(
        motif: MotifRecord,
    ) -> dict[str, Any]:
        """Preserve raw motif annotations in each row."""
        return {
            "protein_id": motif.protein_id,
            "region_type": motif.region_type,
            "start": motif.start,
            "end": motif.end,
            "sequence": motif.sequence,
            "motif_id": motif.motif_id,
            "motif_class": motif.motif_class,
            "pattern": motif.pattern,
            "recognized_by": list(
                motif.recognized_by
            ),
            "interaction_role": (
                motif.interaction_role
            ),
            "evidence_source": (
                motif.evidence_source
            ),
            "evidence_level": (
                motif.evidence_level
            ),
            "confidence": motif.confidence,
        }
