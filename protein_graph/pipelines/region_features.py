"""Build sequence-derived features for arbitrary protein regions."""

from __future__ import annotations

from collections.abc import Iterable
from typing import Any

import pandas as pd

from protein_graph.entities import SequenceRegion
from protein_graph.features.sequence import (
    AminoAcidCompositionExtractor,
    ConsecutiveRunExtractor,
    LocalDensityExtractor,
    PeriodicityExtractor,
    ResidueFractionExtractor,
)


class RegionFeaturePipeline:
    """Calculate sequence features for any SequenceRegion."""

    def __init__(self) -> None:
        self.extractors = [
            AminoAcidCompositionExtractor(),
            ResidueFractionExtractor(),
            LocalDensityExtractor(),
            ConsecutiveRunExtractor(),
            PeriodicityExtractor(),
        ]

    def build(self, region: SequenceRegion) -> dict[str, Any]:
        """Build features for one sequence region."""
        region_input = {
            "protein_id": region.protein_id,
            "sequence": region.sequence,
        }

        features: dict[str, Any] = {
            "protein_id": region.protein_id,
            "region_type": region.region_type,
            "region_start": region.start,
            "region_end": region.end,
            "region_length": region.length,
            "region_description": region.description,
            "region_source": region.source,
        }

        for extractor in self.extractors:
            extracted = extractor.extract(region_input)

            # The composition extractor uses protein_length because it was
            # originally written for complete proteins.
            extracted.pop("protein_length", None)

            duplicated = set(features) & set(extracted)

            if duplicated:
                raise ValueError(
                    f"Region feature collision: {sorted(duplicated)}"
                )

            features.update(extracted)

        return features

    def build_many(
        self,
        regions: Iterable[SequenceRegion],
    ) -> pd.DataFrame:
        """Build one feature row for every sequence region."""
        return pd.DataFrame(
            self.build(region)
            for region in regions
        )

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}("
            f"extractors={len(self.extractors)})"
        )
