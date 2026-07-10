"""Combined sequence property extractor."""

from typing import Any

from protein_graph.features.base import FeatureExtractor
from protein_graph.features.sequence import (
    AminoAcidCompositionExtractor,
    ConsecutiveRunExtractor,
    LocalDensityExtractor,
    PeriodicityExtractor,
    ResidueFractionExtractor,
)


class SequencePropertyExtractor(FeatureExtractor):
    """Run all sequence-derived feature extractors."""

    name = "sequence_properties"

    def __init__(self) -> None:
        self.extractors = [
            AminoAcidCompositionExtractor(),
            ResidueFractionExtractor(),
            LocalDensityExtractor(),
            ConsecutiveRunExtractor(),
            PeriodicityExtractor(),
        ]

    def extract(self, protein: dict[str, Any]) -> dict[str, Any]:
        features: dict[str, Any] = {}

        for extractor in self.extractors:
            features.update(extractor.extract(protein))

        return features
