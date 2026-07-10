"""Feature extraction pipeline."""

from typing import Any

from protein_graph.features.base import FeatureExtractor


class FeaturePipeline:
    """Sequential feature extraction pipeline."""

    def __init__(self) -> None:
        self.extractors: list[FeatureExtractor] = []

    def add(self, extractor: FeatureExtractor) -> None:
        """Register a feature extractor."""
        self.extractors.append(extractor)

    def extract(self, protein: dict[str, Any]) -> dict[str, Any]:
        """Run all registered feature extractors."""

        features: dict[str, Any] = {}

        for extractor in self.extractors:
            features.update(
                extractor.extract(protein)
            )

        return features
