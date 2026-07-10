"""Feature extractor registry."""

from protein_graph.features.base import FeatureExtractor


FEATURE_EXTRACTORS: list[type[FeatureExtractor]] = []
