"""Base classes for ProteinGraph feature extractors."""

from abc import ABC, abstractmethod
from typing import Any


class FeatureExtractor(ABC):
    """Abstract base class for all feature extractors."""

    name: str = "base"

    @abstractmethod
    def extract(self, protein: dict[str, Any]) -> dict[str, Any]:
        """Extract features from a single protein."""
        pass
