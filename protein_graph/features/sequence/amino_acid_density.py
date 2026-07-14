"""Local amino acid density features."""

from typing import Any

from protein_graph.features.base import FeatureExtractor
from protein_graph.features.sequence.constants import LOCAL_GROUPS, WINDOW_SIZES


DENSITY_THRESHOLD = 0.75


class LocalDensityExtractor(FeatureExtractor):
    """Count non-overlapping high-density sequence windows."""

    name = "local_density"

    def extract(self, protein: dict[str, Any]) -> dict[str, int]:
        sequence = protein.get("sequence", "").upper()

        if not sequence:
            raise ValueError("Protein sequence is empty.")

        features: dict[str, int] = {}

        for group_name, residues in LOCAL_GROUPS.items():
            for window_size in WINDOW_SIZES:
                feature_name = (
                    f"high_{group_name}_density_window_count_{window_size}"
                )

                features[feature_name] = self._high_density_window_count(
                    sequence=sequence,
                    residues=residues,
                    window_size=window_size,
                    threshold=DENSITY_THRESHOLD,
                )

        return features

    @staticmethod
    def _high_density_window_count(
        sequence: str,
        residues: set[str],
        window_size: int,
        threshold: float,
    ) -> int:
        if len(sequence) < window_size:
            density = LocalDensityExtractor._fraction(sequence, residues)
            return int(density >= threshold)

        count = 0
        start = 0
        last_start = len(sequence) - window_size

        while start <= last_start:
            window = sequence[start : start + window_size]
            density = LocalDensityExtractor._fraction(window, residues)

            if density >= threshold:
                count += 1
                start += window_size
            else:
                start += 1

        return count

    @staticmethod
    def _fraction(sequence: str, residues: set[str]) -> float:
        residue_count = 0

        for aa in sequence:
            if aa in residues:
                residue_count += 1

        return residue_count / len(sequence)
