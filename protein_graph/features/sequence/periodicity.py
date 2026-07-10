"""Periodic sequence pattern features."""

from typing import Any

from protein_graph.features.base import FeatureExtractor
from protein_graph.features.sequence.constants import CHARGED, HYDROPHOBIC, PERIODS


class PeriodicityExtractor(FeatureExtractor):
    """Detect approximate periodic residue patterns."""

    name = "periodicity"

    def extract(self, protein: dict[str, Any]) -> dict[str, int]:
        sequence = protein.get("sequence", "").upper()

        if not sequence:
            raise ValueError("Protein sequence is empty.")

        features = {}

        groups = {
            "hydrophobic": HYDROPHOBIC,
            "charged": CHARGED,
        }

        for group_name, residues in groups.items():
            for step in PERIODS:
                features[f"longest_{group_name}_periodic_run_step_{step}"] = (
                    self._longest_periodic_run(sequence, residues, step)
                )

        return features

    @staticmethod
    def _longest_periodic_run(sequence: str, residues: set[str], step: int) -> int:
        longest = 0

        for start in range(step):
            current = 0

            for position in range(start, len(sequence), step):
                if sequence[position] in residues:
                    current += 1
                    longest = max(longest, current)
                else:
                    current = 0

        return longest
