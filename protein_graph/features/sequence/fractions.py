"""Global residue group fraction features."""

from typing import Any

from protein_graph.features.base import FeatureExtractor
from protein_graph.features.sequence.constants import GLOBAL_GROUPS


class ResidueFractionExtractor(FeatureExtractor):
    """Compute global residue group fractions."""

    name = "residue_fractions"

    def extract(self, protein: dict[str, Any]) -> dict[str, float]:
        sequence = protein.get("sequence", "").upper()

        if not sequence:
            raise ValueError("Protein sequence is empty.")

        features = {}

        for group_name, residues in GLOBAL_GROUPS.items():
            features[f"{group_name}_fraction"] = self._fraction(sequence, residues)

        return features

    @staticmethod
    def _fraction(sequence: str, residues: set[str]) -> float:
        count = 0

        for aa in sequence:
            if aa in residues:
                count += 1

        return count / len(sequence)
