"""Amino acid composition features."""

from typing import Any

from protein_graph.features.base import FeatureExtractor
from protein_graph.features.sequence.constants import AMINO_ACIDS


class AminoAcidCompositionExtractor(FeatureExtractor):
    """Compute amino acid composition."""

    name = "amino_acid_composition"

    def extract(self, protein: dict[str, Any]) -> dict[str, float | int]:
        sequence = protein.get("sequence", "").upper()

        if not sequence:
            raise ValueError("Protein sequence is empty.")

        length = len(sequence)

        features: dict[str, float | int] = {
            "protein_length": length,
        }

        for aa in AMINO_ACIDS:
            features[f"aa_fraction_{aa}"] = sequence.count(aa) / length

        return features
