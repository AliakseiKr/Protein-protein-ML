"""Consecutive residue run features."""

from typing import Any

from protein_graph.features.base import FeatureExtractor
from protein_graph.features.sequence.constants import LOCAL_GROUPS, RUN_BINS


class ConsecutiveRunExtractor(FeatureExtractor):
    """Compute binned consecutive residue runs."""

    name = "consecutive_runs"

    def extract(self, protein: dict[str, Any]) -> dict[str, int]:
        sequence = protein.get("sequence", "").upper()

        if not sequence:
            raise ValueError("Protein sequence is empty.")

        features = {}

        for group_name, residues in LOCAL_GROUPS.items():
            max_run = self._max_run(sequence, residues)
            features[f"max_{group_name}_run"] = max_run

            for lower, upper in RUN_BINS:
                features[f"{group_name}_run_{lower}_{upper}"] = int(
                    lower <= max_run < upper
                )

            features[f"{group_name}_run_gt_40"] = int(max_run >= 40)

        return features

    @staticmethod
    def _max_run(sequence: str, residues: set[str]) -> int:
        max_run = 0
        current_run = 0

        for aa in sequence:
            if aa in residues:
                current_run += 1
                max_run = max(max_run, current_run)
            else:
                current_run = 0

        return max_run
