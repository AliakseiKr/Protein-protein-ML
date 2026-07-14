from typing import Any

from protein_graph.features.base import FeatureExtractor
from protein_graph.features.sequence.constants import (
    CHARGED,
    HYDROPHOBIC,
    PERIODS,
)


MIN_PERIODIC_RUN = 3

PERIODIC_RUN_BINS = (
    (3, 5),
    (6, 8),
    (9, 12),
    (13, 16),
    (17, 19),
    (20, None),
)


class PeriodicityExtractor(FeatureExtractor):
    """Detect and summarize periodic residue patterns."""

    name = "periodicity"

    def extract(self, protein: dict[str, Any]) -> dict[str, int]:
        sequence = protein.get("sequence", "").upper()

        if not sequence:
            raise ValueError("Protein sequence is empty.")

        features: dict[str, int] = {}

        groups = {
            "hydrophobic": HYDROPHOBIC,
            "charged": CHARGED,
        }

        for group_name, residues in groups.items():
            for step in PERIODS:
                run_lengths = self._periodic_run_lengths(
                    sequence=sequence,
                    residues=residues,
                    step=step,
                    minimum_run=MIN_PERIODIC_RUN,
                )

                prefix = f"{group_name}_periodic_run_step_{step}"

                # Total number of valid periodic regions.
                features[f"{prefix}_count"] = len(run_lengths)

                # Longest valid periodic run.
                features[f"longest_{prefix}"] = (
                    max(run_lengths) if run_lengths else 0
                )

                # Number of runs in every length category.
                for lower, upper in PERIODIC_RUN_BINS:
                    bin_name = self._bin_name(
                        lower=lower,
                        upper=upper,
                    )

                    features[f"{prefix}_count_{bin_name}"] = (
                        self._count_runs_in_bin(
                            run_lengths=run_lengths,
                            lower=lower,
                            upper=upper,
                        )
                    )

        return features

    @staticmethod
    def _periodic_run_lengths(
        sequence: str,
        residues: set[str],
        step: int,
        minimum_run: int = MIN_PERIODIC_RUN,
    ) -> list[int]:
        """Return the lengths of all qualifying periodic runs."""
        if step <= 0:
            raise ValueError("step must be a positive integer.")

        if minimum_run <= 0:
            raise ValueError(
                "minimum_run must be a positive integer."
            )

        run_lengths: list[int] = []

        # Each start value represents one positional phase.
        #
        # For step=3:
        # start=0 checks positions 0, 3, 6, 9, ...
        # start=1 checks positions 1, 4, 7, 10, ...
        # start=2 checks positions 2, 5, 8, 11, ...
        for start in range(step):
            current_run = 0

            for position in range(
                start,
                len(sequence),
                step,
            ):
                if sequence[position] in residues:
                    current_run += 1
                else:
                    if current_run >= minimum_run:
                        run_lengths.append(current_run)

                    current_run = 0

            # Save a run that reaches the end of this phase.
            if current_run >= minimum_run:
                run_lengths.append(current_run)

        return run_lengths

    @staticmethod
    def _count_runs_in_bin(
        run_lengths: list[int],
        lower: int,
        upper: int | None,
    ) -> int:
        """Count periodic runs inside one non-overlapping bin."""
        if upper is None:
            return sum(
                run_length >= lower
                for run_length in run_lengths
            )

        return sum(
            lower <= run_length <= upper
            for run_length in run_lengths
        )

    @staticmethod
    def _bin_name(
        lower: int,
        upper: int | None,
    ) -> str:
        """Create a stable feature-name suffix for a run bin."""
        if upper is None:
            return f"ge_{lower}"

        return f"{lower}_to_{upper}"
