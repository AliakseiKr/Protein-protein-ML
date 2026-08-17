from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

from protein_graph.features.motif_properties_of_certain_domain.sequence_features import (
    AMINO_ACIDS,
    AMINO_ACID_GROUPS,
    MotifSequenceFeatureExtractor,
)


@dataclass(frozen=True, slots=True)
class SequenceRegion:
    start: int
    end: int

    def __post_init__(self) -> None:
        if self.start < 0:
            raise ValueError("start must be >= 0")

        if self.end <= self.start:
            raise ValueError(
                "end must be greater than start"
            )

    def overlaps(
        self,
        other: SequenceRegion,
    ) -> bool:
        return (
            self.start < other.end
            and other.start < self.end
        )


@dataclass(frozen=True, slots=True)
class MotifContextConfig:
    flank_size: int = 5
    terminal_region_fraction: float = 0.1

    def __post_init__(self) -> None:
        if self.flank_size < 0:
            raise ValueError(
                "flank_size must be >= 0"
            )

        if not 0.0 <= self.terminal_region_fraction <= 0.5:
            raise ValueError(
                "terminal_region_fraction must be "
                "between 0 and 0.5"
            )


class MotifContextFeatureExtractor:
    """
    Извлекает признаки положения и окружения кандидатного мотива.
    """

    def __init__(
        self,
        config: MotifContextConfig | None = None,
    ) -> None:
        self._config = config or MotifContextConfig()
        self._sequence_extractor = (
            MotifSequenceFeatureExtractor()
        )

    def extract(
        self,
        *,
        protein_sequence: str,
        motif_start: int,
        motif_end: int,
        domain_regions: Iterable[SequenceRegion] = (),
        transmembrane_regions: Iterable[
            SequenceRegion
        ] = (),
    ) -> dict[str, float | int]:
        protein = self._normalize(protein_sequence)

        if motif_start < 0:
            raise ValueError(
                "motif_start must be >= 0"
            )

        if motif_end <= motif_start:
            raise ValueError(
                "motif_end must be greater than motif_start"
            )

        if motif_end > len(protein):
            raise ValueError(
                "motif_end exceeds protein length"
            )

        motif_region = SequenceRegion(
            motif_start,
            motif_end,
        )

        flank_start = max(
            0,
            motif_start - self._config.flank_size,
        )

        flank_end = min(
            len(protein),
            motif_end + self._config.flank_size,
        )

        left_flank = protein[flank_start:motif_start]
        right_flank = protein[motif_end:flank_end]
        local_context = protein[flank_start:flank_end]
        motif_sequence = protein[motif_start:motif_end]

        protein_length = len(protein)
        motif_length = len(motif_sequence)

        terminal_limit = max(
            1,
            int(
                protein_length
                * self._config.terminal_region_fraction
            ),
        )

        features: dict[str, float | int] = {
            "motif_start": motif_start,
            "motif_end": motif_end,
            "motif_length": motif_length,
            "relative_start": (
                motif_start / protein_length
            ),
            "relative_center": (
                (
                    motif_start
                    + motif_end
                )
                / 2
                / protein_length
            ),
            "distance_to_n_terminus": motif_start,
            "distance_to_c_terminus": (
                protein_length - motif_end
            ),
            "relative_distance_to_n": (
                motif_start / protein_length
            ),
            "relative_distance_to_c": (
                (protein_length - motif_end)
                / protein_length
            ),
            "inside_n_terminal_region": int(
                motif_start < terminal_limit
            ),
            "inside_c_terminal_region": int(
                protein_length - motif_end
                < terminal_limit
            ),
            "left_flank_length": len(left_flank),
            "right_flank_length": len(right_flank),
            "local_context_length": len(local_context),
            "overlaps_known_domain": int(
                any(
                    motif_region.overlaps(region)
                    for region in domain_regions
                )
            ),
            "overlaps_transmembrane_region": int(
                any(
                    motif_region.overlaps(region)
                    for region
                    in transmembrane_regions
                )
            ),
        }

        features.update(
            self._local_context_features(
                motif_sequence=motif_sequence,
                left_flank=left_flank,
                right_flank=right_flank,
                local_context=local_context,
            )
        )

        return features

    def _local_context_features(
        self,
        *,
        motif_sequence: str,
        left_flank: str,
        right_flank: str,
        local_context: str,
    ) -> dict[str, float]:
        features: dict[str, float] = {}

        motif_features = (
            self._sequence_extractor.extract(
                motif_sequence
            )
        )

        context_features = (
            self._sequence_extractor.extract(
                local_context
            )
        )

        for feature_name in (
            "hydrophobic_fraction",
            "polar_fraction",
            "positive_fraction",
            "negative_fraction",
            "flexible_fraction",
            "helix_breaking_fraction",
            "sequence_entropy",
            "charge_balance",
        ):
            motif_value = float(
                motif_features[feature_name]
            )

            context_value = float(
                context_features[feature_name]
            )

            features[
                f"context_{feature_name}"
            ] = context_value

            features[
                f"motif_context_delta_{feature_name}"
            ] = motif_value - context_value

        features.update(
            self._simple_flank_features(
                prefix="left_flank",
                sequence=left_flank,
            )
        )

        features.update(
            self._simple_flank_features(
                prefix="right_flank",
                sequence=right_flank,
            )
        )

        return features

    @staticmethod
    def _simple_flank_features(
        *,
        prefix: str,
        sequence: str,
    ) -> dict[str, float]:
        if not sequence:
            return {
                f"{prefix}_hydrophobic_fraction": 0.0,
                f"{prefix}_charged_fraction": 0.0,
                f"{prefix}_flexible_fraction": 0.0,
            }

        length = len(sequence)

        return {
            f"{prefix}_hydrophobic_fraction": (
                sum(
                    amino_acid
                    in AMINO_ACID_GROUPS["hydrophobic"]
                    for amino_acid in sequence
                )
                / length
            ),
            f"{prefix}_charged_fraction": (
                sum(
                    amino_acid
                    in (
                        AMINO_ACID_GROUPS["positive"]
                        | AMINO_ACID_GROUPS["negative"]
                    )
                    for amino_acid in sequence
                )
                / length
            ),
            f"{prefix}_flexible_fraction": (
                sum(
                    amino_acid
                    in AMINO_ACID_GROUPS["flexible"]
                    for amino_acid in sequence
                )
                / length
            ),
        }

    @staticmethod
    def _normalize(
        sequence: str,
    ) -> str:
        normalized = sequence.strip().upper()

        if not normalized:
            raise ValueError(
                "protein_sequence must not be empty"
            )

        invalid_symbols = set(normalized) - AMINO_ACIDS

        if invalid_symbols:
            raise ValueError(
                "Unsupported amino-acid symbols: "
                f"{sorted(invalid_symbols)}"
            )

        return normalized
