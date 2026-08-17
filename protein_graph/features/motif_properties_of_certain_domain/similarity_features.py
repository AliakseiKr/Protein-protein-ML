from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import Collection

from protein_graph.features.motif_properties_of_certain_domain.sequence_features import (
    AMINO_ACIDS,
    AMINO_ACID_GROUPS,
)


@dataclass(frozen=True, slots=True)
class MotifSimilarityConfig:
    critical_positions: tuple[int, ...] = ()
    kmer_sizes: tuple[int, ...] = (2, 3)
    require_equal_length: bool = True
    use_blosum62: bool = True


class MotifSimilarityFeatureExtractor:
    """
    Сравнивает кандидатный участок с одним известным мотивом.
    """

    def __init__(
        self,
        config: MotifSimilarityConfig | None = None,
    ) -> None:
        self._config = config or MotifSimilarityConfig()

    def extract(
        self,
        *,
        reference_motif: str,
        candidate_sequence: str,
    ) -> dict[str, float | int]:
        reference = self._normalize(reference_motif)
        candidate = self._normalize(candidate_sequence)

        if (
            self._config.require_equal_length
            and len(reference) != len(candidate)
        ):
            raise ValueError(
                "Sequences must have equal length: "
                f"{len(reference)} != {len(candidate)}"
            )

        aligned_length = min(
            len(reference),
            len(candidate),
        )

        aligned_reference = reference[:aligned_length]
        aligned_candidate = candidate[:aligned_length]

        exact_count = sum(
            reference_amino_acid == candidate_amino_acid
            for reference_amino_acid, candidate_amino_acid
            in zip(aligned_reference, aligned_candidate)
        )

        physicochemical_count = sum(
            self._share_group(
                reference_amino_acid,
                candidate_amino_acid,
            )
            for reference_amino_acid, candidate_amino_acid
            in zip(aligned_reference, aligned_candidate)
        )

        features: dict[str, float | int] = {
            "exact_match": int(reference == candidate),
            "sequence_identity": (
                exact_count / aligned_length
            ),
            "mismatch_count": (
                aligned_length - exact_count
            ),
            "physicochemical_match_fraction": (
                physicochemical_count / aligned_length
            ),
            "length_difference": abs(
                len(reference) - len(candidate)
            ),
            "longest_exact_substring_fraction": (
                self._longest_common_substring_length(
                    aligned_reference,
                    aligned_candidate,
                )
                / aligned_length
            ),
            "longest_exact_prefix_fraction": (
                self._common_prefix_length(
                    aligned_reference,
                    aligned_candidate,
                )
                / aligned_length
            ),
            "longest_exact_suffix_fraction": (
                self._common_suffix_length(
                    aligned_reference,
                    aligned_candidate,
                )
                / aligned_length
            ),
            "mismatch_cluster_count": (
                self._mismatch_cluster_count(
                    aligned_reference,
                    aligned_candidate,
                )
            ),
        }

        features.update(
            self._critical_position_features(
                aligned_reference,
                aligned_candidate,
            )
        )

        features.update(
            self._positional_features(
                aligned_reference,
                aligned_candidate,
            )
        )

        features.update(
            self._kmer_features(
                aligned_reference,
                aligned_candidate,
            )
        )

        if self._config.use_blosum62:
            features.update(
                self._blosum_features(
                    aligned_reference,
                    aligned_candidate,
                )
            )

        return features

    def _critical_position_features(
        self,
        reference: str,
        candidate: str,
    ) -> dict[str, float | int]:
        positions = tuple(
            position
            for position in self._config.critical_positions
            if 0 <= position < len(reference)
        )

        if not positions:
            return {
                "critical_position_count": 0,
                "critical_exact_fraction": 0.0,
                "critical_physicochemical_fraction": 0.0,
            }

        exact_count = sum(
            reference[position] == candidate[position]
            for position in positions
        )

        group_count = sum(
            self._share_group(
                reference[position],
                candidate[position],
            )
            for position in positions
        )

        return {
            "critical_position_count": len(positions),
            "critical_exact_fraction": (
                exact_count / len(positions)
            ),
            "critical_physicochemical_fraction": (
                group_count / len(positions)
            ),
        }

    @classmethod
    def _positional_features(
        cls,
        reference: str,
        candidate: str,
    ) -> dict[str, int]:
        features: dict[str, int] = {}

        for position, (
            reference_amino_acid,
            candidate_amino_acid,
        ) in enumerate(zip(reference, candidate)):
            features[
                f"position_{position}_exact"
            ] = int(
                reference_amino_acid
                == candidate_amino_acid
            )

            features[
                f"position_{position}_physicochemical"
            ] = int(
                cls._share_group(
                    reference_amino_acid,
                    candidate_amino_acid,
                )
            )

        return features

    def _kmer_features(
        self,
        reference: str,
        candidate: str,
    ) -> dict[str, float]:
        features: dict[str, float] = {}

        for kmer_size in self._config.kmer_sizes:
            if kmer_size <= 0:
                raise ValueError(
                    "kmer sizes must be positive"
                )

            reference_kmers = self._kmers(
                reference,
                kmer_size,
            )

            candidate_kmers = self._kmers(
                candidate,
                kmer_size,
            )

            intersection = (
                reference_kmers & candidate_kmers
            )

            union = reference_kmers | candidate_kmers

            features[
                f"kmer_{kmer_size}_jaccard"
            ] = (
                len(intersection) / len(union)
                if union
                else 1.0
            )

            features[
                f"kmer_{kmer_size}_reference_coverage"
            ] = (
                len(intersection) / len(reference_kmers)
                if reference_kmers
                else 1.0
            )

        return features

    @staticmethod
    def _blosum_features(
        reference: str,
        candidate: str,
    ) -> dict[str, float | int]:
        matrix = _load_blosum62()

        scores = [
            float(matrix[
                reference_amino_acid,
                candidate_amino_acid,
            ])
            for reference_amino_acid, candidate_amino_acid
            in zip(reference, candidate)
        ]

        return {
            "blosum62_total_score": sum(scores),
            "blosum62_mean_score": (
                sum(scores) / len(scores)
            ),
            "blosum62_positive_fraction": (
                sum(score > 0 for score in scores)
                / len(scores)
            ),
            "blosum62_negative_fraction": (
                sum(score < 0 for score in scores)
                / len(scores)
            ),
            "blosum62_minimum_position_score": min(scores),
        }

    @staticmethod
    def _share_group(
        first: str,
        second: str,
    ) -> bool:
        if first == second:
            return True

        return any(
            first in group and second in group
            for group in AMINO_ACID_GROUPS.values()
        )

    @staticmethod
    def _kmers(
        sequence: str,
        kmer_size: int,
    ) -> set[str]:
        if kmer_size > len(sequence):
            return set()

        return {
            sequence[index:index + kmer_size]
            for index in range(
                len(sequence) - kmer_size + 1
            )
        }

    @staticmethod
    def _longest_common_substring_length(
        first: str,
        second: str,
    ) -> int:
        previous = [0] * (len(second) + 1)
        maximum = 0

        for first_amino_acid in first:
            current = [0]

            for index, second_amino_acid in enumerate(
                second,
                start=1,
            ):
                if first_amino_acid == second_amino_acid:
                    value = previous[index - 1] + 1
                    maximum = max(maximum, value)
                else:
                    value = 0

                current.append(value)

            previous = current

        return maximum

    @staticmethod
    def _common_prefix_length(
        first: str,
        second: str,
    ) -> int:
        count = 0

        for first_amino_acid, second_amino_acid in zip(
            first,
            second,
        ):
            if first_amino_acid != second_amino_acid:
                break

            count += 1

        return count

    @staticmethod
    def _common_suffix_length(
        first: str,
        second: str,
    ) -> int:
        count = 0

        for first_amino_acid, second_amino_acid in zip(
            reversed(first),
            reversed(second),
        ):
            if first_amino_acid != second_amino_acid:
                break

            count += 1

        return count

    @staticmethod
    def _mismatch_cluster_count(
        reference: str,
        candidate: str,
    ) -> int:
        mismatch_clusters = 0
        inside_cluster = False

        for reference_amino_acid, candidate_amino_acid in zip(
            reference,
            candidate,
        ):
            mismatch = (
                reference_amino_acid
                != candidate_amino_acid
            )

            if mismatch and not inside_cluster:
                mismatch_clusters += 1
                inside_cluster = True
            elif not mismatch:
                inside_cluster = False

        return mismatch_clusters

    @staticmethod
    def _normalize(
        sequence: str,
    ) -> str:
        normalized = sequence.strip().upper()

        if not normalized:
            raise ValueError("Sequence must not be empty")

        invalid_symbols = set(normalized) - AMINO_ACIDS

        if invalid_symbols:
            raise ValueError(
                "Unsupported amino-acid symbols: "
                f"{sorted(invalid_symbols)}"
            )

        return normalized


@lru_cache(maxsize=1)
def _load_blosum62():
    try:
        from Bio.Align import substitution_matrices
    except ImportError as error:
        raise RuntimeError(
            "Biopython is required for BLOSUM62 features. "
            "Install it or set use_blosum62=False."
        ) from error

    return substitution_matrices.load("BLOSUM62")
