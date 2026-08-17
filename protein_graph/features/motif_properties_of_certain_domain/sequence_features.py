from __future__ import annotations

from collections import Counter
from math import log2
from typing import Final


AMINO_ACIDS: Final[frozenset[str]] = frozenset(
    "ACDEFGHIKLMNPQRSTVWY"
)

AMINO_ACID_GROUPS: Final[dict[str, frozenset[str]]] = {
    "hydrophobic": frozenset("AILMFWVY"),
    "polar": frozenset("STNQ"),
    "positive": frozenset("KRH"),
    "negative": frozenset("DE"),
    "aromatic": frozenset("FWY"),
    "small": frozenset("AGSTC"),
    "flexible": frozenset("GSP"),
    "turn_promoting": frozenset("GPSND"),
    "helix_breaking": frozenset("PG"),
}


class MotifSequenceFeatureExtractor:
    """
    Извлекает свойства коротких линейных мотивов.

    Экстрактор предназначен специально для мотивов, а не для анализа
    полной белковой последовательности.
    """

    def extract(
        self,
        sequence: str,
    ) -> dict[str, float | int]:
        normalized_sequence = self._normalize(sequence)

        length = len(normalized_sequence)
        counts = Counter(normalized_sequence)

        features: dict[str, float | int] = {
            "length": length,
            "unique_amino_acid_count": len(counts),
            "sequence_entropy": self._entropy(
                normalized_sequence
            ),
            "normalized_entropy": self._normalized_entropy(
                normalized_sequence
            ),
            "dominant_amino_acid_fraction": (
                max(counts.values()) / length
            ),
            "repeated_position_fraction": (
                sum(
                    count
                    for count in counts.values()
                    if count > 1
                )
                / length
            ),
            "adjacent_identity_fraction": (
                self._adjacent_identity_fraction(
                    normalized_sequence
                )
            ),
            "amino_acid_transition_fraction": (
                self._transition_fraction(
                    normalized_sequence
                )
            ),
        }

        features.update(
            self._composition_features(
                sequence=normalized_sequence,
                counts=counts,
            )
        )

        features.update(
            self._group_fraction_features(
                normalized_sequence
            )
        )

        features.update(
            self._run_features(
                normalized_sequence
            )
        )

        features.update(
            self._spacing_features(
                normalized_sequence
            )
        )

        features.update(
            self._terminal_features(
                normalized_sequence
            )
        )

        features.update(
            self._positional_moment_features(
                normalized_sequence
            )
        )

        return features

    @staticmethod
    def _composition_features(
        *,
        sequence: str,
        counts: Counter[str],
    ) -> dict[str, float]:
        """
        Доля каждой аминокислоты в мотиве.
        """

        length = len(sequence)

        return {
            f"fraction_{amino_acid}": (
                counts[amino_acid] / length
            )
            for amino_acid in sorted(AMINO_ACIDS)
        }

    @staticmethod
    def _group_fraction_features(
        sequence: str,
    ) -> dict[str, float]:
        """
        Доли физико-химических групп аминокислот.
        """

        length = len(sequence)

        features = {
            f"{group_name}_fraction": (
                sum(
                    amino_acid in amino_acids
                    for amino_acid in sequence
                )
                / length
            )
            for group_name, amino_acids
            in AMINO_ACID_GROUPS.items()
        }

        positive_count = sum(
            amino_acid in AMINO_ACID_GROUPS["positive"]
            for amino_acid in sequence
        )

        negative_count = sum(
            amino_acid in AMINO_ACID_GROUPS["negative"]
            for amino_acid in sequence
        )

        features["charge_balance"] = (
            positive_count - negative_count
        ) / length

        features["charged_fraction"] = (
            positive_count + negative_count
        ) / length

        return features

    @classmethod
    def _run_features(
        cls,
        sequence: str,
    ) -> dict[str, float | int]:
        """
        Признаки непрерывных серий одинаковых аминокислот
        или аминокислот одной физико-химической группы.
        """

        length = len(sequence)

        maximum_identical_run = (
            cls._maximum_identical_run(sequence)
        )

        features: dict[str, float | int] = {
            "maximum_identical_run": (
                maximum_identical_run
            ),
            "maximum_identical_run_fraction": (
                maximum_identical_run / length
            ),
        }

        for group_name, amino_acids in (
            AMINO_ACID_GROUPS.items()
        ):
            maximum_run = cls._maximum_group_run(
                sequence,
                amino_acids,
            )

            features[
                f"maximum_{group_name}_run"
            ] = maximum_run

            features[
                f"maximum_{group_name}_run_fraction"
            ] = maximum_run / length

        return features

    @classmethod
    def _spacing_features(
        cls,
        sequence: str,
    ) -> dict[str, float | int]:
        """
        Расстояния между аминокислотами, часто значимыми
        для коротких линейных мотивов.
        """

        features: dict[str, float | int] = {}

        for amino_acid in (
            "P",
            "G",
            "R",
            "K",
            "S",
            "T",
        ):
            distances = cls._distances_between(
                sequence,
                frozenset({amino_acid}),
            )

            features[
                f"{amino_acid}_pair_count"
            ] = len(distances)

            features[
                f"{amino_acid}_mean_spacing"
            ] = (
                sum(distances) / len(distances)
                if distances
                else 0.0
            )

            features[
                f"{amino_acid}_minimum_spacing"
            ] = (
                min(distances)
                if distances
                else 0
            )

        # Шаблон P-X-P:
        # между двумя пролинами находится одна аминокислота.
        features["proline_x_proline_count"] = (
            cls._count_spaced_pair(
                sequence,
                first_amino_acid="P",
                second_amino_acid="P",
                residues_between=1,
            )
        )

        # Шаблон P-X-X-P:
        # между двумя пролинами находятся две аминокислоты.
        features["proline_xx_proline_count"] = (
            cls._count_spaced_pair(
                sequence,
                first_amino_acid="P",
                second_amino_acid="P",
                residues_between=2,
            )
        )

        # Соседние пролины PP.
        features["adjacent_proline_pair_count"] = (
            cls._count_spaced_pair(
                sequence,
                first_amino_acid="P",
                second_amino_acid="P",
                residues_between=0,
            )
        )

        return features

    @staticmethod
    def _terminal_features(
        sequence: str,
    ) -> dict[str, int]:
        """
        One-hot признаки первой и последней аминокислот мотива.
        """

        n_terminal_features = {
            f"n_terminal_is_{amino_acid}": int(
                sequence[0] == amino_acid
            )
            for amino_acid in sorted(AMINO_ACIDS)
        }

        c_terminal_features = {
            f"c_terminal_is_{amino_acid}": int(
                sequence[-1] == amino_acid
            )
            for amino_acid in sorted(AMINO_ACIDS)
        }

        return {
            **n_terminal_features,
            **c_terminal_features,
        }

    @staticmethod
    def _positional_moment_features(
        sequence: str,
    ) -> dict[str, float]:
        """
        Показывает, в какой части мотива расположены аминокислоты
        определённой физико-химической группы.

        Значения находятся примерно в диапазоне [-1, 1]:

        - отрицательное значение — смещение к N-концу;
        - положительное значение — смещение к C-концу;
        - около нуля — расположение около центра;
        - 0 также используется при отсутствии группы.
        """

        length = len(sequence)

        if length == 1:
            normalized_positions = [0.0]
        else:
            normalized_positions = [
                -1.0
                + 2.0 * index / (length - 1)
                for index in range(length)
            ]

        features: dict[str, float] = {}

        for group_name, amino_acids in (
            AMINO_ACID_GROUPS.items()
        ):
            positions = [
                normalized_positions[index]
                for index, amino_acid
                in enumerate(sequence)
                if amino_acid in amino_acids
            ]

            features[
                f"{group_name}_position_mean"
            ] = (
                sum(positions) / len(positions)
                if positions
                else 0.0
            )

        return features

    @staticmethod
    def _entropy(
        sequence: str,
    ) -> float:
        """
        Энтропия Шеннона аминокислотного состава.
        """

        counts = Counter(sequence)
        length = len(sequence)

        return -sum(
            (count / length)
            * log2(count / length)
            for count in counts.values()
        )

    @classmethod
    def _normalized_entropy(
        cls,
        sequence: str,
    ) -> float:
        """
        Энтропия, нормированная на количество уникальных
        аминокислот в последовательности.
        """

        unique_count = len(set(sequence))

        if unique_count <= 1:
            return 0.0

        return (
            cls._entropy(sequence)
            / log2(unique_count)
        )

    @staticmethod
    def _adjacent_identity_fraction(
        sequence: str,
    ) -> float:
        """
        Доля соседних пар, состоящих из одинаковых аминокислот.
        """

        if len(sequence) == 1:
            return 0.0

        identical_pairs = sum(
            first == second
            for first, second in zip(
                sequence,
                sequence[1:],
            )
        )

        return identical_pairs / (len(sequence) - 1)

    @staticmethod
    def _transition_fraction(
        sequence: str,
    ) -> float:
        """
        Доля соседних позиций, в которых аминокислота меняется.
        """

        if len(sequence) == 1:
            return 0.0

        transitions = sum(
            first != second
            for first, second in zip(
                sequence,
                sequence[1:],
            )
        )

        return transitions / (len(sequence) - 1)

    @staticmethod
    def _maximum_identical_run(
        sequence: str,
    ) -> int:
        """
        Максимальная длина серии одной и той же аминокислоты.
        """

        maximum = 1
        current = 1

        for previous, current_amino_acid in zip(
            sequence,
            sequence[1:],
        ):
            if previous == current_amino_acid:
                current += 1
                maximum = max(maximum, current)
            else:
                current = 1

        return maximum

    @staticmethod
    def _maximum_group_run(
        sequence: str,
        amino_acids: frozenset[str],
    ) -> int:
        """
        Максимальная непрерывная серия аминокислот,
        принадлежащих одной группе.
        """

        maximum = 0
        current = 0

        for amino_acid in sequence:
            if amino_acid in amino_acids:
                current += 1
                maximum = max(maximum, current)
            else:
                current = 0

        return maximum

    @staticmethod
    def _distances_between(
        sequence: str,
        amino_acids: frozenset[str],
    ) -> list[int]:
        """
        Расстояния между соседними вхождениями указанных
        аминокислот в последовательности.
        """

        positions = [
            index
            for index, amino_acid
            in enumerate(sequence)
            if amino_acid in amino_acids
        ]

        return [
            second - first
            for first, second in zip(
                positions,
                positions[1:],
            )
        ]

    @staticmethod
    def _count_spaced_pair(
        sequence: str,
        *,
        first_amino_acid: str,
        second_amino_acid: str,
        residues_between: int,
    ) -> int:
        """
        Считает число шаблонов вида:

            A-X...X-B

        Например:

        residues_between=0 -> AB
        residues_between=1 -> AXB
        residues_between=2 -> AXXB

        Проверяются все возможные пары, включая перекрывающиеся.
        """

        if first_amino_acid not in AMINO_ACIDS:
            raise ValueError(
                "Unsupported first_amino_acid: "
                f"{first_amino_acid!r}"
            )

        if second_amino_acid not in AMINO_ACIDS:
            raise ValueError(
                "Unsupported second_amino_acid: "
                f"{second_amino_acid!r}"
            )

        if residues_between < 0:
            raise ValueError(
                "residues_between must be >= 0"
            )

        position_difference = residues_between + 1

        if position_difference >= len(sequence):
            return 0

        return sum(
            sequence[index] == first_amino_acid
            and sequence[
                index + position_difference
            ] == second_amino_acid
            for index in range(
                len(sequence) - position_difference
            )
        )

    @staticmethod
    def _normalize(
        sequence: str,
    ) -> str:
        normalized = sequence.strip().upper()

        if not normalized:
            raise ValueError(
                "sequence must not be empty"
            )

        invalid_symbols = (
            set(normalized) - AMINO_ACIDS
        )

        if invalid_symbols:
            raise ValueError(
                "Unsupported amino-acid symbols: "
                f"{sorted(invalid_symbols)}"
            )

        return normalized
