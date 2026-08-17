from __future__ import annotations

from dataclasses import dataclass
from math import isfinite
from statistics import fmean, pstdev
from typing import Any, Iterable, Mapping, Protocol


Numeric = int | float


class MotifSequenceFeatureExtractorProtocol(Protocol):
    """
    Интерфейс экстрактора свойств одной последовательности мотива.
    """

    def extract(self, sequence: str) -> Mapping[str, Any]:
        ...


@dataclass(frozen=True, slots=True)
class KnownMotif:
    """
    Одна экспериментально известная последовательность мотива,
    распознаваемая определённым доменом.
    """

    motif_id: str
    sequence: str
    motif_class: str | None = None
    source_protein_id: str | None = None
    start: int | None = None
    end: int | None = None
    confidence: float | None = None

    def __post_init__(self) -> None:
        motif_id = self.motif_id.strip()
        sequence = self.sequence.strip().upper()

        if not motif_id:
            raise ValueError("motif_id must not be empty")

        if not sequence:
            raise ValueError("sequence must not be empty")

        invalid_symbols = set(sequence) - set(
            "ACDEFGHIKLMNPQRSTVWY"
        )

        if invalid_symbols:
            raise ValueError(
                "Unsupported amino-acid symbols in motif "
                f"{motif_id!r}: {sorted(invalid_symbols)}"
            )

        if (self.start is None) != (self.end is None):
            raise ValueError(
                "start and end must either both be provided "
                "or both be None"
            )

        if self.start is not None and self.end is not None:
            if self.start < 0:
                raise ValueError("start must be >= 0")

            if self.end <= self.start:
                raise ValueError("end must be greater than start")

            if self.end - self.start != len(sequence):
                raise ValueError(
                    "Motif coordinates do not match sequence length"
                )

        if self.confidence is not None:
            if not 0.0 <= self.confidence <= 1.0:
                raise ValueError(
                    "confidence must be between 0 and 1"
                )

        object.__setattr__(self, "motif_id", motif_id)
        object.__setattr__(self, "sequence", sequence)


@dataclass(frozen=True, slots=True)
class KnownMotifFeatures:
    """
    Известный мотив вместе с рассчитанными признаками.
    """

    motif: KnownMotif
    features: Mapping[str, float]


@dataclass(frozen=True, slots=True)
class FeatureDistribution:
    """
    Описательная статистика признака по известным мотивам домена.

    Эти статистики являются дополнительными признаками.
    Они не заменяют отдельные известные мотивы.
    """

    mean: float
    standard_deviation: float
    minimum: float
    maximum: float


@dataclass(frozen=True, slots=True)
class DomainMotifProfile:
    """
    Набор известных мотивов, распознаваемых одним классом домена.
    """

    domain_class: str
    known_motifs: tuple[KnownMotifFeatures, ...]
    feature_distributions: Mapping[
        str,
        FeatureDistribution,
    ]

    @property
    def motif_count(self) -> int:
        return len(self.known_motifs)

    @property
    def motif_sequences(self) -> tuple[str, ...]:
        return tuple(
            known_motif.motif.sequence
            for known_motif in self.known_motifs
        )

    @property
    def motif_ids(self) -> tuple[str, ...]:
        return tuple(
            known_motif.motif.motif_id
            for known_motif in self.known_motifs
        )

    @property
    def motif_lengths(self) -> tuple[int, ...]:
        return tuple(
            sorted({
                len(known_motif.motif.sequence)
                for known_motif in self.known_motifs
            })
        )

    def motifs_of_length(
        self,
        length: int,
    ) -> tuple[KnownMotifFeatures, ...]:
        """
        Возвращает только референсные мотивы заданной длины.
        """

        if length <= 0:
            raise ValueError("length must be positive")

        return tuple(
            known_motif
            for known_motif in self.known_motifs
            if len(known_motif.motif.sequence) == length
        )


class DomainMotifProfileBuilder:
    """
    Создаёт профиль известных мотивов определённого доменного класса.

    Каждый мотив сохраняется отдельно. Дополнительно рассчитываются
    распределения числовых признаков по всему набору.
    """

    def __init__(
        self,
        feature_extractor: MotifSequenceFeatureExtractorProtocol,
    ) -> None:
        self._feature_extractor = feature_extractor

    def build(
        self,
        *,
        domain_class: str,
        motifs: Iterable[KnownMotif],
    ) -> DomainMotifProfile:
        normalized_domain_class = domain_class.strip()

        if not normalized_domain_class:
            raise ValueError("domain_class must not be empty")

        motif_items = tuple(motifs)

        if not motif_items:
            raise ValueError(
                "At least one known motif is required for "
                f"domain class {normalized_domain_class!r}"
            )

        self._validate_unique_motif_ids(motif_items)

        known_motifs = tuple(
            KnownMotifFeatures(
                motif=motif,
                features=self._extract_numeric_features(
                    motif.sequence
                ),
            )
            for motif in motif_items
        )

        feature_distributions = (
            self._calculate_feature_distributions(
                known_motifs
            )
        )

        return DomainMotifProfile(
            domain_class=normalized_domain_class,
            known_motifs=known_motifs,
            feature_distributions=feature_distributions,
        )

    def _extract_numeric_features(
        self,
        sequence: str,
    ) -> dict[str, float]:
        raw_features = self._feature_extractor.extract(sequence)

        numeric_features: dict[str, float] = {}

        for feature_name, feature_value in raw_features.items():
            if isinstance(feature_value, bool):
                numeric_value = float(feature_value)
            elif isinstance(feature_value, (int, float)):
                numeric_value = float(feature_value)
            else:
                continue

            if not isfinite(numeric_value):
                raise ValueError(
                    f"Feature {feature_name!r} returned a "
                    f"non-finite value for sequence {sequence!r}"
                )

            numeric_features[feature_name] = numeric_value

        if not numeric_features:
            raise ValueError(
                "The motif feature extractor did not return "
                "any numeric features"
            )

        return numeric_features

    @staticmethod
    def _validate_unique_motif_ids(
        motifs: tuple[KnownMotif, ...],
    ) -> None:
        motif_ids = [
            motif.motif_id
            for motif in motifs
        ]

        duplicated_ids = {
            motif_id
            for motif_id in motif_ids
            if motif_ids.count(motif_id) > 1
        }

        if duplicated_ids:
            raise ValueError(
                "Duplicate motif IDs: "
                f"{sorted(duplicated_ids)}"
            )

    @staticmethod
    def _calculate_feature_distributions(
        known_motifs: tuple[KnownMotifFeatures, ...],
    ) -> dict[str, FeatureDistribution]:
        """
        Рассчитывает статистики только для признаков,
        присутствующих у каждого известного мотива.
        """

        common_feature_names = set(
            known_motifs[0].features
        )

        for known_motif in known_motifs[1:]:
            common_feature_names.intersection_update(
                known_motif.features
            )

        distributions: dict[str, FeatureDistribution] = {}

        for feature_name in sorted(common_feature_names):
            values = [
                known_motif.features[feature_name]
                for known_motif in known_motifs
            ]

            distributions[feature_name] = FeatureDistribution(
                mean=fmean(values),
                standard_deviation=(
                    pstdev(values)
                    if len(values) > 1
                    else 0.0
                ),
                minimum=min(values),
                maximum=max(values),
            )

        return distributions
