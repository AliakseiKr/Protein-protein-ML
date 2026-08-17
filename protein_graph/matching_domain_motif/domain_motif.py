from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Iterable, Mapping, Protocol, Sequence


@dataclass(frozen=True, slots=True)
class DetectedDomain:
    """Домен, обнаруженный в последовательности белка."""

    domain_class: str
    start: int
    end: int
    confidence: float
    domain_id: str | None = None

    def __post_init__(self) -> None:
        if not self.domain_class:
            raise ValueError("domain_class must not be empty")

        if self.start < 1:
            raise ValueError(
                "domain start must be >= 1 (1-based inclusive)"
            )

        if self.end < self.start:
            raise ValueError(
                "domain end must be greater than or equal to start"
            )

        if not 0.0 <= self.confidence <= 1.0:
            raise ValueError("domain confidence must be between 0 and 1")


@dataclass(frozen=True, slots=True)
class MotifDefinition:
    """
    Известный мотив, который может распознаваться определённым доменным классом.
    """

    motif_id: str
    motif_class: str
    recognized_by: tuple[str, ...]
    regex: str

    def __post_init__(self) -> None:
        if not self.motif_id:
            raise ValueError("motif_id must not be empty")

        if not self.motif_class:
            raise ValueError("motif_class must not be empty")

        if not self.recognized_by:
            raise ValueError("recognized_by must contain at least one domain class")

        if not self.regex:
            raise ValueError("regex must not be empty")

        try:
            re.compile(self.regex)
        except re.error as exc:
            raise ValueError(
                f"Invalid regex for motif {self.motif_id!r}: {self.regex!r}"
            ) from exc


@dataclass(frozen=True, slots=True)
class DomainMotifCandidate:
    """
    Потенциальное межбелковое соответствие:

    домен source-белка → мотив или motif-like участок target-белка.
    """

    source_protein_id: str
    target_protein_id: str

    domain_class: str
    domain_start: int
    domain_end: int
    domain_confidence: float
    domain_id: str | None

    target_start: int
    target_end: int
    target_sequence: str

    motif_score: float
    combined_score: float

    motif_id: str | None = None
    motif_class: str | None = None
    matched_known_regex: bool = False


class MotifScanner(Protocol):
    """Общий интерфейс для regex- и ML-сканеров."""

    def scan(
        self,
        *,
        domain: DetectedDomain,
        source_protein_id: str,
        target_protein_id: str,
        target_sequence: str,
    ) -> list[DomainMotifCandidate]:
        ...


class DomainMotifRegistry:
    """
    Индекс:

        domain_class → список мотивов, распознаваемых этим доменом.
    """

    def __init__(self, motifs: Iterable[MotifDefinition]) -> None:
        by_domain: dict[str, list[MotifDefinition]] = {}

        for motif in motifs:
            for domain_class in motif.recognized_by:
                normalized_domain = self._normalize(domain_class)
                by_domain.setdefault(normalized_domain, []).append(motif)

        self._by_domain = {
            domain_class: tuple(definitions)
            for domain_class, definitions in by_domain.items()
        }

    def motifs_for_domain(
        self,
        domain_class: str,
    ) -> tuple[MotifDefinition, ...]:
        return self._by_domain.get(self._normalize(domain_class), ())

    def supports_domain(self, domain_class: str) -> bool:
        return bool(self.motifs_for_domain(domain_class))

    @staticmethod
    def _normalize(value: str) -> str:
        return value.strip().upper()


class KnownMotifScanner:
    """
    Сканирует партнёрский белок только по мотивам,
    соответствующим найденному доменному классу.
    """

    def __init__(
        self,
        registry: DomainMotifRegistry,
    ) -> None:
        self._registry = registry

    def scan(
        self,
        *,
        domain: DetectedDomain,
        source_protein_id: str,
        target_protein_id: str,
        target_sequence: str,
    ) -> list[DomainMotifCandidate]:
        sequence = target_sequence.strip().upper()

        if not sequence:
            return []

        motif_definitions = self._registry.motifs_for_domain(
            domain.domain_class
        )

        candidates: list[DomainMotifCandidate] = []

        for motif in motif_definitions:
            pattern = re.compile(motif.regex)

            # lookahead позволяет находить перекрывающиеся совпадения
            overlapping_pattern = re.compile(
                f"(?=({pattern.pattern}))"
            )

            for match in overlapping_pattern.finditer(sequence):
                matched_sequence = match.group(1)
                # Convert Python's 0-based half-open match indices to the
                # repository's canonical biological 1-based inclusive
                # coordinates (see entities.region.SequenceRegion).
                start = match.start() + 1
                end = start + len(matched_sequence) - 1

                motif_score = 1.0
                combined_score = domain.confidence * motif_score

                candidates.append(
                    DomainMotifCandidate(
                        source_protein_id=source_protein_id,
                        target_protein_id=target_protein_id,
                        domain_class=domain.domain_class,
                        domain_start=domain.start,
                        domain_end=domain.end,
                        domain_confidence=domain.confidence,
                        domain_id=domain.domain_id,
                        target_start=start,
                        target_end=end,
                        target_sequence=matched_sequence,
                        motif_score=motif_score,
                        combined_score=combined_score,
                        motif_id=motif.motif_id,
                        motif_class=motif.motif_class,
                        matched_known_regex=True,
                    )
                )

        return candidates
