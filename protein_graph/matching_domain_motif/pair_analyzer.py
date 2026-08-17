from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, Sequence

from protein_graph.matching_domain_motif.domain_motif import (
    DetectedDomain,
    DomainMotifCandidate,
    DomainMotifRegistry,
    MotifScanner,
)


class DomainDetector(Protocol):
    """Интерфейс для Pfam/HMMER/InterPro или собственного детектора."""

    def detect(
        self,
        sequence: str,
    ) -> Sequence[DetectedDomain]:
        ...


@dataclass(frozen=True, slots=True)
class ProteinInput:
    protein_id: str
    sequence: str

    def __post_init__(self) -> None:
        if not self.protein_id:
            raise ValueError("protein_id must not be empty")

        if not self.sequence:
            raise ValueError("protein sequence must not be empty")


@dataclass(frozen=True, slots=True)
class PairAnalysisResult:
    protein_a_id: str
    protein_b_id: str

    domains_a: tuple[DetectedDomain, ...]
    domains_b: tuple[DetectedDomain, ...]

    candidates_a_to_b: tuple[DomainMotifCandidate, ...]
    candidates_b_to_a: tuple[DomainMotifCandidate, ...]

    skipped_domains_a: tuple[DetectedDomain, ...]
    skipped_domains_b: tuple[DetectedDomain, ...]

    @property
    def all_candidates(self) -> tuple[DomainMotifCandidate, ...]:
        return self.candidates_a_to_b + self.candidates_b_to_a

    @property
    def num_candidates(self) -> int:
        return len(self.all_candidates)

    @property
    def has_bidirectional_candidates(self) -> bool:
        return bool(
            self.candidates_a_to_b
            and self.candidates_b_to_a
        )

    @property
    def max_combined_score(self) -> float:
        if not self.all_candidates:
            return 0.0

        return max(
            candidate.combined_score
            for candidate in self.all_candidates
        )


class ProteinPairAnalyzer:
    """
    1. Находит домены в обоих белках.
    2. Для доменов A запускает соответствующие сканеры по B.
    3. Для доменов B запускает соответствующие сканеры по A.
    4. Не запускает сканер для неизвестных или слишком слабых доменов.
    """

    def __init__(
        self,
        *,
        domain_detector: DomainDetector,
        motif_scanner: MotifScanner,
        motif_registry: DomainMotifRegistry,
        domain_confidence_threshold: float = 0.5,
    ) -> None:
        if not 0.0 <= domain_confidence_threshold <= 1.0:
            raise ValueError(
                "domain_confidence_threshold must be between 0 and 1"
            )

        self._domain_detector = domain_detector
        self._motif_scanner = motif_scanner
        self._motif_registry = motif_registry
        self._domain_confidence_threshold = (
            domain_confidence_threshold
        )

    def analyze(
        self,
        protein_a: ProteinInput,
        protein_b: ProteinInput,
    ) -> PairAnalysisResult:
        domains_a = tuple(
            self._domain_detector.detect(protein_a.sequence)
        )
        domains_b = tuple(
            self._domain_detector.detect(protein_b.sequence)
        )

        candidates_a_to_b, skipped_a = self._scan_direction(
            source=protein_a,
            target=protein_b,
            domains=domains_a,
        )

        candidates_b_to_a, skipped_b = self._scan_direction(
            source=protein_b,
            target=protein_a,
            domains=domains_b,
        )

        return PairAnalysisResult(
            protein_a_id=protein_a.protein_id,
            protein_b_id=protein_b.protein_id,
            domains_a=domains_a,
            domains_b=domains_b,
            candidates_a_to_b=tuple(candidates_a_to_b),
            candidates_b_to_a=tuple(candidates_b_to_a),
            skipped_domains_a=tuple(skipped_a),
            skipped_domains_b=tuple(skipped_b),
        )

    def _scan_direction(
        self,
        *,
        source: ProteinInput,
        target: ProteinInput,
        domains: Sequence[DetectedDomain],
    ) -> tuple[
        list[DomainMotifCandidate],
        list[DetectedDomain],
    ]:
        candidates: list[DomainMotifCandidate] = []
        skipped_domains: list[DetectedDomain] = []

        for domain in domains:
            should_scan = (
                domain.confidence
                >= self._domain_confidence_threshold
                and self._motif_registry.supports_domain(
                    domain.domain_class
                )
            )

            if not should_scan:
                skipped_domains.append(domain)
                continue

            domain_candidates = self._motif_scanner.scan(
                domain=domain,
                source_protein_id=source.protein_id,
                target_protein_id=target.protein_id,
                target_sequence=target.sequence,
            )

            candidates.extend(domain_candidates)

        return candidates, skipped_domains
