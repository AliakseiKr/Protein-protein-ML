from __future__ import annotations

from collections.abc import Mapping, Sequence

from protein_graph.matching_domain_motif.domain_motif import (
    DetectedDomain,
)


class StaticDomainDetector:
    def __init__(
        self,
        domains_by_sequence: Mapping[
            str,
            Sequence[DetectedDomain],
        ],
    ) -> None:
        self._domains_by_sequence = {
            sequence.strip().upper(): tuple(domains)
            for sequence, domains in domains_by_sequence.items()
        }

    def detect(
        self,
        sequence: str,
    ) -> tuple[DetectedDomain, ...]:
        normalized_sequence = sequence.strip().upper()

        return self._domains_by_sequence.get(
            normalized_sequence,
            (),
        )
