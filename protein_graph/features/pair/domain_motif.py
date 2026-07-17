from collections.abc import Iterable
from dataclasses import dataclass

from protein_graph.entities.domain import DomainRecord
from protein_graph.entities.motif import MotifRecord


@dataclass(frozen=True, slots=True)
class DomainMotifMatch:
    """A compatible domain-motif combination between two proteins."""

    domain_protein_id: str
    motif_protein_id: str

    domain_family: str
    motif_id: str
    motif_sequence: str

    motif_start: int
    motif_end: int

    direction: str


class DomainMotifCompatibilityExtractor:
    """Build pairwise features from compatible domains and motifs."""

    def build(
        self,
        *,
        domains_a: Iterable[DomainRecord],
        motifs_a: Iterable[MotifRecord],
        domains_b: Iterable[DomainRecord],
        motifs_b: Iterable[MotifRecord],
    ) -> dict[str, object]:
        domains_a = tuple(domains_a)
        motifs_a = tuple(motifs_a)
        domains_b = tuple(domains_b)
        motifs_b = tuple(motifs_b)

        matches_ab = self._find_matches(
            domains=domains_a,
            motifs=motifs_b,
            direction="a_domain_to_b_motif",
        )

        matches_ba = self._find_matches(
            domains=domains_b,
            motifs=motifs_a,
            direction="b_domain_to_a_motif",
        )

        matches = matches_ab + matches_ba

        return {
            "has_compatible_domain_motif_pair": int(bool(matches)),
            "compatible_domain_motif_count": len(matches),
            "a_domain_to_b_motif_count": len(matches_ab),
            "b_domain_to_a_motif_count": len(matches_ba),
            "compatible_domain_families": tuple(
                sorted(
                    {
                        match.domain_family
                        for match in matches
                    }
                )
            ),
            "compatible_motif_ids": tuple(
                sorted(
                    {
                        match.motif_id
                        for match in matches
                    }
                )
            ),
            "domain_motif_matches": tuple(matches),
        }

    def _find_matches(
        self,
        *,
        domains: Iterable[DomainRecord],
        motifs: Iterable[MotifRecord],
        direction: str,
    ) -> list[DomainMotifMatch]:
        matches: list[DomainMotifMatch] = []

        for domain in domains:
            domain_labels = self._domain_labels(domain)

            if not domain_labels:
                continue

            for motif in motifs:
                recognized_by = {
                    value.casefold()
                    for value in motif.recognized_by
                }

                if domain_labels.isdisjoint(recognized_by):
                    continue

                domain_family = (
                    domain.family
                    or domain.description
                    or domain.pfam_id
                    or domain.interpro_id
                    or "unknown"
                )

                matches.append(
                    DomainMotifMatch(
                        domain_protein_id=domain.protein_id,
                        motif_protein_id=motif.protein_id,
                        domain_family=domain_family,
                        motif_id=motif.motif_id or "unknown",
                        motif_sequence=motif.sequence,
                        motif_start=motif.start,
                        motif_end=motif.end,
                        direction=direction,
                    )
                )

        return matches

    @staticmethod
    def _domain_labels(
        domain: DomainRecord,
    ) -> set[str]:
        values = (
            domain.family,
            domain.description,
            domain.pfam_id,
            domain.interpro_id,
            domain.smart_id,
            domain.cdd_id,
        )

        return {
            value.casefold()
            for value in values
            if isinstance(value, str) and value
        }

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}()"
