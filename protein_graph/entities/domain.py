"""Protein domain entities."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from protein_graph.entities.region import SequenceRegion


@dataclass(frozen=True)
class DomainRecord(SequenceRegion):
    """A protein domain with identifiers and biological annotations."""

    domain_id: str | None = None

    pfam_id: str | None = None
    interpro_id: str | None = None
    smart_id: str | None = None
    cdd_id: str | None = None

    family: str | None = None
    superfamily: str | None = None
    clan: str | None = None

    functional_roles: tuple[str, ...] = field(default_factory=tuple)
    binding_types: tuple[str, ...] = field(default_factory=tuple)
    recognized_motifs: tuple[str, ...] = field(default_factory=tuple)

    structural_class: str | None = None
    structural_architecture: str | None = None
    structural_fold: str | None = None

    known_ligands: tuple[str, ...] = field(default_factory=tuple)
    known_ions: tuple[str, ...] = field(default_factory=tuple)
    known_cofactors: tuple[str, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        super().__post_init__()

        if self.region_type != "domain":
            raise ValueError(
                "DomainRecord region_type must be 'domain'."
            )

        identifiers = (
            self.domain_id,
            self.pfam_id,
            self.interpro_id,
            self.smart_id,
            self.cdd_id,
        )

        if not any(identifiers):
            raise ValueError(
                "DomainRecord requires at least one domain identifier."
            )

    def to_dict(self) -> dict[str, Any]:
        """Convert the domain to a serializable dictionary."""
        data = super().to_dict()

        data.update(
            {
                "domain_id": self.domain_id,
                "pfam_id": self.pfam_id,
                "interpro_id": self.interpro_id,
                "smart_id": self.smart_id,
                "cdd_id": self.cdd_id,
                "family": self.family,
                "superfamily": self.superfamily,
                "clan": self.clan,
                "functional_roles": list(self.functional_roles),
                "binding_types": list(self.binding_types),
                "recognized_motifs": list(
                    self.recognized_motifs
                ),
                "structural_class": self.structural_class,
                "structural_architecture": (
                    self.structural_architecture
                ),
                "structural_fold": self.structural_fold,
                "known_ligands": list(self.known_ligands),
                "known_ions": list(self.known_ions),
                "known_cofactors": list(self.known_cofactors),
            }
        )

        return data
