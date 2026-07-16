"""Build domain entities from canonical protein records."""

from __future__ import annotations

from typing import Any

from protein_graph.entities import DomainRecord


class DomainBuilder:
    """Create DomainRecord objects from protein domain annotations."""

    def build_many(
        self,
        protein: dict[str, Any],
    ) -> list[DomainRecord]:
        """Build every annotated domain for one protein."""
        protein_id = str(protein.get("protein_id", "")).strip()
        sequence = str(protein.get("sequence", "")).upper()

        if not protein_id:
            raise ValueError("Protein record is missing protein_id.")

        if not sequence:
            raise ValueError(
                f"Protein {protein_id} is missing its sequence."
            )

        annotations = protein.get("domain_annotations", [])

        if annotations is None:
            return []

        if not isinstance(annotations, list):
            raise ValueError(
                "domain_annotations must be a list."
            )

        domains = []

        for annotation in annotations:
            domain = self._build_domain(
                protein_id=protein_id,
                protein_sequence=sequence,
                annotation=annotation,
            )

            domains.append(domain)

        return domains

    def _build_domain(
        self,
        protein_id: str,
        protein_sequence: str,
        annotation: dict[str, Any],
    ) -> DomainRecord:
        """Build one domain entity."""
        start = self._required_position(
            annotation,
            name="start",
        )
        end = self._required_position(
            annotation,
            name="end",
        )

        if end > len(protein_sequence):
            raise ValueError(
                f"Domain end {end} exceeds protein "
                f"{protein_id} length {len(protein_sequence)}."
            )

        domain_sequence = protein_sequence[start - 1 : end]

        return DomainRecord(
            protein_id=protein_id,
            region_type="domain",
            start=start,
            end=end,
            sequence=domain_sequence,
            description=annotation.get("description"),
            source=annotation.get("source", "UniProt"),
            metadata=annotation.get("metadata", {}),
            domain_id=annotation.get("domain_id"),
            pfam_id=annotation.get("pfam_id"),
            interpro_id=annotation.get("interpro_id"),
            smart_id=annotation.get("smart_id"),
            cdd_id=annotation.get("cdd_id"),
            family=annotation.get("family"),
            superfamily=annotation.get("superfamily"),
            clan=annotation.get("clan"),
            functional_roles=self._as_tuple(
                annotation.get("functional_roles")
            ),
            binding_types=self._as_tuple(
                annotation.get("binding_types")
            ),
            recognized_motifs=self._as_tuple(
                annotation.get("recognized_motifs")
            ),
            structural_class=annotation.get(
                "structural_class"
            ),
            structural_architecture=annotation.get(
                "structural_architecture"
            ),
            structural_fold=annotation.get(
                "structural_fold"
            ),
            known_ligands=self._as_tuple(
                annotation.get("known_ligands")
            ),
            known_ions=self._as_tuple(
                annotation.get("known_ions")
            ),
            known_cofactors=self._as_tuple(
                annotation.get("known_cofactors")
            ),
        )

    @staticmethod
    def _required_position(
        annotation: dict[str, Any],
        name: str,
    ) -> int:
        """Return a required positive integer coordinate."""
        value = annotation.get(name)

        if not isinstance(value, int):
            raise ValueError(
                f"Domain annotation {name} must be an integer."
            )

        if value < 1:
            raise ValueError(
                f"Domain annotation {name} must be at least 1."
            )

        return value

    @staticmethod
    def _as_tuple(value: Any) -> tuple[str, ...]:
        """Normalize annotation values into a string tuple."""
        if value is None:
            return ()

        if isinstance(value, str):
            value = [value]

        if not isinstance(value, (list, tuple, set)):
            raise ValueError(
                "Domain annotation collection must be "
                "a string, list, tuple, or set."
            )

        return tuple(
            str(item).strip()
            for item in value
            if str(item).strip()
        )

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}()"
