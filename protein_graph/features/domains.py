"""Feature extraction for protein domains."""

from dataclasses import dataclass

from protein_graph.entities.domain import DomainRecord
from protein_graph.features.base import FeatureMap


@dataclass(frozen=True)
class DomainFeatureContext:
    """Additional information required for domain feature extraction."""

    protein_length: int

    def __post_init__(self) -> None:
        if self.protein_length <= 0:
            raise ValueError(
                "protein_length must be greater than zero."
            )


class DomainGeometryFeatureExtractor:
    """Extract geometric features from a protein domain."""

    def extract(
        self,
        domain: DomainRecord,
        context: DomainFeatureContext,
    ) -> FeatureMap:
        protein_length = context.protein_length

        if domain.end > protein_length:
            raise ValueError(
                "Domain coordinates exceed protein length."
            )

        position_category = domain.position_category(
            protein_length
        )

        return {
            "domain_length": domain.length,
            "domain_relative_position": (
                domain.relative_position(protein_length)
            ),
            "domain_fraction_of_protein": (
                domain.length / protein_length
            ),
            "domain_distance_to_n": domain.start - 1,
            "domain_distance_to_c": (
                protein_length - domain.end
            ),
            "domain_position_n_terminal": int(
                position_category == "n_terminal"
            ),
            "domain_position_middle": int(
                position_category == "middle"
            ),
            "domain_position_c_terminal": int(
                position_category == "c_terminal"
            ),
        }

class DomainFunctionFeatureExtractor:
    """Extract functional annotation features from a protein domain."""

    def extract(
        self,
        domain: DomainRecord,
        context: DomainFeatureContext,
    ) -> FeatureMap:
        return {
            "num_functional_roles": len(
                domain.functional_roles
            ),
            "num_binding_types": len(
                domain.binding_types
            ),
            "num_recognized_motifs": len(
                domain.recognized_motifs
            ),
            "num_known_ligands": len(
                domain.known_ligands
            ),
            "num_known_ions": len(
                domain.known_ions
            ),
            "num_known_cofactors": len(
                domain.known_cofactors
            ),
            "has_functional_annotation": int(
                bool(domain.functional_roles)
            ),
            "has_binding_annotation": int(
                bool(domain.binding_types)
            ),
            "has_recognized_motifs": int(
                bool(domain.recognized_motifs)
            ),
            "has_known_ligands": int(
                bool(domain.known_ligands)
            ),
            "has_known_ions": int(
                bool(domain.known_ions)
            ),
            "has_known_cofactors": int(
                bool(domain.known_cofactors)
            ),
        }

class DomainAnnotationFeatureExtractor:
    """Extract identifier and classification annotation features."""

    def extract(
        self,
        domain: DomainRecord,
        context: DomainFeatureContext,
    ) -> FeatureMap:
        return {
            "has_domain_id": int(
                bool(domain.domain_id)
            ),
            "has_pfam_id": int(
                bool(domain.pfam_id)
            ),
            "has_interpro_id": int(
                bool(domain.interpro_id)
            ),
            "has_smart_id": int(
                bool(domain.smart_id)
            ),
            "has_cdd_id": int(
                bool(domain.cdd_id)
            ),
            "has_family": int(
                bool(domain.family)
            ),
            "has_superfamily": int(
                bool(domain.superfamily)
            ),
            "has_clan": int(
                bool(domain.clan)
            ),
            "num_external_identifiers": sum(
                (
                    bool(domain.pfam_id),
                    bool(domain.interpro_id),
                    bool(domain.smart_id),
                    bool(domain.cdd_id),
                )
            ),
            "num_classification_annotations": sum(
                (
                    bool(domain.family),
                    bool(domain.superfamily),
                    bool(domain.clan),
                )
            ),
        }

class DomainStructureFeatureExtractor:
    """Extract structural annotation features from a protein domain."""

    def extract(
        self,
        domain: DomainRecord,
        context: DomainFeatureContext,
    ) -> FeatureMap:
        has_structural_class = bool(domain.structural_class)
        has_structural_architecture = bool(
            domain.structural_architecture
        )
        has_structural_fold = bool(domain.structural_fold)

        return {
            "has_structural_class": int(
                has_structural_class
            ),
            "has_structural_architecture": int(
                has_structural_architecture
            ),
            "has_structural_fold": int(
                has_structural_fold
            ),
            "has_structure_annotation": int(
                has_structural_class
                or has_structural_architecture
                or has_structural_fold
            ),
            "num_structure_annotations": sum(
                (
                    has_structural_class,
                    has_structural_architecture,
                    has_structural_fold,
                )
            ),
        }
