"""Feature extraction for individual ELM motif occurrences."""

from __future__ import annotations

from dataclasses import dataclass

from protein_graph.entities.motif import MotifRecord
from protein_graph.features.base import FeatureMap


@dataclass(frozen=True)
class MotifFeatureContext:
    """Protein-level context required by motif extractors."""

    protein_length: int

    def __post_init__(self) -> None:
        if self.protein_length <= 0:
            raise ValueError(
                "protein_length must be greater than zero."
            )


class MotifGeometryFeatureExtractor:
    """Extract geometric features from one motif occurrence."""

    def extract(
        self,
        motif: MotifRecord,
        context: MotifFeatureContext,
    ) -> FeatureMap:
        """Extract motif position and length features."""
        protein_length = context.protein_length

        if motif.end > protein_length:
            raise ValueError(
                "Motif coordinates exceed protein length: "
                f"motif end={motif.end}, "
                f"protein length={protein_length}."
            )

        position_category = motif.position_category(
            protein_length
        )

        return {
            "motif_length": motif.length,
            "motif_relative_position": (
                motif.relative_position(protein_length)
            ),
            "motif_fraction_of_protein": (
                motif.length / protein_length
            ),
            "motif_distance_to_n": motif.start - 1,
            "motif_distance_to_c": (
                protein_length - motif.end
            ),
            "motif_position_n_terminal": int(
                position_category == "n_terminal"
            ),
            "motif_position_middle": int(
                position_category == "middle"
            ),
            "motif_position_c_terminal": int(
                position_category == "c_terminal"
            ),
        }

class MotifAnnotationFeatureExtractor:
    """Extract annotation features from one motif occurrence."""

    def extract(
        self,
        motif: MotifRecord,
        context: MotifFeatureContext,
    ) -> FeatureMap:
        """Extract motif annotation availability and counts."""
        del context

        has_motif_id = bool(motif.motif_id)
        has_motif_class = bool(motif.motif_class)
        has_pattern = bool(motif.pattern)
        has_recognized_by = bool(motif.recognized_by)
        has_interaction_role = bool(
            motif.interaction_role
        )
        has_evidence_source = bool(
            motif.evidence_source
        )
        has_evidence_level = bool(
            motif.evidence_level
        )
        has_confidence = (
            motif.confidence is not None
        )

        return {
            "has_motif_id": int(has_motif_id),
            "has_motif_class": int(
                has_motif_class
            ),
            "has_motif_pattern": int(
                has_pattern
            ),
            "has_recognized_by": int(
                has_recognized_by
            ),
            "has_interaction_role": int(
                has_interaction_role
            ),
            "has_evidence_source": int(
                has_evidence_source
            ),
            "has_evidence_level": int(
                has_evidence_level
            ),
            "has_motif_confidence": int(
                has_confidence
            ),
            "num_recognized_domains": len(
                motif.recognized_by
            ),
            "num_motif_annotations": sum(
                (
                    has_motif_id,
                    has_motif_class,
                    has_pattern,
                    has_recognized_by,
                    has_interaction_role,
                    has_evidence_source,
                    has_evidence_level,
                    has_confidence,
                )
            ),
            "motif_confidence": (
                motif.confidence
                if motif.confidence is not None
                else 0.0
            ),
        }
