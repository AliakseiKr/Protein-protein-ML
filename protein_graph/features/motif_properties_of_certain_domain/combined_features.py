from __future__ import annotations

from typing import Iterable

from protein_graph.features.motif_properties_of_certain_domain.context_features import (
    MotifContextFeatureExtractor,
)
from protein_graph.features.motif_properties_of_certain_domain.window import (
    WindowInterval,
)
from protein_graph.features.motif_properties_of_certain_domain.sequence_features import (
    MotifSequenceFeatureExtractor,
)
from protein_graph.features.motif_properties_of_certain_domain.similarity_features import (
    MotifSimilarityFeatureExtractor,
)


class CombinedMotifFeatureExtractor:
    """
    Собирает полный набор признаков пары:

    reference motif ↔ candidate region.
    """

    def __init__(
        self,
        *,
        sequence_extractor: (
            MotifSequenceFeatureExtractor | None
        ) = None,
        similarity_extractor: (
            MotifSimilarityFeatureExtractor | None
        ) = None,
        context_extractor: (
            MotifContextFeatureExtractor | None
        ) = None,
    ) -> None:
        self._sequence_extractor = (
            sequence_extractor
            or MotifSequenceFeatureExtractor()
        )

        self._similarity_extractor = (
            similarity_extractor
            or MotifSimilarityFeatureExtractor()
        )

        self._context_extractor = (
            context_extractor
            or MotifContextFeatureExtractor()
        )

    def extract(
        self,
        *,
        reference_motif: str,
        protein_sequence: str,
        candidate_start: int,
        candidate_end: int,
        domain_regions: Iterable[WindowInterval] = (),
        transmembrane_regions: Iterable[
            WindowInterval
        ] = (),
    ) -> dict[str, float | int]:
        candidate_sequence = protein_sequence[
            candidate_start:candidate_end
        ]

        reference_features = (
            self._sequence_extractor.extract(
                reference_motif
            )
        )

        candidate_features = (
            self._sequence_extractor.extract(
                candidate_sequence
            )
        )

        similarity_features = (
            self._similarity_extractor.extract(
                reference_motif=reference_motif,
                candidate_sequence=candidate_sequence,
            )
        )

        context_features = (
            self._context_extractor.extract(
                protein_sequence=protein_sequence,
                motif_start=candidate_start,
                motif_end=candidate_end,
                domain_regions=domain_regions,
                transmembrane_regions=(
                    transmembrane_regions
                ),
            )
        )

        combined: dict[str, float | int] = {}

        combined.update({
            f"reference_{name}": value
            for name, value
            in reference_features.items()
        })

        combined.update({
            f"candidate_{name}": value
            for name, value
            in candidate_features.items()
        })

        combined.update({
            f"similarity_{name}": value
            for name, value
            in similarity_features.items()
        })

        combined.update({
            f"context_{name}": value
            for name, value
            in context_features.items()
        })

        return combined
