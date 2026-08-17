from protein_graph.features.motif_properties_of_certain_domain.combined_features import (
    CombinedMotifFeatureExtractor,
)
from protein_graph.features.motif_properties_of_certain_domain.similarity_features import (
    MotifSimilarityConfig,
    MotifSimilarityFeatureExtractor,
)


def test_combines_all_feature_groups() -> None:
    extractor = CombinedMotifFeatureExtractor(
        similarity_extractor=(
            MotifSimilarityFeatureExtractor(
                MotifSimilarityConfig(
                    critical_positions=(1, 3, 4, 6),
                    use_blosum62=False,
                )
            )
        )
    )

    protein = "AAAAAKPLPPMPRGGGGG"

    features = extractor.extract(
        reference_motif="RPLPPLPR",
        protein_sequence=protein,
        candidate_start=5,
        candidate_end=13,
    )

    assert "reference_fraction_P" in features
    assert "candidate_fraction_P" in features
    assert "similarity_sequence_identity" in features
    assert "similarity_critical_exact_fraction" in features
    assert "context_relative_start" in features
