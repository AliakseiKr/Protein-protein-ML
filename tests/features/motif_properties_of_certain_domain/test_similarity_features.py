import pytest

from protein_graph.features.motif_properties_of_certain_domain.similarity_features import (
    MotifSimilarityConfig,
    MotifSimilarityFeatureExtractor,
)


def test_identical_sequences_have_maximum_similarity() -> None:
    extractor = MotifSimilarityFeatureExtractor(
        MotifSimilarityConfig(
            critical_positions=(1, 3, 4, 6),
            use_blosum62=False,
        )
    )

    features = extractor.extract(
        reference_motif="RPLPPLPR",
        candidate_sequence="RPLPPLPR",
    )

    assert features["exact_match"] == 1
    assert features["sequence_identity"] == 1.0
    assert features["critical_exact_fraction"] == 1.0
    assert features["kmer_2_jaccard"] == 1.0


def test_preserved_critical_positions_are_detected() -> None:
    extractor = MotifSimilarityFeatureExtractor(
        MotifSimilarityConfig(
            critical_positions=(1, 3, 4, 6),
            use_blosum62=False,
        )
    )

    features = extractor.extract(
        reference_motif="RPLPPLPR",
        candidate_sequence="KPLPPMPR",
    )

    assert features["sequence_identity"] == pytest.approx(
        0.75
    )

    assert features["critical_exact_fraction"] == 1.0


def test_changed_order_reduces_similarity() -> None:
    extractor = MotifSimilarityFeatureExtractor(
        MotifSimilarityConfig(
            use_blosum62=False,
        )
    )

    features = extractor.extract(
        reference_motif="PAPP",
        candidate_sequence="PPAA",
    )

    assert features["sequence_identity"] < 1.0
    assert features["kmer_2_jaccard"] < 1.0
