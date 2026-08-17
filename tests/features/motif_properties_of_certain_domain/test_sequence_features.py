import pytest

from protein_graph.features.motif_properties_of_certain_domain.sequence_features import (
    MotifSequenceFeatureExtractor,
)


def test_extracts_motif_specific_features() -> None:
    extractor = MotifSequenceFeatureExtractor()

    features = extractor.extract("RPLPPLPR")

    assert features["length"] == 8
    assert features["fraction_P"] == pytest.approx(
        0.5
    )
    assert features["positive_fraction"] == pytest.approx(
        0.25
    )
    assert features["proline_xx_proline_count"] == 2
    assert features["sequence_entropy"] > 0.0


def test_counts_overlapping_proline_patterns() -> None:
    extractor = MotifSequenceFeatureExtractor()

    features = extractor.extract("RPLPPLPR")

    # Пролины находятся на позициях 1, 3, 4 и 6.
    #
    # P-X-X-P:
    # позиции 1 -> 4
    # позиции 3 -> 6
    assert features["proline_xx_proline_count"] == 2

    # P-X-P:
    # позиции 1 -> 3
    # позиции 4 -> 6
    assert features["proline_x_proline_count"] == 2

    # PP:
    # позиции 3 -> 4
    assert features["adjacent_proline_pair_count"] == 1


def test_order_changes_positional_features() -> None:
    extractor = MotifSequenceFeatureExtractor()

    first = extractor.extract("RRPP")
    second = extractor.extract("PPRR")

    assert (
        first["positive_position_mean"]
        != second["positive_position_mean"]
    )

    assert first["positive_position_mean"] < 0.0
    assert second["positive_position_mean"] > 0.0


def test_single_amino_acid_sequence() -> None:
    extractor = MotifSequenceFeatureExtractor()

    features = extractor.extract("P")

    assert features["length"] == 1
    assert features["fraction_P"] == 1.0
    assert features["sequence_entropy"] == 0.0
    assert features["normalized_entropy"] == 0.0
    assert features["maximum_identical_run"] == 1
    assert features["proline_x_proline_count"] == 0
    assert features["proline_xx_proline_count"] == 0


def test_sequence_is_normalized() -> None:
    extractor = MotifSequenceFeatureExtractor()

    features = extractor.extract("  rplpplpr  ")

    assert features["length"] == 8
    assert features["fraction_P"] == pytest.approx(
        0.5
    )


def test_empty_sequence_raises() -> None:
    extractor = MotifSequenceFeatureExtractor()

    with pytest.raises(
        ValueError,
        match="sequence must not be empty",
    ):
        extractor.extract("   ")


def test_invalid_amino_acid_raises() -> None:
    extractor = MotifSequenceFeatureExtractor()

    with pytest.raises(
        ValueError,
        match="Unsupported amino-acid",
    ):
        extractor.extract("RPLXPLPR")
