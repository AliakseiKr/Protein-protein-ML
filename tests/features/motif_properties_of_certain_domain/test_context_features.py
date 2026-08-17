from protein_graph.features.motif_properties_of_certain_domain.context_features import (
    MotifContextConfig,
    MotifContextFeatureExtractor,
)
from protein_graph.features.motif_properties_of_certain_domain.window import (
    WindowInterval,
)


def test_extracts_position_and_flank_features() -> None:
    extractor = MotifContextFeatureExtractor(
        MotifContextConfig(
            flank_size=3,
        )
    )

    protein = "AAAAARPLPPLPRGGGGG"

    features = extractor.extract(
        protein_sequence=protein,
        motif_start=5,
        motif_end=13,
    )

    assert features["motif_length"] == 8
    assert features["left_flank_length"] == 3
    assert features["right_flank_length"] == 3
    assert features["distance_to_n_terminus"] == 5


def test_detects_domain_overlap() -> None:
    extractor = MotifContextFeatureExtractor()

    features = extractor.extract(
        protein_sequence="AAAAARPLPPLPRGGGGG",
        motif_start=5,
        motif_end=13,
        domain_regions=(
            WindowInterval(
                start=3,
                end=10,
            ),
        ),
    )

    assert features["overlaps_known_domain"] == 1
