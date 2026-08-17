from protein_graph.features.motif_properties_of_certain_domain.domain_motif_profile import (
    DomainMotifProfileBuilder,
    KnownMotif,
)
from protein_graph.features.motif_properties_of_certain_domain.sequence_features import (
    MotifSequenceFeatureExtractor,
)


def test_builds_profile_with_real_motif_feature_extractor() -> None:
    builder = DomainMotifProfileBuilder(
        feature_extractor=MotifSequenceFeatureExtractor(),
    )

    profile = builder.build(
        domain_class="SH3",
        motifs=(
            KnownMotif(
                motif_id="SH3_MOTIF_1",
                motif_class="LIG_SH3_1",
                sequence="RPLPPLPR",
                confidence=0.95,
            ),
            KnownMotif(
                motif_id="SH3_MOTIF_2",
                motif_class="LIG_SH3_1",
                sequence="KPPVPPPR",
                confidence=0.90,
            ),
            KnownMotif(
                motif_id="SH3_MOTIF_3",
                motif_class="LIG_SH3_1",
                sequence="APPLPPRK",
                confidence=0.85,
            ),
        ),
    )

    assert profile.domain_class == "SH3"
    assert profile.motif_count == 3
    assert profile.motif_lengths == (8,)

    first = profile.known_motifs[0]

    assert first.motif.motif_id == "SH3_MOTIF_1"
    assert first.features["length"] == 8.0
    assert first.features["fraction_P"] == 0.5
    assert first.features["proline_xx_proline_count"] == 2.0

    assert "fraction_P" in profile.feature_distributions
    assert (
        profile.feature_distributions[
            "fraction_P"
        ].standard_deviation
        >= 0.0
    )
