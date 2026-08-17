import pytest

from protein_graph.features.motif_properties_of_certain_domain.domain_motif_profile import (
    DomainMotifProfileBuilder,
    KnownMotif,
)


class FakeMotifSequenceFeatureExtractor:
    """
    Упрощённый экстрактор для тестирования ProfileBuilder
    независимо от настоящих sequence features.
    """

    def extract(
        self,
        sequence: str,
    ) -> dict[str, float]:
        length = len(sequence)

        return {
            "length": float(length),
            "proline_fraction": sequence.count("P") / length,
            "basic_fraction": (
                sequence.count("K")
                + sequence.count("R")
                + sequence.count("H")
            ) / length,
        }


@pytest.fixture
def builder() -> DomainMotifProfileBuilder:
    return DomainMotifProfileBuilder(
        feature_extractor=FakeMotifSequenceFeatureExtractor()
    )


@pytest.fixture
def sh3_motifs() -> tuple[KnownMotif, ...]:
    return (
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
    )


def test_build_preserves_every_known_motif(
    builder: DomainMotifProfileBuilder,
    sh3_motifs: tuple[KnownMotif, ...],
) -> None:
    profile = builder.build(
        domain_class="SH3",
        motifs=sh3_motifs,
    )

    assert profile.domain_class == "SH3"
    assert profile.motif_count == 3

    assert profile.motif_ids == (
        "SH3_MOTIF_1",
        "SH3_MOTIF_2",
        "SH3_MOTIF_3",
    )

    assert profile.motif_sequences == (
        "RPLPPLPR",
        "KPPVPPPR",
        "APPLPPRK",
    )


def test_build_extracts_features_for_each_motif(
    builder: DomainMotifProfileBuilder,
    sh3_motifs: tuple[KnownMotif, ...],
) -> None:
    profile = builder.build(
        domain_class="SH3",
        motifs=sh3_motifs,
    )

    first_motif = profile.known_motifs[0]

    assert first_motif.motif.sequence == "RPLPPLPR"
    assert first_motif.features["length"] == 8.0
    assert first_motif.features["proline_fraction"] == 0.5


def test_build_calculates_additional_distributions(
    builder: DomainMotifProfileBuilder,
    sh3_motifs: tuple[KnownMotif, ...],
) -> None:
    profile = builder.build(
        domain_class="SH3",
        motifs=sh3_motifs,
    )

    distribution = profile.feature_distributions[
        "proline_fraction"
    ]

    assert distribution.minimum <= distribution.mean
    assert distribution.mean <= distribution.maximum
    assert distribution.standard_deviation >= 0.0


def test_profile_returns_available_motif_lengths(
    builder: DomainMotifProfileBuilder,
    sh3_motifs: tuple[KnownMotif, ...],
) -> None:
    profile = builder.build(
        domain_class="SH3",
        motifs=sh3_motifs,
    )

    assert profile.motif_lengths == (8,)


def test_profile_filters_motifs_by_length(
    builder: DomainMotifProfileBuilder,
) -> None:
    profile = builder.build(
        domain_class="SH3",
        motifs=(
            KnownMotif(
                motif_id="M1",
                sequence="PPLP",
            ),
            KnownMotif(
                motif_id="M2",
                sequence="RPLPPLPR",
            ),
        ),
    )

    motifs = profile.motifs_of_length(4)

    assert len(motifs) == 1
    assert motifs[0].motif.motif_id == "M1"


def test_empty_motif_collection_raises(
    builder: DomainMotifProfileBuilder,
) -> None:
    with pytest.raises(
        ValueError,
        match="At least one known motif",
    ):
        builder.build(
            domain_class="SH3",
            motifs=(),
        )


def test_duplicate_motif_ids_raise(
    builder: DomainMotifProfileBuilder,
) -> None:
    motifs = (
        KnownMotif(
            motif_id="DUPLICATED",
            sequence="RPLPPLPR",
        ),
        KnownMotif(
            motif_id="DUPLICATED",
            sequence="KPPVPPPR",
        ),
    )

    with pytest.raises(
        ValueError,
        match="Duplicate motif IDs",
    ):
        builder.build(
            domain_class="SH3",
            motifs=motifs,
        )


def test_invalid_motif_sequence_raises() -> None:
    with pytest.raises(
        ValueError,
        match="Unsupported amino-acid",
    ):
        KnownMotif(
            motif_id="INVALID",
            sequence="RPLXPLPR",
        )
