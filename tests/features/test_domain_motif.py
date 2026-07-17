from protein_graph.entities.domain import DomainRecord
from protein_graph.entities.motif import MotifRecord
from protein_graph.features.pair.domain_motif import (
    DomainMotifCompatibilityExtractor,
)


def test_build_finds_domain_motif_compatibility() -> None:
    domain = DomainRecord(
        protein_id="PROTEIN_A",
        region_type="domain",
        start=10,
        end=60,
        sequence="A" * 51,
        domain_id="DOMAIN_SH3",
        family="SH3",
    )

    motif = MotifRecord(
        protein_id="PROTEIN_B",
        region_type="motif",
        start=100,
        end=103,
        sequence="PLLP",
        motif_id="LIG_SH3_GENERIC",
        motif_class="SH3-binding motif",
        pattern=r"P..P",
        recognized_by=("SH3",),
        evidence_source="ELM",
    )

    extractor = DomainMotifCompatibilityExtractor()

    features = extractor.build(
        domains_a=[domain],
        motifs_a=[],
        domains_b=[],
        motifs_b=[motif],
    )

    assert (
        features[
            "has_compatible_domain_motif_pair"
        ]
        == 1
    )

    assert (
        features[
            "compatible_domain_motif_count"
        ]
        == 1
    )

    assert (
        features[
            "a_domain_to_b_motif_count"
        ]
        == 1
    )

    assert (
        features[
            "b_domain_to_a_motif_count"
        ]
        == 0
    )

    match = features["domain_motif_matches"][0]

    assert match.domain_protein_id == "PROTEIN_A"
    assert match.motif_protein_id == "PROTEIN_B"
    assert match.motif_sequence == "PLLP"


def test_build_checks_both_directions() -> None:
    domain_a = DomainRecord(
        protein_id="PROTEIN_A",
        region_type="domain",
        start=1,
        end=30,
        sequence="A" * 30,
        domain_id="DOMAIN_SH3_A",
        family="SH3",
    )

    motif_a = MotifRecord(
        protein_id="PROTEIN_A",
        region_type="motif",
        start=50,
        end=53,
        sequence="PQQP",
        motif_id="LIG_SH3_GENERIC",
        recognized_by=("SH3",),
    )

    domain_b = DomainRecord(
        protein_id="PROTEIN_B",
        region_type="domain",
        start=1,
        end=30,
        sequence="A" * 30,
        domain_id="DOMAIN_SH3_B",
        family="SH3",
    )

    motif_b = MotifRecord(
        protein_id="PROTEIN_B",
        region_type="motif",
        start=50,
        end=53,
        sequence="PRRP",
        motif_id="LIG_SH3_GENERIC",
        recognized_by=("SH3",),
    )

    extractor = DomainMotifCompatibilityExtractor()

    features = extractor.build(
        domains_a=[domain_a],
        motifs_a=[motif_a],
        domains_b=[domain_b],
        motifs_b=[motif_b],
    )

    assert (
        features[
            "compatible_domain_motif_count"
        ]
        == 2
    )

    assert (
        features[
            "a_domain_to_b_motif_count"
        ]
        == 1
    )

    assert (
        features[
            "b_domain_to_a_motif_count"
        ]
        == 1
    )


def test_build_returns_zero_features_without_matches() -> None:
    domain = DomainRecord(
        protein_id="PROTEIN_A",
        region_type="domain",
        start=1,
        end=30,
        sequence="A" * 30,
        domain_id="DOMAIN_PDZ",
        family="PDZ",
    )

    motif = MotifRecord(
        protein_id="PROTEIN_B",
        region_type="motif",
        start=50,
        end=53,
        sequence="PQQP",
        motif_id="LIG_SH3_GENERIC",
        recognized_by=("SH3",),
    )

    extractor = DomainMotifCompatibilityExtractor()

    features = extractor.build(
        domains_a=[domain],
        motifs_a=[],
        domains_b=[],
        motifs_b=[motif],
    )

    assert (
        features[
            "has_compatible_domain_motif_pair"
        ]
        == 0
    )

    assert (
        features[
            "compatible_domain_motif_count"
        ]
        == 0
    )
