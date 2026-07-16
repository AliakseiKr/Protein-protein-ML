"""Tests for region feature extraction."""

from protein_graph.entities import SequenceRegion
from protein_graph.pipelines import RegionFeaturePipeline


def test_builds_region_features() -> None:
    region = SequenceRegion(
        protein_id="P12345",
        region_type="transmembrane",
        start=10,
        end=19,
        sequence="AVILMFWYAV",
        source="UniProt",
    )

    features = RegionFeaturePipeline().build(region)

    assert features["protein_id"] == "P12345"
    assert features["region_type"] == "transmembrane"
    assert features["region_length"] == 10
    assert "hydrophobic_fraction" in features
    assert "aa_fraction_A" in features
    assert "max_hydrophobic_run" in features
    assert "protein_length" not in features


def test_build_many() -> None:
    regions = [
        SequenceRegion(
            protein_id="P1",
            region_type="domain",
            start=1,
            end=5,
            sequence="MAVIL",
        ),
        SequenceRegion(
            protein_id="P1",
            region_type="signal_peptide",
            start=6,
            end=10,
            sequence="MFWYA",
        ),
    ]

    features = RegionFeaturePipeline().build_many(regions)

    assert len(features) == 2
    assert set(features["region_type"]) == {
        "domain",
        "signal_peptide",
    }


def test_repr() -> None:
    assert "extractors=5" in repr(RegionFeaturePipeline())
