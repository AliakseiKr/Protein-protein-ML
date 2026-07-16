"""Tests for domain feature extraction."""

from protein_graph.entities import DomainRecord
from protein_graph.pipelines import DomainFeaturePipeline


def test_builds_domain_feature_table() -> None:
    protein = {
        "protein_id": "P12345",
        "sequence": "MAVILMFWYKDE",
        "domain_annotations": [
            {
                "domain_id": "SH3",
                "pfam_id": "PF00018",
                "start": 3,
                "end": 7,
                "description": "SH3 domain",
                "functional_roles": ["binding"],
                "binding_types": ["domain_peptide"],
                "recognized_motifs": ["PXXP"],
            }
        ],
    }

    features = DomainFeaturePipeline().build(protein)

    assert len(features) == 1
    assert features.iloc[0]["domain_id"] == "SH3"
    assert features.iloc[0]["pfam_id"] == "PF00018"
    assert features.iloc[0]["region_length"] == 5
    assert features.iloc[0]["position_category"] == "middle"
    assert "hydrophobic_fraction" in features.columns


def test_build_from_domain_records() -> None:
    domain = DomainRecord(
        protein_id="P1",
        region_type="domain",
        start=1,
        end=5,
        sequence="MAVIL",
        domain_id="D1",
        structural_fold="beta sandwich",
    )

    features = DomainFeaturePipeline().build_from_domains(
        domains=[domain],
        protein_length=20,
    )

    assert len(features) == 1
    assert features.iloc[0]["structural_fold"] == "beta sandwich"
    assert features.iloc[0]["position_category"] == "n_terminal"


def test_empty_domains_return_empty_table() -> None:
    protein = {
        "protein_id": "P1",
        "sequence": "MAVIL",
        "domain_annotations": [],
    }

    features = DomainFeaturePipeline().build(protein)

    assert features.empty


def test_repr() -> None:
    assert repr(DomainFeaturePipeline()) == "DomainFeaturePipeline()"
