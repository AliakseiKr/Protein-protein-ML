"""Tests for domain feature extraction."""

from protein_graph.entities import DomainRecord
from protein_graph.pipelines import DomainFeaturePipeline

import pytest

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

def test_pipeline_combines_all_default_domain_features():
    domain = DomainRecord(
        protein_id="P12345",
        start=11,
        end=30,
        sequence="A" * 20,
        region_type="domain",
        domain_id="DOMAIN_1",
        pfam_id="PF00001",
        interpro_id="IPR000001",
        family="Protein kinase family",
        functional_roles={"enzyme", "signaling"},
        binding_types={"ATP"},
        recognized_motifs={"P-loop"},
        structural_class="Alpha and beta proteins",
        structural_architecture="Three-layer sandwich",
        structural_fold="Rossmann-like fold",
        known_ligands={"ATP", "ADP"},
        known_ions={"Mg2+"},
        known_cofactors={"FAD"},
    )

    pipeline = DomainFeaturePipeline()

    table = pipeline.build_from_domains(
        domains=[domain],
        protein_length=100,
    )

    assert len(table) == 1

    row = table.iloc[0]

    # Geometry features
    assert row["domain_length"] == 20
    assert row["domain_relative_position"] == pytest.approx(
        0.205
    )
    assert row["domain_fraction_of_protein"] == pytest.approx(
        0.2
    )
    assert row["domain_distance_to_n"] == 10
    assert row["domain_distance_to_c"] == 70
    assert row["domain_position_n_terminal"] == 1
    assert row["domain_position_middle"] == 0
    assert row["domain_position_c_terminal"] == 0

    # Functional features
    assert row["num_functional_roles"] == 2
    assert row["num_binding_types"] == 1
    assert row["num_recognized_motifs"] == 1
    assert row["num_known_ligands"] == 2
    assert row["num_known_ions"] == 1
    assert row["num_known_cofactors"] == 1

    assert row["has_functional_annotation"] == 1
    assert row["has_binding_annotation"] == 1
    assert row["has_recognized_motifs"] == 1
    assert row["has_known_ligands"] == 1
    assert row["has_known_ions"] == 1
    assert row["has_known_cofactors"] == 1

    # Identifier and classification features
    assert row["has_domain_id"] == 1
    assert row["has_pfam_id"] == 1
    assert row["has_interpro_id"] == 1
    assert row["has_smart_id"] == 0
    assert row["has_cdd_id"] == 0
    assert row["has_family"] == 1
    assert row["has_superfamily"] == 0
    assert row["has_clan"] == 0
    assert row["num_external_identifiers"] == 2
    assert row["num_classification_annotations"] == 1

    # Structural features
    assert row["has_structural_class"] == 1
    assert row["has_structural_architecture"] == 1
    assert row["has_structural_fold"] == 1
    assert row["has_structure_annotation"] == 1
    assert row["num_structure_annotations"] == 3

    # Raw annotations are still preserved
    assert row["domain_id"] == "DOMAIN_1"
    assert row["pfam_id"] == "PF00001"
    assert row["family"] == "Protein kinase family"
    assert row["functional_roles"] == [
        "enzyme",
        "signaling",
    ] or row["functional_roles"] == [
        "signaling",
        "enzyme",
    ]
