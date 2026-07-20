import pytest

from protein_graph.entities.domain import DomainRecord
from protein_graph.features.domains import (
    DomainFeatureContext,
    DomainFunctionFeatureExtractor,
)


def make_domain() -> DomainRecord:
    return DomainRecord(
        protein_id="P12345",
        start=10,
        end=30,
        sequence="A" * 21,
        region_type="domain",
        domain_id="DOMAIN_1",
        functional_roles={"kinase", "enzyme"},
        binding_types={"ATP"},
        recognized_motifs={"P-loop"},
        known_ligands={"ATP", "ADP"},
        known_ions={"Mg2+"},
        known_cofactors={"FAD"},
    )


def test_extracts_function_features():
    extractor = DomainFunctionFeatureExtractor()

    features = extractor.extract(
        make_domain(),
        DomainFeatureContext(protein_length=100),
    )

    assert features["num_functional_roles"] == 2
    assert features["num_binding_types"] == 1
    assert features["num_recognized_motifs"] == 1
    assert features["num_known_ligands"] == 2
    assert features["num_known_ions"] == 1
    assert features["num_known_cofactors"] == 1

    assert features["has_functional_annotation"] == 1
    assert features["has_binding_annotation"] == 1
    assert features["has_recognized_motifs"] == 1
    assert features["has_known_ligands"] == 1
    assert features["has_known_ions"] == 1
    assert features["has_known_cofactors"] == 1


def test_extracts_empty_function_features():
    domain = DomainRecord(
        protein_id="P12345",
        start=10,
        end=30,
        sequence="A" * 21,
        region_type="domain",
        domain_id="DOMAIN_1",
    )

    extractor = DomainFunctionFeatureExtractor()

    features = extractor.extract(
        domain,
        DomainFeatureContext(protein_length=100),
    )

    assert features["num_functional_roles"] == 0
    assert features["num_binding_types"] == 0
    assert features["num_recognized_motifs"] == 0
    assert features["num_known_ligands"] == 0
    assert features["num_known_ions"] == 0
    assert features["num_known_cofactors"] == 0

    assert features["has_functional_annotation"] == 0
    assert features["has_binding_annotation"] == 0
    assert features["has_recognized_motifs"] == 0
    assert features["has_known_ligands"] == 0
    assert features["has_known_ions"] == 0
    assert features["has_known_cofactors"] == 0
