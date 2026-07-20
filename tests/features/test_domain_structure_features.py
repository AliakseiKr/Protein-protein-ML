"""Tests for domain structure feature extraction."""

from protein_graph.entities.domain import DomainRecord
from protein_graph.features.domains import (
    DomainFeatureContext,
    DomainStructureFeatureExtractor,
)


def make_domain(
    **overrides: object,
) -> DomainRecord:
    values = {
        "protein_id": "P12345",
        "start": 10,
        "end": 30,
        "sequence": "A" * 21,
        "region_type": "domain",
        "domain_id": "DOMAIN_1",
    }

    values.update(overrides)

    return DomainRecord(**values)


def test_extracts_complete_structure_features():
    domain = make_domain(
        structural_class="Alpha and beta proteins",
        structural_architecture="Three-layer sandwich",
        structural_fold="Rossmann-like fold",
    )

    features = DomainStructureFeatureExtractor().extract(
        domain,
        DomainFeatureContext(protein_length=100),
    )

    assert features["has_structural_class"] == 1
    assert features["has_structural_architecture"] == 1
    assert features["has_structural_fold"] == 1
    assert features["has_structure_annotation"] == 1
    assert features["num_structure_annotations"] == 3


def test_extracts_partial_structure_features():
    domain = make_domain(
        structural_fold="Rossmann-like fold",
    )

    features = DomainStructureFeatureExtractor().extract(
        domain,
        DomainFeatureContext(protein_length=100),
    )

    assert features["has_structural_class"] == 0
    assert features["has_structural_architecture"] == 0
    assert features["has_structural_fold"] == 1
    assert features["has_structure_annotation"] == 1
    assert features["num_structure_annotations"] == 1


def test_extracts_missing_structure_features():
    domain = make_domain()

    features = DomainStructureFeatureExtractor().extract(
        domain,
        DomainFeatureContext(protein_length=100),
    )

    assert features["has_structural_class"] == 0
    assert features["has_structural_architecture"] == 0
    assert features["has_structural_fold"] == 0
    assert features["has_structure_annotation"] == 0
    assert features["num_structure_annotations"] == 0
