"""Tests for domain annotation feature extraction."""

from protein_graph.entities.domain import DomainRecord
from protein_graph.features.domains import (
    DomainAnnotationFeatureExtractor,
    DomainFeatureContext,
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


def test_extracts_identifier_features():
    domain = make_domain(
        pfam_id="PF00001",
        interpro_id="IPR000001",
        smart_id="SM00001",
        cdd_id="cd00001",
    )

    features = DomainAnnotationFeatureExtractor().extract(
        domain,
        DomainFeatureContext(protein_length=100),
    )

    assert features["has_domain_id"] == 1
    assert features["has_pfam_id"] == 1
    assert features["has_interpro_id"] == 1
    assert features["has_smart_id"] == 1
    assert features["has_cdd_id"] == 1
    assert features["num_external_identifiers"] == 4


def test_extracts_classification_features():
    domain = make_domain(
        family="Protein kinase family",
        superfamily="Kinase superfamily",
        clan="CL0016",
    )

    features = DomainAnnotationFeatureExtractor().extract(
        domain,
        DomainFeatureContext(protein_length=100),
    )

    assert features["has_family"] == 1
    assert features["has_superfamily"] == 1
    assert features["has_clan"] == 1
    assert features["num_classification_annotations"] == 3


def test_extracts_missing_annotation_features():
    domain = make_domain()

    features = DomainAnnotationFeatureExtractor().extract(
        domain,
        DomainFeatureContext(protein_length=100),
    )

    assert features["has_domain_id"] == 1
    assert features["has_pfam_id"] == 0
    assert features["has_interpro_id"] == 0
    assert features["has_smart_id"] == 0
    assert features["has_cdd_id"] == 0
    assert features["has_family"] == 0
    assert features["has_superfamily"] == 0
    assert features["has_clan"] == 0
    assert features["num_external_identifiers"] == 0
    assert features["num_classification_annotations"] == 0


def test_counts_partial_annotations():
    domain = make_domain(
        pfam_id="PF00001",
        cdd_id="cd00001",
        family="Protein kinase family",
    )

    features = DomainAnnotationFeatureExtractor().extract(
        domain,
        DomainFeatureContext(protein_length=100),
    )

    assert features["num_external_identifiers"] == 2
    assert features["num_classification_annotations"] == 1
