"""Tests for protein domain feature extractors."""

import pytest

from protein_graph.entities.domain import DomainRecord
from protein_graph.features.domains import (
    DomainFeatureContext,
    DomainGeometryFeatureExtractor,
)


def make_domain(
    *,
    start: int = 11,
    end: int = 30,
) -> DomainRecord:
    """Create a domain record for feature extractor tests."""
    return DomainRecord(
        protein_id="P12345",
        start=start,
        end=end,
        sequence="A" * (end - start + 1),
        region_type="domain",
        domain_id="DOMAIN_1",
    )


@pytest.mark.parametrize("protein_length", [0, -1, -100])
def test_domain_feature_context_rejects_non_positive_length(
    protein_length: int,
) -> None:
    with pytest.raises(
        ValueError,
        match="protein_length must be greater than zero",
    ):
        DomainFeatureContext(protein_length=protein_length)


def test_extracts_domain_geometry_features() -> None:
    domain = make_domain(start=11, end=30)
    context = DomainFeatureContext(protein_length=100)
    extractor = DomainGeometryFeatureExtractor()

    features = extractor.extract(domain, context)

    assert features["domain_length"] == 20
    assert features["domain_relative_position"] == pytest.approx(
        0.205
    )
    assert features["domain_fraction_of_protein"] == pytest.approx(
        0.2
    )
    assert features["domain_distance_to_n"] == 10
    assert features["domain_distance_to_c"] == 70

    assert features["domain_position_n_terminal"] == 1
    assert features["domain_position_middle"] == 0
    assert features["domain_position_c_terminal"] == 0


def test_rejects_domain_outside_protein() -> None:
    domain = make_domain(start=90, end=110)
    context = DomainFeatureContext(protein_length=100)
    extractor = DomainGeometryFeatureExtractor()

    with pytest.raises(
        ValueError,
        match="Domain coordinates exceed protein length",
    ):
        extractor.extract(domain, context)


def test_domain_covering_entire_protein() -> None:
    domain = make_domain(start=1, end=100)
    context = DomainFeatureContext(protein_length=100)
    extractor = DomainGeometryFeatureExtractor()

    features = extractor.extract(domain, context)

    assert features["domain_length"] == 100
    assert features["domain_relative_position"] == pytest.approx(
        0.505
    )
    assert features["domain_fraction_of_protein"] == pytest.approx(
        1.0
    )
    assert features["domain_distance_to_n"] == 0
    assert features["domain_distance_to_c"] == 0

    assert features["domain_position_n_terminal"] == 0
    assert features["domain_position_middle"] == 1
    assert features["domain_position_c_terminal"] == 0


def test_extracts_n_terminal_position_category() -> None:
    domain = make_domain(start=1, end=20)
    context = DomainFeatureContext(protein_length=100)
    extractor = DomainGeometryFeatureExtractor()

    features = extractor.extract(domain, context)

    assert features["domain_position_n_terminal"] == 1
    assert features["domain_position_middle"] == 0
    assert features["domain_position_c_terminal"] == 0


def test_extracts_middle_position_category() -> None:
    domain = make_domain(start=41, end=60)
    context = DomainFeatureContext(protein_length=100)
    extractor = DomainGeometryFeatureExtractor()

    features = extractor.extract(domain, context)

    assert features["domain_position_n_terminal"] == 0
    assert features["domain_position_middle"] == 1
    assert features["domain_position_c_terminal"] == 0


def test_extracts_c_terminal_position_category() -> None:
    domain = make_domain(start=81, end=100)
    context = DomainFeatureContext(protein_length=100)
    extractor = DomainGeometryFeatureExtractor()

    features = extractor.extract(domain, context)

    assert features["domain_position_n_terminal"] == 0
    assert features["domain_position_middle"] == 0
    assert features["domain_position_c_terminal"] == 1


def test_position_category_is_one_hot_encoded() -> None:
    domain = make_domain(start=41, end=60)
    context = DomainFeatureContext(protein_length=100)
    extractor = DomainGeometryFeatureExtractor()

    features = extractor.extract(domain, context)

    category_values = (
        features["domain_position_n_terminal"],
        features["domain_position_middle"],
        features["domain_position_c_terminal"],
    )

    assert sum(category_values) == 1
