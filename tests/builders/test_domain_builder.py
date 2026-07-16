"""Tests for the domain builder."""

import pytest

from protein_graph.builders import DomainBuilder


def test_builds_domain_sequence() -> None:
    protein = {
        "protein_id": "P12345",
        "sequence": "MAVILMFWYKDE",
        "domain_annotations": [
            {
                "domain_id": "PF00018",
                "pfam_id": "PF00018",
                "description": "SH3 domain",
                "start": 3,
                "end": 7,
                "source": "UniProt",
                "functional_roles": ["binding"],
                "binding_types": ["domain_peptide"],
                "recognized_motifs": ["PXXP"],
            }
        ],
    }

    domains = DomainBuilder().build_many(protein)

    assert len(domains) == 1

    domain = domains[0]

    assert domain.start == 3
    assert domain.end == 7
    assert domain.sequence == "VILMF"
    assert domain.pfam_id == "PF00018"
    assert domain.functional_roles == ("binding",)


def test_builds_multiple_domains() -> None:
    protein = {
        "protein_id": "P12345",
        "sequence": "MAVILMFWYKDE",
        "domain_annotations": [
            {
                "domain_id": "D1",
                "start": 1,
                "end": 3,
            },
            {
                "domain_id": "D2",
                "start": 8,
                "end": 12,
            },
        ],
    }

    domains = DomainBuilder().build_many(protein)

    assert len(domains) == 2
    assert domains[0].sequence == "MAV"
    assert domains[1].sequence == "WYKDE"


def test_missing_annotations_returns_empty_list() -> None:
    protein = {
        "protein_id": "P12345",
        "sequence": "MAVIL",
    }

    assert DomainBuilder().build_many(protein) == []


def test_rejects_domain_outside_protein() -> None:
    protein = {
        "protein_id": "P12345",
        "sequence": "MAVIL",
        "domain_annotations": [
            {
                "domain_id": "D1",
                "start": 2,
                "end": 10,
            }
        ],
    }

    with pytest.raises(
        ValueError,
        match="exceeds protein",
    ):
        DomainBuilder().build_many(protein)


def test_rejects_invalid_coordinates() -> None:
    protein = {
        "protein_id": "P12345",
        "sequence": "MAVIL",
        "domain_annotations": [
            {
                "domain_id": "D1",
                "start": "1",
                "end": 5,
            }
        ],
    }

    with pytest.raises(
        ValueError,
        match="start must be an integer",
    ):
        DomainBuilder().build_many(protein)


def test_repr() -> None:
    assert repr(DomainBuilder()) == "DomainBuilder()"
