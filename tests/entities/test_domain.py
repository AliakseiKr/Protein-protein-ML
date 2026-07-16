"""Tests for protein domain entities."""

import pytest

from protein_graph.entities import DomainRecord


def test_domain_record() -> None:
    domain = DomainRecord(
        protein_id="P12345",
        region_type="domain",
        start=10,
        end=19,
        sequence="AVILMFWYKD",
        pfam_id="PF00018",
        interpro_id="IPR001452",
        family="SH3 domain family",
        functional_roles=("binding",),
        binding_types=("domain_peptide",),
        recognized_motifs=("PXXP",),
    )

    assert domain.length == 10
    assert domain.pfam_id == "PF00018"
    assert domain.recognized_motifs == ("PXXP",)


def test_domain_requires_identifier() -> None:
    with pytest.raises(
        ValueError,
        match="at least one domain identifier",
    ):
        DomainRecord(
            protein_id="P12345",
            region_type="domain",
            start=1,
            end=5,
            sequence="MAVIL",
        )


def test_domain_requires_domain_region_type() -> None:
    with pytest.raises(
        ValueError,
        match="region_type must be 'domain'",
    ):
        DomainRecord(
            protein_id="P12345",
            region_type="motif",
            start=1,
            end=5,
            sequence="MAVIL",
            pfam_id="PF00018",
        )


def test_domain_to_dict() -> None:
    domain = DomainRecord(
        protein_id="P12345",
        region_type="domain",
        start=1,
        end=5,
        sequence="MAVIL",
        pfam_id="PF00018",
        binding_types=("domain_peptide",),
        known_ligands=("peptide",),
    )

    data = domain.to_dict()

    assert data["pfam_id"] == "PF00018"
    assert data["binding_types"] == ["domain_peptide"]
    assert data["known_ligands"] == ["peptide"]
