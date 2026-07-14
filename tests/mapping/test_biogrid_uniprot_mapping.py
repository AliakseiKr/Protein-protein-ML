"""Tests for BioGRID to UniProt mapping."""

import pandas as pd
import pytest

from protein_graph.mapping import BioGRIDUniProtMapper


def test_prefers_swissprot_over_trembl() -> None:
    interactions = pd.DataFrame(
        {
            "swissprot_accessions_a": ["P12345"],
            "swissprot_accessions_b": ["Q11111"],
            "trembl_accessions_a": ["A0A000"],
            "trembl_accessions_b": ["A0A111"],
        }
    )

    mapper = BioGRIDUniProtMapper()
    mapped = mapper.map(interactions)

    assert mapped.iloc[0]["protein_a_id"] == "P12345"
    assert mapped.iloc[0]["protein_b_id"] == "Q11111"
    assert mapped.iloc[0]["protein_a_mapping_source"] == "Swiss-Prot"


def test_falls_back_to_trembl() -> None:
    interactions = pd.DataFrame(
        {
            "swissprot_accessions_a": [None],
            "swissprot_accessions_b": ["-"],
            "trembl_accessions_a": ["A0A000"],
            "trembl_accessions_b": ["A0A111"],
        }
    )

    mapped = BioGRIDUniProtMapper().map(interactions)

    assert mapped.iloc[0]["protein_a_id"] == "A0A000"
    assert mapped.iloc[0]["protein_b_id"] == "A0A111"
    assert mapped.iloc[0]["mapping_complete"]


def test_uses_first_accession() -> None:
    interactions = pd.DataFrame(
        {
            "swissprot_accessions_a": ["P12345|Q99999"],
            "swissprot_accessions_b": ["Q11111;Q22222"],
            "trembl_accessions_a": [None],
            "trembl_accessions_b": [None],
        }
    )

    mapped = BioGRIDUniProtMapper().map(interactions)

    assert mapped.iloc[0]["protein_a_id"] == "P12345"
    assert mapped.iloc[0]["protein_b_id"] == "Q11111"


def test_builds_undirected_pair_id() -> None:
    interactions = pd.DataFrame(
        {
            "swissprot_accessions_a": ["Q11111", "P12345"],
            "swissprot_accessions_b": ["P12345", "Q11111"],
            "trembl_accessions_a": [None, None],
            "trembl_accessions_b": [None, None],
        }
    )

    mapped = BioGRIDUniProtMapper().map(interactions)

    assert (
        mapped.iloc[0]["uniprot_pair_id"]
        == mapped.iloc[1]["uniprot_pair_id"]
    )


def test_coverage_statistics() -> None:
    interactions = pd.DataFrame(
        {
            "swissprot_accessions_a": ["P1", "P2", None],
            "swissprot_accessions_b": ["Q1", None, None],
            "trembl_accessions_a": [None, None, None],
            "trembl_accessions_b": [None, None, None],
        }
    )

    mapper = BioGRIDUniProtMapper()
    mapped = mapper.map(interactions)
    coverage = mapper.coverage(mapped)

    assert coverage["total_rows"] == 3
    assert coverage["fully_mapped_rows"] == 1
    assert coverage["partially_mapped_rows"] == 1
    assert coverage["unmapped_rows"] == 1


def test_missing_columns_are_rejected() -> None:
    with pytest.raises(ValueError, match="Missing BioGRID mapping columns"):
        BioGRIDUniProtMapper().map(pd.DataFrame())
