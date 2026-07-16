"""Tests for sequence region entities."""

import pytest

from protein_graph.entities import SequenceRegion


def test_region_length() -> None:
    region = SequenceRegion(
        protein_id="P12345",
        region_type="signal_peptide",
        start=1,
        end=5,
        sequence="MAVIL",
    )

    assert region.length == 5


def test_region_position_categories() -> None:
    n_terminal = SequenceRegion(
        protein_id="P1",
        region_type="domain",
        start=1,
        end=20,
        sequence="A" * 20,
    )

    middle = SequenceRegion(
        protein_id="P1",
        region_type="domain",
        start=41,
        end=60,
        sequence="A" * 20,
    )

    c_terminal = SequenceRegion(
        protein_id="P1",
        region_type="domain",
        start=81,
        end=100,
        sequence="A" * 20,
    )

    assert n_terminal.position_category(100) == "n_terminal"
    assert middle.position_category(100) == "middle"
    assert c_terminal.position_category(100) == "c_terminal"


def test_region_rejects_invalid_coordinates() -> None:
    with pytest.raises(
        ValueError,
        match="start must be at least 1",
    ):
        SequenceRegion(
            protein_id="P1",
            region_type="domain",
            start=0,
            end=5,
            sequence="AAAAAA",
        )


def test_region_rejects_sequence_length_mismatch() -> None:
    with pytest.raises(
        ValueError,
        match="sequence length does not match",
    ):
        SequenceRegion(
            protein_id="P1",
            region_type="domain",
            start=1,
            end=5,
            sequence="AAA",
        )


def test_region_to_dict() -> None:
    region = SequenceRegion(
        protein_id="P1",
        region_type="transmembrane",
        start=10,
        end=14,
        sequence="AVILM",
        source="UniProt",
    )

    data = region.to_dict()

    assert data["length"] == 5
    assert data["source"] == "UniProt"
