"""Tests for the canonical protein builder."""

import pandas as pd
import pytest

from protein_graph.protein_builder import ProteinBuilder


def test_build_returns_one_record(monkeypatch, tmp_path) -> None:
    expected = pd.DataFrame(
        [
            {
                "protein_id": "P12345",
                "sequence": "MAVIL",
                "protein_length": 5,
                "organism": "Rattus norvegicus",
                "taxon_id": 10116,
            }
        ]
    )

    def fake_build_many(
        self,
        accessions,
        force_download=False,
    ):
        assert accessions == ["P12345"]
        assert force_download is True
        return expected

    monkeypatch.setattr(
        ProteinBuilder,
        "build_many",
        fake_build_many,
    )

    builder = ProteinBuilder(raw_dir=tmp_path)

    protein = builder.build(
        "P12345",
        force_download=True,
    )

    assert protein["protein_id"] == "P12345"
    assert protein["protein_length"] == 5


def test_validate_accepts_valid_table() -> None:
    proteins = pd.DataFrame(
        [
            {
                "protein_id": "P12345",
                "sequence": "MAVIL",
                "protein_length": 5,
                "organism": "Rattus norvegicus",
                "taxon_id": 10116,
            }
        ]
    )

    ProteinBuilder._validate_table(proteins)


def test_validate_rejects_missing_columns() -> None:
    proteins = pd.DataFrame(
        [
            {
                "protein_id": "P12345",
            }
        ]
    )

    with pytest.raises(
        ValueError,
        match="missing columns",
    ):
        ProteinBuilder._validate_table(proteins)


def test_validate_rejects_duplicate_ids() -> None:
    proteins = pd.DataFrame(
        [
            {
                "protein_id": "P12345",
                "sequence": "MAVIL",
                "protein_length": 5,
                "organism": "Rattus norvegicus",
                "taxon_id": 10116,
            },
            {
                "protein_id": "P12345",
                "sequence": "MAVIL",
                "protein_length": 5,
                "organism": "Rattus norvegicus",
                "taxon_id": 10116,
            },
        ]
    )

    with pytest.raises(
        ValueError,
        match="duplicate protein IDs",
    ):
        ProteinBuilder._validate_table(proteins)


def test_validate_rejects_missing_sequences() -> None:
    proteins = pd.DataFrame(
        [
            {
                "protein_id": "P12345",
                "sequence": None,
                "protein_length": 0,
                "organism": "Rattus norvegicus",
                "taxon_id": 10116,
            }
        ]
    )

    with pytest.raises(
        ValueError,
        match="missing sequences",
    ):
        ProteinBuilder._validate_table(proteins)


def test_repr(tmp_path) -> None:
    builder = ProteinBuilder(
        raw_dir=tmp_path,
        timeout=45,
    )

    representation = repr(builder)

    assert "ProteinBuilder" in representation
    assert "timeout=45" in representation
