from pathlib import Path

import pandas as pd
import pytest

from protein_graph.data_sources.biogrid import BioGRIDDataSource


@pytest.fixture
def source(tmp_path: Path) -> BioGRIDDataSource:
    return BioGRIDDataSource(
        species="Rattus norvegicus",
        taxonomy_id=10116,
        interaction_type="physical",
        raw_dir=tmp_path,
    )


@pytest.fixture
def biogrid_dataframe() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "#BioGRID Interaction ID": [
                "1001",
                "1002",
                "1003",
                "1004",
            ],
            "Experimental System Type": [
                "physical",
                "physical",
                "genetic",
                "physical",
            ],
            "Experimental System": [
                "Two-hybrid",
                "Affinity Capture-MS",
                "Synthetic Lethality",
                "Two-hybrid",
            ],
            "Official Symbol Interactor A": [
                "GeneA",
                "GeneB",
                "GeneA",
                "GeneA",
            ],
            "Official Symbol Interactor B": [
                "GeneB",
                "GeneA",
                "GeneC",
                "HumanGene",
            ],
            "BioGRID ID Interactor A": [
                "10",
                "20",
                "10",
                "10",
            ],
            "BioGRID ID Interactor B": [
                "20",
                "10",
                "30",
                "40",
            ],
            "Entrez Gene Interactor A": [
                "101",
                "102",
                "101",
                "101",
            ],
            "Entrez Gene Interactor B": [
                "102",
                "101",
                "103",
                "999",
            ],
            "Organism ID Interactor A": [
                "10116",
                "10116",
                "10116",
                "10116",
            ],
            "Organism ID Interactor B": [
                "10116",
                "10116",
                "10116",
                "9606",
            ],
            "Publication Source": [
                "PUBMED:1",
                "PUBMED:2",
                "PUBMED:3",
                "PUBMED:4",
            ],
            "Source Database": [
                "BIOGRID",
                "BIOGRID",
                "BIOGRID",
                "BIOGRID",
            ],
            "Throughput": [
                "Low Throughput",
                "High Throughput",
                "Low Throughput",
                "Low Throughput",
            ],
            "Score": [
                "-",
                "0.9",
                "-",
                "-",
            ],
            "SWISS-PROT Accessions Interactor A": [
                "P00001",
                "P00002",
                "P00001",
                "P00001",
            ],
            "SWISS-PROT Accessions Interactor B": [
                "P00002",
                "P00001",
                "P00003",
                "Q99999",
            ],
            "TREMBL Accessions Interactor A": [
                "-",
                "-",
                "-",
                "-",
            ],
            "TREMBL Accessions Interactor B": [
                "-",
                "-",
                "-",
                "-",
            ],
            "REFSEQ Accessions Interactor A": [
                "NP_001",
                "NP_002",
                "NP_001",
                "NP_001",
            ],
            "REFSEQ Accessions Interactor B": [
                "NP_002",
                "NP_001",
                "NP_003",
                "NP_999",
            ],
        },
        dtype="string",
    )


def test_normalize_filters_species_and_interaction_type(
    source: BioGRIDDataSource,
    biogrid_dataframe: pd.DataFrame,
) -> None:
    result = source._normalize(biogrid_dataframe)

    assert len(result) == 2
    assert set(result["biogrid_interaction_id"]) == {
        "1001",
        "1002",
    }


def test_same_pair_has_same_biogrid_pair_id(
    source: BioGRIDDataSource,
    biogrid_dataframe: pd.DataFrame,
) -> None:
    result = source._normalize(biogrid_dataframe)

    assert result.iloc[0]["biogrid_pair_id"] == "10__20"
    assert result.iloc[1]["biogrid_pair_id"] == "10__20"


def test_evidence_rows_are_not_collapsed(
    source: BioGRIDDataSource,
    biogrid_dataframe: pd.DataFrame,
) -> None:
    result = source._normalize(biogrid_dataframe)

    assert len(result) == 2
    assert result["biogrid_pair_id"].nunique() == 1


def test_optional_columns_can_be_missing(
    source: BioGRIDDataSource,
    biogrid_dataframe: pd.DataFrame,
) -> None:
    dataframe = biogrid_dataframe.drop(
        columns=[
            "TREMBL Accessions Interactor A",
            "TREMBL Accessions Interactor B",
        ]
    )

    result = source._normalize(dataframe)

    assert result["trembl_accessions_a"].isna().all()
    assert result["trembl_accessions_b"].isna().all()


def test_validate_required_columns(
    source: BioGRIDDataSource,
    biogrid_dataframe: pd.DataFrame,
) -> None:
    dataframe = biogrid_dataframe.drop(
        columns=["Experimental System Type"]
    )

    with pytest.raises(
        KeyError,
        match="Experimental System Type",
    ):
        source._validate_required_columns(dataframe)


def test_repr(source: BioGRIDDataSource) -> None:
    representation = repr(source)

    assert "BioGRIDDataSource" in representation
    assert "Rattus norvegicus" in representation
    assert "10116" in representation
    assert "physical" in representation
