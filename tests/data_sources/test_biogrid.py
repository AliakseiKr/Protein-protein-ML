from pathlib import Path
from zipfile import ZipFile

import pandas as pd
import pytest

from protein_graph.data_sources.biogrid import BioGRIDDataSource


@pytest.fixture
def source(tmp_path: Path) -> BioGRIDDataSource:
    return BioGRIDDataSource(
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
            "Organism Name Interactor A": [
                "Rattus norvegicus",
                "Rattus norvegicus",
                "Rattus norvegicus",
                "Rattus norvegicus",
            ],
            "Organism Name Interactor B": [
                "Rattus norvegicus",
                "Rattus norvegicus",
                "Rattus norvegicus",
                "Homo sapiens",
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


def test_normalize_filters_genetic_interactions_only(
    source: BioGRIDDataSource,
    biogrid_dataframe: pd.DataFrame,
) -> None:
    """Genetic-type evidence must never enter the physical evidence table."""
    result = source._normalize(biogrid_dataframe)

    assert len(result) == 3
    assert set(result["biogrid_interaction_id"]) == {
        "1001",
        "1002",
        "1004",
    }
    assert set(result["interaction_type"]) == {"physical"}


def test_global_dataset_retains_cross_taxonomy_evidence(
    source: BioGRIDDataSource,
    biogrid_dataframe: pd.DataFrame,
) -> None:
    """Production (no taxonomy_id set) must not reject differing taxonomy."""
    result = source._normalize(biogrid_dataframe)

    cross_species_row = result.loc[
        result["biogrid_interaction_id"] == "1004"
    ].iloc[0]

    assert cross_species_row["organism_id_a"] == 10116
    assert cross_species_row["organism_id_b"] == 9606


def test_taxonomy_id_filter_is_opt_in_and_inclusive_of_partners(
    tmp_path: Path,
    biogrid_dataframe: pd.DataFrame,
) -> None:
    """An explicit taxonomy_id keeps rows where EITHER side matches -
    cross-species pairs involving that organism are not excluded."""
    scoped_source = BioGRIDDataSource(
        interaction_type="physical",
        taxonomy_id=9606,
        raw_dir=tmp_path,
    )

    result = scoped_source._normalize(biogrid_dataframe)

    assert set(result["biogrid_interaction_id"]) == {"1004"}


def test_same_pair_has_same_biogrid_pair_id(
    source: BioGRIDDataSource,
    biogrid_dataframe: pd.DataFrame,
) -> None:
    result = source._normalize(biogrid_dataframe)

    row_1001 = result.loc[result["biogrid_interaction_id"] == "1001"].iloc[0]
    row_1002 = result.loc[result["biogrid_interaction_id"] == "1002"].iloc[0]

    assert row_1001["biogrid_pair_id"] == "10__20"
    assert row_1002["biogrid_pair_id"] == "10__20"


def test_evidence_rows_are_not_collapsed(
    source: BioGRIDDataSource,
    biogrid_dataframe: pd.DataFrame,
) -> None:
    result = source._normalize(biogrid_dataframe)

    same_pair = result.loc[result["biogrid_pair_id"] == "10__20"]

    assert len(same_pair) == 2


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


def test_extract_release_version_from_archive_member_name() -> None:
    version = BioGRIDDataSource._extract_release_version(
        "BIOGRID-ALL-5.0.260.tab3.txt"
    )

    assert version == "5.0.260"


def test_extract_release_version_rejects_unexpected_name() -> None:
    with pytest.raises(ValueError, match="release version"):
        BioGRIDDataSource._extract_release_version(
            "BIOGRID-ORGANISM-Rat-5.0.260.tab3.txt"
        )


def test_sample_limit_must_be_positive(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="sample_limit"):
        BioGRIDDataSource(raw_dir=tmp_path, sample_limit=0)


def test_default_download_url_is_latest_release(tmp_path: Path) -> None:
    source = BioGRIDDataSource(raw_dir=tmp_path)

    assert source.release is None
    assert source.download_url == BioGRIDDataSource.LATEST_RELEASE_URL
    assert source.archive_path.name == "BIOGRID-ALL-LATEST.tab3.zip"


def test_pinned_release_resolves_to_release_archive_url(
    tmp_path: Path,
) -> None:
    source = BioGRIDDataSource(raw_dir=tmp_path, release="5.0.255")

    assert source.download_url == (
        "https://downloads.thebiogrid.org/Download/BioGRID/"
        "Release-Archive/BIOGRID-5.0.255/BIOGRID-ALL-5.0.255.tab3.zip"
    )
    assert source.archive_path.name == "BIOGRID-ALL-5.0.255.tab3.zip"


def test_pinned_release_uses_a_distinct_cache_file_from_latest(
    tmp_path: Path,
) -> None:
    """Pinning a different release must never silently reuse a cached
    LATEST archive from a previous build."""
    latest = BioGRIDDataSource(raw_dir=tmp_path)
    pinned = BioGRIDDataSource(raw_dir=tmp_path, release="5.0.255")

    assert latest.archive_path != pinned.archive_path


def test_release_must_not_be_empty(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="release"):
        BioGRIDDataSource(raw_dir=tmp_path, release="  ")


@pytest.mark.parametrize(
    "malicious_release",
    [
        "../../../etc/passwd",
        "5.0.260/../../evil",
        "a/b",
        "5.0.260\\evil",
        "5.0.260; rm -rf",
        "LATEST",
        "..",
        "/etc/passwd",
    ],
)
def test_release_rejects_path_traversal_and_unsafe_characters(
    tmp_path: Path,
    malicious_release: str,
) -> None:
    """Regression test: a release string is interpolated directly into a
    local cache file path (archive_path = raw_dir / f"BIOGRID-ALL-
    {release}...") and into a download URL. Without validation, a value
    like "../../../etc/passwd" makes archive_path resolve OUTSIDE raw_dir,
    turning a routine download into an arbitrary local file write."""
    with pytest.raises(ValueError, match="release"):
        BioGRIDDataSource(raw_dir=tmp_path, release=malicious_release)


def test_release_path_traversal_cannot_escape_raw_dir(
    tmp_path: Path,
) -> None:
    """Even if validation regressed, archive_path must never resolve
    outside raw_dir for any release string that passes validation."""
    for valid_release in ["5.0.260", "5.0.255", "4.2.192"]:
        source = BioGRIDDataSource(raw_dir=tmp_path, release=valid_release)
        resolved = source.archive_path.resolve()
        assert str(resolved).startswith(str(tmp_path.resolve()))


def test_valid_release_formats_are_accepted(tmp_path: Path) -> None:
    for valid_release in ["5.0.260", "5.0.255", "4.2.192", "5"]:
        source = BioGRIDDataSource(raw_dir=tmp_path, release=valid_release)
        assert source.release == valid_release


def test_parse_rejects_archive_that_does_not_match_pinned_release(
    tmp_path: Path,
    biogrid_dataframe: pd.DataFrame,
) -> None:
    """If a pinned release is requested, the archive that was actually
    downloaded must resolve to that same release - never silently fall
    back to whatever version the file turned out to contain."""
    archive_path = tmp_path / "BIOGRID-ALL-5.0.255.tab3.zip"

    with ZipFile(archive_path, "w") as archive:
        archive.writestr(
            "BIOGRID-ALL-5.0.260.tab3.txt",
            biogrid_dataframe.to_csv(sep="\t", index=False),
        )

    source = BioGRIDDataSource(raw_dir=tmp_path, release="5.0.255")
    source.archive_path = archive_path

    with pytest.raises(RuntimeError, match="Pinned BioGRID release"):
        source.load()


def test_repr_includes_release(tmp_path: Path) -> None:
    source = BioGRIDDataSource(raw_dir=tmp_path, release="5.0.255")

    assert "5.0.255" in repr(source)


def test_sample_mode_and_full_mode_share_the_same_normalization_path(
    tmp_path: Path,
    biogrid_dataframe: pd.DataFrame,
) -> None:
    """Sample mode must only reduce the row count read from disk - the
    normalization/qualification code path must be identical, so a sample
    build is a strict prefix of the full build, never a different
    transformation."""
    archive_path = tmp_path / "BIOGRID-ALL-LATEST.tab3.zip"

    with ZipFile(archive_path, "w") as archive:
        archive.writestr(
            "BIOGRID-ALL-5.0.260.tab3.txt",
            biogrid_dataframe.to_csv(sep="\t", index=False),
        )

    full_source = BioGRIDDataSource(raw_dir=tmp_path)
    full_source.archive_path = archive_path
    full_result = full_source.load()

    sample_source = BioGRIDDataSource(raw_dir=tmp_path, sample_limit=2)
    sample_source.archive_path = archive_path
    sample_result = sample_source.load()

    assert len(sample_result) < len(full_result)
    pd.testing.assert_frame_equal(
        sample_result,
        full_result.head(len(sample_result)),
    )


def test_repr(source: BioGRIDDataSource) -> None:
    representation = repr(source)

    assert "BioGRIDDataSource" in representation
    assert "physical" in representation
