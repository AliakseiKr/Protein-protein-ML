"""Tests for the UniProt data source."""

import json

import pytest

from protein_graph.data_sources.uniprot import UniProtDataSource


@pytest.fixture
def uniprot_record() -> dict:
    return {
        "primaryAccession": "P12345",
        "uniProtkbId": "TEST_RAT",
        "entryType": "UniProtKB reviewed (Swiss-Prot)",
        "proteinDescription": {
            "recommendedName": {
                "fullName": {
                    "value": "Test protein",
                }
            }
        },
        "genes": [
            {
                "geneName": {
                    "value": "Test1",
                }
            }
        ],
        "organism": {
            "scientificName": "Rattus norvegicus",
            "taxonId": 10116,
        },
        "sequence": {
            "value": "MAVILMFWYKDE",
        },
        "comments": [
            {
                "commentType": "SUBCELLULAR LOCATION",
                "subcellularLocations": [
                    {
                        "location": {
                            "value": "Cell membrane",
                        }
                    }
                ],
            }
        ],
        "features": [
            {
                "type": "Signal",
                "location": {
                    "start": {"value": 1},
                    "end": {"value": 3},
                },
            },
            {
                "type": "Transmembrane",
                "location": {
                    "start": {"value": 4},
                    "end": {"value": 8},
                },
            },
            {
                "type": "Domain",
                "description": "Example domain",
                "location": {
                    "start": {"value": 4},
                    "end": {"value": 10},
                },
            },
        ],
        "uniProtKBCrossReferences": [
            {
                "database": "InterPro",
                "id": "IPR000001",
            },
            {
                "database": "Pfam",
                "id": "PF00001",
            },
            {
                "database": "RefSeq",
                "id": "NP_000001.1",
            },
        ],
    }


def test_normalize_record(
    tmp_path,
    uniprot_record,
) -> None:
    source = UniProtDataSource(
        accessions=["P12345"],
        raw_dir=tmp_path,
    )

    protein = source._normalize_record(uniprot_record)

    assert protein["protein_id"] == "P12345"
    assert protein["gene_name"] == "Test1"
    assert protein["protein_name"] == "Test protein"
    assert protein["taxon_id"] == 10116
    assert protein["protein_length"] == 12
    assert protein["reviewed"] is True

    assert protein["localization_terms"] == ["Cell membrane"]
    assert protein["domains"] == ["Example domain"]
    assert protein["domain_count"] == 1

    assert protein["has_signal_peptide"] is True
    assert protein["signal_peptide_length"] == 3

    assert protein["has_transmembrane_region"] is True
    assert protein["transmembrane_region_count"] == 1

    assert protein["interpro_ids"] == ["IPR000001"]
    assert protein["pfam_ids"] == ["PF00001"]
    assert protein["refseq_ids"] == ["NP_000001.1"]


def test_load_from_cached_record(
    tmp_path,
    uniprot_record,
) -> None:
    path = tmp_path / "P12345.json"
    path.write_text(
        json.dumps(uniprot_record),
        encoding="utf-8",
    )

    source = UniProtDataSource(
        accessions=["P12345"],
        raw_dir=tmp_path,
    )

    proteins = source.load()

    assert len(proteins) == 1
    assert proteins.iloc[0]["protein_id"] == "P12345"


def test_accessions_are_normalized_and_deduplicated(
    tmp_path,
) -> None:
    source = UniProtDataSource(
        accessions=[" p12345 ", "P12345", "Q99999"],
        raw_dir=tmp_path,
    )

    assert source.accessions == ["P12345", "Q99999"]


def test_empty_accessions_are_rejected(tmp_path) -> None:
    with pytest.raises(
        ValueError,
        match="At least one UniProt accession",
    ):
        UniProtDataSource(
            accessions=[],
            raw_dir=tmp_path,
        )


def test_missing_accession_is_rejected(tmp_path) -> None:
    source = UniProtDataSource(
        accessions=["P12345"],
        raw_dir=tmp_path,
    )

    with pytest.raises(
        ValueError,
        match="missing primaryAccession",
    ):
        source._normalize_record({
            "sequence": {
                "value": "MAVIL",
            }
        })


def test_missing_sequence_is_rejected(tmp_path) -> None:
    source = UniProtDataSource(
        accessions=["P12345"],
        raw_dir=tmp_path,
    )

    with pytest.raises(
        ValueError,
        match="has no sequence",
    ):
        source._normalize_record({
            "primaryAccession": "P12345",
        })


def test_repr(tmp_path) -> None:
    source = UniProtDataSource(
        accessions=["P12345", "Q99999"],
        raw_dir=tmp_path,
    )

    representation = repr(source)

    assert "UniProtDataSource" in representation
    assert "accessions=2" in representation

def test_iter_batches() -> None:
    batches = list(
        UniProtDataSource._iter_batches(
            ["P1", "P2", "P3", "P4", "P5"],
            size=2,
        )
    )

    assert batches == [
        ["P1", "P2"],
        ["P3", "P4"],
        ["P5"],
    ]


def test_download_uses_batch_request(
    monkeypatch,
    tmp_path,
    uniprot_record,
) -> None:
    source = UniProtDataSource(
        accessions=["P12345"],
        raw_dir=tmp_path,
    )

    calls = []

    def fake_download_batch(accessions):
        calls.append(accessions)
        return [uniprot_record]

    monkeypatch.setattr(
        source,
        "_download_batch",
        fake_download_batch,
    )

    source._download()

    assert calls == [["P12345"]]
    assert (tmp_path / "P12345.json").exists()


def test_download_skips_cached_records(
    monkeypatch,
    tmp_path,
    uniprot_record,
) -> None:
    cached_path = tmp_path / "P12345.json"
    cached_path.write_text(
        json.dumps(uniprot_record),
        encoding="utf-8",
    )

    source = UniProtDataSource(
        accessions=["P12345"],
        raw_dir=tmp_path,
    )

    def fail_if_called(accessions):
        raise AssertionError("Batch download should not be called.")

    monkeypatch.setattr(
        source,
        "_download_batch",
        fail_if_called,
    )

    source._download()




def test_download_records_unresolved_accessions(
    monkeypatch,
    tmp_path,
) -> None:
    source = UniProtDataSource(
        accessions=["P12345", "Q99999"],
        raw_dir=tmp_path,
    )

    monkeypatch.setattr(
        source,
        "_download_batch",
        lambda accessions: [
            {
                "primaryAccession": "P12345",
            }
        ],
    )

    source._download()

    assert source.unresolved_accessions == ["Q99999"]
    assert (tmp_path / "P12345.json").exists()
    assert not (tmp_path / "Q99999.json").exists()


def test_parse_skips_unresolved_accessions(
    tmp_path,
    uniprot_record,
) -> None:
    path = tmp_path / "P12345.json"
    path.write_text(
        json.dumps(uniprot_record),
        encoding="utf-8",
    )

    source = UniProtDataSource(
        accessions=["P12345", "Q99999"],
        raw_dir=tmp_path,
    )

    source.unresolved_accessions = ["Q99999"]

    records = source._parse()

    assert len(records) == 1
    assert records[0]["primaryAccession"] == "P12345"
