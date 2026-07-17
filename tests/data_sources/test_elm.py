from pathlib import Path
from unittest.mock import Mock

import pytest
import requests

from protein_graph.data_sources.elm import ELMDataSource


ELM_TSV = (
    "# ELM motif classes\n"
    "Accession\tFunctionalSiteName\tRegex\tEvidence\n"
    "LIG_SH3_1\tSH3-binding motif\tP..P\tcurated\n"
    "LIG_PDZ_1\tPDZ-binding motif\t"
    "...[ST].[VIL]$\t\n"
)


def write_cached_file(
    source: ELMDataSource,
    content: str = ELM_TSV,
) -> None:
    source.raw_dir.mkdir(
        parents=True,
        exist_ok=True,
    )
    source.raw_path.write_text(
        content,
        encoding="utf-8",
    )


def test_loads_cached_elm_classes(
    tmp_path: Path,
) -> None:
    source = ELMDataSource(
        raw_dir=tmp_path / "elm",
    )
    write_cached_file(source)

    records = source.load()

    assert records == [
        {
            "motif_id": "LIG_SH3_1",
            "motif_class": "SH3-binding motif",
            "pattern": "P..P",
            "recognized_by": (),
            "interaction_role": None,
            "source": "ELM",
            "evidence_level": "curated",
        },
        {
            "motif_id": "LIG_PDZ_1",
            "motif_class": "PDZ-binding motif",
            "pattern": "...[ST].[VIL]$",
            "recognized_by": (),
            "interaction_role": None,
            "source": "ELM",
            "evidence_level": None,
        },
    ]


def test_cached_file_skips_download(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source = ELMDataSource(
        raw_dir=tmp_path / "elm",
    )
    write_cached_file(source)

    get = Mock()
    monkeypatch.setattr(requests, "get", get)

    source.load()

    get.assert_not_called()


def test_missing_file_without_url_is_rejected(
    tmp_path: Path,
) -> None:
    source = ELMDataSource(
        raw_dir=tmp_path / "elm",
    )

    with pytest.raises(
        RuntimeError,
        match="ELM URL is not configured",
    ):
        source.load()


def test_downloads_configured_elm_file(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    response = Mock()
    response.text = ELM_TSV
    response.raise_for_status.return_value = None

    get = Mock(return_value=response)
    monkeypatch.setattr(requests, "get", get)

    source = ELMDataSource(
        raw_dir=tmp_path / "elm",
        url="https://example.test/elm.tsv",
    )

    records = source.load()

    assert len(records) == 2
    assert source.raw_path.read_text(
        encoding="utf-8",
    ) == ELM_TSV

    get.assert_called_once_with(
        "https://example.test/elm.tsv",
        timeout=30,
        headers={
            "Accept": (
                "text/tab-separated-values, "
                "text/plain;q=0.9, */*;q=0.1"
            ),
            "User-Agent": "ProteinGraph/0.1",
        },
    )


def test_force_download_replaces_cached_file(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    response = Mock()
    response.text = ELM_TSV
    response.raise_for_status.return_value = None

    get = Mock(return_value=response)
    monkeypatch.setattr(requests, "get", get)

    source = ELMDataSource(
        raw_dir=tmp_path / "elm",
        url="https://example.test/elm.tsv",
    )
    write_cached_file(
        source,
        content=(
            "Accession\tRegex\n"
            "OLD\tAAA\n"
        ),
    )

    records = source.load(
        force_download=True,
    )

    assert len(records) == 2
    assert records[0]["motif_id"] == "LIG_SH3_1"
    get.assert_called_once()


def test_rejects_unknown_header(
    tmp_path: Path,
) -> None:
    source = ELMDataSource(
        raw_dir=tmp_path / "elm",
    )
    write_cached_file(
        source,
        content=(
            "Foo\tBar\tBaz\n"
            "a\tb\tc\n"
        ),
    )

    with pytest.raises(
        ValueError,
        match="Unsupported ELM TSV header",
    ):
        source.load()


def test_motif_class_falls_back_to_id(
    tmp_path: Path,
) -> None:
    source = ELMDataSource(
        raw_dir=tmp_path / "elm",
    )
    write_cached_file(
        source,
        content=(
            "Accession\tRegex\n"
            "LIG_TEST\tA..A\n"
        ),
    )

    records = source.load()

    assert records[0]["motif_class"] == "LIG_TEST"


def test_rejects_empty_required_value(
    tmp_path: Path,
) -> None:
    source = ELMDataSource(
        raw_dir=tmp_path / "elm",
    )
    write_cached_file(
        source,
        content=(
            "Accession\tRegex\n"
            "LIG_TEST\t\n"
        ),
    )

    with pytest.raises(
        ValueError,
        match="empty 'pattern'",
    ):
        source.load()


def test_rejects_duplicate_motif_ids(
    tmp_path: Path,
) -> None:
    source = ELMDataSource(
        raw_dir=tmp_path / "elm",
    )
    write_cached_file(
        source,
        content=(
            "Accession\tRegex\n"
            "LIG_TEST\tA..A\n"
            "LIG_TEST\tB..B\n"
        ),
    )

    with pytest.raises(
        ValueError,
        match="duplicate motif ID",
    ):
        source.load()


def test_download_error_is_wrapped(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    get = Mock(
        side_effect=requests.RequestException(
            "network error"
        )
    )
    monkeypatch.setattr(requests, "get", get)

    source = ELMDataSource(
        raw_dir=tmp_path / "elm",
        url="https://example.test/elm.tsv",
    )

    with pytest.raises(
        RuntimeError,
        match="Failed to download ELM",
    ):
        source.load()


def test_rejects_empty_download(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    response = Mock()
    response.text = " \n"
    response.raise_for_status.return_value = None

    monkeypatch.setattr(
        requests,
        "get",
        Mock(return_value=response),
    )

    source = ELMDataSource(
        raw_dir=tmp_path / "elm",
        url="https://example.test/elm.tsv",
    )

    with pytest.raises(
        ValueError,
        match="empty motif definition file",
    ):
        source.load()


def test_rejects_non_positive_timeout(
    tmp_path: Path,
) -> None:
    with pytest.raises(
        ValueError,
        match="timeout must be greater than zero",
    ):
        ELMDataSource(
            raw_dir=tmp_path,
            timeout=0,
        )


def test_repr(tmp_path: Path) -> None:
    source = ELMDataSource(raw_dir=tmp_path)

    representation = repr(source)

    assert "ELMDataSource" in representation
    assert "url=None" in representation
    assert "timeout=30" in representation
