from pathlib import Path
from typing import Any

from protein_graph.data_sources.base import BaseDataSource


class DummyDataSource(BaseDataSource):
    """Dummy implementation used to test BaseDataSource."""

    def __init__(self, raw_dir: str | Path) -> None:
        super().__init__(raw_dir=raw_dir)
        self.download_calls: list[bool] = []

    def _download(self, force: bool = False) -> None:
        self.download_calls.append(force)

    def _parse(self) -> list[int]:
        return [1, 2, 3]

    def _normalize(self, data: Any) -> list[int]:
        return [value * 2 for value in data]


def test_load_runs_complete_lifecycle(tmp_path: Path) -> None:
    raw_dir = tmp_path / "raw"
    source = DummyDataSource(raw_dir=raw_dir)

    result = source.load()

    assert result == [2, 4, 6]
    assert source.download_calls == [False]


def test_load_passes_force_download_flag(tmp_path: Path) -> None:
    raw_dir = tmp_path / "raw"
    source = DummyDataSource(raw_dir=raw_dir)

    source.load(force_download=True)

    assert source.download_calls == [True]


def test_load_creates_raw_directory(tmp_path: Path) -> None:
    raw_dir = tmp_path / "nested" / "raw"
    source = DummyDataSource(raw_dir=raw_dir)

    assert not raw_dir.exists()

    source.load()

    assert raw_dir.is_dir()
