from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any

import pandas as pd


class BaseDataSource(ABC):
    """Base class for all external biological data sources.

    A data source is responsible for fetching data from one external
    biological resource and returning a normalized pandas DataFrame.

    Public API:
        fetch()

    Internal lifecycle:
        _download()
        _parse()
        _normalize()
    """

    def __init__(
        self,
        raw_dir: str | Path = "data/raw",
    ) -> None:
        self.raw_dir = Path(raw_dir)
        self.raw_dir.mkdir(parents=True, exist_ok=True)

    def fetch(self) -> pd.DataFrame:
        """Fetch, parse, normalize, and return the data source table."""
        self._download()
        raw_data = self._parse()
        normalized_data = self._normalize(raw_data)
        return normalized_data

    def _download(self) -> None:
        """Download raw data if needed.

        This method is optional. Data sources that do not need downloading
        can use the default no-op implementation.
        """
        return None

    @abstractmethod
    def _parse(self) -> pd.DataFrame:
        """Parse raw source data into a pandas DataFrame."""
        raise NotImplementedError

    @abstractmethod
    def _normalize(self, data: pd.DataFrame) -> pd.DataFrame:
        """Normalize source-specific data into ProteinGraph format."""
        raise NotImplementedError

    def __repr__(self) -> str:
        params = {
            "raw_dir": str(self.raw_dir),
        }
        formatted = ", ".join(f"{key}={value!r}" for key, value in params.items())
        return f"{self.__class__.__name__}({formatted})"
