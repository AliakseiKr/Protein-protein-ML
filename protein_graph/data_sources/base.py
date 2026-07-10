from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any


class BaseDataSource(ABC):
    """Base class for external biological data sources.

    A data source follows a three-stage lifecycle:

    1. Download or verify raw source data.
    2. Parse source-specific raw data.
    3. Normalize parsed data into a project-friendly format.

    The public entry point is :meth:`load`.
    """

    def __init__(self, raw_dir: str | Path) -> None:
        self.raw_dir = Path(raw_dir)

    @abstractmethod
    def _download(self, force: bool = False) -> None:
        """Ensure that the required raw data is available locally.

        Parameters
        ----------
        force:
            Download the source again even when a cached raw file exists.
        """

    @abstractmethod
    def _parse(self) -> Any:
        """Parse locally stored source-specific data."""

    @abstractmethod
    def _normalize(self, data: Any) -> Any:
        """Normalize parsed source-specific data."""

    def load(self, force_download: bool = False) -> Any:
        """Download, parse, and normalize data.

        Parameters
        ----------
        force_download:
            Force a fresh download even when cached raw data exists.

        Returns
        -------
        Any
            Normalized source data.
        """
        self.raw_dir.mkdir(parents=True, exist_ok=True)

        self._download(force=force_download)
        parsed_data = self._parse()

        return self._normalize(parsed_data)
