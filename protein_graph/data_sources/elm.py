"""ELM motif definition data source."""

from __future__ import annotations

import csv
from collections.abc import Iterable
from pathlib import Path
from typing import Any, TextIO

import requests

from protein_graph.data_sources.base import BaseDataSource


class ELMDataSource(BaseDataSource):
    """Load and normalize ELM motif-class definitions.

    The data source can operate in two modes:

    1. Read an existing TSV file from ``raw_dir``.
    2. Download the TSV file from an explicitly configured URL.

    No default remote URL is fixed until an official stable export
    endpoint and its real header have been verified.
    """

    RAW_FILENAME = "elm_classes.tsv"

    _COLUMN_ALIASES: dict[str, tuple[str, ...]] = {
        "motif_id": (
            "ELMIdentifier",
            "ELM identifier",
            "elm_identifier",
            "Accession",
            "accession",
            "motif_id",
        ),
        "motif_class": (
            "FunctionalSiteName",
            "Functional Site Name",
            "functional_site_name",
            "Description",
            "description",
            "motif_class",
        ),
        "pattern": (
            "Regex",
            "regex",
            "Pattern",
            "pattern",
        ),
        "evidence_level": (
            "Evidence",
            "evidence",
            "EvidenceLevel",
            "Evidence Level",
            "evidence_level",
        ),
    }

    _REQUIRED_FIELDS = (
        "motif_id",
        "pattern",
    )

    def __init__(
        self,
        raw_dir: str | Path = "data/raw/elm",
        *,
        url: str | None = None,
        timeout: int = 30,
    ) -> None:
        super().__init__(raw_dir=raw_dir)

        if url is not None:
            url = url.strip()

            if not url:
                raise ValueError(
                    "ELM source URL cannot be empty."
                )

        if timeout <= 0:
            raise ValueError(
                "timeout must be greater than zero."
            )

        self.url = url
        self.timeout = timeout

    @property
    def raw_path(self) -> Path:
        """Return the local ELM TSV path."""
        return self.raw_dir / self.RAW_FILENAME

    def _download(self, force: bool = False) -> None:
        """Ensure that the ELM TSV file is available locally."""
        if self.raw_path.exists() and not force:
            return

        if self.url is None:
            raise RuntimeError(
                "ELM URL is not configured and no cached "
                f"file is available at {self.raw_path}."
            )

        try:
            response = requests.get(
                self.url,
                timeout=self.timeout,
                headers={
                    "Accept": (
                        "text/tab-separated-values, "
                        "text/plain;q=0.9, */*;q=0.1"
                    ),
                    "User-Agent": "ProteinGraph/0.1",
                },
            )
            response.raise_for_status()
        except requests.RequestException as error:
            raise RuntimeError(
                "Failed to download ELM motif definitions "
                f"from {self.url!r}."
            ) from error

        content = response.text

        if not content.strip():
            raise ValueError(
                "ELM returned an empty motif definition file."
            )

        self.raw_path.write_text(
            content,
            encoding="utf-8",
        )

    def _parse(self) -> list[dict[str, str]]:
        """Parse the locally cached ELM TSV file."""
        if not self.raw_path.exists():
            raise FileNotFoundError(
                f"ELM raw file does not exist: "
                f"{self.raw_path}"
            )

        with self.raw_path.open(
            "r",
            encoding="utf-8-sig",
            newline="",
        ) as handle:
            return self._parse_tsv(handle)

    def _normalize(
        self,
        records: list[dict[str, str]],
    ) -> list[dict[str, Any]]:
        """Convert ELM rows into canonical motif records."""
        normalized: list[dict[str, Any]] = []
        seen_ids: set[str] = set()

        for row_number, record in enumerate(
            records,
            start=2,
        ):
            motif_id = self._required_value(
                record,
                canonical_name="motif_id",
                row_number=row_number,
            )

            pattern = self._required_value(
                record,
                canonical_name="pattern",
                row_number=row_number,
            )

            motif_class = self._optional_value(
                record,
                canonical_name="motif_class",
            )

            evidence_level = self._optional_value(
                record,
                canonical_name="evidence_level",
            )

            if motif_id in seen_ids:
                raise ValueError(
                    "ELM data contains duplicate motif ID "
                    f"{motif_id!r}."
                )

            seen_ids.add(motif_id)

            normalized.append(
                {
                    "motif_id": motif_id,
                    "motif_class": motif_class or motif_id,
                    "pattern": pattern,
                    "recognized_by": (),
                    "interaction_role": None,
                    "source": "ELM",
                    "evidence_level": evidence_level,
                }
            )

        return normalized

    def _parse_tsv(
        self,
        handle: TextIO,
    ) -> list[dict[str, str]]:
        """Parse a TSV stream and ignore empty/comment lines."""
        reader = csv.DictReader(
            self._non_comment_lines(handle),
            delimiter="\t",
        )

        if reader.fieldnames is None:
            raise ValueError(
                "ELM TSV file does not contain a header."
            )

        fieldnames = {
            field.strip()
            for field in reader.fieldnames
            if field is not None and field.strip()
        }

        if not fieldnames:
            raise ValueError(
                "ELM TSV header is empty."
            )

        self._validate_header(fieldnames)

        records: list[dict[str, str]] = []

        for row in reader:
            cleaned = {
                str(key).strip(): (
                    value.strip()
                    if isinstance(value, str)
                    else ""
                )
                for key, value in row.items()
                if key is not None
            }

            if any(cleaned.values()):
                records.append(cleaned)

        return records

    @staticmethod
    def _non_comment_lines(
        handle: TextIO,
    ) -> Iterable[str]:
        """Yield non-empty lines that are not comments."""
        for line in handle:
            if not line.strip():
                continue

            if line.lstrip().startswith("#"):
                continue

            yield line

    def _validate_header(
        self,
        fieldnames: set[str],
    ) -> None:
        """Check that the TSV header contains required fields."""
        missing_fields = [
            canonical_name
            for canonical_name in self._REQUIRED_FIELDS
            if not self._has_alias(
                fieldnames,
                canonical_name,
            )
        ]

        if not missing_fields:
            return

        expected = {
            name: self._COLUMN_ALIASES[name]
            for name in missing_fields
        }

        raise ValueError(
            "Unsupported ELM TSV header. Missing canonical "
            f"fields {missing_fields!r}. Accepted aliases: "
            f"{expected!r}. Actual columns: "
            f"{sorted(fieldnames)!r}."
        )

    def _has_alias(
        self,
        fieldnames: set[str],
        canonical_name: str,
    ) -> bool:
        """Return whether a canonical field has a known column."""
        return any(
            alias in fieldnames
            for alias in self._COLUMN_ALIASES[
                canonical_name
            ]
        )

    def _required_value(
        self,
        record: dict[str, str],
        *,
        canonical_name: str,
        row_number: int,
    ) -> str:
        """Read a required normalized field from one ELM row."""
        value = self._optional_value(
            record,
            canonical_name=canonical_name,
        )

        if value is None:
            raise ValueError(
                f"ELM row {row_number} has an empty "
                f"{canonical_name!r} value."
            )

        return value

    def _optional_value(
        self,
        record: dict[str, str],
        *,
        canonical_name: str,
    ) -> str | None:
        """Read an optional value through its known aliases."""
        aliases = self._COLUMN_ALIASES[
            canonical_name
        ]

        for alias in aliases:
            value = record.get(alias)

            if isinstance(value, str) and value.strip():
                return value.strip()

        return None

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}("
            f"raw_dir={str(self.raw_dir)!r}, "
            f"url={self.url!r}, "
            f"timeout={self.timeout}"
            f")"
        )
