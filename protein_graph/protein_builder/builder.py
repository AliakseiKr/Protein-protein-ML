"""Canonical protein record builder."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from protein_graph.data_sources import UniProtDataSource


class ProteinBuilder:
    """Build canonical protein records from external data sources."""

    def __init__(
        self,
        raw_dir: str | Path = "data/raw/uniprot",
        timeout: int = 30,
    ) -> None:
        self.raw_dir = Path(raw_dir)
        self.timeout = timeout

    def build(
        self,
        accession: str,
        force_download: bool = False,
    ) -> dict[str, Any]:
        """Build one canonical protein record."""
        proteins = self.build_many(
            accessions=[accession],
            force_download=force_download,
        )

        if len(proteins) != 1:
            raise RuntimeError(
                f"Expected one protein record, received {len(proteins)}."
            )

        return proteins.iloc[0].to_dict()

    def build_many(
        self,
        accessions: list[str],
        force_download: bool = False,
    ) -> pd.DataFrame:
        """Build canonical records for multiple proteins."""
        source = UniProtDataSource(
            accessions=accessions,
            raw_dir=self.raw_dir,
            timeout=self.timeout,
        )

        proteins = source.load(
            force_download=force_download,
        )

        self._validate_table(proteins)

        return proteins

    @staticmethod
    def _validate_table(proteins: pd.DataFrame) -> None:
        """Validate the canonical protein table."""
        required_columns = {
            "protein_id",
            "sequence",
            "protein_length",
            "organism",
            "taxon_id",
        }

        missing_columns = required_columns - set(proteins.columns)

        if missing_columns:
            raise ValueError(
                "Canonical protein table is missing columns: "
                f"{sorted(missing_columns)}"
            )

        if proteins["protein_id"].isna().any():
            raise ValueError(
                "Canonical protein table contains missing protein IDs."
            )

        if proteins["protein_id"].duplicated().any():
            raise ValueError(
                "Canonical protein table contains duplicate protein IDs."
            )

        if proteins["sequence"].isna().any():
            raise ValueError(
                "Canonical protein table contains missing sequences."
            )

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}("
            f"raw_dir='{self.raw_dir}', "
            f"timeout={self.timeout})"
        )
