"""Map BioGRID interaction records to UniProt accessions."""

from __future__ import annotations

from typing import Any

import pandas as pd


class BioGRIDUniProtMapper:
    """Resolve BioGRID interactors to canonical UniProt accessions."""

    def map(self, interactions: pd.DataFrame) -> pd.DataFrame:
        """Add canonical UniProt identifiers to BioGRID records."""
        required_columns = {
            "swissprot_accessions_a",
            "swissprot_accessions_b",
            "trembl_accessions_a",
            "trembl_accessions_b",
        }

        missing_columns = required_columns - set(interactions.columns)

        if missing_columns:
            raise ValueError(
                f"Missing BioGRID mapping columns: {sorted(missing_columns)}"
            )

        mapped = interactions.copy()

        protein_a = mapped.apply(
            lambda row: self._resolve_interactor(
                swissprot=row["swissprot_accessions_a"],
                trembl=row["trembl_accessions_a"],
            ),
            axis=1,
        )

        protein_b = mapped.apply(
            lambda row: self._resolve_interactor(
                swissprot=row["swissprot_accessions_b"],
                trembl=row["trembl_accessions_b"],
            ),
            axis=1,
        )

        mapped["protein_a_id"] = [
            result["accession"]
            for result in protein_a
        ]
        mapped["protein_b_id"] = [
            result["accession"]
            for result in protein_b
        ]

        mapped["protein_a_mapping_source"] = [
            result["source"]
            for result in protein_a
        ]
        mapped["protein_b_mapping_source"] = [
            result["source"]
            for result in protein_b
        ]

        mapped["mapping_complete"] = (
            mapped["protein_a_id"].notna()
            & mapped["protein_b_id"].notna()
        )

        mapped["uniprot_pair_id"] = mapped.apply(
            self._build_pair_id,
            axis=1,
        )

        return mapped

    def coverage(self, mapped: pd.DataFrame) -> dict[str, int | float]:
        """Calculate UniProt mapping coverage statistics."""
        required_columns = {
            "protein_a_id",
            "protein_b_id",
            "mapping_complete",
        }

        missing_columns = required_columns - set(mapped.columns)

        if missing_columns:
            raise ValueError(
                f"Missing mapped columns: {sorted(missing_columns)}"
            )

        total = len(mapped)

        fully_mapped = int(mapped["mapping_complete"].sum())

        partially_mapped = int(
            (
                mapped["protein_a_id"].notna()
                ^ mapped["protein_b_id"].notna()
            ).sum()
        )

        unmapped = total - fully_mapped - partially_mapped

        coverage_percentage = (
            fully_mapped / total * 100
            if total
            else 0.0
        )

        return {
            "total_rows": total,
            "fully_mapped_rows": fully_mapped,
            "partially_mapped_rows": partially_mapped,
            "unmapped_rows": unmapped,
            "fully_mapped_percentage": coverage_percentage,
        }

    @classmethod
    def _resolve_interactor(
        cls,
        swissprot: Any,
        trembl: Any,
    ) -> dict[str, str | None]:
        swissprot_accession = cls._first_accession(swissprot)

        if swissprot_accession:
            return {
                "accession": swissprot_accession,
                "source": "Swiss-Prot",
            }

        trembl_accession = cls._first_accession(trembl)

        if trembl_accession:
            return {
                "accession": trembl_accession,
                "source": "TrEMBL",
            }

        return {
            "accession": None,
            "source": None,
        }

    @staticmethod
    def _first_accession(value: Any) -> str | None:
        if value is None or pd.isna(value):
            return None

        text = str(value).strip()

        if not text or text == "-":
            return None

        for separator in ("|", ";", ","):
            text = text.replace(separator, " ")

        accessions = [
            item.strip().upper()
            for item in text.split()
            if item.strip()
        ]

        return accessions[0] if accessions else None

    @staticmethod
    def _build_pair_id(row: pd.Series) -> str | None:
        protein_a = row["protein_a_id"]
        protein_b = row["protein_b_id"]

        if pd.isna(protein_a) or pd.isna(protein_b):
            return None

        first, second = sorted(
            [str(protein_a), str(protein_b)]
        )

        return f"{first}--{second}"
