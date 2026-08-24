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

        (
            mapped["protein_a_id"],
            mapped["protein_a_mapping_source"],
        ) = self._resolve_column(
            swissprot=mapped["swissprot_accessions_a"],
            trembl=mapped["trembl_accessions_a"],
        )
        (
            mapped["protein_b_id"],
            mapped["protein_b_mapping_source"],
        ) = self._resolve_column(
            swissprot=mapped["swissprot_accessions_b"],
            trembl=mapped["trembl_accessions_b"],
        )

        mapped["mapping_complete"] = (
            mapped["protein_a_id"].notna()
            & mapped["protein_b_id"].notna()
        )

        mapped["uniprot_pair_id"] = self._build_pair_id_vectorized(
            mapped["protein_a_id"],
            mapped["protein_b_id"],
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
    def _resolve_column(
        cls,
        swissprot: pd.Series,
        trembl: pd.Series,
    ) -> tuple[pd.Series, pd.Series]:
        """Resolve one interactor side for an entire evidence table at once.

        Uses ``Series.map`` (one scalar call per cell via a fast C-level
        loop) rather than ``DataFrame.apply(axis=1)`` (which reconstructs a
        full ``Series`` for every row). This matters once evidence tables
        cover BioGRID-ALL rather than a single organism.
        """
        swissprot_accession = swissprot.map(cls._first_accession)
        trembl_accession = trembl.map(cls._first_accession)

        has_swissprot = swissprot_accession.notna()
        has_trembl_only = ~has_swissprot & trembl_accession.notna()

        accession = swissprot_accession.where(
            has_swissprot,
            trembl_accession,
        )

        source = pd.Series(
            pd.NA,
            index=swissprot.index,
            dtype="object",
        )
        source[has_swissprot] = "Swiss-Prot"
        source[has_trembl_only] = "TrEMBL"

        return accession, source

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
    def _build_pair_id_vectorized(
        protein_a_id: pd.Series,
        protein_b_id: pd.Series,
    ) -> pd.Series:
        """Build the canonical undirected UniProt pair id for a whole table.

        ``A--B`` and ``B--A`` must produce the same id, so each pair is
        ordered lexicographically before joining.
        """
        a = protein_a_id.astype("string")
        b = protein_b_id.astype("string")

        missing = a.isna() | b.isna()

        first = a.where(a <= b, b)
        second = b.where(a <= b, a)

        pair_id = first + "--" + second

        return pair_id.mask(missing, other=None).astype("object")
