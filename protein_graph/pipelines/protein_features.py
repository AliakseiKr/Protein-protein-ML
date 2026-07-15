"""Build protein-level feature tables."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from protein_graph.features.sequence import (
    AminoAcidCompositionExtractor,
    ConsecutiveRunExtractor,
    LocalDensityExtractor,
    PeriodicityExtractor,
    ResidueFractionExtractor,
)
from protein_graph.protein_builder import ProteinBuilder


class ProteinFeaturePipeline:
    """Build ML-ready features for UniProt proteins."""

    def __init__(
        self,
        protein_builder: ProteinBuilder | None = None,
    ) -> None:
        self.protein_builder = protein_builder or ProteinBuilder()

        self.extractors = [
            AminoAcidCompositionExtractor(),
            ResidueFractionExtractor(),
            LocalDensityExtractor(),
            ConsecutiveRunExtractor(),
            PeriodicityExtractor(),
        ]

    def build(
        self,
        accessions: list[str],
        force_download: bool = False,
    ) -> pd.DataFrame:
        """Download proteins and calculate sequence features."""
        proteins = self.protein_builder.build_many(
            accessions=accessions,
            force_download=force_download,
        )

        feature_rows = []

        for protein in proteins.to_dict(orient="records"):
            feature_rows.append(self._extract_protein_features(protein))

        features = pd.DataFrame(feature_rows)

        self._validate_table(features)

        return features

    def save(
        self,
        features: pd.DataFrame,
        path: str | Path = "data/processed/protein_features.parquet",
    ) -> Path:
        """Save the feature table as Parquet."""
        destination = Path(path)
        destination.parent.mkdir(parents=True, exist_ok=True)

        features.to_parquet(destination, index=False)

        return destination

    def build_and_save(
        self,
        accessions: list[str],
        path: str | Path = "data/processed/protein_features.parquet",
        force_download: bool = False,
    ) -> pd.DataFrame:
        """Build and save protein features."""
        features = self.build(
            accessions=accessions,
            force_download=force_download,
        )

        self.save(features, path)

        return features

    def _extract_protein_features(
        self,
        protein: dict[str, Any],
    ) -> dict[str, Any]:
        features: dict[str, Any] = {
            "protein_id": protein["protein_id"],
            "gene_name": protein.get("gene_name"),
            "organism": protein.get("organism"),
            "taxon_id": protein.get("taxon_id"),
        }

        for extractor in self.extractors:
            extracted = extractor.extract(protein)

            duplicated = set(features) & set(extracted)

            if duplicated:
                raise ValueError(
                    "Feature name collision: "
                    f"{sorted(duplicated)}"
                )

            features.update(extracted)

        return features

    @staticmethod
    def _validate_table(features: pd.DataFrame) -> None:
        required_columns = {
            "protein_id",
            "protein_length",
        }

        missing = required_columns - set(features.columns)

        if missing:
            raise ValueError(
                f"Protein feature table is missing columns: {sorted(missing)}"
            )

        if features["protein_id"].duplicated().any():
            raise ValueError(
                "Protein feature table contains duplicate protein IDs."
            )

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}("
            f"extractors={len(self.extractors)})"
        )
