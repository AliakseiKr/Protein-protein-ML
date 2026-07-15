"""Build the complete ProteinGraph training dataset."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from protein_graph.mapping import BioGRIDUniProtMapper
from protein_graph.pipelines.pair_features import PairFeaturePipeline
from protein_graph.pipelines.protein_features import ProteinFeaturePipeline


class TrainingDatasetPipeline:
    """Build a pair-level dataset from normalized BioGRID interactions."""

    def __init__(
        self,
        protein_feature_pipeline: ProteinFeaturePipeline | None = None,
        pair_feature_pipeline: PairFeaturePipeline | None = None,
        mapper: BioGRIDUniProtMapper | None = None,
    ) -> None:
        self.protein_feature_pipeline = (
            protein_feature_pipeline or ProteinFeaturePipeline()
        )
        self.pair_feature_pipeline = (
            pair_feature_pipeline or PairFeaturePipeline()
        )
        self.mapper = mapper or BioGRIDUniProtMapper()

    def build(
        self,
        interactions: pd.DataFrame,
        force_download: bool = False,
    ) -> tuple[pd.DataFrame, dict[str, Any]]:
        """Build a pair-level feature dataset."""
        mapped_interactions = self.mapper.map(interactions)
        coverage = self.mapper.coverage(mapped_interactions)

        complete_interactions = mapped_interactions.loc[
            mapped_interactions["mapping_complete"]
        ].copy()

        accessions = self._collect_accessions(complete_interactions)

        if not accessions:
            raise ValueError(
                "No completely mapped UniProt interaction pairs were found."
            )

        protein_features = self.protein_feature_pipeline.build(
            accessions=accessions,
            force_download=force_download,
        )

        pair_features = self.pair_feature_pipeline.build(
            interactions=complete_interactions,
            protein_features=protein_features,
        )

        report = {
            **coverage,
            "unique_uniprot_accessions": len(accessions),
            "protein_feature_rows": len(protein_features),
            "pair_feature_rows": len(pair_features),
        }

        return pair_features, report

    def build_and_save(
        self,
        interactions: pd.DataFrame,
        output_path: str | Path = (
            "data/processed/biogrid_rat_pair_features.parquet"
        ),
        report_path: str | Path = (
            "data/processed/biogrid_rat_pair_features_report.json"
        ),
        force_download: bool = False,
    ) -> tuple[pd.DataFrame, dict[str, Any]]:
        """Build and save the pair-level feature dataset."""
        import json

        pair_features, report = self.build(
            interactions=interactions,
            force_download=force_download,
        )

        output = Path(output_path)
        output.parent.mkdir(parents=True, exist_ok=True)
        pair_features.to_parquet(output, index=False)

        report_destination = Path(report_path)
        report_destination.parent.mkdir(parents=True, exist_ok=True)
        report_destination.write_text(
            json.dumps(report, indent=2),
            encoding="utf-8",
        )

        return pair_features, report

    @staticmethod
    def _collect_accessions(
        interactions: pd.DataFrame,
    ) -> list[str]:
        """Collect unique UniProt accessions from both pair sides."""
        accessions = pd.concat(
            [
                interactions["protein_a_id"],
                interactions["protein_b_id"],
            ],
            ignore_index=True,
        )

        accessions = (
            accessions
            .dropna()
            .astype(str)
            .str.strip()
        )

        accessions = accessions[accessions.ne("")]

        return sorted(accessions.unique().tolist())

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}()"
