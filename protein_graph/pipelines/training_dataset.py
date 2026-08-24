"""Build the complete ProteinGraph training dataset.

Production flow (see also
:class:`protein_graph.pipelines.biogrid_pair_acceptance.BioGRIDPairAcceptancePipeline`)::

    BIOGRID-ALL evidence
            |
    evidence-type filtering (RNA rejected)
            |
    UniProt normalization
            |
    canonical undirected grouping
            |
    strong/direct evidence qualification
            |
    one accepted positive pair record
            |
    protein feature extraction
            |
    pair feature extraction
            |
    explicit model-feature selection (PairFeaturePipeline.select_model_features)

Raw BioGRID evidence rows are never independent training examples: the
acceptance pipeline collapses each canonical pair's evidence into exactly
one row before any feature extraction happens here.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from protein_graph.pipelines.biogrid_pair_acceptance import (
    BioGRIDPairAcceptancePipeline,
)
from protein_graph.pipelines.pair_features import PairFeaturePipeline
from protein_graph.pipelines.protein_features import ProteinFeaturePipeline


class TrainingDatasetPipeline:
    """Build a pair-level dataset from raw global BioGRID evidence."""

    def __init__(
        self,
        protein_feature_pipeline: ProteinFeaturePipeline | None = None,
        pair_feature_pipeline: PairFeaturePipeline | None = None,
        pair_acceptance_pipeline: (
            BioGRIDPairAcceptancePipeline | None
        ) = None,
    ) -> None:
        self.protein_feature_pipeline = (
            protein_feature_pipeline or ProteinFeaturePipeline()
        )
        self.pair_feature_pipeline = (
            pair_feature_pipeline or PairFeaturePipeline()
        )
        self.pair_acceptance_pipeline = (
            pair_acceptance_pipeline or BioGRIDPairAcceptancePipeline()
        )

    def build(
        self,
        evidence: pd.DataFrame,
        force_download: bool = False,
    ) -> tuple[pd.DataFrame, dict[str, Any]]:
        """Build a pair-level feature dataset from raw BioGRID evidence."""
        accepted_pairs, acceptance_report = (
            self.pair_acceptance_pipeline.build(evidence)
        )

        if accepted_pairs.empty:
            raise ValueError(
                "No BioGRID pairs satisfied the strong/direct evidence "
                "policy (see biogrid_evidence_policy). This does not mean "
                "the pairs are biologically negative - only that no row "
                "qualifies for the current strong-positive training "
                "criterion."
            )

        accessions = self._collect_accessions(accepted_pairs)

        protein_features = self.protein_feature_pipeline.build(
            accessions=accessions,
            force_download=force_download,
        )

        resolved_accessions = set(
            protein_features["protein_id"].astype(str)
        )

        unresolved_accessions = sorted(
            set(accessions) - resolved_accessions
        )

        available_pairs = accepted_pairs.loc[
            accepted_pairs["protein_a_id"].isin(
                resolved_accessions
            )
            & accepted_pairs["protein_b_id"].isin(
                resolved_accessions
            )
        ].copy()

        pairs_dropped_missing_features = (
            len(accepted_pairs)
            - len(available_pairs)
        )

        pair_features = self.pair_feature_pipeline.build(
            interactions=available_pairs,
            protein_features=protein_features,
        )

        report = {
            **acceptance_report,
            "requested_uniprot_accessions": len(accessions),
            "resolved_uniprot_accessions": len(
                resolved_accessions
            ),
            "unresolved_uniprot_accessions": (
                unresolved_accessions
            ),
            "protein_feature_rows": len(protein_features),
            "pairs_dropped_missing_protein_features": (
                pairs_dropped_missing_features
            ),
            "pair_feature_rows": len(pair_features),
        }

        return pair_features, report

    def build_model_matrix(
        self,
        pair_features: pd.DataFrame,
    ) -> pd.DataFrame:
        """Build the production model feature matrix.

        This is the explicit fail-closed boundary between REFERENCE/PAIR
        DATA (canonical IDs, taxonomy/context, provenance, extracted
        biological features - everything ``build`` returns) and the MODEL
        FEATURE MATRIX (only reviewed, approved columns). A new or
        unreviewed column present in ``pair_features`` never leaks into the
        matrix just because it exists - that filtering is delegated to
        :meth:`PairFeaturePipeline.select_model_features`, which stays
        permissive for lower-level/non-production callers.

        This method additionally enforces SCHEMA COMPLETENESS: every
        column in ``PairFeaturePipeline.MODEL_FEATURE_COLUMNS`` must be
        present in ``pair_features``, or this raises with the exact
        missing column name(s) rather than silently returning a smaller
        matrix. In the real production path every approved column is
        unconditionally produced (the a_*/b_* sequence-feature columns by
        ProteinFeaturePipeline's fixed extractor set, the pair-level
        columns by ``PairFeaturePipeline._add_pair_features``), so a
        missing column here means the upstream pipeline is broken or was
        called with a hand-built table - never a case to paper over.

        On success, columns are returned in a fixed, deterministic order
        (``sorted(MODEL_FEATURE_COLUMNS)``) rather than whichever order
        they happened to appear in ``pair_features`` - the resulting
        matrix's shape and column order are therefore stable across calls
        and independent of upstream column ordering.
        """
        approved_columns = PairFeaturePipeline.MODEL_FEATURE_COLUMNS

        missing_columns = sorted(
            approved_columns - set(pair_features.columns)
        )

        if missing_columns:
            raise ValueError(
                "Model matrix is missing required approved feature "
                f"column(s): {missing_columns}. The production "
                "model-matrix boundary (TrainingDatasetPipeline."
                "build_model_matrix) requires the complete reviewed "
                "MODEL_FEATURE_COLUMNS schema - fix the upstream feature "
                "pipeline rather than building an incomplete matrix."
            )

        model_features = PairFeaturePipeline.select_model_features(
            pair_features
        )

        return model_features[sorted(approved_columns)]

    def build_and_save(
        self,
        evidence: pd.DataFrame,
        output_path: str | Path = (
            "data/processed/biogrid_pair_features.parquet"
        ),
        model_features_output_path: str | Path = (
            "data/processed/biogrid_pair_model_features.parquet"
        ),
        report_path: str | Path = (
            "data/processed/biogrid_pair_features_report.json"
        ),
        force_download: bool = False,
    ) -> tuple[pd.DataFrame, dict[str, Any]]:
        """Build and save the pair-level dataset and its model matrix.

        Saves two separate artifacts: the full reference/pair table
        (``output_path``, includes IDs and provenance) and the fail-closed
        model feature matrix (``model_features_output_path``, only
        ``PairFeaturePipeline.MODEL_FEATURE_COLUMNS``). Keeping them as
        separate files keeps provenance from ever being mistaken for a
        model input.
        """
        pair_features, report = self.build(
            evidence=evidence,
            force_download=force_download,
        )

        model_features = self.build_model_matrix(pair_features)

        output = Path(output_path)
        output.parent.mkdir(parents=True, exist_ok=True)
        pair_features.to_parquet(output, index=False)

        model_output = Path(model_features_output_path)
        model_output.parent.mkdir(parents=True, exist_ok=True)
        model_features.to_parquet(model_output, index=False)

        report["model_feature_columns"] = len(model_features.columns)
        report["model_feature_rows"] = len(model_features)

        report_destination = Path(report_path)
        report_destination.parent.mkdir(
            parents=True,
            exist_ok=True,
        )
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
