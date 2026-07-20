"""Tests for the motif feature pipeline."""

from typing import Any

import pandas as pd
import pytest

from protein_graph.entities.motif import MotifRecord
from protein_graph.features.base import FeatureMap
from protein_graph.features.elm_motifs import (
    MotifFeatureContext,
)
from protein_graph.features.motif_pipeline import (
    MotifFeaturePipeline,
)


def make_motif(**overrides: Any) -> MotifRecord:
    values: dict[str, Any] = {
        "protein_id": "P12345",
        "region_type": "motif",
        "start": 11,
        "end": 14,
        "sequence": "PQQP",
        "motif_id": "LIG_SH3_GENERIC",
        "motif_class": "SH3-binding motif",
        "pattern": r"P..P",
        "recognized_by": (
            "SH3",
            "WW",
        ),
        "interaction_role": "ligand",
        "evidence_source": "ELM",
        "evidence_level": "curated",
        "confidence": 0.8,
    }

    values.update(overrides)

    return MotifRecord(**values)


def test_pipeline_combines_default_motif_features() -> None:
    motif = make_motif()

    pipeline = MotifFeaturePipeline()

    table = pipeline.build_from_motifs(
        motifs=[motif],
        protein_length=100,
    )

    assert isinstance(table, pd.DataFrame)
    assert len(table) == 1

    row = table.iloc[0]

    # Raw metadata
    assert row["protein_id"] == "P12345"
    assert row["motif_id"] == "LIG_SH3_GENERIC"
    assert row["motif_class"] == (
        "SH3-binding motif"
    )
    assert set(row["recognized_by"]) == {
        "SH3",
        "WW",
    }

    # Geometry
    assert row["motif_length"] == 4
    assert row[
        "motif_relative_position"
    ] == pytest.approx(0.125)
    assert row[
        "motif_fraction_of_protein"
    ] == pytest.approx(0.04)
    assert row["motif_distance_to_n"] == 10
    assert row["motif_distance_to_c"] == 86

    # Position category
    assert row[
        "motif_position_n_terminal"
    ] == 1
    assert row["motif_position_middle"] == 0
    assert row[
        "motif_position_c_terminal"
    ] == 0

    # Annotations
    assert row["has_motif_id"] == 1
    assert row["has_motif_class"] == 1
    assert row["has_motif_pattern"] == 1
    assert row["has_recognized_by"] == 1
    assert row["has_interaction_role"] == 1
    assert row["has_evidence_source"] == 1
    assert row["has_evidence_level"] == 1
    assert row["has_motif_confidence"] == 1

    assert row["num_recognized_domains"] == 2
    assert row["num_motif_annotations"] == 8

    assert row[
        "motif_confidence"
    ] == pytest.approx(0.8)


def test_pipeline_builds_one_row_per_motif() -> None:
    motifs = [
        make_motif(
            motif_id="MOTIF_A",
            start=1,
            end=4,
        ),
        make_motif(
            motif_id="MOTIF_B",
            start=51,
            end=54,
        ),
    ]

    table = MotifFeaturePipeline().build_from_motifs(
        motifs=motifs,
        protein_length=100,
    )

    assert len(table) == 2

    assert list(table["motif_id"]) == [
        "MOTIF_A",
        "MOTIF_B",
    ]


def test_pipeline_returns_empty_table_for_no_motifs() -> None:
    table = MotifFeaturePipeline().build_from_motifs(
        motifs=[],
        protein_length=100,
    )

    assert isinstance(table, pd.DataFrame)
    assert table.empty


def test_pipeline_accepts_custom_extractors() -> None:
    class ConstantExtractor:
        def extract(
            self,
            motif: MotifRecord,
            context: MotifFeatureContext,
        ) -> FeatureMap:
            del motif
            del context

            return {
                "custom_feature": 42,
            }

    pipeline = MotifFeaturePipeline(
        extractors=[
            ConstantExtractor(),
        ]
    )

    table = pipeline.build_from_motifs(
        motifs=[make_motif()],
        protein_length=100,
    )

    row = table.iloc[0]

    assert row["custom_feature"] == 42

    assert "motif_length" not in row
    assert "has_motif_id" not in row


def test_pipeline_rejects_duplicate_feature_keys() -> None:
    class FirstExtractor:
        def extract(
            self,
            motif: MotifRecord,
            context: MotifFeatureContext,
        ) -> FeatureMap:
            del motif
            del context

            return {
                "shared_feature": 1,
            }

    class SecondExtractor:
        def extract(
            self,
            motif: MotifRecord,
            context: MotifFeatureContext,
        ) -> FeatureMap:
            del motif
            del context

            return {
                "shared_feature": 2,
            }

    pipeline = MotifFeaturePipeline(
        extractors=[
            FirstExtractor(),
            SecondExtractor(),
        ]
    )

    with pytest.raises(
        ValueError,
        match="Duplicate motif feature keys",
    ):
        pipeline.build_from_motifs(
            motifs=[make_motif()],
            protein_length=100,
        )


def test_pipeline_rejects_feature_overwriting_metadata() -> None:
    class InvalidExtractor:
        def extract(
            self,
            motif: MotifRecord,
            context: MotifFeatureContext,
        ) -> FeatureMap:
            del motif
            del context

            return {
                "motif_id": "OVERWRITTEN",
            }

    pipeline = MotifFeaturePipeline(
        extractors=[
            InvalidExtractor(),
        ]
    )

    with pytest.raises(
        ValueError,
        match="Duplicate motif feature keys",
    ):
        pipeline.build_from_motifs(
            motifs=[make_motif()],
            protein_length=100,
        )
