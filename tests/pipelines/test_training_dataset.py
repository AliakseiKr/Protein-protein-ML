"""Tests for the complete training dataset pipeline."""

import pandas as pd
import pytest

from protein_graph.pipelines.training_dataset import TrainingDatasetPipeline


class FakeMapper:
    """Return already mapped interaction records."""

    def map(self, interactions: pd.DataFrame) -> pd.DataFrame:
        mapped = interactions.copy()

        mapped["mapping_complete"] = (
            mapped["protein_a_id"].notna()
            & mapped["protein_b_id"].notna()
        )

        return mapped

    def coverage(self, mapped: pd.DataFrame) -> dict:
        fully_mapped = int(mapped["mapping_complete"].sum())

        partially_mapped = int(
            (
                mapped["protein_a_id"].notna()
                ^ mapped["protein_b_id"].notna()
            ).sum()
        )

        unmapped = len(mapped) - fully_mapped - partially_mapped

        return {
            "total_rows": len(mapped),
            "fully_mapped_rows": fully_mapped,
            "partially_mapped_rows": partially_mapped,
            "unmapped_rows": unmapped,
            "fully_mapped_percentage": (
                fully_mapped / len(mapped) * 100
                if len(mapped)
                else 0.0
            ),
        }


class FakeProteinFeaturePipeline:
    """Return features for every requested protein."""

    def build(
        self,
        accessions,
        force_download=False,
    ) -> pd.DataFrame:
        rows = []

        for accession in accessions:
            rows.append(
                {
                    "protein_id": accession,
                    "protein_length": 100,
                    "hydrophobic_fraction": 0.4,
                    "charged_fraction": 0.2,
                    "positive_fraction": 0.1,
                    "negative_fraction": 0.1,
                    "taxon_id": 10116,
                }
            )

        return pd.DataFrame(rows)


class IncompleteProteinFeaturePipeline:
    """Return features for only part of the requested proteins."""

    def build(
        self,
        accessions,
        force_download=False,
    ) -> pd.DataFrame:
        return pd.DataFrame(
            [
                {
                    "protein_id": "P1",
                    "protein_length": 100,
                    "hydrophobic_fraction": 0.4,
                    "charged_fraction": 0.2,
                    "positive_fraction": 0.1,
                    "negative_fraction": 0.1,
                    "taxon_id": 10116,
                },
                {
                    "protein_id": "P2",
                    "protein_length": 120,
                    "hydrophobic_fraction": 0.5,
                    "charged_fraction": 0.3,
                    "positive_fraction": 0.2,
                    "negative_fraction": 0.1,
                    "taxon_id": 10116,
                },
            ]
        )


class FakePairFeaturePipeline:
    """Join feature rows for both proteins in each pair."""

    def build(
        self,
        interactions,
        protein_features,
    ) -> pd.DataFrame:
        pairs = interactions.merge(
            protein_features.add_prefix("a_"),
            left_on="protein_a_id",
            right_on="a_protein_id",
            how="inner",
        )

        return pairs.merge(
            protein_features.add_prefix("b_"),
            left_on="protein_b_id",
            right_on="b_protein_id",
            how="inner",
        )


def test_build_creates_complete_dataset() -> None:
    interactions = pd.DataFrame(
        [
            {
                "protein_a_id": "P1",
                "protein_b_id": "P2",
            },
            {
                "protein_a_id": "P2",
                "protein_b_id": "P3",
            },
        ]
    )

    pipeline = TrainingDatasetPipeline(
        mapper=FakeMapper(),
        protein_feature_pipeline=FakeProteinFeaturePipeline(),
        pair_feature_pipeline=FakePairFeaturePipeline(),
    )

    dataset, report = pipeline.build(interactions)

    assert len(dataset) == 2
    assert report["requested_uniprot_accessions"] == 3
    assert report["resolved_uniprot_accessions"] == 3
    assert report["unresolved_uniprot_accessions"] == []
    assert report["pairs_dropped_missing_protein_features"] == 0
    assert report["pair_feature_rows"] == 2


def test_collect_accessions_deduplicates_ids() -> None:
    interactions = pd.DataFrame(
        {
            "protein_a_id": ["P1", "P2"],
            "protein_b_id": ["P2", "P1"],
        }
    )

    accessions = TrainingDatasetPipeline._collect_accessions(
        interactions
    )

    assert accessions == ["P1", "P2"]


def test_build_rejects_no_complete_pairs() -> None:
    interactions = pd.DataFrame(
        [
            {
                "protein_a_id": None,
                "protein_b_id": None,
            }
        ]
    )

    pipeline = TrainingDatasetPipeline(
        mapper=FakeMapper(),
        protein_feature_pipeline=FakeProteinFeaturePipeline(),
        pair_feature_pipeline=FakePairFeaturePipeline(),
    )

    with pytest.raises(
        ValueError,
        match="No completely mapped",
    ):
        pipeline.build(interactions)


def test_build_drops_pairs_with_missing_protein_features() -> None:
    interactions = pd.DataFrame(
        [
            {
                "protein_a_id": "P1",
                "protein_b_id": "P2",
            },
            {
                "protein_a_id": "P2",
                "protein_b_id": "P3",
            },
        ]
    )

    pipeline = TrainingDatasetPipeline(
        mapper=FakeMapper(),
        protein_feature_pipeline=IncompleteProteinFeaturePipeline(),
        pair_feature_pipeline=FakePairFeaturePipeline(),
    )

    dataset, report = pipeline.build(interactions)

    assert len(dataset) == 1
    assert report["requested_uniprot_accessions"] == 3
    assert report["resolved_uniprot_accessions"] == 2
    assert report["unresolved_uniprot_accessions"] == ["P3"]
    assert report["pairs_dropped_missing_protein_features"] == 1
    assert report["pair_feature_rows"] == 1


def test_repr() -> None:
    pipeline = TrainingDatasetPipeline(
        mapper=FakeMapper(),
        protein_feature_pipeline=FakeProteinFeaturePipeline(),
        pair_feature_pipeline=FakePairFeaturePipeline(),
    )

    assert repr(pipeline) == "TrainingDatasetPipeline()"
