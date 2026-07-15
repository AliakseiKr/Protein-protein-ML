"""Tests for the protein feature pipeline."""

import pandas as pd
import pytest

from protein_graph.pipelines import ProteinFeaturePipeline


class FakeProteinBuilder:
    """Return protein records without network access."""

    def build_many(
        self,
        accessions,
        force_download=False,
    ) -> pd.DataFrame:
        return pd.DataFrame(
            [
                {
                    "protein_id": "P12345",
                    "gene_name": "Test1",
                    "organism": "Rattus norvegicus",
                    "taxon_id": 10116,
                    "sequence": "MKKKKKAVILMFWYAVILMFWYDEDEDE",
                },
                {
                    "protein_id": "Q99999",
                    "gene_name": "Test2",
                    "organism": "Rattus norvegicus",
                    "taxon_id": 10116,
                    "sequence": "MSTNQCYLQQLQQLQQLAVILMFWY",
                },
            ]
        )


def test_build_creates_feature_table() -> None:
    pipeline = ProteinFeaturePipeline(
        protein_builder=FakeProteinBuilder(),
    )

    features = pipeline.build(["P12345", "Q99999"])

    assert len(features) == 2
    assert set(features["protein_id"]) == {"P12345", "Q99999"}
    assert "hydrophobic_fraction" in features.columns
    assert "max_hydrophobic_run" in features.columns
    assert (
        "high_hydrophobic_density_window_count_5"
        in features.columns
    )


def test_build_rejects_duplicate_proteins() -> None:
    features = pd.DataFrame(
        {
            "protein_id": ["P12345", "P12345"],
            "protein_length": [10, 10],
        }
    )

    with pytest.raises(
        ValueError,
        match="duplicate protein IDs",
    ):
        ProteinFeaturePipeline._validate_table(features)


def test_save_creates_parquet(tmp_path) -> None:
    pipeline = ProteinFeaturePipeline(
        protein_builder=FakeProteinBuilder(),
    )

    features = pipeline.build(["P12345", "Q99999"])
    destination = tmp_path / "protein_features.parquet"

    result = pipeline.save(features, destination)

    loaded = pd.read_parquet(result)

    assert result == destination
    assert len(loaded) == 2
    assert set(loaded["protein_id"]) == {"P12345", "Q99999"}


def test_repr() -> None:
    pipeline = ProteinFeaturePipeline(
        protein_builder=FakeProteinBuilder(),
    )

    assert "extractors=5" in repr(pipeline)
