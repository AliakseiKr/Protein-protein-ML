"""Tests for pair feature pipeline."""

import pandas as pd
import pytest

from protein_graph.pipelines import PairFeaturePipeline


def make_interactions() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "protein_a_id": "P1",
                "protein_b_id": "P2",
                "uniprot_pair_id": "P1--P2",
            }
        ]
    )


def make_protein_features() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "protein_id": "P1",
                "protein_length": 100,
                "hydrophobic_fraction": 0.40,
                "charged_fraction": 0.20,
                "positive_fraction": 0.12,
                "negative_fraction": 0.08,
                "taxon_id": 10116,
            },
            {
                "protein_id": "P2",
                "protein_length": 140,
                "hydrophobic_fraction": 0.55,
                "charged_fraction": 0.25,
                "positive_fraction": 0.10,
                "negative_fraction": 0.15,
                "taxon_id": 10116,
            },
        ]
    )


def test_build_creates_pair_features() -> None:
    pipeline = PairFeaturePipeline()

    pairs = pipeline.build(
        interactions=make_interactions(),
        protein_features=make_protein_features(),
    )

    assert len(pairs) == 1
    assert pairs.iloc[0]["protein_length_difference"] == 40
    assert pairs.iloc[0]["same_taxonomy"] == 1
    assert "a_hydrophobic_fraction" in pairs.columns
    assert "b_hydrophobic_fraction" in pairs.columns


def test_missing_interaction_columns_are_rejected() -> None:
    with pytest.raises(
        ValueError,
        match="Interaction table is missing columns",
    ):
        PairFeaturePipeline().build(
            interactions=pd.DataFrame(),
            protein_features=make_protein_features(),
        )


def test_duplicate_protein_ids_are_rejected() -> None:
    features = pd.concat(
        [
            make_protein_features(),
            make_protein_features().iloc[[0]],
        ],
        ignore_index=True,
    )

    with pytest.raises(
        ValueError,
        match="duplicate protein IDs",
    ):
        PairFeaturePipeline().build(
            interactions=make_interactions(),
            protein_features=features,
        )


def test_unmapped_proteins_are_removed() -> None:
    interactions = pd.DataFrame(
        [
            {
                "protein_a_id": "P1",
                "protein_b_id": "UNKNOWN",
            }
        ]
    )

    pairs = PairFeaturePipeline().build(
        interactions=interactions,
        protein_features=make_protein_features(),
    )

    assert pairs.empty


def test_repr() -> None:
    assert repr(PairFeaturePipeline()) == "PairFeaturePipeline()"
