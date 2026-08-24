"""Tests for pair feature pipeline."""

import pandas as pd
import pytest

from protein_graph.features.sequence.constants import AMINO_ACIDS
from protein_graph.pipelines import PairFeaturePipeline, ProteinFeaturePipeline


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


def test_same_taxonomy_is_unknown_when_taxon_id_missing() -> None:
    features = make_protein_features()
    features.loc[
        features["protein_id"] == "P1", "taxon_id"
    ] = None

    pairs = PairFeaturePipeline().build(
        interactions=make_interactions(),
        protein_features=features,
    )

    assert pd.isna(pairs.iloc[0]["same_taxonomy"])


def test_same_taxonomy_is_zero_for_different_known_taxa() -> None:
    features = make_protein_features()
    features.loc[
        features["protein_id"] == "P2", "taxon_id"
    ] = 9606

    pairs = PairFeaturePipeline().build(
        interactions=make_interactions(),
        protein_features=features,
    )

    assert pairs.iloc[0]["same_taxonomy"] == 0


def test_same_taxonomy_is_unknown_when_b_taxon_id_missing() -> None:
    features = make_protein_features()
    features.loc[
        features["protein_id"] == "P2", "taxon_id"
    ] = None

    pairs = PairFeaturePipeline().build(
        interactions=make_interactions(),
        protein_features=features,
    )

    assert pd.isna(pairs.iloc[0]["same_taxonomy"])


def test_same_taxonomy_is_unknown_when_both_taxon_ids_missing() -> None:
    """Unknown must never collapse to 0 ("definitely different taxa")."""
    features = make_protein_features()
    features["taxon_id"] = None

    pairs = PairFeaturePipeline().build(
        interactions=make_interactions(),
        protein_features=features,
    )

    assert pd.isna(pairs.iloc[0]["same_taxonomy"])


def test_cross_species_host_pathogen_pair_is_retained_not_filtered() -> None:
    """Taxonomy must not become a hidden eligibility filter: a host-pathogen
    pair (human vs. SARS-CoV-2) gets same_taxonomy=0 as a feature value, but
    the pair itself is never dropped because of it."""
    features = make_protein_features()
    features.loc[features["protein_id"] == "P1", "taxon_id"] = 9606
    features.loc[features["protein_id"] == "P2", "taxon_id"] = 2697049

    pairs = PairFeaturePipeline().build(
        interactions=make_interactions(),
        protein_features=features,
    )

    assert len(pairs) == 1
    assert pairs.iloc[0]["same_taxonomy"] == 0


def test_model_features_retains_approved_columns() -> None:
    pairs = PairFeaturePipeline().build(
        interactions=make_interactions(),
        protein_features=make_protein_features(),
    )

    model_features = PairFeaturePipeline.select_model_features(pairs)

    assert "same_taxonomy" in model_features.columns
    assert "protein_length_difference" in model_features.columns
    assert "hydrophobic_fraction_difference" in model_features.columns
    assert "charged_fraction_difference" in model_features.columns
    assert "positive_negative_balance" in model_features.columns
    assert "a_hydrophobic_fraction" in model_features.columns
    assert "b_hydrophobic_fraction" in model_features.columns
    assert "a_protein_length" in model_features.columns
    assert "b_protein_length" in model_features.columns


def test_model_features_excludes_unknown_new_columns() -> None:
    pairs = PairFeaturePipeline().build(
        interactions=make_interactions(),
        protein_features=make_protein_features(),
    )
    pairs = pairs.assign(
        totally_new_column_nobody_reviewed="unreviewed"
    )

    model_features = PairFeaturePipeline.select_model_features(pairs)

    assert (
        "totally_new_column_nobody_reviewed"
        not in model_features.columns
    )


def test_missing_approved_columns_are_absent_not_fabricated() -> None:
    """If an approved MODEL_FEATURE_COLUMNS entry is simply absent from
    the input (e.g. an extractor was removed, or a caller passes a
    hand-built table), select_model_features must not silently invent a
    placeholder value (0/NaN) for it - the column must just be missing
    from the output, not present-but-wrong."""
    pairs = pd.DataFrame(
        [{"protein_a_id": "P1", "protein_b_id": "P2", "same_taxonomy": 1}]
    )

    model_features = PairFeaturePipeline.select_model_features(pairs)

    assert list(model_features.columns) == ["same_taxonomy"]
    assert "protein_length_difference" not in model_features.columns
    assert "a_protein_length" not in model_features.columns


def test_model_features_excludes_provenance_columns() -> None:
    pairs = PairFeaturePipeline().build(
        interactions=make_interactions(),
        protein_features=make_protein_features(),
    )
    pairs = pairs.assign(
        score=9.9,
        publication_id="PMID:12345",
        biogrid_interaction_id="999",
        experimental_system="Two-hybrid",
        throughput="High Throughput",
    )

    model_features = PairFeaturePipeline.select_model_features(pairs)

    assert "score" not in model_features.columns
    assert "publication_id" not in model_features.columns
    assert "biogrid_interaction_id" not in model_features.columns
    assert "experimental_system" not in model_features.columns
    assert "throughput" not in model_features.columns


def test_model_features_excludes_biogrid_pair_acceptance_provenance() -> None:
    """Minimal provenance emitted by BioGRIDPairAcceptancePipeline
    (representative_* / biogrid_release_version) answers "why was this
    pair admitted" - it must never leak into the model feature matrix."""
    pairs = PairFeaturePipeline().build(
        interactions=make_interactions(),
        protein_features=make_protein_features(),
    )
    pairs = pairs.assign(
        representative_biogrid_interaction_id="999",
        representative_experimental_system="Two-hybrid",
        representative_publication_id="PMID:12345",
        biogrid_release_version="5.0.260",
    )

    model_features = PairFeaturePipeline.select_model_features(pairs)

    assert "representative_biogrid_interaction_id" not in model_features.columns
    assert "representative_experimental_system" not in model_features.columns
    assert "representative_publication_id" not in model_features.columns
    assert "biogrid_release_version" not in model_features.columns


def test_model_features_excludes_identity_columns() -> None:
    pairs = PairFeaturePipeline().build(
        interactions=make_interactions(),
        protein_features=make_protein_features(),
    )

    model_features = PairFeaturePipeline.select_model_features(pairs)

    assert "a_protein_id" not in model_features.columns
    assert "b_protein_id" not in model_features.columns
    assert "protein_a_id" not in model_features.columns
    assert "protein_b_id" not in model_features.columns
    assert "uniprot_pair_id" not in model_features.columns
    assert "a_taxon_id" not in model_features.columns
    assert "b_taxon_id" not in model_features.columns


def test_model_feature_columns_matches_live_extractor_output() -> None:
    """MODEL_FEATURE_COLUMNS' a_/b_ sequence-feature portion must exactly
    match what ProteinFeaturePipeline's currently configured extractors
    actually produce.

    This guards against drift in both directions: a new extractor
    feature that was never added to the allowlist (would silently be
    dropped from any model input, not a leakage risk but a coverage
    regression), and a stale allowlist entry no longer produced by any
    extractor (dead entry, or worse - a sign the allowlist was hand-
    edited out of sync with the real pipeline).

    Reads ProteinFeaturePipeline().extractors directly, so it tracks
    whatever extractors are actually configured today rather than a
    second hardcoded extractor list that could itself drift.
    """
    live_columns: set[str] = set()

    for extractor in ProteinFeaturePipeline().extractors:
        live_columns.update(
            extractor.extract({"sequence": AMINO_ACIDS}).keys()
        )

    allowlist = PairFeaturePipeline.MODEL_FEATURE_COLUMNS

    a_prefixed = {
        name[len("a_"):]
        for name in allowlist
        if name.startswith("a_")
    }
    b_prefixed = {
        name[len("b_"):]
        for name in allowlist
        if name.startswith("b_")
    }

    assert a_prefixed == live_columns, (
        "a_-prefixed portion of MODEL_FEATURE_COLUMNS diverges from "
        "ProteinFeaturePipeline's live extractor output.\n"
        f"Missing from allowlist: {sorted(live_columns - a_prefixed)}\n"
        f"Stale in allowlist: {sorted(a_prefixed - live_columns)}"
    )
    assert b_prefixed == live_columns, (
        "b_-prefixed portion of MODEL_FEATURE_COLUMNS diverges from "
        "ProteinFeaturePipeline's live extractor output.\n"
        f"Missing from allowlist: {sorted(live_columns - b_prefixed)}\n"
        f"Stale in allowlist: {sorted(b_prefixed - live_columns)}"
    )


def test_repr() -> None:
    assert repr(PairFeaturePipeline()) == "PairFeaturePipeline()"
