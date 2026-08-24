"""Tests for the complete training dataset pipeline."""

import pandas as pd
import pytest

from protein_graph.pipelines.pair_features import PairFeaturePipeline
from protein_graph.pipelines.training_dataset import TrainingDatasetPipeline


class FakePairAcceptancePipeline:
    """Return already-accepted pairs, bypassing real BioGRID evidence."""

    def __init__(self, accepted: pd.DataFrame) -> None:
        self._accepted = accepted

    def build(self, evidence: pd.DataFrame):
        return self._accepted, {
            "total_evidence_rows": len(evidence),
            "accepted_pairs": len(self._accepted),
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


class CompletePairFeaturePipeline:
    """Like FakePairFeaturePipeline, but fills in every real
    PairFeaturePipeline.MODEL_FEATURE_COLUMNS entry with a dummy value.

    build_model_matrix() enforces that the complete reviewed schema is
    present; tests exercising that success path need a fixture producing
    the full (~300-column) allowlist without hand-enumerating every
    feature name, and without going stale if the real allowlist changes.
    """

    def build(
        self,
        interactions: pd.DataFrame,
        protein_features: pd.DataFrame,
    ) -> pd.DataFrame:
        pairs = FakePairFeaturePipeline().build(
            interactions,
            protein_features,
        )

        missing_columns = [
            column
            for column in PairFeaturePipeline.MODEL_FEATURE_COLUMNS
            if column not in pairs.columns
        ]

        filler = pd.DataFrame(
            0.0,
            index=pairs.index,
            columns=missing_columns,
        )

        completed = pd.concat([pairs, filler], axis=1)

        return completed


def _accepted_pairs(rows: list[dict]) -> pd.DataFrame:
    return pd.DataFrame(rows)


def test_build_creates_complete_dataset() -> None:
    accepted = _accepted_pairs(
        [
            {"protein_a_id": "P1", "protein_b_id": "P2"},
            {"protein_a_id": "P2", "protein_b_id": "P3"},
        ]
    )

    pipeline = TrainingDatasetPipeline(
        pair_acceptance_pipeline=FakePairAcceptancePipeline(accepted),
        protein_feature_pipeline=FakeProteinFeaturePipeline(),
        pair_feature_pipeline=FakePairFeaturePipeline(),
    )

    dataset, report = pipeline.build(pd.DataFrame({"experimental_system": []}))

    assert len(dataset) == 2
    assert report["requested_uniprot_accessions"] == 3
    assert report["resolved_uniprot_accessions"] == 3
    assert report["unresolved_uniprot_accessions"] == []
    assert report["pairs_dropped_missing_protein_features"] == 0
    assert report["pair_feature_rows"] == 2
    assert report["accepted_pairs"] == 2


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


def test_build_rejects_when_no_pair_qualifies() -> None:
    pipeline = TrainingDatasetPipeline(
        pair_acceptance_pipeline=FakePairAcceptancePipeline(
            _accepted_pairs([])
        ),
        protein_feature_pipeline=FakeProteinFeaturePipeline(),
        pair_feature_pipeline=FakePairFeaturePipeline(),
    )

    with pytest.raises(
        ValueError,
        match="No BioGRID pairs satisfied",
    ):
        pipeline.build(pd.DataFrame({"experimental_system": []}))


def test_build_drops_pairs_with_missing_protein_features() -> None:
    accepted = _accepted_pairs(
        [
            {"protein_a_id": "P1", "protein_b_id": "P2"},
            {"protein_a_id": "P2", "protein_b_id": "P3"},
        ]
    )

    pipeline = TrainingDatasetPipeline(
        pair_acceptance_pipeline=FakePairAcceptancePipeline(accepted),
        protein_feature_pipeline=IncompleteProteinFeaturePipeline(),
        pair_feature_pipeline=FakePairFeaturePipeline(),
    )

    dataset, report = pipeline.build(pd.DataFrame({"experimental_system": []}))

    assert len(dataset) == 1
    assert report["requested_uniprot_accessions"] == 3
    assert report["resolved_uniprot_accessions"] == 2
    assert report["unresolved_uniprot_accessions"] == ["P3"]
    assert report["pairs_dropped_missing_protein_features"] == 1
    assert report["pair_feature_rows"] == 1


def test_build_model_matrix_excludes_identity_and_provenance() -> None:
    accepted = _accepted_pairs(
        [{"protein_a_id": "P1", "protein_b_id": "P2"}]
    )

    pipeline = TrainingDatasetPipeline(
        pair_acceptance_pipeline=FakePairAcceptancePipeline(accepted),
        protein_feature_pipeline=FakeProteinFeaturePipeline(),
        pair_feature_pipeline=CompletePairFeaturePipeline(),
    )

    pair_features, _ = pipeline.build(pd.DataFrame({"experimental_system": []}))
    model_features = pipeline.build_model_matrix(pair_features)

    assert "protein_a_id" not in model_features.columns
    assert "protein_b_id" not in model_features.columns
    assert "a_protein_id" not in model_features.columns


def test_build_model_matrix_succeeds_with_the_full_approved_schema() -> None:
    """When every PairFeaturePipeline.MODEL_FEATURE_COLUMNS entry is
    present, build_model_matrix returns exactly that column set - the
    complete reviewed schema, no more, no less."""
    accepted = _accepted_pairs(
        [{"protein_a_id": "P1", "protein_b_id": "P2"}]
    )

    pipeline = TrainingDatasetPipeline(
        pair_acceptance_pipeline=FakePairAcceptancePipeline(accepted),
        protein_feature_pipeline=FakeProteinFeaturePipeline(),
        pair_feature_pipeline=CompletePairFeaturePipeline(),
    )

    pair_features, _ = pipeline.build(pd.DataFrame({"experimental_system": []}))
    model_features = pipeline.build_model_matrix(pair_features)

    assert set(model_features.columns) == PairFeaturePipeline.MODEL_FEATURE_COLUMNS
    assert len(model_features) == 1


def test_build_model_matrix_column_order_is_deterministic() -> None:
    """The output column order must be a fixed function of the approved
    allowlist (sorted(MODEL_FEATURE_COLUMNS)), not incidentally whatever
    order the upstream pair-feature table happened to produce."""
    accepted = _accepted_pairs(
        [{"protein_a_id": "P1", "protein_b_id": "P2"}]
    )

    pipeline = TrainingDatasetPipeline(
        pair_acceptance_pipeline=FakePairAcceptancePipeline(accepted),
        protein_feature_pipeline=FakeProteinFeaturePipeline(),
        pair_feature_pipeline=CompletePairFeaturePipeline(),
    )

    pair_features, _ = pipeline.build(pd.DataFrame({"experimental_system": []}))
    model_features = pipeline.build_model_matrix(pair_features)

    assert list(model_features.columns) == sorted(
        PairFeaturePipeline.MODEL_FEATURE_COLUMNS
    )

    # Column order must be stable across repeated calls too, not just
    # incidentally correct once.
    model_features_again = pipeline.build_model_matrix(pair_features)
    assert list(model_features.columns) == list(model_features_again.columns)


def test_build_model_matrix_fails_on_one_missing_approved_column() -> None:
    """A single missing approved column must fail clearly and name it -
    never silently return a smaller matrix or fabricate the value."""
    accepted = _accepted_pairs(
        [{"protein_a_id": "P1", "protein_b_id": "P2"}]
    )

    pipeline = TrainingDatasetPipeline(
        pair_acceptance_pipeline=FakePairAcceptancePipeline(accepted),
        protein_feature_pipeline=FakeProteinFeaturePipeline(),
        pair_feature_pipeline=CompletePairFeaturePipeline(),
    )

    pair_features, _ = pipeline.build(pd.DataFrame({"experimental_system": []}))
    incomplete = pair_features.drop(columns=["same_taxonomy"])

    with pytest.raises(ValueError, match="same_taxonomy"):
        pipeline.build_model_matrix(incomplete)


def test_build_model_matrix_reports_all_missing_approved_columns() -> None:
    """Several missing approved columns must all be named in the error,
    not just the first one encountered."""
    accepted = _accepted_pairs(
        [{"protein_a_id": "P1", "protein_b_id": "P2"}]
    )

    pipeline = TrainingDatasetPipeline(
        pair_acceptance_pipeline=FakePairAcceptancePipeline(accepted),
        protein_feature_pipeline=FakeProteinFeaturePipeline(),
        pair_feature_pipeline=CompletePairFeaturePipeline(),
    )

    pair_features, _ = pipeline.build(pd.DataFrame({"experimental_system": []}))
    incomplete = pair_features.drop(
        columns=[
            "same_taxonomy",
            "protein_length_difference",
            "a_protein_length",
        ]
    )

    with pytest.raises(ValueError) as excinfo:
        pipeline.build_model_matrix(incomplete)

    message = str(excinfo.value)
    assert "same_taxonomy" in message
    assert "protein_length_difference" in message
    assert "a_protein_length" in message


def test_build_model_matrix_never_synthesizes_a_missing_feature() -> None:
    """A missing approved column must never be silently backfilled with a
    fabricated value (0/NaN/otherwise) - the only acceptable outcomes are
    "complete matrix returned" or "clear failure", never a matrix that
    quietly contains an invented column."""
    accepted = _accepted_pairs(
        [{"protein_a_id": "P1", "protein_b_id": "P2"}]
    )

    pipeline = TrainingDatasetPipeline(
        pair_acceptance_pipeline=FakePairAcceptancePipeline(accepted),
        protein_feature_pipeline=FakeProteinFeaturePipeline(),
        pair_feature_pipeline=CompletePairFeaturePipeline(),
    )

    pair_features, _ = pipeline.build(pd.DataFrame({"experimental_system": []}))
    incomplete = pair_features.drop(columns=["same_taxonomy"])

    try:
        pipeline.build_model_matrix(incomplete)
        raised = False
    except ValueError:
        raised = True

    assert raised, (
        "build_model_matrix must fail rather than silently synthesize "
        "the missing 'same_taxonomy' column"
    )


def test_build_model_matrix_excludes_extra_provenance_and_numeric_columns() -> (
    None
):
    """Extra unreviewed columns (provenance or otherwise) must stay
    excluded even when the schema is otherwise complete."""
    accepted = _accepted_pairs(
        [{"protein_a_id": "P1", "protein_b_id": "P2"}]
    )

    pipeline = TrainingDatasetPipeline(
        pair_acceptance_pipeline=FakePairAcceptancePipeline(accepted),
        protein_feature_pipeline=FakeProteinFeaturePipeline(),
        pair_feature_pipeline=CompletePairFeaturePipeline(),
    )

    pair_features, _ = pipeline.build(pd.DataFrame({"experimental_system": []}))
    pair_features = pair_features.assign(
        representative_experimental_system="Two-hybrid",
        biogrid_release_version="5.0.260",
        some_future_numeric_provenance_score=0.99,
    )

    model_features = pipeline.build_model_matrix(pair_features)

    assert "representative_experimental_system" not in model_features.columns
    assert "biogrid_release_version" not in model_features.columns
    assert "some_future_numeric_provenance_score" not in model_features.columns
    assert set(model_features.columns) == PairFeaturePipeline.MODEL_FEATURE_COLUMNS


def test_build_and_save_writes_two_separate_row_aligned_files(
    tmp_path,
) -> None:
    """The reference table and the model matrix must be genuinely separate
    artifacts that cannot be accidentally swapped: distinct output paths,
    but the same number of rows and the same row order (the model matrix
    is a strict column-subset of the reference table, never re-sorted or
    re-joined), so row N in one file always corresponds to row N in the
    other."""
    accepted = _accepted_pairs(
        [
            {"protein_a_id": "P1", "protein_b_id": "P2"},
            {"protein_a_id": "P2", "protein_b_id": "P3"},
            {"protein_a_id": "P1", "protein_b_id": "P3"},
        ]
    )

    pipeline = TrainingDatasetPipeline(
        pair_acceptance_pipeline=FakePairAcceptancePipeline(accepted),
        protein_feature_pipeline=FakeProteinFeaturePipeline(),
        pair_feature_pipeline=CompletePairFeaturePipeline(),
    )

    output_path = tmp_path / "reference.parquet"
    model_path = tmp_path / "model.parquet"
    report_path = tmp_path / "report.json"

    pair_features, report = pipeline.build_and_save(
        evidence=pd.DataFrame({"experimental_system": []}),
        output_path=output_path,
        model_features_output_path=model_path,
        report_path=report_path,
    )

    assert output_path.exists()
    assert model_path.exists()
    assert output_path != model_path

    reloaded_reference = pd.read_parquet(output_path)
    reloaded_model = pd.read_parquet(model_path)

    assert len(reloaded_reference) == len(pair_features)
    assert len(reloaded_model) == len(pair_features)
    assert len(reloaded_reference) == len(reloaded_model)

    # The model matrix is a column-subset of the same rows, in the same
    # order - every model-matrix column's values must match the
    # corresponding reference-table column position-by-position.
    for column in reloaded_model.columns:
        assert (
            reloaded_model[column].tolist()
            == reloaded_reference[column].tolist()
        )

    # Identity/provenance columns present in the reference table must be
    # absent from the model matrix file on disk (not just in memory).
    assert "protein_a_id" not in reloaded_model.columns
    assert "protein_a_id" in reloaded_reference.columns

    assert report["model_feature_rows"] == len(reloaded_model)
    assert report["model_feature_columns"] == len(reloaded_model.columns)


def test_repr() -> None:
    pipeline = TrainingDatasetPipeline(
        pair_acceptance_pipeline=FakePairAcceptancePipeline(
            _accepted_pairs([])
        ),
        protein_feature_pipeline=FakeProteinFeaturePipeline(),
        pair_feature_pipeline=FakePairFeaturePipeline(),
    )

    assert repr(pipeline) == "TrainingDatasetPipeline()"
