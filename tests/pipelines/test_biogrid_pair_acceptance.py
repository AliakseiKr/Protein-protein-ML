"""Tests for BioGRIDPairAcceptancePipeline."""

import pandas as pd
import pytest

from protein_graph.pipelines.biogrid_pair_acceptance import (
    BioGRIDPairAcceptancePipeline,
)


def _evidence_row(
    *,
    biogrid_interaction_id: str,
    experimental_system: str,
    swissprot_a: str,
    swissprot_b: str,
    publication_id: str = "PUBMED:1",
    biogrid_release_version: str = "5.0.260",
) -> dict:
    return {
        "biogrid_interaction_id": biogrid_interaction_id,
        "experimental_system": experimental_system,
        "publication_id": publication_id,
        "biogrid_release_version": biogrid_release_version,
        "swissprot_accessions_a": swissprot_a,
        "swissprot_accessions_b": swissprot_b,
        "trembl_accessions_a": None,
        "trembl_accessions_b": None,
    }


def test_pair_with_single_strong_direct_experiment_is_accepted() -> None:
    """A single reviewed strong/direct experiment is sufficient - no
    corroboration is required."""
    evidence = pd.DataFrame(
        [
            _evidence_row(
                biogrid_interaction_id="1",
                experimental_system="Co-crystal Structure",
                swissprot_a="P00001",
                swissprot_b="P00002",
            )
        ]
    )

    accepted, report = BioGRIDPairAcceptancePipeline().build(evidence)

    assert len(accepted) == 1
    assert report["accepted_pairs"] == 1
    assert set(accepted["uniprot_pair_id"]) == {"P00001--P00002"}


def test_pair_with_only_weak_evidence_is_not_accepted() -> None:
    evidence = pd.DataFrame(
        [
            _evidence_row(
                biogrid_interaction_id="1",
                experimental_system="Affinity Capture-MS",
                swissprot_a="P00001",
                swissprot_b="P00002",
            )
        ]
    )

    accepted, report = BioGRIDPairAcceptancePipeline().build(evidence)

    assert accepted.empty
    assert report["accepted_pairs"] == 0


def test_genetic_style_or_unknown_method_fails_closed() -> None:
    evidence = pd.DataFrame(
        [
            _evidence_row(
                biogrid_interaction_id="1",
                experimental_system="Some Brand New Method",
                swissprot_a="P00001",
                swissprot_b="P00002",
            )
        ]
    )

    accepted, report = BioGRIDPairAcceptancePipeline().build(evidence)

    assert accepted.empty
    assert report["unknown_experimental_systems"] == [
        "Some Brand New Method"
    ]


def test_unknown_experimental_systems_deduplicated_across_many_rows() -> None:
    """unknown_experimental_systems is computed from the distinct values
    (.unique()) rather than iterating every evidence row - verify repeated
    occurrences of the same unknown method across many pairs are still
    reported exactly once, and every occurrence is still caught."""
    evidence = pd.DataFrame(
        [
            _evidence_row(
                biogrid_interaction_id=str(index),
                experimental_system="Some Brand New Method",
                swissprot_a=f"P{index:05d}",
                swissprot_b=f"Q{index:05d}",
            )
            for index in range(10)
        ]
    )

    accepted, report = BioGRIDPairAcceptancePipeline().build(evidence)

    assert accepted.empty
    assert report["unknown_experimental_systems"] == [
        "Some Brand New Method"
    ]


def test_rna_evidence_never_becomes_a_ppi_pair() -> None:
    """Affinity Capture-RNA / Protein-RNA must never enter the accepted
    protein-protein pair set, no matter how strong the RNA-binding
    evidence looks."""
    evidence = pd.DataFrame(
        [
            _evidence_row(
                biogrid_interaction_id="1",
                experimental_system="Affinity Capture-RNA",
                swissprot_a="P00001",
                swissprot_b="P00002",
            ),
            _evidence_row(
                biogrid_interaction_id="2",
                experimental_system="Protein-RNA",
                swissprot_a="P00001",
                swissprot_b="P00002",
            ),
        ]
    )

    accepted, report = BioGRIDPairAcceptancePipeline().build(evidence)

    assert accepted.empty
    assert report["rna_evidence_rows_rejected"] == 2


def test_reverse_orientation_collapses_to_one_pair() -> None:
    """A-B and B-A rows for the same biological pair must not create two
    training examples."""
    evidence = pd.DataFrame(
        [
            _evidence_row(
                biogrid_interaction_id="1",
                experimental_system="Two-hybrid",
                swissprot_a="P00001",
                swissprot_b="P00002",
            ),
            _evidence_row(
                biogrid_interaction_id="2",
                experimental_system="Co-crystal Structure",
                swissprot_a="P00002",
                swissprot_b="P00001",
            ),
        ]
    )

    accepted, report = BioGRIDPairAcceptancePipeline().build(evidence)

    assert len(accepted) == 1
    assert report["accepted_pairs"] == 1
    assert accepted.iloc[0]["protein_a_id"] == "P00001"
    assert accepted.iloc[0]["protein_b_id"] == "P00002"


def test_multiple_evidence_rows_yield_one_pair_row_not_twenty() -> None:
    """A pair with many BioGRID records (weak and strong mixed) must still
    emit exactly one row - no groupby(...).agg('first') on raw evidence,
    but also no per-evidence-row duplication of the accepted pair."""
    evidence = pd.DataFrame(
        [
            _evidence_row(
                biogrid_interaction_id=str(index),
                experimental_system="Affinity Capture-MS",
                swissprot_a="P00001",
                swissprot_b="P00002",
            )
            for index in range(4)
        ]
        + [
            _evidence_row(
                biogrid_interaction_id="99",
                experimental_system="Co-crystal Structure",
                swissprot_a="P00001",
                swissprot_b="P00002",
            )
        ]
    )

    accepted, report = BioGRIDPairAcceptancePipeline().build(evidence)

    assert len(accepted) == 1
    assert (
        accepted.iloc[0]["representative_experimental_system"]
        == "Co-crystal Structure"
    )
    assert accepted.iloc[0]["representative_biogrid_interaction_id"] == "99"


def test_weak_evidence_does_not_hide_a_qualifying_row_in_the_same_group() -> None:
    """Weak rows sorted before the strong row must not suppress acceptance -
    qualification is decided per row, not by which row is encountered
    first in the group."""
    evidence = pd.DataFrame(
        [
            _evidence_row(
                biogrid_interaction_id="1",
                experimental_system="Co-localization",
                swissprot_a="P00001",
                swissprot_b="P00002",
            ),
            _evidence_row(
                biogrid_interaction_id="2",
                experimental_system="Co-purification",
                swissprot_a="P00001",
                swissprot_b="P00002",
            ),
            _evidence_row(
                biogrid_interaction_id="3",
                experimental_system="Thermal Shift Assay",
                swissprot_a="P00001",
                swissprot_b="P00002",
            ),
        ]
    )

    accepted, _ = BioGRIDPairAcceptancePipeline().build(evidence)

    assert len(accepted) == 1
    assert (
        accepted.iloc[0]["representative_experimental_system"]
        == "Thermal Shift Assay"
    )


def test_cross_taxonomy_pair_is_accepted_like_any_other() -> None:
    """Taxonomy is not an eligibility filter: this pipeline does not even
    see organism columns, so cross-species pairs are naturally eligible."""
    evidence = pd.DataFrame(
        [
            _evidence_row(
                biogrid_interaction_id="1",
                experimental_system="Two-hybrid",
                swissprot_a="P00001",
                swissprot_b="Q99999",
            )
        ]
    )

    accepted, report = BioGRIDPairAcceptancePipeline().build(evidence)

    assert len(accepted) == 1
    assert report["accepted_pairs"] == 1


def test_self_interaction_homodimer_is_accepted_documented_behavior() -> None:
    """Documents current, intentional behavior: a homodimer (BioGRID
    Interactor A == Interactor B, a real and common biological record -
    e.g. many proteins self-associate) is accepted like any other pair
    when it has qualifying evidence. It is not silently dropped."""
    evidence = pd.DataFrame(
        [
            _evidence_row(
                biogrid_interaction_id="1",
                experimental_system="Co-crystal Structure",
                swissprot_a="P00001",
                swissprot_b="P00001",
            )
        ]
    )

    accepted, report = BioGRIDPairAcceptancePipeline().build(evidence)

    assert len(accepted) == 1
    assert accepted.iloc[0]["protein_a_id"] == "P00001"
    assert accepted.iloc[0]["protein_b_id"] == "P00001"
    assert accepted.iloc[0]["uniprot_pair_id"] == "P00001--P00001"


def test_provenance_columns_are_present_but_minimal() -> None:
    evidence = pd.DataFrame(
        [
            _evidence_row(
                biogrid_interaction_id="1",
                experimental_system="Co-crystal Structure",
                swissprot_a="P00001",
                swissprot_b="P00002",
                publication_id="PUBMED:42",
            )
        ]
    )

    accepted, _ = BioGRIDPairAcceptancePipeline().build(evidence)

    row = accepted.iloc[0]

    assert row["representative_biogrid_interaction_id"] == "1"
    assert row["representative_experimental_system"] == "Co-crystal Structure"
    assert row["representative_publication_id"] == "PUBMED:42"
    assert row["biogrid_release_version"] == "5.0.260"

    assert "evidence_count" not in accepted.columns
    assert "publication_count" not in accepted.columns
    assert "experimental_system_count" not in accepted.columns


def test_representative_provenance_never_mixes_rows_via_null_columns() -> None:
    """Regression test: pandas' ``groupby(...).first()`` picks the first
    *non-null* value independently PER COLUMN, not the values of one
    consistent row. If a qualifying row has a null in one provenance
    column while a sibling qualifying row has a non-null value there, a
    naive ``.first()`` can silently stitch together provenance from two
    different BioGRID records (e.g. biogrid_interaction_id from row A
    combined with publication_id from row B). Every representative
    provenance field must come from exactly one physical evidence row."""
    evidence = pd.DataFrame(
        [
            _evidence_row(
                biogrid_interaction_id="100",
                experimental_system="Two-hybrid",
                swissprot_a="P00001",
                swissprot_b="P00002",
                publication_id=None,
            ),
            _evidence_row(
                biogrid_interaction_id="200",
                experimental_system="PCA",
                swissprot_a="P00001",
                swissprot_b="P00002",
                publication_id="PUBMED:99",
            ),
        ]
    )

    accepted, _ = BioGRIDPairAcceptancePipeline().build(evidence)

    assert len(accepted) == 1
    row = accepted.iloc[0]

    # "100" sorts before "200", so row "100" (Two-hybrid, null pub) is the
    # deterministic representative - its own null publication_id must be
    # preserved, not silently backfilled from row "200".
    assert row["representative_biogrid_interaction_id"] == "100"
    assert row["representative_experimental_system"] == "Two-hybrid"
    assert pd.isna(row["representative_publication_id"])


def test_missing_required_columns_are_rejected() -> None:
    with pytest.raises(ValueError, match="missing columns"):
        BioGRIDPairAcceptancePipeline().build(pd.DataFrame())


def test_repr() -> None:
    assert repr(BioGRIDPairAcceptancePipeline()) == "BioGRIDPairAcceptancePipeline()"
