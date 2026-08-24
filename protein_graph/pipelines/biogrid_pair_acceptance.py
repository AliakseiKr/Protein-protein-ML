"""Turn raw BioGRID physical evidence rows into accepted strong/direct pairs.

BioGRID publishes individual experimental EVIDENCE ROWS. ProteinGraph needs
biological PROTEIN PAIRS: ``A-B`` and ``B-A`` are the same pair, and a pair
heavily studied across twenty BioGRID records must not implicitly receive
twenty times the training weight of a pair supported by a single strong
direct experiment.

This module implements that reduction as an explicit, order-sensitive
pipeline (see ``BioGRIDPairAcceptancePipeline.build``):

1. reject evidence rows that are not protein-protein (protein-RNA);
2. resolve interactors to canonical UniProt accessions;
3. build the canonical undirected ``uniprot_pair_id`` (``A-B == B-A``);
4. group every remaining evidence row by that pair id;
5. accept the pair iff at least one row in its group satisfies the
   reviewed strong/direct evidence policy
   (:mod:`protein_graph.data_sources.biogrid_evidence_policy`) - a single
   qualifying experiment is sufficient, no corroboration required;
6. emit exactly one row per accepted pair, carrying only the minimal
   provenance needed to answer "why was this pair admitted".

Evidence is never collapsed before qualification is decided (no
``groupby(...).agg("first")`` as the acceptance rule): a pair with four weak
co-complex rows and one strong direct row must still be accepted, and the
decision must be driven by inspecting each row's own experimental system,
not by which row happens to sort first.

A pair that fails to qualify is not "biologically negative" - it simply does
not satisfy the current strong-positive training criterion. Such pairs are
not persisted anywhere by this pipeline (see the project's absence-is-not-
evidence-of-non-interaction policy).
"""

from __future__ import annotations

from typing import Any

import pandas as pd

from protein_graph.data_sources.biogrid_evidence_policy import (
    EVIDENCE_POLICY_VERSION,
    RNA_EXPERIMENTAL_SYSTEMS,
    is_strong_direct_evidence,
    unknown_experimental_systems,
)
from protein_graph.mapping import BioGRIDUniProtMapper

#: Minimal provenance retained per accepted pair, answering only "why was
#: this pair admitted" - not biological/model features (see
#: PairFeaturePipeline.MODEL_FEATURE_COLUMNS, which never includes these).
_REPRESENTATIVE_EVIDENCE_COLUMNS: dict[str, str] = {
    "biogrid_interaction_id": "representative_biogrid_interaction_id",
    "experimental_system": "representative_experimental_system",
    "publication_id": "representative_publication_id",
    "biogrid_release_version": "biogrid_release_version",
}


class BioGRIDPairAcceptancePipeline:
    """Reduce BioGRID evidence rows to accepted strong/direct positive pairs."""

    def __init__(self, mapper: BioGRIDUniProtMapper | None = None) -> None:
        self.mapper = mapper or BioGRIDUniProtMapper()

    def build(
        self,
        evidence: pd.DataFrame,
    ) -> tuple[pd.DataFrame, dict[str, Any]]:
        """Build the accepted-pair table and a reproducibility report."""
        self._validate_inputs(evidence)

        total_evidence_rows = len(evidence)

        non_rna = evidence.loc[
            ~evidence["experimental_system"].isin(RNA_EXPERIMENTAL_SYSTEMS)
        ].copy()

        rna_evidence_rows_rejected = total_evidence_rows - len(non_rna)

        mapped = self.mapper.map(non_rna)
        mapping_coverage = self.mapper.coverage(mapped)

        complete = mapped.loc[mapped["mapping_complete"]].copy()

        # .unique() first: only the ~20 distinct Experimental System
        # values need checking against the policy table, not every one of
        # up to millions of evidence rows individually.
        unknown_systems = sorted(
            unknown_experimental_systems(
                complete["experimental_system"].unique()
            )
        )

        # Categorical, per-row qualification decision. Every row for a pair
        # is inspected here; nothing is discarded before this point.
        complete["is_strong_direct_evidence"] = complete[
            "experimental_system"
        ].map(is_strong_direct_evidence)

        candidate_pairs = complete["uniprot_pair_id"].nunique()

        accepted = self._accept_pairs(complete)

        report = {
            "evidence_policy_version": EVIDENCE_POLICY_VERSION,
            "biogrid_release_version": (
                complete["biogrid_release_version"].dropna().iloc[0]
                if not complete["biogrid_release_version"].dropna().empty
                else None
            ),
            "total_evidence_rows": total_evidence_rows,
            "rna_evidence_rows_rejected": rna_evidence_rows_rejected,
            **mapping_coverage,
            "evidence_rows_with_complete_mapping": len(complete),
            "unknown_experimental_systems": unknown_systems,
            "candidate_pairs_grouped": int(candidate_pairs),
            "accepted_pairs": len(accepted),
        }

        return accepted, report

    @staticmethod
    def _accept_pairs(complete: pd.DataFrame) -> pd.DataFrame:
        """Emit one row per pair with >=1 strong/direct evidence row.

        Filtering to qualifying rows *before* grouping is equivalent to
        "group everything, then check whether any row qualifies": row-level
        qualification never depends on sibling rows, so nothing about the
        rest of a pair's evidence is lost by this ordering. The subsequent
        row-selection step only breaks ties among rows that already
        qualify, to pick one deterministic representative for provenance -
        it is not the acceptance rule itself.

        Deliberately uses ``drop_duplicates(keep="first")`` on the sorted
        table rather than ``groupby(...).first()``: pandas' groupby
        ``first()`` picks the first *non-null value independently per
        column*, which can silently stitch together provenance from two
        *different* evidence rows within the same group whenever one row
        has a null in one provenance column and another row has a null in
        a different one (e.g. row A's biogrid_interaction_id paired with
        row B's publication_id). ``drop_duplicates`` instead keeps one
        single physical row intact, guaranteeing every representative
        provenance field originates from the same qualifying record.
        """
        qualifying = complete.loc[
            complete["is_strong_direct_evidence"]
        ].copy()

        if qualifying.empty:
            return pd.DataFrame(
                columns=[
                    "protein_a_id",
                    "protein_b_id",
                    "uniprot_pair_id",
                    *_REPRESENTATIVE_EVIDENCE_COLUMNS.values(),
                ]
            )

        qualifying = qualifying.sort_values(
            ["uniprot_pair_id", "biogrid_interaction_id"]
        )

        representative = qualifying.drop_duplicates(
            subset="uniprot_pair_id",
            keep="first",
        ).reset_index(drop=True)

        # Reorient to the canonical sorted pair order regardless of which
        # BioGRID row (and therefore which A/B orientation) was picked as
        # the representative, so protein_a_id/protein_b_id always agree
        # with uniprot_pair_id's own sorted "A--B" ordering.
        split = representative["uniprot_pair_id"].str.split(
            "--",
            n=1,
            expand=True,
        )

        accepted = pd.DataFrame(
            {
                "protein_a_id": split[0],
                "protein_b_id": split[1],
                "uniprot_pair_id": representative["uniprot_pair_id"],
            }
        )

        for source_column, output_column in (
            _REPRESENTATIVE_EVIDENCE_COLUMNS.items()
        ):
            accepted[output_column] = representative[source_column]

        return accepted.reset_index(drop=True)

    @staticmethod
    def _validate_inputs(evidence: pd.DataFrame) -> None:
        required_columns = {
            "experimental_system",
            "biogrid_interaction_id",
            "publication_id",
            "biogrid_release_version",
            "swissprot_accessions_a",
            "swissprot_accessions_b",
            "trembl_accessions_a",
            "trembl_accessions_b",
        }

        missing = required_columns - set(evidence.columns)

        if missing:
            raise ValueError(
                f"BioGRID evidence table is missing columns: {sorted(missing)}"
            )

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}()"
