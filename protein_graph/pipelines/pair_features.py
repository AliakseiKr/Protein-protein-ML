"""Build interaction pair feature tables."""

from __future__ import annotations

import pandas as pd


class PairFeaturePipeline:
    """Combine two protein feature vectors into one pair row."""

    def build(
        self,
        interactions: pd.DataFrame,
        protein_features: pd.DataFrame,
    ) -> pd.DataFrame:
        """Join protein features to mapped interaction pairs."""
        self._validate_inputs(interactions, protein_features)

        protein_a_features = protein_features.add_prefix("a_")
        protein_b_features = protein_features.add_prefix("b_")

        pairs = interactions.merge(
            protein_a_features,
            left_on="protein_a_id",
            right_on="a_protein_id",
            how="inner",
        )

        pairs = pairs.merge(
            protein_b_features,
            left_on="protein_b_id",
            right_on="b_protein_id",
            how="inner",
        )

        pairs = self._add_pair_features(pairs)

        return pairs.reset_index(drop=True)

    @staticmethod
    def _add_pair_features(pairs: pd.DataFrame) -> pd.DataFrame:
        """Add numerical compatibility features."""
        result = pairs.copy()

        result["protein_length_difference"] = (
            result["a_protein_length"]
            - result["b_protein_length"]
        ).abs()

        result["hydrophobic_fraction_difference"] = (
            result["a_hydrophobic_fraction"]
            - result["b_hydrophobic_fraction"]
        ).abs()

        result["charged_fraction_difference"] = (
            result["a_charged_fraction"]
            - result["b_charged_fraction"]
        ).abs()

        result["positive_negative_balance"] = (
            result["a_positive_fraction"]
            * result["b_negative_fraction"]
            + result["b_positive_fraction"]
            * result["a_negative_fraction"]
        )

        taxonomy_known = (
            result["a_taxon_id"].notna()
            & result["b_taxon_id"].notna()
        )

        result["same_taxonomy"] = (
            (result["a_taxon_id"] == result["b_taxon_id"])
            .where(taxonomy_known)
            .astype("Int64")
        )

        return result

    #: Unprefixed per-protein sequence-feature column names produced by
    #: the extractors wired into ProteinFeaturePipeline
    #: (protein_graph.features.sequence: AminoAcidCompositionExtractor,
    #: ResidueFractionExtractor, LocalDensityExtractor,
    #: ConsecutiveRunExtractor, PeriodicityExtractor). Captured by
    #: executing each extractor against a full-alphabet sample sequence
    #: and reading back its real output keys - not hand-transcribed -
    #: so this list cannot contain a typo'd or invented name.
    _SEQUENCE_FEATURE_COLUMNS: frozenset[str] = frozenset(
        {
            "protein_length",
            "aa_fraction_A",
            "aa_fraction_C",
            "aa_fraction_D",
            "aa_fraction_E",
            "aa_fraction_F",
            "aa_fraction_G",
            "aa_fraction_H",
            "aa_fraction_I",
            "aa_fraction_K",
            "aa_fraction_L",
            "aa_fraction_M",
            "aa_fraction_N",
            "aa_fraction_P",
            "aa_fraction_Q",
            "aa_fraction_R",
            "aa_fraction_S",
            "aa_fraction_T",
            "aa_fraction_V",
            "aa_fraction_W",
            "aa_fraction_Y",
            "hydrophobic_fraction",
            "hydrophilic_fraction",
            "charged_fraction",
            "positive_fraction",
            "negative_fraction",
            "polar_fraction",
            "aromatic_fraction",
            "high_hydrophobic_density_window_count_5",
            "high_hydrophobic_density_window_count_10",
            "high_hydrophobic_density_window_count_20",
            "high_hydrophobic_density_window_count_30",
            "high_charged_density_window_count_5",
            "high_charged_density_window_count_10",
            "high_charged_density_window_count_20",
            "high_charged_density_window_count_30",
            "high_positive_density_window_count_5",
            "high_positive_density_window_count_10",
            "high_positive_density_window_count_20",
            "high_positive_density_window_count_30",
            "high_negative_density_window_count_5",
            "high_negative_density_window_count_10",
            "high_negative_density_window_count_20",
            "high_negative_density_window_count_30",
            "max_hydrophobic_run",
            "hydrophobic_run_0_5",
            "hydrophobic_run_5_10",
            "hydrophobic_run_10_15",
            "hydrophobic_run_15_20",
            "hydrophobic_run_20_25",
            "hydrophobic_run_25_30",
            "hydrophobic_run_30_35",
            "hydrophobic_run_35_40",
            "hydrophobic_run_gt_40",
            "max_charged_run",
            "charged_run_0_5",
            "charged_run_5_10",
            "charged_run_10_15",
            "charged_run_15_20",
            "charged_run_20_25",
            "charged_run_25_30",
            "charged_run_30_35",
            "charged_run_35_40",
            "charged_run_gt_40",
            "max_positive_run",
            "positive_run_0_5",
            "positive_run_5_10",
            "positive_run_10_15",
            "positive_run_15_20",
            "positive_run_20_25",
            "positive_run_25_30",
            "positive_run_30_35",
            "positive_run_35_40",
            "positive_run_gt_40",
            "max_negative_run",
            "negative_run_0_5",
            "negative_run_5_10",
            "negative_run_10_15",
            "negative_run_15_20",
            "negative_run_20_25",
            "negative_run_25_30",
            "negative_run_30_35",
            "negative_run_35_40",
            "negative_run_gt_40",
            "hydrophobic_periodic_run_step_2_count",
            "longest_hydrophobic_periodic_run_step_2",
            "hydrophobic_periodic_run_step_2_count_3_to_5",
            "hydrophobic_periodic_run_step_2_count_6_to_8",
            "hydrophobic_periodic_run_step_2_count_9_to_12",
            "hydrophobic_periodic_run_step_2_count_13_to_16",
            "hydrophobic_periodic_run_step_2_count_17_to_19",
            "hydrophobic_periodic_run_step_2_count_ge_20",
            "hydrophobic_periodic_run_step_3_count",
            "longest_hydrophobic_periodic_run_step_3",
            "hydrophobic_periodic_run_step_3_count_3_to_5",
            "hydrophobic_periodic_run_step_3_count_6_to_8",
            "hydrophobic_periodic_run_step_3_count_9_to_12",
            "hydrophobic_periodic_run_step_3_count_13_to_16",
            "hydrophobic_periodic_run_step_3_count_17_to_19",
            "hydrophobic_periodic_run_step_3_count_ge_20",
            "hydrophobic_periodic_run_step_4_count",
            "longest_hydrophobic_periodic_run_step_4",
            "hydrophobic_periodic_run_step_4_count_3_to_5",
            "hydrophobic_periodic_run_step_4_count_6_to_8",
            "hydrophobic_periodic_run_step_4_count_9_to_12",
            "hydrophobic_periodic_run_step_4_count_13_to_16",
            "hydrophobic_periodic_run_step_4_count_17_to_19",
            "hydrophobic_periodic_run_step_4_count_ge_20",
            "hydrophobic_periodic_run_step_7_count",
            "longest_hydrophobic_periodic_run_step_7",
            "hydrophobic_periodic_run_step_7_count_3_to_5",
            "hydrophobic_periodic_run_step_7_count_6_to_8",
            "hydrophobic_periodic_run_step_7_count_9_to_12",
            "hydrophobic_periodic_run_step_7_count_13_to_16",
            "hydrophobic_periodic_run_step_7_count_17_to_19",
            "hydrophobic_periodic_run_step_7_count_ge_20",
            "charged_periodic_run_step_2_count",
            "longest_charged_periodic_run_step_2",
            "charged_periodic_run_step_2_count_3_to_5",
            "charged_periodic_run_step_2_count_6_to_8",
            "charged_periodic_run_step_2_count_9_to_12",
            "charged_periodic_run_step_2_count_13_to_16",
            "charged_periodic_run_step_2_count_17_to_19",
            "charged_periodic_run_step_2_count_ge_20",
            "charged_periodic_run_step_3_count",
            "longest_charged_periodic_run_step_3",
            "charged_periodic_run_step_3_count_3_to_5",
            "charged_periodic_run_step_3_count_6_to_8",
            "charged_periodic_run_step_3_count_9_to_12",
            "charged_periodic_run_step_3_count_13_to_16",
            "charged_periodic_run_step_3_count_17_to_19",
            "charged_periodic_run_step_3_count_ge_20",
            "charged_periodic_run_step_4_count",
            "longest_charged_periodic_run_step_4",
            "charged_periodic_run_step_4_count_3_to_5",
            "charged_periodic_run_step_4_count_6_to_8",
            "charged_periodic_run_step_4_count_9_to_12",
            "charged_periodic_run_step_4_count_13_to_16",
            "charged_periodic_run_step_4_count_17_to_19",
            "charged_periodic_run_step_4_count_ge_20",
            "charged_periodic_run_step_7_count",
            "longest_charged_periodic_run_step_7",
            "charged_periodic_run_step_7_count_3_to_5",
            "charged_periodic_run_step_7_count_6_to_8",
            "charged_periodic_run_step_7_count_9_to_12",
            "charged_periodic_run_step_7_count_13_to_16",
            "charged_periodic_run_step_7_count_17_to_19",
            "charged_periodic_run_step_7_count_ge_20",
        }
    )

    #: Explicit allowlist of columns approved as ML model input: every
    #: a_*/b_* sequence-feature column derived above, plus the pair-level
    #: features computed in `_add_pair_features`. Any column NOT listed
    #: here - whether a raw BioGRID/provenance column, an identity
    #: column, or a brand-new column nobody has reviewed yet - is
    #: excluded by default. A new legitimate feature must be added here
    #: deliberately; nothing is included automatically.
    MODEL_FEATURE_COLUMNS: frozenset[str] = frozenset(
        {f"a_{name}" for name in _SEQUENCE_FEATURE_COLUMNS}
        | {f"b_{name}" for name in _SEQUENCE_FEATURE_COLUMNS}
        | {
            "protein_length_difference",
            "hydrophobic_fraction_difference",
            "charged_fraction_difference",
            "positive_negative_balance",
            "same_taxonomy",
        }
    )

    @classmethod
    def select_model_features(
        cls,
        pairs: pd.DataFrame,
    ) -> pd.DataFrame:
        """Return only columns approved as ML model input.

        Uses MODEL_FEATURE_COLUMNS as an explicit allowlist: any column
        not listed - regardless of its source - is excluded by default,
        including any future raw interaction/provenance/identity column
        that has not been deliberately reviewed and added.
        """
        feature_columns = [
            column
            for column in pairs.columns
            if column in cls.MODEL_FEATURE_COLUMNS
        ]

        return pairs[feature_columns]

    @staticmethod
    def _validate_inputs(
        interactions: pd.DataFrame,
        protein_features: pd.DataFrame,
    ) -> None:
        interaction_columns = {
            "protein_a_id",
            "protein_b_id",
        }

        missing_interaction_columns = (
            interaction_columns - set(interactions.columns)
        )

        if missing_interaction_columns:
            raise ValueError(
                "Interaction table is missing columns: "
                f"{sorted(missing_interaction_columns)}"
            )

        protein_columns = {
            "protein_id",
            "protein_length",
            "hydrophobic_fraction",
            "charged_fraction",
            "positive_fraction",
            "negative_fraction",
            "taxon_id",
        }

        missing_protein_columns = (
            protein_columns - set(protein_features.columns)
        )

        if missing_protein_columns:
            raise ValueError(
                "Protein feature table is missing columns: "
                f"{sorted(missing_protein_columns)}"
            )

        if protein_features["protein_id"].duplicated().any():
            raise ValueError(
                "Protein feature table contains duplicate protein IDs."
            )

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}()"
