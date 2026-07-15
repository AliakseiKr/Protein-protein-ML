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

        result["same_taxonomy"] = (
            result["a_taxon_id"] == result["b_taxon_id"]
        ).astype(int)

        return result

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
