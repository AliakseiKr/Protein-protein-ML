"""Build the first real ProteinGraph training dataset."""

from __future__ import annotations

from protein_graph.data_sources import BioGRIDDataSource
from protein_graph.pipelines import TrainingDatasetPipeline


def main() -> None:
    biogrid = BioGRIDDataSource(
        raw_dir="data/raw/biogrid",
        species="Rattus norvegicus",
        taxonomy_id=10116,
        interaction_type="physical",
    )

    print("Loading BioGRID interactions...")

    interactions = biogrid.load().head(20)

    print(f"Normalized BioGRID rows: {len(interactions):,}")

    pipeline = TrainingDatasetPipeline()

    print("Mapping proteins and downloading UniProt records...")

    pair_features, report = pipeline.build_and_save(
        interactions=interactions,
        output_path=(
            "data/processed/"
            "biogrid_rat_pair_features.parquet"
        ),
        report_path=(
            "data/processed/"
            "biogrid_rat_pair_features_report.json"
        ),
        force_download=False,
    )

    print("\nDataset report")

    for name, value in report.items():
        print(f"{name}: {value}")

    print("\nResult")
    print(f"Shape: {pair_features.shape}")
    print(f"Columns: {len(pair_features.columns)}")

    print(
        "\nSaved dataset:"
        "\ndata/processed/biogrid_rat_pair_features.parquet"
    )

    print(
        "\nSaved report:"
        "\ndata/processed/biogrid_rat_pair_features_report.json"
    )


if __name__ == "__main__":
    main()
