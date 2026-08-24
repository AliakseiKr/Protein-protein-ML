"""Build the ProteinGraph global training dataset from BioGRID-ALL.

Usage:
    python scripts/build_training_dataset.py                     # full production build (latest release)
    python scripts/build_training_dataset.py --sample             # quick dev/debug run
    python scripts/build_training_dataset.py --release 5.0.255    # reproduce a build from a pinned release
"""

from __future__ import annotations

import argparse

from protein_graph.data_sources import BioGRIDDataSource
from protein_graph.pipelines import TrainingDatasetPipeline


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sample",
        action="store_true",
        help=(
            "Development/debug mode: only read the first 500 raw BioGRID "
            "rows. Production runs must omit this flag so the full "
            "selected BioGRID-ALL dataset is processed."
        ),
    )
    parser.add_argument(
        "--force-download",
        action="store_true",
        help="Re-download BioGRID even if a cached archive exists.",
    )
    parser.add_argument(
        "--release",
        default=None,
        help=(
            "Pin a specific BioGRID release (e.g. 5.0.255) instead of "
            "always tracking Latest-Release. Use this to reproduce an "
            "earlier build exactly."
        ),
    )
    args = parser.parse_args()

    sample_limit = 500 if args.sample else None

    biogrid = BioGRIDDataSource(
        raw_dir="data/raw/biogrid",
        interaction_type="physical",
        sample_limit=sample_limit,
        release=args.release,
    )

    print(f"Loading BioGRID-ALL interactions (sample_limit={sample_limit})...")

    evidence = biogrid.load(force_download=args.force_download)

    print(f"Normalized BioGRID physical evidence rows: {len(evidence):,}")
    print(f"BioGRID release version: {biogrid.release_version}")

    pipeline = TrainingDatasetPipeline()

    print("Qualifying pairs and downloading UniProt records...")

    pair_features, report = pipeline.build_and_save(
        evidence=evidence,
        output_path="data/processed/biogrid_pair_features.parquet",
        report_path="data/processed/biogrid_pair_features_report.json",
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
        "\ndata/processed/biogrid_pair_features.parquet"
    )

    print(
        "\nSaved report:"
        "\ndata/processed/biogrid_pair_features_report.json"
    )


if __name__ == "__main__":
    main()
