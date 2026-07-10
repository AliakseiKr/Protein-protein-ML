from pathlib import Path

from protein_graph.data_sources.biogrid import BioGRIDDataSource


def main() -> None:
    source = BioGRIDDataSource(
        species="Rattus norvegicus",
        taxonomy_id=10116,
        interaction_type="physical",
        raw_dir=Path("data/raw/biogrid"),
    )

    print("=" * 80)
    print("Loading BioGRID...")
    print("=" * 80)

    data = source.load()

    print()
    print("=" * 80)
    print("SUCCESS")
    print("=" * 80)

    print(f"Rows: {len(data):,}")
    print(f"Columns: {len(data.columns)}")

    print("\nColumns:")
    for column in data.columns:
        print(f"  - {column}")

    print("\nFirst five rows:")
    print(data.head())

    print("\nData types:")
    print(data.dtypes)

    print("\nArchive location:")
    print(source.archive_path)

    print("\nArchive exists:")
    print(source.archive_path.exists())


if __name__ == "__main__":
    main()
