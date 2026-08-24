from protein_graph.data_sources import BioGRIDDataSource


def main() -> None:
    loader = BioGRIDDataSource(sample_limit=500)

    print(loader)

    interactions = loader.load()

    print()
    print(interactions.head())
    print()
    print(f"Shape: {interactions.shape}")


if __name__ == "__main__":
    main()
