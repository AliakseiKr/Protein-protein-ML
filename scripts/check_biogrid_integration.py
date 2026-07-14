from protein_graph.data_sources import BioGRIDDataSource


def main() -> None:
    loader = BioGRIDDataSource()

    print(loader)

    interactions = loader.fetch()

    print()
    print(interactions.head())
    print()
    print(f"Shape: {interactions.shape}")


if __name__ == "__main__":
    main()
