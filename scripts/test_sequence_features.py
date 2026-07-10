"""Integration test for sequence-derived features."""

from protein_graph.features.sequence import (
    AminoAcidCompositionExtractor,
    ConsecutiveRunExtractor,
    LocalDensityExtractor,
    PeriodicityExtractor,
    ResidueFractionExtractor,
)


def main() -> None:
    protein = {
        "protein_id": "TEST_PROTEIN",
        "sequence": (
            "MKKKKK"
            "AVILMFWYAVILMFWYAVIL"
            "DEDEDE"
            "LXXLXXLXXL"
            "STNQCY"
        ),
    }

    extractors = [
        AminoAcidCompositionExtractor(),
        ResidueFractionExtractor(),
        LocalDensityExtractor(),
        ConsecutiveRunExtractor(),
        PeriodicityExtractor(),
    ]

    features = {
        "protein_id": protein["protein_id"],
    }

    for extractor in extractors:
        features.update(extractor.extract(protein))

    for name, value in sorted(features.items()):
        print(f"{name}: {value}")

    assert features["protein_length"] == len(protein["sequence"])
    assert 0.0 <= features["hydrophobic_fraction"] <= 1.0
    assert features["max_hydrophobic_run"] > 0
    assert features["high_hydrophobic_density_window_count_5"] > 0
    assert features["longest_hydrophobic_periodic_run_step_3"] > 0

    print("\nAll sequence feature checks passed.")


if __name__ == "__main__":
    main()
