"""Tests for protein-level ELM motif profile extraction."""

import pandas as pd
import pytest

from protein_graph.features.elm_profile import (
    ELMMotifExtractor,
)


def test_extracts_aggregated_elm_profile() -> None:
    motifs = pd.DataFrame(
        [
            {
                "elm_id": "LIG_SH3_1",
                "motif_class": "LIG_SH3",
                "regex": r"P..P",
            },
            {
                "elm_id": "MOD_CK2_1",
                "motif_class": "MOD_CK2",
                "regex": r"S..E",
            },
            {
                "elm_id": "TRG_NLS_1",
                "motif_class": "TRG_NLS",
                "regex": r"KKKK",
            },
        ]
    )

    extractor = ELMMotifExtractor(
        motifs=motifs
    )

    features = extractor.extract(
        {
            "sequence": (
                "AAPQQPAA"
                "SAAE"
                "AAPRRPA"
                "SQQE"
                "SYYE"
            ),
        }
    )

    assert features == {
        "elm_total_motif_hits": 5,
        "elm_unique_motif_classes": 2,
        "elm_matched_motif_ids": [
            "LIG_SH3_1",
            "MOD_CK2_1",
        ],
        "elm_matched_motif_classes": [
            "LIG_SH3",
            "MOD_CK2",
        ],
        "elm_motif_class_counts": {
            "LIG_SH3": 2,
            "MOD_CK2": 3,
        },
    }


def test_counts_multiple_occurrences_of_same_motif() -> None:
    motifs = pd.DataFrame(
        [
            {
                "elm_id": "MOD_CK2_1",
                "motif_class": "MOD_CK2",
                "regex": r"S..E",
            },
        ]
    )

    extractor = ELMMotifExtractor(
        motifs=motifs
    )

    features = extractor.extract(
        {
            "sequence": "SAAESQQESYYE",
        }
    )

    assert features == {
        "elm_total_motif_hits": 3,
        "elm_unique_motif_classes": 1,
        "elm_matched_motif_ids": [
            "MOD_CK2_1",
        ],
        "elm_matched_motif_classes": [
            "MOD_CK2",
        ],
        "elm_motif_class_counts": {
            "MOD_CK2": 3,
        },
    }


def test_combines_ids_from_same_motif_class() -> None:
    motifs = pd.DataFrame(
        [
            {
                "elm_id": "LIG_SH3_1",
                "motif_class": "LIG_SH3",
                "regex": r"P..P",
            },
            {
                "elm_id": "LIG_SH3_2",
                "motif_class": "LIG_SH3",
                "regex": r"R..P",
            },
        ]
    )

    extractor = ELMMotifExtractor(
        motifs=motifs
    )

    features = extractor.extract(
        {
            "sequence": "PQQPAARSSP",
        }
    )

    assert features == {
        "elm_total_motif_hits": 2,
        "elm_unique_motif_classes": 1,
        "elm_matched_motif_ids": [
            "LIG_SH3_1",
            "LIG_SH3_2",
        ],
        "elm_matched_motif_classes": [
            "LIG_SH3",
        ],
        "elm_motif_class_counts": {
            "LIG_SH3": 2,
        },
    }


def test_returns_empty_profile_when_no_motifs_match() -> None:
    motifs = pd.DataFrame(
        [
            {
                "elm_id": "LIG_SH3_1",
                "motif_class": "LIG_SH3",
                "regex": r"P..P",
            },
        ]
    )

    extractor = ELMMotifExtractor(
        motifs=motifs
    )

    features = extractor.extract(
        {
            "sequence": "AAAAAAAAAA",
        }
    )

    assert features == {
        "elm_total_motif_hits": 0,
        "elm_unique_motif_classes": 0,
        "elm_matched_motif_ids": [],
        "elm_matched_motif_classes": [],
        "elm_motif_class_counts": {},
    }


def test_matches_lowercase_protein_sequence() -> None:
    motifs = pd.DataFrame(
        [
            {
                "elm_id": "LIG_SH3_1",
                "motif_class": "LIG_SH3",
                "regex": r"P..P",
            },
        ]
    )

    extractor = ELMMotifExtractor(
        motifs=motifs
    )

    features = extractor.extract(
        {
            "sequence": "aapqqpaa",
        }
    )

    assert features[
        "elm_matched_motif_ids"
    ] == [
        "LIG_SH3_1",
    ]

    assert features[
        "elm_total_motif_hits"
    ] == 1


def test_derives_motif_class_from_elm_id() -> None:
    motifs = pd.DataFrame(
        [
            {
                "elm_id": "DOC_MAPK_gen_1",
                "regex": r"K..K",
            },
            {
                "elm_id": "TRG_NLS_MonoExtN_4",
                "regex": r"KKKK",
            },
        ]
    )

    extractor = ELMMotifExtractor(
        motifs=motifs
    )

    features = extractor.extract(
        {
            "sequence": "KAAKAAKKKK",
        }
    )

    assert features[
        "elm_matched_motif_classes"
    ] == [
        "DOC_MAPK",
        "TRG_NLS",
    ]


def test_prefers_explicit_motif_class_column() -> None:
    motifs = pd.DataFrame(
        [
            {
                "elm_id": "CUSTOM_IDENTIFIER_1",
                "motif_class": "CUSTOM_CLASS",
                "regex": r"P..P",
            },
        ]
    )

    extractor = ELMMotifExtractor(
        motifs=motifs
    )

    features = extractor.extract(
        {
            "sequence": "PQQP",
        }
    )

    assert features[
        "elm_matched_motif_classes"
    ] == [
        "CUSTOM_CLASS",
    ]


def test_skips_invalid_regex() -> None:
    motifs = pd.DataFrame(
        [
            {
                "elm_id": "INVALID_MOTIF",
                "motif_class": "INVALID",
                "regex": r"[",
            },
            {
                "elm_id": "LIG_SH3_1",
                "motif_class": "LIG_SH3",
                "regex": r"P..P",
            },
        ]
    )

    extractor = ELMMotifExtractor(
        motifs=motifs
    )

    features = extractor.extract(
        {
            "sequence": "PQQP",
        }
    )

    assert features[
        "elm_matched_motif_ids"
    ] == [
        "LIG_SH3_1",
    ]

    assert features[
        "elm_total_motif_hits"
    ] == 1


def test_skips_rows_with_missing_values() -> None:
    motifs = pd.DataFrame(
        [
            {
                "elm_id": None,
                "motif_class": None,
                "regex": r"P..P",
            },
            {
                "elm_id": "EMPTY_REGEX",
                "motif_class": "EMPTY",
                "regex": None,
            },
            {
                "elm_id": "VALID_MOTIF",
                "motif_class": "VALID",
                "regex": r"K..K",
            },
        ]
    )

    extractor = ELMMotifExtractor(
        motifs=motifs
    )

    features = extractor.extract(
        {
            "sequence": "PQQPKAAK",
        }
    )

    assert features[
        "elm_matched_motif_ids"
    ] == [
        "VALID_MOTIF",
    ]


def test_rejects_empty_sequence() -> None:
    motifs = pd.DataFrame(
        [
            {
                "elm_id": "LIG_SH3_1",
                "regex": r"P..P",
            },
        ]
    )

    extractor = ELMMotifExtractor(
        motifs=motifs
    )

    with pytest.raises(
        ValueError,
        match="Protein sequence is empty",
    ):
        extractor.extract(
            {
                "sequence": "",
            }
        )


def test_rejects_missing_sequence() -> None:
    motifs = pd.DataFrame(
        [
            {
                "elm_id": "LIG_SH3_1",
                "regex": r"P..P",
            },
        ]
    )

    extractor = ELMMotifExtractor(
        motifs=motifs
    )

    with pytest.raises(
        ValueError,
        match="Protein sequence is empty",
    ):
        extractor.extract({})


def test_rejects_missing_required_columns() -> None:
    motifs = pd.DataFrame(
        [
            {
                "elm_id": "LIG_SH3_1",
            },
        ]
    )

    with pytest.raises(
        ValueError,
        match="missing required columns: regex",
    ):
        ELMMotifExtractor(
            motifs=motifs
        )
