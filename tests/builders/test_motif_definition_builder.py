import pytest

from protein_graph.builders.motif_definition_builder import (
    MotifDefinitionBuilder,
)


def test_build_creates_motif_definition() -> None:
    builder = MotifDefinitionBuilder()

    definition = builder.build(
        {
            "motif_id": "LIG_SH3_1",
            "motif_class": "SH3-binding motif",
            "pattern": r"P..P",
            "recognized_by": ("SH3",),
            "interaction_role": "ligand",
            "source": "ELM",
            "evidence_level": "curated",
        }
    )

    assert definition.motif_id == "LIG_SH3_1"
    assert definition.motif_class == "SH3-binding motif"
    assert definition.pattern == r"P..P"
    assert definition.recognized_by == ("SH3",)
    assert definition.interaction_role == "ligand"
    assert definition.source == "ELM"
    assert definition.evidence_level == "curated"


def test_build_uses_unknown_source_by_default() -> None:
    builder = MotifDefinitionBuilder()

    definition = builder.build(
        {
            "motif_id": "TEST",
            "motif_class": "Test motif",
            "pattern": "AAA",
        }
    )

    assert definition.source == "unknown"


def test_build_normalizes_comma_separated_recognized_by() -> None:
    builder = MotifDefinitionBuilder()

    definition = builder.build(
        {
            "motif_id": "TEST",
            "motif_class": "Test motif",
            "pattern": "AAA",
            "recognized_by": "SH3, WW, SH3",
        }
    )

    assert definition.recognized_by == (
        "SH3",
        "WW",
    )


def test_build_many_creates_definitions() -> None:
    builder = MotifDefinitionBuilder()

    definitions = builder.build_many(
        [
            {
                "motif_id": "MOTIF_A",
                "motif_class": "Motif A",
                "pattern": "A..A",
            },
            {
                "motif_id": "MOTIF_B",
                "motif_class": "Motif B",
                "pattern": "B..B",
            },
        ]
    )

    assert [
        definition.motif_id
        for definition in definitions
    ] == [
        "MOTIF_A",
        "MOTIF_B",
    ]


def test_build_rejects_missing_pattern() -> None:
    builder = MotifDefinitionBuilder()

    with pytest.raises(
        ValueError,
        match="pattern",
    ):
        builder.build(
            {
                "motif_id": "TEST",
                "motif_class": "Test motif",
            }
        )


def test_build_rejects_invalid_regex() -> None:
    builder = MotifDefinitionBuilder()

    with pytest.raises(
        ValueError,
        match="Invalid motif regex",
    ):
        builder.build(
            {
                "motif_id": "TEST",
                "motif_class": "Test motif",
                "pattern": "(",
            }
        )


def test_build_many_rejects_duplicate_ids() -> None:
    builder = MotifDefinitionBuilder()

    records = [
        {
            "motif_id": "TEST",
            "motif_class": "First",
            "pattern": "AAA",
        },
        {
            "motif_id": "TEST",
            "motif_class": "Second",
            "pattern": "BBB",
        },
    ]

    with pytest.raises(
        ValueError,
        match="Duplicate motif definition ID",
    ):
        builder.build_many(records)


def test_rejects_invalid_recognized_by_type() -> None:
    builder = MotifDefinitionBuilder()

    with pytest.raises(
        TypeError,
        match="recognized_by",
    ):
        builder.build(
            {
                "motif_id": "TEST",
                "motif_class": "Test motif",
                "pattern": "AAA",
                "recognized_by": 42,
            }
        )


def test_repr() -> None:
    builder = MotifDefinitionBuilder()

    assert repr(builder) == "MotifDefinitionBuilder()"
