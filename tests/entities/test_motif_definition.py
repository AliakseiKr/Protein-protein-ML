import pytest

from protein_graph.entities.motif_definition import (
    MotifDefinition,
)


def test_creates_valid_motif_definition() -> None:
    definition = MotifDefinition(
        motif_id="LIG_SH3_1",
        motif_class="SH3-binding motif",
        pattern=r"P..P",
        recognized_by=("SH3",),
        source="ELM",
    )

    assert definition.motif_id == "LIG_SH3_1"
    assert definition.motif_class == "SH3-binding motif"
    assert definition.pattern == r"P..P"
    assert definition.recognized_by == ("SH3",)
    assert definition.source == "ELM"


def test_strips_string_fields() -> None:
    definition = MotifDefinition(
        motif_id=" TEST ",
        motif_class=" Test motif ",
        pattern=" A..A ",
        source=" local ",
        interaction_role=" ligand ",
        evidence_level=" curated ",
    )

    assert definition.motif_id == "TEST"
    assert definition.motif_class == "Test motif"
    assert definition.pattern == "A..A"
    assert definition.source == "local"
    assert definition.interaction_role == "ligand"
    assert definition.evidence_level == "curated"


def test_rejects_invalid_regex() -> None:
    with pytest.raises(
        ValueError,
        match="Invalid motif regex",
    ):
        MotifDefinition(
            motif_id="BROKEN",
            motif_class="Broken motif",
            pattern="(",
        )


def test_rejects_empty_motif_class() -> None:
    with pytest.raises(
        ValueError,
        match="motif_class cannot be empty",
    ):
        MotifDefinition(
            motif_id="TEST",
            motif_class=" ",
            pattern="AAA",
        )


def test_normalizes_recognized_by() -> None:
    definition = MotifDefinition(
        motif_id="TEST",
        motif_class="Test motif",
        pattern="AAA",
        recognized_by=(
            " SH3 ",
            "",
            "SH3",
            " WW ",
        ),
    )

    assert definition.recognized_by == (
        "SH3",
        "WW",
    )


def test_rejects_string_recognized_by() -> None:
    with pytest.raises(
        TypeError,
        match="tuple of strings",
    ):
        MotifDefinition(
            motif_id="TEST",
            motif_class="Test motif",
            pattern="AAA",
            recognized_by="SH3",  # type: ignore[arg-type]
        )
