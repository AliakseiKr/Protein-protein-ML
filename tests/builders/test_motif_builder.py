from dataclasses import dataclass

import pytest

from protein_graph.builders.motif_builder import MotifBuilder
from protein_graph.entities.motif_definition import MotifDefinition


@dataclass
class ProteinStub:
    protein_id: str
    sequence: str


@pytest.fixture
def sh3_definition() -> MotifDefinition:
    return MotifDefinition(
        motif_id="LIG_SH3_GENERIC",
        motif_class="SH3-binding motif",
        pattern=r"P..P",
        recognized_by=("SH3",),
        interaction_role="ligand",
        source="ELM",
    )


def test_build_finds_concrete_motif_sequence(
    sh3_definition,
) -> None:
    protein = ProteinStub(
        protein_id="P12345",
        sequence="AAAPLLPQQQ",
    )

    builder = MotifBuilder([sh3_definition])

    motifs = builder.build(protein)

    assert len(motifs) == 1

    motif = motifs[0]

    assert motif.protein_id == "P12345"
    assert motif.motif_id == "LIG_SH3_GENERIC"
    assert motif.sequence == "PLLP"
    assert motif.matched_sequence == "PLLP"
    assert motif.start == 4
    assert motif.end == 7
    assert motif.recognized_by == ("SH3",)


def test_build_returns_empty_list_when_no_match(
    sh3_definition,
) -> None:
    protein = ProteinStub(
        protein_id="P12345",
        sequence="AAAAAAAAAA",
    )

    builder = MotifBuilder([sh3_definition])

    assert builder.build(protein) == []


def test_build_finds_overlapping_matches() -> None:
    definition = MotifDefinition(
        motif_id="TEST",
        motif_class="Test motif",
        pattern=r"AAA",
        source="test",
    )

    protein = ProteinStub(
        protein_id="P12345",
        sequence="AAAA",
    )

    builder = MotifBuilder(
        [definition],
        allow_overlapping=True,
    )

    motifs = builder.build(protein)

    assert [motif.start for motif in motifs] == [1, 2]
    assert [motif.end for motif in motifs] == [3, 4]


def test_build_can_disable_overlapping_matches() -> None:
    definition = MotifDefinition(
        motif_id="TEST",
        motif_class="Test motif",
        pattern=r"AAA",
        source="test",
    )

    protein = ProteinStub(
        protein_id="P12345",
        sequence="AAAA",
    )

    builder = MotifBuilder(
        [definition],
        allow_overlapping=False,
    )

    motifs = builder.build(protein)

    assert len(motifs) == 1
    assert motifs[0].start == 1
    assert motifs[0].end == 3


def test_invalid_regex_cannot_reach_motif_builder() -> None:
    """Invalid motif definitions must fail before MotifBuilder."""
    with pytest.raises(
        ValueError,
        match="Invalid motif regex",
    ):
        MotifDefinition(
            motif_id="BROKEN",
            motif_class="Broken motif",
            pattern=r"[",
        )



def test_missing_protein_sequence_is_rejected(
    sh3_definition,
) -> None:
    protein = ProteinStub(
        protein_id="P12345",
        sequence="",
    )

    builder = MotifBuilder([sh3_definition])

    with pytest.raises(
        ValueError,
        match="sequence",
    ):
        builder.build(protein)


def test_repr(
    sh3_definition,
) -> None:
    builder = MotifBuilder([sh3_definition])

    assert repr(builder) == (
        "MotifBuilder("
        "definitions=1, "
        "allow_overlapping=True"
        ")"
    )
