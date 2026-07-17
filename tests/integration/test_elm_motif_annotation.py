from pathlib import Path

from protein_graph.builders.motif_definition_builder import (
    MotifDefinitionBuilder,
)
from protein_graph.builders.motif_builder import MotifBuilder
from protein_graph.data_sources.elm import ELMDataSource


def test_elm_definition_can_annotate_protein(
    tmp_path: Path,
) -> None:
    raw_dir = tmp_path / "elm"
    raw_dir.mkdir()

    source = ELMDataSource(raw_dir=raw_dir)

    source.raw_path.write_text(
        (
            "Accession\tFunctionalSiteName\tRegex\n"
            "LIG_SH3_1\tSH3-binding motif\tP..P\n"
        ),
        encoding="utf-8",
    )

    records = source.load()

    definitions = (
        MotifDefinitionBuilder().build_many(records)
    )

    motifs = MotifBuilder(definitions).build(
        {
            "protein_id": "P12345",
            "sequence": "AAAPLLPQQQ",
        }
    )

    assert len(motifs) == 1

    motif = motifs[0]

    assert motif.protein_id == "P12345"
    assert motif.motif_id == "LIG_SH3_1"
    assert motif.sequence == "PLLP"
    assert motif.start == 4
    assert motif.end == 7
