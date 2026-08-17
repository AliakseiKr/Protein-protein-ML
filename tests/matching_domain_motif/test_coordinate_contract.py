from protein_graph.builders.motif_builder import MotifBuilder
from protein_graph.entities.motif_definition import (
    MotifDefinition as CanonicalMotif,
)
from protein_graph.matching_domain_motif.domain_motif import (
    DetectedDomain,
    DomainMotifRegistry,
    KnownMotifScanner,
    MotifDefinition as MatchingMotif,
)


def test_scanner_and_builder_agree_on_coordinates() -> None:
    sequence = "AAPXXPAA"

    canonical = CanonicalMotif(
        motif_id="M1",
        motif_class="TEST",
        pattern="P..P",
        recognized_by=("SH3",),
    )
    built = MotifBuilder([canonical]).build(
        {"protein_id": "TARGET", "sequence": sequence}
    )

    matching = MatchingMotif(
        motif_id="M1",
        motif_class="TEST",
        recognized_by=("SH3",),
        regex="P..P",
    )
    scanner = KnownMotifScanner(DomainMotifRegistry([matching]))
    candidates = scanner.scan(
        domain=DetectedDomain(
            domain_class="SH3", start=1, end=2, confidence=1.0,
        ),
        source_protein_id="SOURCE",
        target_protein_id="TARGET",
        target_sequence=sequence,
    )

    assert built[0].start == 3
    assert built[0].end == 6

    assert candidates[0].target_start == 3
    assert candidates[0].target_end == 6

    assert (built[0].start, built[0].end) == (
        candidates[0].target_start,
        candidates[0].target_end,
    )
