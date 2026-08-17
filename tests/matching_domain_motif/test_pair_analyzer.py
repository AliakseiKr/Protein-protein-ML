from protein_graph.domain_detector.static_detector import (
    StaticDomainDetector,
)
from protein_graph.matching_domain_motif.domain_motif import (
    DetectedDomain,
    DomainMotifRegistry,
    KnownMotifScanner,
    MotifDefinition,
)
from protein_graph.matching_domain_motif.pair_analyzer import (
    ProteinInput,
    ProteinPairAnalyzer,
)

def build_analyzer(
    domains_by_sequence,
    *,
    threshold: float = 0.5,
):
    motifs = [
        MotifDefinition(
            motif_id="LIG_SH3_TEST",
            motif_class="LIG_SH3",
            recognized_by=("SH3",),
            regex=r"P..P",
        ),
        MotifDefinition(
            motif_id="LIG_PDZ_TEST",
            motif_class="LIG_PDZ",
            recognized_by=("PDZ",),
            regex=r"[ST].[VIL]$",
        ),
    ]

    registry = DomainMotifRegistry(motifs)

    return ProteinPairAnalyzer(
        domain_detector=StaticDomainDetector(
            domains_by_sequence
        ),
        motif_scanner=KnownMotifScanner(registry),
        motif_registry=registry,
        domain_confidence_threshold=threshold,
    )


def test_scans_partner_for_recognized_domain():
    sequence_a = "MAAAAAAAAAAAAAAAAAAAA"
    sequence_b = "MSSSRPLPPLPRGG"

    analyzer = build_analyzer(
        {
            sequence_a: [
                DetectedDomain(
                    domain_class="SH3",
                    start=2,
                    end=18,
                    confidence=0.9,
                )
            ],
            sequence_b: [],
        }
    )

    result = analyzer.analyze(
        ProteinInput("A", sequence_a),
        ProteinInput("B", sequence_b),
    )

    assert result.candidates_a_to_b
    assert not result.candidates_b_to_a

    candidate = result.candidates_a_to_b[0]

    assert candidate.source_protein_id == "A"
    assert candidate.target_protein_id == "B"
    assert candidate.domain_class == "SH3"
    assert candidate.motif_class == "LIG_SH3"
    assert candidate.matched_known_regex is True
    assert candidate.combined_score == 0.9


def test_does_not_scan_for_unknown_domain():
    sequence_a = "MAAAAAAAAAAAAAAAAAAAA"
    sequence_b = "MSSSRPLPPLPRGG"

    analyzer = build_analyzer(
        {
            sequence_a: [
                DetectedDomain(
                    domain_class="UNKNOWN_DOMAIN",
                    start=2,
                    end=18,
                    confidence=0.99,
                )
            ],
            sequence_b: [],
        }
    )

    result = analyzer.analyze(
        ProteinInput("A", sequence_a),
        ProteinInput("B", sequence_b),
    )

    assert not result.candidates_a_to_b
    assert len(result.skipped_domains_a) == 1


def test_does_not_scan_low_confidence_domain():
    sequence_a = "MAAAAAAAAAAAAAAAAAAAA"
    sequence_b = "MSSSRPLPPLPRGG"

    analyzer = build_analyzer(
        {
            sequence_a: [
                DetectedDomain(
                    domain_class="SH3",
                    start=2,
                    end=18,
                    confidence=0.3,
                )
            ],
            sequence_b: [],
        },
        threshold=0.5,
    )

    result = analyzer.analyze(
        ProteinInput("A", sequence_a),
        ProteinInput("B", sequence_b),
    )

    assert not result.candidates_a_to_b
    assert len(result.skipped_domains_a) == 1


def test_scans_in_both_directions():
    # PLAP соответствует SH3 regex P..P.
    # STAV на C-конце соответствует PDZ regex [ST].[VIL]$.
    sequence_a = "MAAARPLAPAAAASTAV"

    # STLV на C-конце соответствует PDZ regex [ST].[VIL]$.
    sequence_b = "MBBBRPLPBBBBSTLV"

    analyzer = build_analyzer(
        {
            sequence_a: [
                DetectedDomain(
                    domain_class="PDZ",
                    start=1,
                    end=10,
                    confidence=0.8,
                )
            ],
            sequence_b: [
                DetectedDomain(
                    domain_class="SH3",
                    start=2,
                    end=12,
                    confidence=0.9,
                )
            ],
        }
    )

    result = analyzer.analyze(
        ProteinInput("A", sequence_a),
        ProteinInput("B", sequence_b),
    )

    assert result.candidates_a_to_b
    assert result.candidates_b_to_a
    assert result.has_bidirectional_candidates is True

    assert result.candidates_a_to_b[0].domain_class == "PDZ"
    assert result.candidates_a_to_b[0].motif_class == "LIG_PDZ"

    assert result.candidates_b_to_a[0].domain_class == "SH3"
    assert result.candidates_b_to_a[0].motif_class == "LIG_SH3"
    assert result.candidates_b_to_a[0].target_sequence == "PLAP"


def test_no_domains_does_not_mean_negative_interaction():
    sequence_a = "MAAAAAAAAA"
    sequence_b = "MBBBBBBBBB"

    analyzer = build_analyzer(
        {
            sequence_a: [],
            sequence_b: [],
        }
    )

    result = analyzer.analyze(
        ProteinInput("A", sequence_a),
        ProteinInput("B", sequence_b),
    )

    assert result.domains_a == ()
    assert result.domains_b == ()
    assert result.num_candidates == 0
    assert result.max_combined_score == 0.0
