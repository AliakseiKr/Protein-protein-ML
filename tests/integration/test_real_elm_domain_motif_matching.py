"""End-to-end test that real ELM data drives real bidirectional scanning.

Uses the actual cached ELM exports under data/raw/elm/ (elm_classes.tsv and
elm_interaction_domains.tsv, both fetched from the official ELM database -
see protein_graph.data_sources.elm and
protein_graph.data_sources.elm_interaction_domains) to prove the full chain
that used to be dead due to the hardcoded ``recognized_by = ()`` now
actually works with real domain<->motif knowledge:

    ELM motif classes + ELM interaction domains
        -> MotifDefinitionBuilder(recognized_by_overrides=...)
        -> matching_domain_motif adapter (entities.MotifDefinition
           -> matching_domain_motif.MotifDefinition scanner rule)
        -> DomainMotifRegistry
        -> KnownMotifScanner
        -> ProteinPairAnalyzer (bidirectional A<->B scanning)
"""

from pathlib import Path

import pytest

from protein_graph.data_sources.elm import ELMDataSource
from protein_graph.data_sources.elm_interaction_domains import (
    ELMInteractionDomainsDataSource,
    build_recognized_by_map,
)
from protein_graph.builders.motif_definition_builder import (
    MotifDefinitionBuilder,
)
from protein_graph.domain_detector.static_detector import (
    StaticDomainDetector,
)
from protein_graph.matching_domain_motif import (
    DetectedDomain,
    DomainMotifRegistry,
    KnownMotifScanner,
    ProteinInput,
    ProteinPairAnalyzer,
    to_scanner_rules,
)

REAL_ELM_RAW_DIR = Path("data/raw/elm")


def _require_cached_elm_data() -> None:
    classes_path = REAL_ELM_RAW_DIR / ELMDataSource.RAW_FILENAME
    domains_path = (
        REAL_ELM_RAW_DIR / ELMInteractionDomainsDataSource.RAW_FILENAME
    )

    if not classes_path.exists() or not domains_path.exists():
        pytest.skip(
            "Cached ELM raw files are not available at "
            f"{REAL_ELM_RAW_DIR} in this environment."
        )


def test_real_elm_domain_motif_data_drives_bidirectional_exact_scan() -> None:
    _require_cached_elm_data()

    motif_classes = ELMDataSource(raw_dir=REAL_ELM_RAW_DIR).load()
    interaction_domains = ELMInteractionDomainsDataSource(
        raw_dir=REAL_ELM_RAW_DIR
    ).load()

    recognized_by_map = build_recognized_by_map(interaction_domains)

    # Two real ELM motif classes with real, curated recognizing domains:
    # - CLV_NRD_NRD_1 (NRD cleavage site) <-> Pfam PF00675 (Peptidase_M16)
    # - CLV_PCSK_FUR_1 (Furin cleavage site) <-> Pfam PF00082 (Peptidase_S8)
    assert recognized_by_map["CLV_NRD_NRD_1"] == ("PF00675",)
    assert "PF00082" in recognized_by_map["CLV_PCSK_FUR_1"]

    general_definitions = MotifDefinitionBuilder().build_many(
        motif_classes,
        recognized_by_overrides=recognized_by_map,
    )

    definitions_by_id = {
        definition.motif_id: definition
        for definition in general_definitions
    }

    nrd_motif = definitions_by_id["CLV_NRD_NRD_1"]
    fur_motif = definitions_by_id["CLV_PCSK_FUR_1"]

    assert nrd_motif.pattern == "(.RK)|(RR[^KR])"
    assert nrd_motif.recognized_by == ("PF00675",)
    assert fur_motif.pattern == "R.[RK]R."
    assert "PF00082" in fur_motif.recognized_by

    scanner_rules = to_scanner_rules(general_definitions)
    registry = DomainMotifRegistry(scanner_rules)

    assert registry.supports_domain("PF00675")
    assert registry.supports_domain("PF00082")

    # Protein A carries domain PF00675 (recognizes CLV_NRD_NRD_1) AND its
    # own sequence contains a real CLV_PCSK_FUR_1 match ("RARRA"). Protein
    # B carries domain PF00082 (recognizes CLV_PCSK_FUR_1) AND its own
    # sequence contains a real CLV_NRD_NRD_1 match ("ARK"). One single
    # analyze() call therefore proves BOTH directions independently
    # produce a real, DIFFERENT candidate from real ELM data - not one
    # direction finding something while the other trivially finds nothing.
    protein_a_sequence = "GGGGRARRAGGGG"  # contains "RARRA" -> CLV_PCSK_FUR_1
    protein_b_sequence = "GGGGARKGGGG"  # contains "ARK" -> CLV_NRD_NRD_1

    domain_detector = StaticDomainDetector(
        {
            protein_a_sequence: (
                DetectedDomain(
                    domain_class="PF00675",
                    start=1,
                    end=4,
                    confidence=1.0,
                    domain_id="PF00675",
                ),
            ),
            protein_b_sequence: (
                DetectedDomain(
                    domain_class="PF00082",
                    start=1,
                    end=4,
                    confidence=1.0,
                    domain_id="PF00082",
                ),
            ),
        }
    )

    analyzer = ProteinPairAnalyzer(
        domain_detector=domain_detector,
        motif_scanner=KnownMotifScanner(registry),
        motif_registry=registry,
    )

    protein_a = ProteinInput(protein_id="PROT_A", sequence=protein_a_sequence)
    protein_b = ProteinInput(protein_id="PROT_B", sequence=protein_b_sequence)

    result = analyzer.analyze(protein_a, protein_b)

    # A -> B: A's domain PF00675 only searches for CLV_NRD_NRD_1, which B
    # does contain ("ARK").
    assert len(result.candidates_a_to_b) == 1
    candidate_a_to_b = result.candidates_a_to_b[0]
    assert candidate_a_to_b.motif_id == "CLV_NRD_NRD_1"
    assert candidate_a_to_b.target_sequence == "ARK"
    assert candidate_a_to_b.source_protein_id == "PROT_A"
    assert candidate_a_to_b.target_protein_id == "PROT_B"

    # B -> A: B's domain PF00082 searches only for the motifs it recognizes
    # (several related PCSK-family cleavage sites - PF00082/Peptidase_S8
    # is a real, documented case of "one domain recognizes multiple
    # motifs"). A's sequence ("RARRA") is a real match for CLV_PCSK_FUR_1,
    # proving B->A independently found a genuinely different real motif
    # than the A->B direction did (not a mirrored duplicate) - other real
    # PCSK-family motifs recognized by the same domain may also match the
    # same short sequence, which is expected, not a bug.
    candidate_motif_ids_b_to_a = {
        candidate.motif_id for candidate in result.candidates_b_to_a
    }
    assert "CLV_PCSK_FUR_1" in candidate_motif_ids_b_to_a

    fur_candidate = next(
        candidate
        for candidate in result.candidates_b_to_a
        if candidate.motif_id == "CLV_PCSK_FUR_1"
    )
    assert fur_candidate.target_sequence == "RARRA"
    assert fur_candidate.source_protein_id == "PROT_B"
    assert fur_candidate.target_protein_id == "PROT_A"

    assert result.has_bidirectional_candidates is True
