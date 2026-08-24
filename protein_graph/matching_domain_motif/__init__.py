"""Domain -> motif scanning (exact-match matching layer).

See ``protein_graph.matching_domain_motif.motif_definition_adapter`` for
how general motif knowledge (``entities.motif_definition.MotifDefinition``)
becomes an executable scanner rule here.
"""

from protein_graph.matching_domain_motif.domain_motif import (
    DetectedDomain,
    DomainMotifCandidate,
    DomainMotifRegistry,
    KnownMotifScanner,
    MotifDefinition,
    MotifScanner,
)
from protein_graph.matching_domain_motif.motif_definition_adapter import (
    to_scanner_rule,
    to_scanner_rules,
)
from protein_graph.matching_domain_motif.pair_analyzer import (
    PairAnalysisResult,
    ProteinInput,
    ProteinPairAnalyzer,
)

__all__ = [
    "DetectedDomain",
    "DomainMotifCandidate",
    "DomainMotifRegistry",
    "KnownMotifScanner",
    "MotifDefinition",
    "MotifScanner",
    "PairAnalysisResult",
    "ProteinInput",
    "ProteinPairAnalyzer",
    "to_scanner_rule",
    "to_scanner_rules",
]
