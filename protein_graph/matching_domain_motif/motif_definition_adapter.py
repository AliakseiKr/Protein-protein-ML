"""Adapter from general motif knowledge to scanner-ready motif rules.

Bridges :class:`protein_graph.entities.motif_definition.MotifDefinition`
(general biological motif-class knowledge, used for protein-level motif
annotation - see ``features.elm_motifs`` / ``features.motif_pipeline``) to
:class:`protein_graph.matching_domain_motif.domain_motif.MotifDefinition`
(an executable scanner rule, used by ``DomainMotifRegistry`` /
``KnownMotifScanner`` to scan a specific partner protein for a motif known
to be recognized by a specific domain class).

These are deliberately two different dataclasses, not accidental
duplicates (see the Phase 10 investigation in the project history):

- ``entities.motif_definition.MotifDefinition`` describes a motif on its
  own (``pattern``, ``source``, ``evidence_level``); ``recognized_by`` may
  legitimately be empty when no recognizing domain is known yet - the
  motif is still valid general knowledge.
- ``matching_domain_motif.domain_motif.MotifDefinition`` additionally
  *requires* a non-empty ``recognized_by``, because a scanner rule with no
  domain to trigger it is not executable. Its ``regex`` field is the same
  underlying string as ``pattern``, renamed to signal "this is the
  compiled, matching-ready form", not a second, independently-derived
  pattern.

A general motif definition converts to a scanner rule only when it already
has a non-empty ``recognized_by`` (populated from real authoritative
domain<->motif knowledge - see
``protein_graph.data_sources.elm_interaction_domains``). Motifs without a
known recognizing domain are simply not convertible: the general
definition remains valid, it just cannot become a scanner rule yet.
Fabricating a placeholder domain to force a conversion would misrepresent
unreviewed biology as known, which this adapter never does.
"""

from __future__ import annotations

from collections.abc import Iterable

from protein_graph.entities.motif_definition import (
    MotifDefinition as GeneralMotifDefinition,
)
from protein_graph.matching_domain_motif.domain_motif import (
    MotifDefinition as ScannerMotifDefinition,
)


def to_scanner_rule(
    definition: GeneralMotifDefinition,
) -> ScannerMotifDefinition | None:
    """Convert one general motif definition into a scanner rule.

    Returns ``None`` when ``definition.recognized_by`` is empty: there is
    no recognizing domain to key a scanner rule on.
    """
    if not definition.recognized_by:
        return None

    return ScannerMotifDefinition(
        motif_id=definition.motif_id,
        motif_class=definition.motif_class,
        recognized_by=definition.recognized_by,
        regex=definition.pattern,
    )


def to_scanner_rules(
    definitions: Iterable[GeneralMotifDefinition],
) -> list[ScannerMotifDefinition]:
    """Convert many general motif definitions, dropping unconvertible ones."""
    rules: list[ScannerMotifDefinition] = []

    for definition in definitions:
        rule = to_scanner_rule(definition)

        if rule is not None:
            rules.append(rule)

    return rules
