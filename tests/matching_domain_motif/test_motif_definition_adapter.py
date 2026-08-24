import dataclasses

import pytest

from protein_graph.entities.motif_definition import (
    MotifDefinition as GeneralMotifDefinition,
)
from protein_graph.matching_domain_motif.motif_definition_adapter import (
    to_scanner_rule,
    to_scanner_rules,
)


def test_converts_definition_with_recognized_domain() -> None:
    definition = GeneralMotifDefinition(
        motif_id="CLV_NRD_NRD_1",
        motif_class="Cleavage motif",
        pattern="R.R",
        recognized_by=("PF00675",),
        source="ELM",
    )

    rule = to_scanner_rule(definition)

    assert rule is not None
    assert rule.motif_id == "CLV_NRD_NRD_1"
    assert rule.motif_class == "Cleavage motif"
    assert rule.regex == "R.R"
    assert rule.recognized_by == ("PF00675",)


def test_definition_without_recognized_domain_is_not_convertible() -> None:
    """No fabricated placeholder domain is invented to force a rule."""
    definition = GeneralMotifDefinition(
        motif_id="TEST",
        motif_class="Test motif",
        pattern="AAA",
        source="ELM",
    )

    assert to_scanner_rule(definition) is None


def test_to_scanner_rules_drops_unconvertible_definitions() -> None:
    convertible = GeneralMotifDefinition(
        motif_id="M1",
        motif_class="Motif 1",
        pattern="AAA",
        recognized_by=("PF001",),
        source="ELM",
    )

    unconvertible = GeneralMotifDefinition(
        motif_id="M2",
        motif_class="Motif 2",
        pattern="BBB",
        source="ELM",
    )

    rules = to_scanner_rules([convertible, unconvertible])

    assert [rule.motif_id for rule in rules] == ["M1"]


def test_adapter_does_not_mutate_the_source_entity() -> None:
    """entities.MotifDefinition is a frozen dataclass, so mutation is
    structurally impossible - verify this directly (conversion leaves the
    original object's fields identical, and attempting to assign to it
    still raises)."""
    definition = GeneralMotifDefinition(
        motif_id="CLV_NRD_NRD_1",
        motif_class="Cleavage motif",
        pattern="R.R",
        recognized_by=("PF00675",),
        source="ELM",
    )

    to_scanner_rule(definition)

    assert definition.motif_id == "CLV_NRD_NRD_1"
    assert definition.motif_class == "Cleavage motif"
    assert definition.pattern == "R.R"
    assert definition.recognized_by == ("PF00675",)

    with pytest.raises(dataclasses.FrozenInstanceError):
        definition.recognized_by = ("PF99999",)  # type: ignore[misc]


def test_multiple_recognizing_domains_are_preserved_exactly() -> None:
    """One motif recognized by several domains must carry ALL of them
    through to the scanner rule, in order, none dropped or deduplicated
    incorrectly."""
    definition = GeneralMotifDefinition(
        motif_id="LIG_BIR_III_1",
        motif_class="IAP-binding motif",
        pattern="AAA",
        recognized_by=("PF00653", "PF99999", "PF00001"),
        source="ELM",
    )

    rule = to_scanner_rule(definition)

    assert rule is not None
    assert rule.recognized_by == ("PF00653", "PF99999", "PF00001")


def test_invalid_regex_is_rejected_at_source_construction() -> None:
    """A general MotifDefinition can never hold an invalid regex in the
    first place (entities.MotifDefinition validates re.compile(pattern) at
    construction) - the adapter can therefore never receive one. Documents
    that this failure happens predictably, at the source, not silently
    inside the adapter."""
    with pytest.raises(ValueError, match="Invalid motif regex"):
        GeneralMotifDefinition(
            motif_id="BAD",
            motif_class="Bad motif",
            pattern="(unclosed",
            recognized_by=("PF00001",),
            source="ELM",
        )
