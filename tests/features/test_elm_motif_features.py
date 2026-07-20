from typing import Any

import pytest

from protein_graph.entities.motif import MotifRecord
from protein_graph.features.elm_motifs import (
    MotifAnnotationFeatureExtractor,
    MotifFeatureContext,
)


def make_motif(**overrides: Any) -> MotifRecord:
    values: dict[str, Any] = {
        "protein_id": "P12345",
        "region_type": "motif",
        "start": 11,
        "end": 14,
        "sequence": "PQQP",
        "motif_id": "LIG_SH3_GENERIC",
        "motif_class": "SH3-binding motif",
        "pattern": r"P..P",
        "recognized_by": ("SH3",),
        "interaction_role": None,
        "evidence_source": "ELM",
        "evidence_level": None,
        "confidence": None,
    }

    values.update(overrides)

    return MotifRecord(**values)

from protein_graph.features.elm_motifs import (
    MotifAnnotationFeatureExtractor,
    MotifFeatureContext,
    MotifGeometryFeatureExtractor,
)


def test_extracts_complete_motif_annotations() -> None:
    motif = make_motif(
        recognized_by=(
            "SH3",
            "WW",
        ),
        interaction_role="ligand",
        evidence_source="ELM",
        evidence_level="curated",
        confidence=0.8,
    )

    features = MotifAnnotationFeatureExtractor().extract(
        motif=motif,
        context=MotifFeatureContext(
            protein_length=100
        ),
    )

    assert features["has_motif_id"] == 1
    assert features["has_motif_class"] == 1
    assert features["has_motif_pattern"] == 1

    assert features["has_recognized_by"] == 1
    assert features["has_interaction_role"] == 1

    assert features["has_evidence_source"] == 1
    assert features["has_evidence_level"] == 1
    assert features["has_motif_confidence"] == 1

    assert features["num_recognized_domains"] == 2
    assert features["num_motif_annotations"] == 8

    assert features[
        "motif_confidence"
    ] == pytest.approx(0.8)


def test_extracts_partial_motif_annotations() -> None:
    motif = make_motif(
        motif_class=None,
        pattern=None,
        recognized_by=(),
        interaction_role=None,
        evidence_source="ELM",
        evidence_level=None,
        confidence=None,
    )

    features = MotifAnnotationFeatureExtractor().extract(
        motif=motif,
        context=MotifFeatureContext(
            protein_length=100
        ),
    )

    assert features["has_motif_id"] == 1
    assert features["has_motif_class"] == 0
    assert features["has_motif_pattern"] == 0

    assert features["has_recognized_by"] == 0
    assert features["has_interaction_role"] == 0

    assert features["has_evidence_source"] == 1
    assert features["has_evidence_level"] == 0
    assert features["has_motif_confidence"] == 0

    assert features["num_recognized_domains"] == 0
    assert features["num_motif_annotations"] == 2
    assert features["motif_confidence"] == 0.0


def test_zero_confidence_is_still_present() -> None:
    motif = make_motif(
        confidence=0.0,
    )

    features = MotifAnnotationFeatureExtractor().extract(
        motif=motif,
        context=MotifFeatureContext(
            protein_length=100
        ),
    )

    assert features["has_motif_confidence"] == 1
    assert features["motif_confidence"] == 0.0
