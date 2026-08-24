from pathlib import Path

import pytest

from protein_graph.data_sources.elm_interaction_domains import (
    ELMInteractionDomainsDataSource,
    build_recognized_by_map,
)


INTERACTION_DOMAINS_TSV = (
    '"ELM identifier"\t"Interaction Domain Id"\t'
    '"Interaction Domain Description"\t"Interaction Domain Name"\n'
    '"CLV_NRD_NRD_1"\t"PF00675"\t"Peptidase_M16"\t'
    '"Insulinase (Peptidase family M16)"\n'
    '"CLV_PCSK_FUR_1"\t"PF00082"\t"Peptidase_S8"\t"Subtilase family"\n'
    '"CLV_PCSK_PC1ET2_1"\t"PF00082"\t"Peptidase_S8"\t"Subtilase family"\n'
)


def write_cached_file(
    source: ELMInteractionDomainsDataSource,
    content: str = INTERACTION_DOMAINS_TSV,
) -> None:
    source.raw_dir.mkdir(parents=True, exist_ok=True)
    source.raw_path.write_text(content, encoding="utf-8")


def test_loads_cached_interaction_domains(tmp_path: Path) -> None:
    source = ELMInteractionDomainsDataSource(raw_dir=tmp_path / "elm")
    write_cached_file(source)

    records = source.load()

    assert records == [
        {
            "motif_id": "CLV_NRD_NRD_1",
            "domain_id": "PF00675",
            "domain_namespace": "Pfam",
            "domain_description": "Peptidase_M16",
            "domain_name": "Insulinase (Peptidase family M16)",
            "source": "ELM",
        },
        {
            "motif_id": "CLV_PCSK_FUR_1",
            "domain_id": "PF00082",
            "domain_namespace": "Pfam",
            "domain_description": "Peptidase_S8",
            "domain_name": "Subtilase family",
            "source": "ELM",
        },
        {
            "motif_id": "CLV_PCSK_PC1ET2_1",
            "domain_id": "PF00082",
            "domain_namespace": "Pfam",
            "domain_description": "Peptidase_S8",
            "domain_name": "Subtilase family",
            "source": "ELM",
        },
    ]


def test_duplicate_motif_domain_pairs_are_deduplicated(
    tmp_path: Path,
) -> None:
    source = ELMInteractionDomainsDataSource(raw_dir=tmp_path / "elm")
    write_cached_file(
        source,
        INTERACTION_DOMAINS_TSV
        + '"CLV_NRD_NRD_1"\t"PF00675"\t"Peptidase_M16"\t'
        '"Insulinase (Peptidase family M16)"\n',
    )

    records = source.load()

    assert len(records) == 3


def test_missing_local_file_and_no_url_raises(tmp_path: Path) -> None:
    source = ELMInteractionDomainsDataSource(
        raw_dir=tmp_path / "elm",
        url=None,
    )

    with pytest.raises(RuntimeError, match="URL is not configured"):
        source.load()


def test_build_recognized_by_map_preserves_multiplicity() -> None:
    records = [
        {"motif_id": "CLV_NRD_NRD_1", "domain_id": "PF00675"},
        {"motif_id": "CLV_PCSK_FUR_1", "domain_id": "PF00082"},
        {"motif_id": "CLV_PCSK_FUR_1", "domain_id": "PF99999"},
    ]

    recognized_by = build_recognized_by_map(records)

    assert recognized_by == {
        "CLV_NRD_NRD_1": ("PF00675",),
        "CLV_PCSK_FUR_1": ("PF00082", "PF99999"),
    }


def test_build_recognized_by_map_deduplicates_repeat_domain() -> None:
    records = [
        {"motif_id": "M1", "domain_id": "PF1"},
        {"motif_id": "M1", "domain_id": "PF1"},
    ]

    assert build_recognized_by_map(records) == {"M1": ("PF1",)}


def test_repr(tmp_path: Path) -> None:
    source = ELMInteractionDomainsDataSource(raw_dir=tmp_path / "elm")

    assert "ELMInteractionDomainsDataSource" in repr(source)


def test_domain_namespace_detects_pfam_smart_and_interpro(
    tmp_path: Path,
) -> None:
    """Regression test: real ELM interaction-domains data is NOT
    exclusively Pfam despite the column name - 66/409 real rows are SMART
    or InterPro accessions. The namespace must be detected from the
    identifier's own prefix, not hardcoded to "Pfam" for every row."""
    source = ELMInteractionDomainsDataSource(raw_dir=tmp_path / "elm")
    write_cached_file(
        source,
        '"ELM identifier"\t"Interaction Domain Id"\t'
        '"Interaction Domain Description"\t"Interaction Domain Name"\n'
        '"DEG_APCC_TPR_1"\t"SM00028"\t"TPR"\t"Tetratricopeptide repeats"\n'
        '"LIG_WD40_WDR5_WIN_1"\t"IPR017986"\t"IPR017986"\t"WD40 repeat"\n'
        '"CLV_NRD_NRD_1"\t"PF00675"\t"Peptidase_M16"\t"Insulinase"\n'
        '"DEG_TEST_UNKNOWN"\t"XX99999"\t"Unknown"\t"Unknown scheme"\n',
    )

    records = source.load()
    by_motif = {record["motif_id"]: record for record in records}

    assert by_motif["DEG_APCC_TPR_1"]["domain_namespace"] == "SMART"
    assert by_motif["LIG_WD40_WDR5_WIN_1"]["domain_namespace"] == "InterPro"
    assert by_motif["CLV_NRD_NRD_1"]["domain_namespace"] == "Pfam"
    assert by_motif["DEG_TEST_UNKNOWN"]["domain_namespace"] == "unknown"


def test_real_elm_data_domain_namespaces_are_all_classified(
    tmp_path: Path,
) -> None:
    """Every real ELM interaction-domains row must resolve to a known
    namespace (Pfam/SMART/InterPro), not silently fall through to
    "unknown", confirming the detection covers the real data fully."""
    real_dir = Path("data/raw/elm")

    if not (real_dir / ELMInteractionDomainsDataSource.RAW_FILENAME).exists():
        pytest.skip("Cached real ELM interaction-domains file not available.")

    records = ELMInteractionDomainsDataSource(raw_dir=real_dir).load()

    namespaces = {record["domain_namespace"] for record in records}
    assert namespaces == {"Pfam", "SMART", "InterPro"}
    assert "unknown" not in namespaces
