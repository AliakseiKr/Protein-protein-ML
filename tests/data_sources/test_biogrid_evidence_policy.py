"""Tests for the reviewed BioGRID strong/direct evidence policy."""

from protein_graph.data_sources.biogrid_evidence_policy import (
    RNA_EXPERIMENTAL_SYSTEMS,
    STRONG_DIRECT_EXPERIMENTAL_SYSTEMS,
    WEAK_PHYSICAL_EXPERIMENTAL_SYSTEMS,
    is_rna_evidence,
    is_strong_direct_evidence,
    unknown_experimental_systems,
)


def test_buckets_are_disjoint() -> None:
    assert not (
        STRONG_DIRECT_EXPERIMENTAL_SYSTEMS & WEAK_PHYSICAL_EXPERIMENTAL_SYSTEMS
    )
    assert not (STRONG_DIRECT_EXPERIMENTAL_SYSTEMS & RNA_EXPERIMENTAL_SYSTEMS)
    assert not (WEAK_PHYSICAL_EXPERIMENTAL_SYSTEMS & RNA_EXPERIMENTAL_SYSTEMS)


def test_reviewed_systems_match_official_biogrid_physical_list() -> None:
    """Every official BioGRID physical Experimental System from the task
    specification must be classified into exactly one reviewed bucket."""
    official_physical_systems = {
        "Affinity Capture-Luminescence",
        "Affinity Capture-MS",
        "Affinity Capture-RNA",
        "Affinity Capture-Western",
        "Biochemical Activity",
        "Co-crystal Structure",
        "Co-fractionation",
        "Co-localization",
        "Co-purification",
        "Cross-Linking-MS (XL-MS)",
        "Far Western",
        "FRET",
        "PCA",
        "Protein-peptide",
        "Protein-RNA",
        "Proximity Label-MS",
        "Reconstituted Complex",
        "Two-hybrid",
        "Surface Display",
        "Thermal Shift Assay",
    }

    reviewed = (
        RNA_EXPERIMENTAL_SYSTEMS
        | STRONG_DIRECT_EXPERIMENTAL_SYSTEMS
        | WEAK_PHYSICAL_EXPERIMENTAL_SYSTEMS
    )

    assert official_physical_systems == reviewed


def test_co_crystal_structure_qualifies_as_strong_direct() -> None:
    assert is_strong_direct_evidence("Co-crystal Structure") is True


def test_two_hybrid_qualifies_as_strong_direct() -> None:
    assert is_strong_direct_evidence("Two-hybrid") is True


def test_affinity_capture_ms_does_not_qualify() -> None:
    """A pulldown identified by MS is co-complex evidence, not direct
    binding: it must not automatically qualify as strong/direct."""
    assert is_strong_direct_evidence("Affinity Capture-MS") is False


def test_proximity_label_ms_does_not_qualify() -> None:
    """Proximity does not automatically equal direct binding."""
    assert is_strong_direct_evidence("Proximity Label-MS") is False


def test_co_localization_does_not_qualify() -> None:
    assert is_strong_direct_evidence("Co-localization") is False


def test_far_western_qualifies_as_strong_direct() -> None:
    assert is_strong_direct_evidence("Far Western") is True


def test_thermal_shift_assay_qualifies_as_strong_direct() -> None:
    assert is_strong_direct_evidence("Thermal Shift Assay") is True


def test_pca_qualifies_as_strong_direct() -> None:
    assert is_strong_direct_evidence("PCA") is True


def test_protein_peptide_qualifies_as_strong_direct() -> None:
    assert is_strong_direct_evidence("Protein-peptide") is True


def test_biochemical_activity_qualifies_as_strong_direct() -> None:
    assert is_strong_direct_evidence("Biochemical Activity") is True


def test_surface_display_qualifies_as_strong_direct() -> None:
    assert is_strong_direct_evidence("Surface Display") is True


def test_fret_does_not_qualify() -> None:
    """Policy v2 correction: BioGRID's own wording is "close proximity ...
    in vivo", not "direct contact" - unlike Two-hybrid/PCA, a positive FRET
    signal does not require the two full proteins to physically join for
    the readout to occur, only the fluorophores to be within the Forster
    radius, which two proteins docked on a shared scaffold can satisfy
    without touching. Fails closed for the strict training set."""
    assert is_strong_direct_evidence("FRET") is False


def test_reconstituted_complex_does_not_qualify() -> None:
    """Policy v2 correction: BioGRID's own definition names a GST pull-down
    against a whole cell lysate as an example of "Reconstituted Complex" -
    the same bridging risk as Affinity Capture-MS/-Western - alongside
    genuinely clean two-purified-component methods (SPR/BLI/gel-shift), with
    no structured field distinguishing the two. Fails closed for the strict
    training set."""
    assert is_strong_direct_evidence("Reconstituted Complex") is False


def test_rna_evidence_never_qualifies_as_strong_direct() -> None:
    for system in RNA_EXPERIMENTAL_SYSTEMS:
        assert is_strong_direct_evidence(system) is False


def test_affinity_capture_rna_is_rna_evidence() -> None:
    assert is_rna_evidence("Affinity Capture-RNA") is True
    assert is_rna_evidence("Protein-RNA") is True


def test_protein_protein_evidence_is_not_rna_evidence() -> None:
    assert is_rna_evidence("Two-hybrid") is False


def test_unknown_system_fails_closed() -> None:
    """A BioGRID Experimental System value not yet reviewed must not
    qualify as strong/direct - it fails closed rather than being guessed."""
    assert is_strong_direct_evidence("Some Brand New BioGRID Method") is False


def test_unknown_experimental_systems_reports_unreviewed_values() -> None:
    systems = [
        "Two-hybrid",
        "Affinity Capture-MS",
        "Some Brand New BioGRID Method",
    ]

    assert unknown_experimental_systems(systems) == {
        "Some Brand New BioGRID Method"
    }


def test_unknown_experimental_systems_empty_when_all_reviewed() -> None:
    assert unknown_experimental_systems(["Two-hybrid", "Co-purification"]) == set()


def test_none_experimental_system_is_not_strong_direct() -> None:
    assert is_strong_direct_evidence(None) is False
