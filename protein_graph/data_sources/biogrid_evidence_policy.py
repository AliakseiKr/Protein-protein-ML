"""Reviewed BioGRID physical Experimental System evidence policy.

ProteinGraph's strong/direct positive-pair training set must contain only
protein-protein interaction evidence strong enough to support a direct or
binary interaction claim. BioGRID's own "physical" Experimental System Type
mixes methods with very different evidentiary strength: a co-crystal
structure and a bait pulldown identified by mass spectrometry are both
"physical", but only one of them demonstrates direct binding.

This module is the single reviewed source of truth for that distinction. It
is a categorical policy (qualifies / does not qualify), not a confidence
score: BioGRID does not publish per-method reliability numbers, so inventing
one (e.g. "Affinity Capture-MS = 0.6") would not be a reviewed biological
fact, it would be a guess.

Reviewed against the official, VERBATIM BioGRID Experimental Evidence Codes
documentation (https://wiki.thebiogrid.org/doku.php/experimental_systems,
page last modified 2026-06-01), re-audited on 2026-08-24 (policy v2). The
2026-08-24 (v1) policy conflated "physical", "binary", "proximity" and
"direct molecular contact" for two methods; see the FRET and Reconstituted
Complex entries below for the correction and its reasoning.

The acceptance decision below is this project's own policy judgment applied
to BioGRID's definitions, not a BioGRID-endorsed classification. BioGRID
itself explicitly warns against assuming its records are direct/physical:
"It should not be assumed that the interaction reported in BioGRID is direct
and physical in nature; the experimental system definitions ... indicate the
nature of the supporting evidence." BioGRID also uses "binary relationship"
to mean "the database record links exactly two entities" - not "direct
biophysical contact" - and this policy does not conflate the two meanings.

Every BioGRID physical Experimental System value falls into exactly one of
three reviewed buckets:

- ``RNA_EXPERIMENTAL_SYSTEMS``: the interactor pair is protein-RNA, not
  protein-protein. These rows must never enter the PPI training set
  regardless of how strong the RNA-binding evidence is.
- ``STRONG_DIRECT_EXPERIMENTAL_SYSTEMS``: methods whose BioGRID definition
  structurally guarantees (by the mechanics of the assay itself, not
  merely by being labeled "physical") that a positive result reflects a
  direct contact/binding event between exactly the two recorded proteins,
  with no plausible route for an uninvolved third protein to produce the
  same signal. A method is placed here only when EVERY realistic
  execution of the assay, as BioGRID itself describes it, meets that bar -
  a category that can cover both direct and indirect executions is not
  placed here even if some fraction of its real-world use is direct,
  unless a separate structured BioGRID field reliably distinguishes the
  cases (none does, for any method reviewed below).
- ``WEAK_PHYSICAL_EXPERIMENTAL_SYSTEMS``: methods whose BioGRID definition
  demonstrates co-complex membership, spatial/photophysical proximity, or
  an indirect pulldown/mixture-based assay, without structurally
  guaranteeing that the two proteins contact each other directly. Real,
  valid physical-interaction evidence - just not strong enough for this
  project's current strong/direct training criterion. Ambiguous methods
  (capable of both direct and indirect execution, per BioGRID's own text)
  are deliberately placed here rather than in the strong bucket: this is a
  fail-closed choice for the current strict training goal, not a claim
  that the method is biologically weak evidence in general.

Any Experimental System value BioGRID adds later and not yet reviewed here
is neither of the above: it fails closed (``is_strong_direct_evidence``
returns False) and is reported by ``unknown_experimental_systems`` so it can
be reviewed and explicitly classified, rather than silently guessed at.
"""

from __future__ import annotations

from collections.abc import Iterable

#: Version of this reviewed evidence policy, bumped whenever a bucket
#: assignment changes. Recorded in build reports (see
#: BioGRIDPairAcceptancePipeline) so a dataset build stays reproducible
#: even if the policy itself is later revised.
#:
#: v1 (2026-08-24): initial policy.
#: v2 (2026-08-24): re-audit. Reclassified FRET and Reconstituted Complex
#: from strong/direct to weak/ambiguous - see module docstring and the
#: WEAK_PHYSICAL_EXPERIMENTAL_SYSTEMS entries below for the reasoning.
EVIDENCE_POLICY_VERSION = "2026-08-24.v2"

#: Protein-RNA evidence. BioGRID lists these under "physical", but the
#: interaction is between a protein and an RNA molecule, not two proteins.
#: They must never become protein-protein training examples.
#:
#: - Affinity Capture-RNA: in vivo protein-RNA binding evidence.
#: - Protein-RNA: in vitro protein-RNA binding evidence.
RNA_EXPERIMENTAL_SYSTEMS: frozenset[str] = frozenset(
    {
        "Affinity Capture-RNA",
        "Protein-RNA",
    }
)

#: Methods whose BioGRID definition structurally guarantees a direct/binary
#: protein-protein contact strong enough for the current strong-positive
#: training/reference set. A single row using one of these methods is
#: sufficient to admit a pair - no corroboration is required.
#:
#: - Co-crystal Structure: BioGRID's own text says the interaction is
#:   "directly demonstrated at the atomic level by X-ray crystallography"
#:   (also used for NMR/EM). Atomic-resolution contact is the strongest
#:   possible direct-binding evidence; no mixture/lysate scenario applies.
#: - Two-hybrid: reporter-gene activation requires the DBD-bait and
#:   TAD-prey fusions to be close enough to functionally reconstitute a
#:   transcription factor. This structural requirement - not merely the
#:   word "physical" - is why the field treats Y2H as the canonical
#:   "direct, binary" PPI assay, in contrast to co-complex methods like
#:   AP-MS. (Rare bridging by an endogenous yeast protein is a recognized
#:   Y2H artifact, but it is not how BioGRID's own definition describes
#:   the assay, and the same "reconstitution requires proximity" logic
#:   applies to PCA below.)
#: - PCA (protein-fragment complementation assay): same reconstitution
#:   logic as Two-hybrid - complementary reporter fragments must be
#:   brought together by the bait/prey fusions to fold into an active
#:   reporter (e.g. split-ubiquitin, BiFC).
#: - Far Western: a bait is immobilized at one specific position (after
#:   SDS-PAGE transfer) and probed with a defined, detectable prey
#:   protein; the signal is read at that specific position, not "whatever
#:   sticks to the membrane" - there is no third-protein mixture that
#:   could bridge the recorded pair.
#: - Protein-peptide: binding is tested against a defined peptide, not a
#:   whole-protein mixture, so there is no route for an uninvolved third
#:   protein to bridge the interaction the way a lysate pulldown could.
#: - Thermal Shift Assay: BioGRID's text describes a shift in a purified
#:   target protein's thermal stability upon ligand binding (Nano-DSF,
#:   DSC); a two-component biophysical measurement, no lysate/mixture
#:   context described.
#: - Surface Display: a displayed protein is directly assayed for binding
#:   to a specific candidate partner/ligand; the recorded interaction
#:   reflects a targeted binding test, not incidental co-occurrence.
#: - Biochemical Activity: BioGRID's text is "the biochemical effect of
#:   one protein upon another in vitro ... the bait protein executes the
#:   activity on the substrate hit protein" (e.g. kinase phosphorylation,
#:   GEF-mediated GTP exchange) - mechanisms that require direct
#:   enzyme-substrate contact to occur, and (unlike Reconstituted Complex
#:   below) BioGRID's definition does not describe fishing a substrate out
#:   of a cell lysate/extract.
STRONG_DIRECT_EXPERIMENTAL_SYSTEMS: frozenset[str] = frozenset(
    {
        "Co-crystal Structure",
        "Two-hybrid",
        "Far Western",
        "Protein-peptide",
        "PCA",
        "Thermal Shift Assay",
        "Surface Display",
        "Biochemical Activity",
    }
)

#: Methods whose BioGRID definition demonstrates co-complex membership,
#: spatial/photophysical proximity, or an indirect pulldown rather than a
#: structurally guaranteed direct binary contact. Real physical-interaction
#: evidence, but not accepted into the current strong/direct training set.
#:
#: - Affinity Capture-Luminescence / -MS / -Western: a bait is pulled down
#:   from a cell extract and any co-captured protein is reported as an
#:   interactor; the readout does not distinguish direct binding from
#:   membership in the same complex (a bridging third protein in the
#:   extract cannot be ruled out).
#: - Co-fractionation / Co-purification: proteins are observed together
#:   after one or more purification/fractionation steps - complex-level,
#:   not binary, evidence.
#: - Co-localization: spatial co-occurrence in the cell, which does not by
#:   itself demonstrate physical contact (BioGRID itself requires an
#:   additional mislocalization-on-deletion criterion just to accept this
#:   as evidence at all, underscoring how weak plain co-occurrence is).
#: - Cross-Linking-MS (XL-MS): BioGRID's text explicitly allows this to be
#:   "carried out with live cells or cell lysates in which ALL proteins
#:   are expressed at endogenous levels"; the cross-link shows two
#:   residues were within the crosslinker's reach in a complex cellular
#:   mixture, which can include residues on different subunits of a larger
#:   assembly - proximity, not confirmed pairwise direct contact.
#: - Proximity Label-MS: a bait-enzyme fusion labels any protein within
#:   diffusion range (e.g. BioID), which explicitly includes indirect
#:   neighbors, not just direct binding partners.
#: - FRET (reclassified in policy v2, was strong/direct in v1): BioGRID's
#:   own wording is "close proximity ... detected by fluorescence
#:   resonance energy transfer ... in vivo", not "direct contact" or
#:   "directly demonstrated" (contrast with Co-crystal Structure's
#:   wording). FRET/BRET report photophysical proximity within roughly the
#:   Forster radius (~1-10 nm) between fluorophore attachment points -
#:   comparable to or larger than a folded protein domain - which two
#:   proteins can satisfy while both docked on a shared scaffold/complex
#:   in a crowded live-cell compartment without directly touching each
#:   other; this is a well-documented FRET confound (proximity-induced
#:   signal from co-recruitment or overexpression). Unlike Two-hybrid/PCA,
#:   a positive FRET signal does not require the two full proteins to
#:   physically join for the assay's readout to occur - only the two
#:   fluorophores need to be within range - so FRET lacks the structural
#:   "must touch to work" guarantee that justifies the strong bucket, and
#:   BioGRID exposes no separate field (FRET efficiency, method variant,
#:   controls used) to distinguish rigorous from marginal reports. Fails
#:   closed for the current strict training set.
#: - Reconstituted Complex (reclassified in policy v2, was strong/direct
#:   in v1): BioGRID's own definition explicitly includes, as a named
#:   example, "GST pull-down assays where a GST-tagged protein is first
#:   isolated and then used to fish interactors from CELL LYSATES" as a
#:   Reconstituted Complex. That specific case has exactly the same
#:   bridging risk as Affinity Capture-MS/-Western (a purified bait
#:   exposed to an undefined mixture of cellular proteins, any of which
#:   could mediate an indirect association) - it differs only in whether
#:   the bait was purified before or after being expressed in a cell, not
#:   in whether other proteins were present during the binding step. The
#:   same Experimental System value also covers genuinely clean
#:   two-purified-component methods (gel-shift, SPR, BLI). BioGRID
#:   provides no structured field (Modification/Qualifications/Tags/
#:   Throughput/Score) that reliably separates the lysate-pulldown case
#:   from the two-purified-component case across the TAB3 export, and
#:   inferring it from free-text Qualifications would be exactly the
#:   fuzzy text-matching this policy is required to avoid. Since the
#:   category demonstrably includes indirect-capable executions with no
#:   reliable machine-readable discriminator, it fails closed for the
#:   current strict training set.
WEAK_PHYSICAL_EXPERIMENTAL_SYSTEMS: frozenset[str] = frozenset(
    {
        "Affinity Capture-Luminescence",
        "Affinity Capture-MS",
        "Affinity Capture-Western",
        "Co-fractionation",
        "Co-localization",
        "Co-purification",
        "Cross-Linking-MS (XL-MS)",
        "Proximity Label-MS",
        "FRET",
        "Reconstituted Complex",
    }
)

#: Every BioGRID physical Experimental System value reviewed by this
#: policy, regardless of bucket. A value outside this set is unreviewed and
#: therefore unknown to the policy (fails closed).
REVIEWED_PHYSICAL_EXPERIMENTAL_SYSTEMS: frozenset[str] = (
    RNA_EXPERIMENTAL_SYSTEMS
    | STRONG_DIRECT_EXPERIMENTAL_SYSTEMS
    | WEAK_PHYSICAL_EXPERIMENTAL_SYSTEMS
)


def is_strong_direct_evidence(experimental_system: str | None) -> bool:
    """Return whether one evidence row qualifies as strong/direct PPI evidence.

    Categorical decision only. An unreviewed/unknown Experimental System
    value fails closed (returns False) rather than being guessed at.
    """
    if experimental_system is None:
        return False

    return experimental_system in STRONG_DIRECT_EXPERIMENTAL_SYSTEMS


def is_rna_evidence(experimental_system: str | None) -> bool:
    """Return whether one evidence row is protein-RNA, not protein-protein."""
    if experimental_system is None:
        return False

    return experimental_system in RNA_EXPERIMENTAL_SYSTEMS


def unknown_experimental_systems(
    experimental_systems: Iterable[str],
) -> set[str]:
    """Return distinct values absent from the reviewed policy table.

    Used to surface newly introduced BioGRID Experimental System values for
    manual review rather than silently admitting or discarding them.
    """
    return {
        str(value)
        for value in experimental_systems
        if value is not None
        and str(value) not in REVIEWED_PHYSICAL_EXPERIMENTAL_SYSTEMS
    }
