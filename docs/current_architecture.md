# ProteinGraph — Complete Architecture Map

Analysis-only, based on reading every `.py` file under `protein_graph/` at the current checkpoint (`1bb4ecd`), cross-checked against `tests/`, `scripts/`, and `data/`. No files were modified.

---

## protein_graph/entities/

**Role of the package:** pure data — immutable dataclasses describing biological objects (a region of a sequence, a domain, a motif occurrence, a motif's abstract definition). No I/O, no external dependencies except each other.

### protein_graph/entities/region.py

**Purpose of the file:** defines the single coordinate/sequence contract every "piece of a protein" in this codebase must satisfy — 1-based, inclusive biological coordinates plus the literal subsequence.

**Classes:**
- `SequenceRegion`
  - Purpose: base dataclass for any bounded region inside a protein sequence.
  - Important fields/state: `protein_id`, `region_type`, `start`, `end`, `sequence`, `description`, `source`, `metadata`.
  - Important methods: `length` (property), `relative_position(protein_length)`, `position_category(protein_length)` (n_terminal/middle/c_terminal), `to_dict()`.
  - Input: raw coordinate/sequence values at construction; validates `start≥1`, `end≥start`, and that `len(sequence)` matches the coordinate span.
  - Output: a validated, frozen region object; `to_dict()` for serialization.
  - Used by: `DomainRecord`, `MotifRecord` (subclasses); `RegionFeaturePipeline`, `DomainFeaturePipeline`.

**Place in the pipeline:** foundational entity — everything downstream that talks about "a piece of a protein" (domain or motif) is a `SequenceRegion`.

**Current status:** PRODUCTION PATH (base class actively used by both production and dormant subsystems).

---

### protein_graph/entities/domain.py

**Purpose of the file:** the canonical representation of a protein domain (Pfam/InterPro/SMART/CDD-identified region) with biological annotations.

**Classes:**
- `DomainRecord(SequenceRegion)`
  - Purpose: represents one domain occurrence in one protein, with identifiers (`pfam_id`, `interpro_id`, `smart_id`, `cdd_id`), classification (`family`/`superfamily`/`clan`), function (`functional_roles`, `binding_types`, `recognized_motifs`), structure (`structural_class/architecture/fold`), and known ligand/ion/cofactor annotations.
  - Important fields/state: all of the above, all optional except that at least one domain identifier must be present.
  - Important methods: `__post_init__` (enforces `region_type == "domain"` and requires ≥1 identifier), `to_dict()`.
  - Input: constructed by `DomainBuilder` from UniProt domain-annotation dicts.
  - Output: immutable `DomainRecord`; consumed by `features/domains.py` extractors and `features/pair/domain_motif.py`.
  - Used by: `builders/domain_builder.py`, `pipelines/domain_features.py`, `features/domains.py`, `features/pair/domain_motif.py`.

**Dependencies:** `entities.region.SequenceRegion`.

**Place in the pipeline:**

    UniProtDataSource (domain_annotations)
        ↓
    DomainBuilder
        ↓
    DomainRecord
        ↓
    DomainFeaturePipeline / DomainMotifCompatibilityExtractor

**Current status:** SUPPORTING COMPONENT — fully implemented and tested, but the pipeline that consumes it (`DomainFeaturePipeline`) is not wired into `TrainingDatasetPipeline` (see pipelines section).

---

### protein_graph/entities/motif.py

**Purpose of the file:** the canonical representation of one concrete motif *occurrence* found inside a specific protein sequence (as opposed to the abstract motif *definition* — see `motif_definition.py`).

**Classes:**
- `MotifRecord(SequenceRegion)`
  - Purpose: one matched motif instance, carrying the matched text, which motif class/definition produced it, which domains are known to recognize it (`recognized_by`), and evidence metadata (`evidence_source`, `evidence_level`, `confidence`).
  - Important methods: `__post_init__` (enforces `region_type == "motif"`, requires `motif_id`, validates `matched_sequence == sequence`, validates `confidence ∈ [0,1]`), `to_dict()`.
  - Input: built by `MotifBuilder` when regex-scanning a full protein sequence against a set of `MotifDefinition`s.
  - Output: consumed by `features/elm_motifs.py`, `features/motif_pipeline.py`, `features/pair/domain_motif.py`.
  - Used by: `builders/motif_builder.py`, `features/motif_pipeline.py`, `features/pair/domain_motif.py`.

**Dependencies:** `entities.region.SequenceRegion`.

**Current status:** SUPPORTING COMPONENT — tested and functional, but the motif-occurrence extraction path (`MotifBuilder` → `MotifRecord` → `MotifFeaturePipeline`) is not called by `TrainingDatasetPipeline`.

---

### protein_graph/entities/motif_definition.py — **Special case #1 (of two `MotifDefinition` classes)**

**Purpose of the file:** the general, source-agnostic description of a motif *class* (e.g. "SH3-binding motif", regex `P..P`) — knowledge about a motif independent of any specific occurrence.

**Classes:**
- `MotifDefinition`
  - Purpose: validated immutable record of one motif class: `motif_id`, `motif_class`, `pattern` (regex, compiled once to check validity), optional `recognized_by` (tuple of recognizing domain identifiers — **may legitimately be empty** when no recognizing domain is curated yet), `interaction_role`, `source`, `evidence_level`.
  - Important methods: `__post_init__` normalizes/validates all fields, including compiling `pattern` to catch malformed regex early.
  - Input: built by `MotifDefinitionBuilder` from normalized ELM (or other future source) records.
  - Output: consumed by `MotifBuilder` (to scan whole proteins) and by the `matching_domain_motif` adapter (to become an executable scanner rule).
  - Used by: `builders/motif_builder.py`, `builders/motif_definition_builder.py`, `matching_domain_motif/motif_definition_adapter.py`.

**Current status:** SUPPORTING COMPONENT — real, tested, general-purpose knowledge representation; not currently fed into the production BioGRID training pipeline.

---

### protein_graph/entities/__init__.py
Re-exports `DomainRecord`, `MotifDefinition`, `MotifRecord`, `SequenceRegion`. No behavior of its own.

---

## protein_graph/data_sources/

**Role of the package:** all external-data ingestion. Every class here implements the same three-stage lifecycle (`_download` → `_parse` → `_normalize`) defined by `base.py`, so every data source behaves identically from the outside (`load()`).

### protein_graph/data_sources/base.py

**Classes:**
- `BaseDataSource` (ABC)
  - Purpose: enforces the download → parse → normalize lifecycle for every external data source in the project.
  - Important methods: abstract `_download`, `_parse`, `_normalize`; concrete `load(force_download=False)` orchestrates the three steps and creates `raw_dir`.
  - Used by: `BioGRIDDataSource`, `ELMDataSource`, `ELMInteractionDomainsDataSource`, `UniProtDataSource` (all inherit it).

**Current status:** PRODUCTION PATH — foundational contract for every data source actually used in the production build.

---

### protein_graph/data_sources/biogrid.py

**Purpose of the file:** faithful ingestion of the **global** BioGRID-ALL physical-interaction dataset (not the single-organism or multi-validated-only product — this distinction is explicit and deliberate, per the module docstring). This class's *only* job is downloading, parsing, and normalizing individual evidence rows — it does **not** decide evidence strength or resolve identifiers.

**Classes:**
- `BioGRIDDataSource(BaseDataSource)`
  - Purpose: downloads the `BIOGRID-ALL-<version>.tab3.zip` archive (either the "Latest-Release" alias or a pinned specific `release`), validates the ZIP, extracts the single global TAB3 member, filters to a configured `interaction_type` (default `"physical"`), optionally narrows by `taxonomy_id`/`sample_limit` for dev use, and normalizes columns into a flat evidence table (one row per BioGRID interaction record) including a source-specific undirected `biogrid_pair_id`.
  - Important fields/state: `interaction_type`, `taxonomy_id`, `species`, `sample_limit`, `release`, `release_version` (resolved concrete version, populated after `_parse`), `archive_path`, `download_url` (property).
  - Important methods: `_download` (streamed, atomic via `.part` temp file + ZIP validation), `_parse` (finds the single `.tab3.txt` member, extracts/validates the resolved release version, reads via pandas), `_normalize` (filters by interaction type/taxonomy, builds the flat evidence table, dedups on `biogrid_interaction_id`).
  - Input: none required except config; downloads from `downloads.thebiogrid.org`.
  - Output: `pd.DataFrame` of individual physical-evidence rows (columns like `organism_id_a/b`, `experimental_system`, `swissprot_accessions_a/b`, `biogrid_pair_id`, etc.) — **not yet** UniProt-mapped, **not yet** aggregated into pairs, **not yet** evidence-qualified.
  - Used by: `scripts/build_training_dataset.py` (real production entry point), `pipelines/biogrid_pair_acceptance.py` (consumes its output shape).

**Dependencies:** `data_sources.base.BaseDataSource`, `pandas`, `requests`.

**Current status:** PRODUCTION PATH.

---

### protein_graph/data_sources/biogrid_evidence_policy.py

**Purpose of the file:** the single reviewed source of truth distinguishing which BioGRID "physical" Experimental System values are strong/direct enough for the strict training set versus merely "physical but not direct enough" versus protein-RNA (must be excluded entirely). This is a *policy* module (constants + pure functions), not a class-based component — extensively documented with the biological reasoning for each bucket (e.g. why FRET and Reconstituted Complex were reclassified from strong to weak in policy v2).

**Important functions/constants:**
- `RNA_EXPERIMENTAL_SYSTEMS`, `STRONG_DIRECT_EXPERIMENTAL_SYSTEMS`, `WEAK_PHYSICAL_EXPERIMENTAL_SYSTEMS`, `REVIEWED_PHYSICAL_EXPERIMENTAL_SYSTEMS` (frozensets).
- `is_strong_direct_evidence(experimental_system)` — categorical, fails closed on unreviewed values.
- `is_rna_evidence(experimental_system)`.
- `unknown_experimental_systems(values)` — surfaces new/unreviewed BioGRID method names for manual review rather than silently guessing.
- Used by: `pipelines/biogrid_pair_acceptance.py`.

**Current status:** PRODUCTION PATH — this is the evidence-quality gatekeeper for the entire positive-reference dataset.

---

### protein_graph/data_sources/uniprot.py

**Purpose of the file:** downloads and normalizes UniProtKB records (batched, cached per-accession as JSON) into a flat canonical protein table.

**Classes:**
- `UniProtDataSource(BaseDataSource)`
  - Purpose: given a list of accessions, batch-fetch UniProtKB JSON via the REST search API, cache each record to its own JSON file, and normalize into a table with sequence, organism, domain annotations, signal-peptide/transmembrane ranges, localization terms, and cross-references (InterPro/Pfam/RefSeq).
  - Important fields/state: `accessions`, `unresolved_accessions` (populated for inactive/undownloadable entries — treated as expected churn, not an error).
  - Important methods: `_download` (batches of 100), `_download_batch`, `_parse` (reads cached JSON), `_normalize`/`_normalize_record` (one record → one flat dict; returns `None` for well-formed "Inactive" entries rather than raising), various `_extract_*` helpers for features/cross-refs.
  - Output: `pd.DataFrame` — one row per resolved protein, with `domain_annotations` (list of dicts consumed by `DomainBuilder`).
  - Used by: `protein_builder/builder.py`.

**Current status:** PRODUCTION PATH.

---

### protein_graph/data_sources/elm.py

**Purpose of the file:** loads ELM's motif-*class* definitions (regex patterns) from a TSV file (cached locally, or downloaded if a URL is configured — no default remote URL is hardcoded).

**Classes:**
- `ELMDataSource(BaseDataSource)`
  - Purpose: parse an ELM classes TSV (with column-alias tolerance for varying ELM export header names), normalize into canonical dicts with `motif_id`, `motif_class`, `pattern`, `recognized_by=()` (always empty at this stage — see `elm_interaction_domains.py` for how this gets filled), `source="ELM"`, `evidence_level`.
  - Output: `list[dict]` ready for `MotifDefinitionBuilder`.
  - Used by: `builders/motif_definition_builder.py` (via tests/integration), `scripts` (implicitly, cached data present under `data/raw/elm/`).

**Current status:** SUPPORTING COMPONENT (real, working, cached data exists on disk) — feeds the *domain/motif* subsystem, which is not wired into the BioGRID training pipeline.

---

### protein_graph/data_sources/elm_interaction_domains.py

**Purpose of the file:** loads ELM's curated "Interaction Domains" resource — for each ELM motif class, which Pfam/SMART/InterPro domain(s) recognize it. This is the *authoritative* source that fills in the `recognized_by` field that `ELMDataSource` always leaves empty.

**Important functions:**
- `build_recognized_by_map(records)` — groups normalized `(motif_id, domain_id)` records into `motif_id → tuple[domain_id, ...]`, preserving multiplicity (one motif can be recognized by several domain families).

**Classes:**
- `ELMInteractionDomainsDataSource(BaseDataSource)`
  - Purpose: parses the ELM interaction-domains TSV (default URL: `elm.eu.org/infos/browse_elm_interactiondomains.tsv`), detects each domain ID's namespace (Pfam/SMART/InterPro/CDD) from the identifier's own prefix pattern (never from name/description text), normalizes into `motif_id, domain_id, domain_namespace, domain_description, domain_name, source` records.
  - Used by: `builders/motif_definition_builder.py` (via `recognized_by_overrides`), tested end-to-end in `tests/integration/test_real_elm_domain_motif_matching.py` against real cached data.

**Current status:** SUPPORTING COMPONENT — real, verified against live data, exercised by a genuine integration test, but not reached by the BioGRID training pipeline.

---

### protein_graph/data_sources/__init__.py
Re-exports `BaseDataSource`, `BioGRIDDataSource`, `ELMDataSource`, `ELMInteractionDomainsDataSource`, `build_recognized_by_map`, `UniProtDataSource`.

---

## protein_graph/mapping/

### protein_graph/mapping/uniprot.py

**Purpose of the file:** resolves BioGRID's raw cross-reference columns (`swissprot_accessions_a/b`, `trembl_accessions_a/b`) to a single canonical UniProt accession per interactor, and builds the canonical undirected pair key.

**Classes:**
- `BioGRIDUniProtMapper`
  - Purpose: for each BioGRID evidence row, pick one UniProt accession per interactor — Swiss-Prot preferred over TrEMBL, first listed accession if multiple — and record the mapping source; flag rows as `mapping_complete` only if both sides resolved; build `uniprot_pair_id` as a lexicographically-ordered `"A--B"` string so `A-B` and `B-A` collapse to the same key.
  - Important methods: `map(interactions)` (vectorized via `Series.map`, not row-wise `apply`, for BioGRID-ALL scale), `coverage(mapped)` (mapping-quality statistics).
  - Input: BioGRID evidence `pd.DataFrame` (from `BioGRIDDataSource`).
  - Output: same table + `protein_a_id`, `protein_b_id`, `protein_a_mapping_source`, `protein_b_mapping_source`, `mapping_complete`, `uniprot_pair_id`.
  - Used by: `pipelines/biogrid_pair_acceptance.py` (`BioGRIDPairAcceptancePipeline.__init__` defaults to one).

**Place in the pipeline:**

    BioGRIDDataSource
        ↓
    BioGRIDUniProtMapper
        ↓
    BioGRIDPairAcceptancePipeline (evidence qualification + pair reduction)

**Current status:** PRODUCTION PATH.

---

## protein_graph/builders/

**Role of the package:** turn raw normalized dicts/records from data sources into validated entity objects (`DomainRecord`, `MotifRecord`, `MotifDefinition`).

### protein_graph/builders/domain_builder.py

**Classes:**
- `DomainBuilder`
  - Purpose: converts a protein's `domain_annotations` (list of dicts, produced by `UniProtDataSource`) into `DomainRecord` objects, slicing the domain subsequence out of the full protein sequence.
  - Important methods: `build_many(protein)`, `_build_domain(...)`.
  - Used by: `pipelines/domain_features.py`.

**Current status:** SUPPORTING COMPONENT (used only by the non-wired `DomainFeaturePipeline`).

### protein_graph/builders/motif_builder.py

**Classes:**
- `MotifBuilder`
  - Purpose: given a set of `MotifDefinition`s, regex-scans one protein's full sequence and produces one `MotifRecord` per match (supports overlapping matches via lookahead regex).
  - Important methods: `build(protein)` / `build_many(...)` (alias), `_find_matches`, `_compile_pattern`.
  - Used by: `tests/integration/test_elm_motif_annotation.py` (real usage demonstrated end-to-end); not called by any pipeline.

**Current status:** SUPPORTING COMPONENT / CURRENTLY DORMANT relative to production — works and is proven correct by an integration test, but no pipeline invokes it.

### protein_graph/builders/motif_definition_builder.py

**Classes:**
- `MotifDefinitionBuilder`
  - Purpose: converts a canonical motif record dict (already normalized by a data source) into a validated `entities.motif_definition.MotifDefinition`, optionally filling `recognized_by` from a `recognized_by_overrides` map (i.e. from `build_recognized_by_map`) only when the record itself didn't already specify one.
  - Important methods: `build(record, recognized_by_overrides=None)`, `build_many(records, ...)`.
  - Used by: `tests/integration/test_real_elm_domain_motif_matching.py`, `tests/integration/test_elm_motif_annotation.py`.

**Current status:** SUPPORTING COMPONENT — this is the exact place where ELM's regex definitions and its authoritative domain↔motif knowledge are joined into one object; real and tested but not called by any pipeline.

---

## protein_graph/protein_builder/

### protein_graph/protein_builder/builder.py

**Classes:**
- `ProteinBuilder`
  - Purpose: thin orchestration wrapper around `UniProtDataSource` — builds one or many canonical protein records (as a validated `pd.DataFrame`) from a list of accessions.
  - Important methods: `build(accession, ...)`, `build_many(accessions, ...)`, `_validate_table` (requires `protein_id`, `sequence`, `protein_length`, `organism`, `taxon_id`; no nulls/dupes in `protein_id`, no null sequences).
  - Used by: `pipelines/protein_features.py` (`ProteinFeaturePipeline` defaults to one).

**Place in the pipeline:**

    UniProtDataSource
        ↓
    ProteinBuilder
        ↓
    ProteinFeaturePipeline

**Current status:** PRODUCTION PATH.

---

## protein_graph/domain_detector/

### protein_graph/domain_detector/static_detector.py

**Classes:**
- `StaticDomainDetector`
  - Purpose: a fixed-lookup-table implementation of the `DomainDetector` protocol from `matching_domain_motif.pair_analyzer` — given a pre-populated `{sequence: [DetectedDomain, ...]}` mapping, returns the domains for an exact-matching sequence. This is a **test/dev fixture**, not a real domain predictor (no Pfam/HMMER/InterPro scanning is implemented anywhere in the repo).
  - Used by: `tests/matching_domain_motif/test_pair_analyzer.py`, `tests/integration/test_real_elm_domain_motif_matching.py`.

**Dependencies:** `matching_domain_motif.domain_motif.DetectedDomain`.

**Current status:** TEST/DEVELOPMENT UTILITY — this is explicitly a stand-in until a real domain detector (Pfam/HMMER/InterPro) is implemented; there is no production `DomainDetector` implementation anywhere yet.

`protein_graph/domain_detector/__init__.py` is **empty** (no re-exports).

---

## protein_graph/matching_domain_motif/ — **Special case #2: the second `MotifDefinition`, and the domain↔motif matching layer**

**Role of the package:** an executable "given a domain found in protein A, scan protein B for the motifs that domain is known to recognize (and vice versa)" engine. This is the exact-regex matching layer, distinct from the *feature computation* layer in `features/motif_properties_of_certain_domain/` and from the *label-overlap* layer in `features/pair/domain_motif.py`.

Note: this module's source comments/docstrings are written partly in **Russian** (Cyrillic) — a stylistic inconsistency with the rest of the codebase (all English), worth knowing if editing it, but doesn't affect behavior.

### protein_graph/matching_domain_motif/domain_motif.py

**Classes:**
- `DetectedDomain` — a domain detected in a sequence (`domain_class`, `start`, `end`, `confidence`, `domain_id`); 1-based inclusive coordinates.
- `MotifDefinition` — **the second `MotifDefinition` class in the repo** (see comparison below).
- `DomainMotifCandidate` — one candidate cross-protein domain→motif match: source/target protein IDs, the domain that triggered the scan, the target region matched, `motif_score`, `combined_score` (`domain.confidence * motif_score`), and `matched_known_regex`.
- `MotifScanner` (Protocol) — interface for regex- or ML-based scanners.
- `DomainMotifRegistry`
  - Purpose: builds a `domain_class → tuple[MotifDefinition, ...]` index from a list of scanner-rule `MotifDefinition`s (each motif can be recognized by multiple domain classes, so it's indexed under each).
  - Important methods: `motifs_for_domain(domain_class)`, `supports_domain(domain_class)`.
- `KnownMotifScanner`
  - Purpose: the concrete `MotifScanner` implementation — given a detected domain and a target sequence, looks up only the motifs that domain recognizes (via the registry) and regex-scans the target sequence for them (overlapping matches supported via lookahead), converting Python 0-based match indices to the project's 1-based biological coordinate convention.

**Current status:** SUPPORTING COMPONENT — fully implemented, exact-match domain→motif scanning; proven correct against real ELM data in `tests/integration/test_real_elm_domain_motif_matching.py`, but never invoked from a pipeline.

---

#### Special case: `entities.motif_definition.MotifDefinition` vs `matching_domain_motif.domain_motif.MotifDefinition`

**Verdict: related concepts with different responsibilities, not duplicates** — and the codebase itself documents this explicitly in `matching_domain_motif/motif_definition_adapter.py`'s module docstring.

| | `entities.motif_definition.MotifDefinition` | `matching_domain_motif.domain_motif.MotifDefinition` |
|---|---|---|
| Represents | General knowledge about a motif class | An executable scanner rule |
| `recognized_by` | Optional — may be empty; the motif is still valid knowledge | **Required non-empty** — a scanner rule with no triggering domain isn't executable, so construction raises if empty |
| Field name for the regex | `pattern` | `regex` (same underlying string, renamed to signal "matching-ready") |
| Extra fields | `interaction_role`, `source`, `evidence_level` (provenance) | none — stripped down to only what `KnownMotifScanner` needs |
| Used by | `MotifBuilder` (whole-protein annotation), `MotifDefinitionBuilder` | `DomainMotifRegistry`, `KnownMotifScanner` |

The **bridge** between them is `to_scanner_rule` / `to_scanner_rules` in `motif_definition_adapter.py`: converts a general definition to a scanner rule only when `recognized_by` is already non-empty; returns `None` (dropped by `to_scanner_rules`) otherwise. No domain is ever fabricated to force a conversion.

---

### protein_graph/matching_domain_motif/motif_definition_adapter.py

**Important functions:**
- `to_scanner_rule(definition)` → `ScannerMotifDefinition | None`
- `to_scanner_rules(definitions)` → `list[ScannerMotifDefinition]`
- Used by: `tests/integration/test_real_elm_domain_motif_matching.py`, `tests/matching_domain_motif/test_motif_definition_adapter.py`.

**Current status:** SUPPORTING COMPONENT.

### protein_graph/matching_domain_motif/pair_analyzer.py

**Classes:**
- `DomainDetector` (Protocol) — interface any real domain detector (Pfam/HMMER/InterPro/custom) must satisfy; only implementation today is the test-only `StaticDomainDetector`.
- `ProteinInput` — simple `(protein_id, sequence)` validated pair.
- `PairAnalysisResult` — bidirectional analysis output: domains found in A and B, candidates found in each direction, skipped (low-confidence or unsupported) domains, plus convenience properties `all_candidates`, `num_candidates`, `has_bidirectional_candidates`, `max_combined_score`.
- `ProteinPairAnalyzer`
  - Purpose: orchestrates the full domain↔motif matching for one protein pair — detects domains in both proteins, then for each domain above `domain_confidence_threshold` that the registry supports, scans the *other* protein for the domain's recognized motifs, in both directions.
  - Used by: `tests/matching_domain_motif/test_pair_analyzer.py`, `tests/integration/test_real_elm_domain_motif_matching.py`.

**Place in the pipeline (currently isolated, not connected to production):**

    ELMDataSource + ELMInteractionDomainsDataSource
        ↓
    MotifDefinitionBuilder (recognized_by_overrides = build_recognized_by_map(...))
        ↓
    to_scanner_rules  (entities.MotifDefinition → matching_domain_motif.MotifDefinition)
        ↓
    DomainMotifRegistry
        ↓
    KnownMotifScanner  +  a DomainDetector implementation (test-only stand-in for a real detector)
        ↓
    ProteinPairAnalyzer.analyze(protein_a, protein_b) → PairAnalysisResult

This chain is proven correct end-to-end against **real cached ELM data** in `tests/integration/test_real_elm_domain_motif_matching.py`, but:
- No production `DomainDetector` exists — only the test-only `StaticDomainDetector` lookup table.
- Nothing in `pipelines/` calls `ProteinPairAnalyzer` or consumes `PairAnalysisResult`.
- The separate `MotifBuilder`→`MotifRecord` path (whole-protein motif annotation, not domain-triggered) and `features/motif_properties_of_certain_domain/` (motif similarity scoring) are also fully working but likewise unconnected to each other and to the pipeline above.

**Current status:** SUPPORTING COMPONENT, proven correct end-to-end against real ELM data, but **CURRENTLY DORMANT / NOT INTEGRATED** into any pipeline. This matches the checkpoint's own stated boundary: exact domain-motif matching is validated; motif-instance/fuzzy candidate-generation work has not been started, and none of this subsystem feeds the BioGRID training pipeline yet.

---

## protein_graph/features/

**Role of the package:** feature extractors, at four different granularities: whole-protein (`features/sequence/`), single-domain (`domains.py`), single-motif-occurrence (`elm_motifs.py`), and motif-vs-known-motif comparison (`motif_properties_of_certain_domain/`). Only the whole-protein sequence extractors are wired into the production pipeline.

### protein_graph/features/base.py

**Classes:**
- `FeatureExtractor` (ABC) — one abstract method `extract(protein) -> dict`. `name` class attribute.
  - Used by: every extractor in `features/sequence/`, plus `SequencePropertyExtractor`, `ELMMotifExtractor`.

**Current status:** PRODUCTION PATH (base contract for the extractors actually used).

### protein_graph/features/registry.py & protein_graph/features/pipeline.py

- `registry.py` defines `FEATURE_EXTRACTORS: list[type[FeatureExtractor]] = []` — declared but **never populated anywhere** in the codebase.
- `pipeline.py` defines `FeaturePipeline` (generic `add()`/`extract()` sequential runner) — **never imported outside itself**.

**Current status:** CURRENTLY DORMANT / NOT INTEGRATED — generic scaffolding for a plugin-style extractor registry that was never adopted; the real production pipelines (`ProteinFeaturePipeline`, etc.) each hardcode their own extractor list instead of using this registry/pipeline.

### protein_graph/features/sequence/ (composition.py, fractions.py, amino_acid_density.py, runs.py, periodicity.py, constants.py)

**Purpose of the package:** whole-protein sequence-derived numeric features — the actual feature set used by the production model.

**Classes** (each a `FeatureExtractor`):
- `AminoAcidCompositionExtractor` (composition.py) — per-residue fraction of each of the 20 amino acids + `protein_length`.
- `ResidueFractionExtractor` (fractions.py) — fraction of residues in each physicochemical group (hydrophobic/hydrophilic/charged/positive/negative/polar/aromatic).
- `LocalDensityExtractor` (amino_acid_density.py) — count of non-overlapping windows (sizes 5/10/20/30) exceeding a 0.75 density threshold for hydrophobic/charged/positive/negative groups.
- `ConsecutiveRunExtractor` (runs.py) — max run length + binned run-length indicators per group.
- `PeriodicityExtractor` (periodicity.py) — periodic-pattern run detection (steps 2/3/4/7) for hydrophobic/charged residues.
- `constants.py` — shared amino-acid group definitions (`AMINO_ACIDS`, `HYDROPHOBIC`, `GLOBAL_GROUPS`, `LOCAL_GROUPS`, `WINDOW_SIZES`, `RUN_BINS`, `PERIODS`); not a class, just shared constants.

**Used by:** `features/sequence_properties.py` (`SequencePropertyExtractor`, a convenience bundler), `pipelines/protein_features.py` (`ProteinFeaturePipeline`, directly instantiates the same 5 extractors — **not** via `SequencePropertyExtractor`), `pipelines/region_features.py`.

**Current status:** PRODUCTION PATH — these are the only per-protein feature extractors that actually feed `PairFeaturePipeline.MODEL_FEATURE_COLUMNS`.

### protein_graph/features/sequence_properties.py

**Classes:**
- `SequencePropertyExtractor(FeatureExtractor)` — bundles all 5 sequence extractors above into one `extract()` call.
  - Used by: nothing in production code (`ProteinFeaturePipeline` reimplements the same extractor list independently rather than using this bundler) — only reachable via its own tests, if any exist (none found under `tests/features/` specifically for this file).

**Current status:** SUPPORTING COMPONENT / near-duplicate of the extractor wiring inside `ProteinFeaturePipeline` — not itself a duplicate *implementation* (it delegates to the same underlying extractor classes), but it is a second, unused place where "run all sequence extractors" is assembled.

### protein_graph/features/domains.py

**Classes:**
- `DomainFeatureContext` — carries `protein_length` for domain feature extraction.
- `DomainGeometryFeatureExtractor` — position/length/fraction-of-protein features for one `DomainRecord`.
- `DomainFunctionFeatureExtractor` — counts/flags for functional roles, binding types, recognized motifs, ligands/ions/cofactors.
- `DomainAnnotationFeatureExtractor` — presence flags/counts for identifiers (Pfam/InterPro/SMART/CDD) and classification (family/superfamily/clan).
- `DomainStructureFeatureExtractor` — presence flags/counts for structural class/architecture/fold annotations.

**Used by:** `pipelines/domain_features.py` (`DomainFeaturePipeline`).

**Current status:** SUPPORTING COMPONENT — real, tested, but `DomainFeaturePipeline` is not called from `TrainingDatasetPipeline`.

### protein_graph/features/elm_motifs.py

**Classes:**
- `MotifFeatureContext` — carries `protein_length`.
- `MotifGeometryFeatureExtractor` — position/length/fraction features for one `MotifRecord`.
- `MotifAnnotationFeatureExtractor` — presence flags/counts for motif annotation fields.

**Used by:** `features/motif_pipeline.py` (`MotifFeaturePipeline`).

**Current status:** SUPPORTING COMPONENT — not reached by any pipeline in `pipelines/`.

### protein_graph/features/elm_profile.py

**Classes:**
- `ELMMotifExtractor(FeatureExtractor)` — a *different, coarser* way of getting ELM motif information onto a protein: instead of building `MotifRecord`s via `MotifBuilder`+`MotifDefinition`, it scans a protein sequence directly against a raw ELM motif `pd.DataFrame` (columns `elm_id`, `regex`, optional `motif_class`) and returns an aggregated profile (total hits, unique classes, matched IDs/classes, per-class counts) rather than per-occurrence records.
  - Not related by inheritance/duplication to `MotifBuilder`/`MotifRecord` — it's a self-contained alternate implementation of "count ELM motif hits in a protein," used nowhere else and not consumed by `MotifBuilder`.
  - Used by: only its own test (`tests/features/test_elm_profile.py`).

**Current status:** CURRENTLY DORMANT / NOT INTEGRATED — an independent, working, but unconnected protein-level ELM summarizer; overlaps in *purpose* (not implementation) with the `MotifBuilder`→`MotifRecord`→`MotifFeaturePipeline` chain, since both ultimately answer "what ELM motifs are in this protein," but neither calls the other.

### protein_graph/features/motif_pipeline.py

**Classes:**
- `MotifFeatureExtractor` (Protocol) — interface for motif-level extractors (`extract(motif, context)`).
- `MotifFeaturePipeline` — combines `MotifGeometryFeatureExtractor` + `MotifAnnotationFeatureExtractor` (default) into one feature table, one row per `MotifRecord`, with collision detection between extractors.
  - Used by: only its own test (`tests/features/test_motif_pipeline.py`).

**Current status:** CURRENTLY DORMANT / NOT INTEGRATED.

### protein_graph/features/interface.py, localization.py, embeddings.py, evolution.py, structure.py, surface.py

All six files are **completely empty** — no classes, no functions, just placeholder modules (presumably named for planned future feature categories: protein-protein interface features, subcellular localization features, learned embeddings, evolutionary/conservation features, structural features, surface features).

**Current status:** CURRENTLY DORMANT / NOT INTEGRATED (unimplemented placeholders).

---

### protein_graph/features/motif_properties_of_certain_domain/ — **Special case: `WindowInterval` and the motif-comparison feature layer**

**Role of the package:** compute rich numeric features describing a candidate motif region (composition, similarity to a known reference motif, positional context) — a *feature-engineering* layer, distinct from the domain-motif *matching* layer in `matching_domain_motif/`.

#### protein_graph/features/motif_properties_of_certain_domain/window.py — `WindowInterval`

**Classes:**
- `WindowInterval`
  - Purpose: a **0-based, half-open** Python-slice-style interval — explicitly *not* a biological entity: no sequence, no identity, no protein_id. Exists purely to express local windows and overlap tests for this package's context/combined extractors.
  - Important methods: `overlaps(other)`, `from_biological(start, end)` classmethod — explicit converter from `SequenceRegion`'s 1-based inclusive convention.
  - Used by: `context_features.py`, `combined_features.py`.

**Special case verdict — `WindowInterval` vs `SequenceRegion`/`DomainRecord`/`MotifRecord`:** **not duplicates, different responsibilities by design** — the module's own docstring states this explicitly. `SequenceRegion` (and subclasses) carry biological identity + 1-based coordinates + sequence text; `WindowInterval` is a bare, identity-free, 0-based Python-index helper used only for local overlap math inside this one feature package. `from_biological` is the documented bridge between the two coordinate systems.

**Current status:** SUPPORTING COMPONENT (used only within this feature package, which is itself not wired into production).

#### protein_graph/features/motif_properties_of_certain_domain/sequence_features.py

**Classes:**
- `MotifSequenceFeatureExtractor`
  - Purpose: computes intrinsic composition/entropy/run/spacing/terminal/positional-moment features for a *single* short motif sequence (not a whole protein) — entropy, dominant-AA fraction, physicochemical group fractions, max runs per group, proline-spacing patterns (`P-X-P`, `P-X-X-P`), N/C-terminal one-hot identity, positional moments per group.
  - Used by: `context_features.py`, `combined_features.py`, `domain_motif_profile.py`.

**Current status:** SUPPORTING COMPONENT.

#### protein_graph/features/motif_properties_of_certain_domain/context_features.py

**Classes:**
- `MotifContextConfig` — `flank_size`, `terminal_region_fraction`.
- `MotifContextFeatureExtractor` — given a protein sequence and a candidate motif span, computes position/flank/local-context features, including overlap with supplied domain/transmembrane `WindowInterval`s, and delegates composition comparisons to `MotifSequenceFeatureExtractor`.

**Current status:** SUPPORTING COMPONENT.

#### protein_graph/features/motif_properties_of_certain_domain/similarity_features.py

**Classes:**
- `MotifSimilarityConfig` — `critical_positions`, `kmer_sizes`, `require_equal_length`, `use_blosum62`.
- `MotifSimilarityFeatureExtractor` — compares a candidate sequence against one reference motif: exact/physicochemical identity, k-mer Jaccard overlap, longest common substring/prefix/suffix, mismatch clustering, optional BLOSUM62 scoring (lazy Biopython import via `_load_blosum62`, raises a clear error if Biopython isn't installed).

**Current status:** SUPPORTING COMPONENT.

#### protein_graph/features/motif_properties_of_certain_domain/domain_motif_profile.py

**Classes:**
- `KnownMotif` — one experimentally-known motif sequence recognized by a domain (validated AA alphabet, optional coordinates/confidence).
- `KnownMotifFeatures` — a `KnownMotif` + its computed features.
- `FeatureDistribution` — mean/std/min/max of one feature across a domain's known motifs.
- `DomainMotifProfile` — the full set of known motifs (with features) for one domain class, plus per-feature distributions; convenience properties `motif_count`, `motif_sequences`, `motif_ids`, `motif_lengths`, `motifs_of_length(n)`.
- `DomainMotifProfileBuilder` — builds a `DomainMotifProfile` from `domain_class` + an iterable of `KnownMotif`, extracting numeric features via any `MotifSequenceFeatureExtractorProtocol`-conforming extractor and computing cross-motif distributions for features common to all known motifs.

**Current status:** SUPPORTING COMPONENT — this is the "characterize what a domain's known motifs typically look like" reference-building layer, used to later compare a *candidate* against a domain's known motif population (not yet wired to anything that generates candidates in production).

#### protein_graph/features/motif_properties_of_certain_domain/combined_features.py

**Classes:**
- `CombinedMotifFeatureExtractor` — the top-level assembler for this feature package: given a `reference_motif` string and a candidate region (`protein_sequence`, `candidate_start/end`), runs `MotifSequenceFeatureExtractor` on both reference and candidate, `MotifSimilarityFeatureExtractor` between them, and `MotifContextFeatureExtractor` on the candidate's context, merging everything under `reference_*`/`candidate_*`/`similarity_*`/`context_*` prefixes.

**Place in the pipeline (self-contained, not connected elsewhere):**

    reference motif + candidate region
        ↓
    MotifSequenceFeatureExtractor (×2)  +  MotifSimilarityFeatureExtractor  +  MotifContextFeatureExtractor
        ↓
    CombinedMotifFeatureExtractor → one flat feature dict

**Current status:** CURRENTLY DORMANT / NOT INTEGRATED as a whole package — internally complete and well-tested (7 test files under `tests/features/motif_properties_of_certain_domain/`), but nothing outside this package calls into it, and it isn't connected to `matching_domain_motif` (which finds *candidates*) — this package assumes a candidate region is already given.

---

### protein_graph/features/pair/domain_motif.py — **Special case: third domain↔motif approach**

**Purpose of the file:** pairwise (protein A × protein B) feature extraction based on *label overlap* between a protein's domains and its partner's motifs — a third, independent way the codebase expresses "domain recognizes motif," different from both `matching_domain_motif` (regex scanning by domain class) and `motif_properties_of_certain_domain` (similarity scoring against known motifs).

**Classes:**
- `DomainMotifMatch` — one compatible domain(protein X)–motif(protein Y) pairing, with `direction` (`"a_domain_to_b_motif"` / `"b_domain_to_a_motif"`).
- `DomainMotifCompatibilityExtractor`
  - Purpose: given `domains_a`, `motifs_a`, `domains_b`, `motifs_b` (already-built `DomainRecord`/`MotifRecord` collections for a pair), finds matches in both directions by **string-label overlap** — a domain's `family`/`description`/`pfam_id`/`interpro_id`/`smart_id`/`cdd_id` values, case-folded, checked against a motif's `recognized_by` set — and returns aggregate counts/flags plus the individual match tuples.
  - Input: pre-built domain/motif records for both proteins (no scanning of raw sequences — this operates purely on labels already attached to `DomainRecord`/`MotifRecord`).
  - Used by: only its own test (`tests/features/test_domain_motif.py`).

**Comparison of the three domain↔motif "matching" components:**

| Component | Mechanism | Input | Output |
|---|---|---|---|
| `matching_domain_motif.domain_motif.KnownMotifScanner` | Regex-scans the **raw target sequence** for motifs the detected domain recognizes | `DetectedDomain` + target sequence | `DomainMotifCandidate` (with matched substring, coordinates) |
| `features.pair.domain_motif.DomainMotifCompatibilityExtractor` | **String/label overlap** between already-annotated `DomainRecord.family/pfam_id/...` and `MotifRecord.recognized_by` | Pre-built `DomainRecord`/`MotifRecord` collections | `DomainMotifMatch` (labels only, no new sequence scanning) |
| `features.motif_properties_of_certain_domain.*` | **Similarity scoring** of a given candidate region against a domain's known reference motifs | reference motif string + candidate region | numeric feature dict (no matching decision at all) |

These are **not duplicates** — they operate at different stages (candidate generation vs. label cross-check vs. feature scoring) — but they are three parallel, unconnected implementations of related domain↔motif reasoning, none of which currently feed each other or a shared pipeline.

**Current status:** CURRENTLY DORMANT / NOT INTEGRATED.

### protein_graph/features/__init__.py
**Empty** — no re-exports at all (unusual for an `__init__.py` in this codebase; every other package re-exports its public classes, this one doesn't).

---

## protein_graph/knowledge/

`go.py`, `kegg.py`, `reactome.py`, `string.py`, `syngo.py`, and `__init__.py` are **all completely empty**. These are placeholder modules for planned external knowledge-base integrations (Gene Ontology, KEGG, Reactome, STRING, SynGO) — no code exists yet.

**Current status:** CURRENTLY DORMANT / NOT INTEGRATED (unimplemented placeholders).

---

## protein_graph/pipelines/

**Role of the package:** orchestration layer — combines entities/builders/features/data_sources into the actual end-to-end dataset build. This is where "production path" vs. "supporting/dormant" gets decided for everything upstream.

### protein_graph/pipelines/biogrid_pair_acceptance.py

**Classes:**
- `BioGRIDPairAcceptancePipeline`
  - Purpose: reduces raw BioGRID evidence rows into accepted strong/direct **protein pairs** (not evidence rows). Order-sensitive pipeline: (1) reject protein-RNA rows, (2) UniProt-map via `BioGRIDUniProtMapper`, (3) keep only fully-mapped rows, (4) per-row strong/direct qualification check (never pre-aggregated — a pair with 4 weak rows and 1 strong row must still be accepted), (5) group by `uniprot_pair_id`, keep pairs with ≥1 qualifying row, (6) emit exactly one representative row per accepted pair (via `drop_duplicates`, deliberately not `groupby().first()`, to avoid stitching provenance from two different physical rows).
  - Important methods: `build(evidence) -> (accepted_pairs_df, report_dict)`, `_accept_pairs`, `_validate_inputs`.
  - Input: `BioGRIDDataSource`'s normalized evidence table.
  - Output: one row per accepted pair (`protein_a_id`, `protein_b_id`, `uniprot_pair_id`, representative provenance columns) + a reproducibility report dict.
  - Used by: `pipelines/training_dataset.py` (`TrainingDatasetPipeline`).

**Dependencies:** `data_sources.biogrid_evidence_policy`, `mapping.BioGRIDUniProtMapper`.

**Current status:** PRODUCTION PATH — this is the evidence→pair reduction step that turns raw BioGRID rows into the positive-reference dataset.

### protein_graph/pipelines/protein_features.py

**Classes:**
- `ProteinFeaturePipeline`
  - Purpose: downloads proteins (via `ProteinBuilder`) and computes the 5 whole-protein sequence extractors (composition, fractions, density, runs, periodicity) for each, producing a validated feature table.
  - Important methods: `build(accessions, ...)`, `save(...)` (Parquet), `build_and_save(...)`, `_extract_protein_features`, `_validate_table` (requires `protein_id`, `protein_length`; no dupe IDs).
  - Used by: `pipelines/training_dataset.py`.

**Current status:** PRODUCTION PATH.

### protein_graph/pipelines/pair_features.py

**Classes:**
- `PairFeaturePipeline`
  - Purpose: joins two protein feature vectors (`a_*`/`b_*` prefixed) onto each accepted interaction pair (inner merge on `protein_a_id`/`protein_b_id`), and adds pair-level derived features (`protein_length_difference`, `hydrophobic_fraction_difference`, `charged_fraction_difference`, `positive_negative_balance`, `same_taxonomy`).
  - Important class attributes: `MODEL_FEATURE_COLUMNS` — an **explicit allowlist** frozenset of every column approved as ML model input (every `a_*`/`b_*` sequence-feature column, verified by actually executing the extractors against a sample sequence rather than hand-transcribed, plus the 5 pair-level derived columns). Any column not in this set — raw BioGRID/provenance/identity columns, or any future unreviewed column — is excluded by default.
  - Important methods: `build(interactions, protein_features)`, `select_model_features(pairs)` (applies the allowlist).
  - Used by: `pipelines/training_dataset.py`.

**Current status:** PRODUCTION PATH — this is the fail-closed boundary the checkpoint refers to as "fail-closed MODEL_FEATURE_COLUMNS boundary."

### protein_graph/pipelines/training_dataset.py

**Classes:**
- `TrainingDatasetPipeline`
  - Purpose: the top-level orchestrator — wires `BioGRIDPairAcceptancePipeline` → collect accessions → `ProteinFeaturePipeline` → drop pairs missing protein features → `PairFeaturePipeline` → (optionally) `build_model_matrix` (enforces schema completeness against `PairFeaturePipeline.MODEL_FEATURE_COLUMNS`, raises if any approved column is missing, returns columns in a fixed sorted order).
  - Important methods: `build(evidence, ...) -> (pair_features_df, report)`, `build_model_matrix(pair_features)`, `build_and_save(...)` (writes both the full reference table and the model-only matrix as separate Parquet files, plus a JSON report).
  - Used by: `scripts/build_training_dataset.py` (the real CLI entry point).

**Dependencies:** `BioGRIDPairAcceptancePipeline`, `PairFeaturePipeline`, `ProteinFeaturePipeline`.

**Current status:** PRODUCTION PATH — this is the top of the production call graph, confirmed by `scripts/build_training_dataset.py` and by real output artifacts on disk (`data/processed/biogrid_rat_pair_features.parquet` + its report).

### protein_graph/pipelines/domain_features.py

**Classes:**
- `DomainFeatureExtractor` (Protocol).
- `DomainFeaturePipeline` — builds `DomainBuilder` + all 4 domain feature extractors (`features/domains.py`) + `RegionFeaturePipeline` sequence features into one per-domain feature table.
  - Used by: only its own test (`tests/pipelines/test_domain_features.py`); **not** called from `TrainingDatasetPipeline`.

**Current status:** SUPPORTING COMPONENT / CURRENTLY DORMANT relative to the production training flow.

### protein_graph/pipelines/region_features.py

**Classes:**
- `RegionFeaturePipeline` — the generic version of the sequence-feature bundling: runs all 5 `features/sequence/` extractors against *any* `SequenceRegion` (not just whole proteins), dropping the `protein_length` key the composition extractor emits (since it was written for whole proteins).
  - Used by: `pipelines/domain_features.py`.

**Current status:** SUPPORTING COMPONENT (only consumer is the dormant `DomainFeaturePipeline`).

### protein_graph/pipelines/__init__.py
Re-exports all six pipeline classes.

---

# End-to-End Architecture Map

## 1. Production PPI reference/training pipeline (verified against real call graph and real output on disk)

```
BioGRIDDataSource                         (data_sources/biogrid.py)
    ↓  raw evidence rows (physical, global, all organisms)
BioGRIDPairAcceptancePipeline             (pipelines/biogrid_pair_acceptance.py)
    │
    ├─ reject protein-RNA rows          (data_sources/biogrid_evidence_policy.py: RNA_EXPERIMENTAL_SYSTEMS)
    ├─ BioGRIDUniProtMapper.map(...)    (mapping/uniprot.py) → protein_a_id/protein_b_id/uniprot_pair_id
    ├─ keep mapping_complete rows only
    ├─ is_strong_direct_evidence(...)   (data_sources/biogrid_evidence_policy.py)
    ├─ group by uniprot_pair_id, keep pairs with ≥1 qualifying row
    └─ emit one representative row per accepted pair
    ↓  accepted strong/direct positive pairs
TrainingDatasetPipeline.build(...)        (pipelines/training_dataset.py)
    │
    ├─ collect unique accessions from both pair sides
    ├─ ProteinFeaturePipeline.build(accessions)      (pipelines/protein_features.py)
    │       └─ ProteinBuilder → UniProtDataSource      (protein_builder/builder.py, data_sources/uniprot.py)
    │       └─ 5 sequence extractors                   (features/sequence/*)
    ├─ drop pairs missing protein features
    ├─ PairFeaturePipeline.build(pairs, protein_features)   (pipelines/pair_features.py)
    │       └─ merge a_*/b_* features + 5 pair-level derived features
    ↓  pair_features (full reference table, includes IDs/provenance/features)
TrainingDatasetPipeline.build_model_matrix(pair_features)
    └─ PairFeaturePipeline.select_model_features(...)  → filtered to MODEL_FEATURE_COLUMNS allowlist only
    ↓
model matrix  (ML-ready, IDs/provenance stripped out)
```

Driven in practice by `scripts/build_training_dataset.py`, and real artifacts exist on disk confirming this has actually run: `data/processed/biogrid_rat_pair_features.parquet` + `biogrid_rat_pair_features_report.json`.

**Not part of this flow:** `DomainRecord`/`DomainBuilder`/`DomainFeaturePipeline`, `MotifRecord`/`MotifBuilder`/`MotifFeaturePipeline`, `matching_domain_motif`, `features/pair/domain_motif.py`, `features/motif_properties_of_certain_domain/`, `RegionFeaturePipeline`, `domain_detector`. The production model matrix today is built **only** from whole-protein sequence-composition features (`features/sequence/`) plus a handful of pair-level derived numbers — no domain or motif information reaches the model yet.

## 2. Domain/motif knowledge flow (exact-match layer — implemented, tested, NOT connected to the production pipeline above)

```
ELMDataSource                              (data_sources/elm.py)          — motif class + regex
ELMInteractionDomainsDataSource            (data_sources/elm_interaction_domains.py)  — motif ↔ Pfam/SMART/InterPro domain
    ↓  build_recognized_by_map(...)
MotifDefinitionBuilder(recognized_by_overrides=...)     (builders/motif_definition_builder.py)
    ↓  entities.motif_definition.MotifDefinition  (general knowledge, recognized_by now populated)
to_scanner_rules(...)                       (matching_domain_motif/motif_definition_adapter.py)
    ↓  matching_domain_motif.domain_motif.MotifDefinition  (executable scanner rule)
DomainMotifRegistry                         (matching_domain_motif/domain_motif.py)
    ↓
KnownMotifScanner  +  a DomainDetector implementation
    (production: none exists yet — only StaticDomainDetector, a test fixture)
    ↓
ProteinPairAnalyzer.analyze(protein_a, protein_b)     (matching_domain_motif/pair_analyzer.py)
    ↓
PairAnalysisResult  (bidirectional DomainMotifCandidate list)
```

This chain is proven correct end-to-end against **real cached ELM data** in `tests/integration/test_real_elm_domain_motif_matching.py`, but:
- No production `DomainDetector` exists — only the test-only `StaticDomainDetector` lookup table.
- Nothing in `pipelines/` calls `ProteinPairAnalyzer` or consumes `PairAnalysisResult`.
- The separate `MotifBuilder`→`MotifRecord` path (whole-protein motif annotation, not domain-triggered) and `features/motif_properties_of_certain_domain/` (motif similarity scoring) are also fully working but likewise unconnected to each other and to the pipeline above.

This matches the checkpoint's own boundary statement: exact domain-motif matching is validated; motif-instance/fuzzy candidate-generation work has not been started, and none of this subsystem feeds the BioGRID training pipeline yet.

---

# Final Summary

## 1. CORE FILES (the working production pipeline)
- `protein_graph/data_sources/biogrid.py` — `BioGRIDDataSource`
- `protein_graph/data_sources/biogrid_evidence_policy.py` — evidence policy
- `protein_graph/data_sources/uniprot.py` — `UniProtDataSource`
- `protein_graph/mapping/uniprot.py` — `BioGRIDUniProtMapper`
- `protein_graph/protein_builder/builder.py` — `ProteinBuilder`
- `protein_graph/features/sequence/*.py` — the 5 production feature extractors
- `protein_graph/pipelines/biogrid_pair_acceptance.py` — `BioGRIDPairAcceptancePipeline`
- `protein_graph/pipelines/protein_features.py` — `ProteinFeaturePipeline`
- `protein_graph/pipelines/pair_features.py` — `PairFeaturePipeline` (+ `MODEL_FEATURE_COLUMNS`)
- `protein_graph/pipelines/training_dataset.py` — `TrainingDatasetPipeline`
- `scripts/build_training_dataset.py` — production entry point

## 2. SUPPORTING FILES (reusable, tested, working — but not part of the core call chain above)
- `entities/*` (`SequenceRegion`, `DomainRecord`, `MotifRecord`, `MotifDefinition`)
- `builders/*` (`DomainBuilder`, `MotifBuilder`, `MotifDefinitionBuilder`)
- `data_sources/elm.py`, `data_sources/elm_interaction_domains.py`
- `features/domains.py`, `features/elm_motifs.py`, `features/motif_pipeline.py`
- `pipelines/domain_features.py`, `pipelines/region_features.py`
- `matching_domain_motif/*` (scanner engine, proven correct, unconnected)
- `domain_detector/static_detector.py` (test-only stand-in for a real detector)
- `features/pair/domain_motif.py` (label-overlap matching)
- `features/motif_properties_of_certain_domain/*` (motif similarity/context feature engineering)

## 3. DORMANT OR NOT-YET-INTEGRATED FILES
- `features/registry.py`, `features/pipeline.py` — unused generic plugin scaffolding
- `features/elm_profile.py` (`ELMMotifExtractor`) — independent, unconnected ELM summarizer
- `features/sequence_properties.py` (`SequencePropertyExtractor`) — unused bundler, duplicated by `ProteinFeaturePipeline`'s own inline extractor list
- `features/interface.py`, `localization.py`, `embeddings.py`, `evolution.py`, `structure.py`, `surface.py` — completely empty placeholders
- `knowledge/go.py`, `kegg.py`, `reactome.py`, `string.py`, `syngo.py` — completely empty placeholders
- The entire domain/motif matching subsystem's connection *to the training pipeline* (the components themselves are tested and working, but no pipeline calls them — see section 2)

## 4. Architectural Responsibility Map

| File | Main classes/functions | Responsibility | Used by | Status |
|---|---|---|---|---|
| `entities/region.py` | `SequenceRegion` | Coordinate/sequence base contract | Domain/Motif records, RegionFeaturePipeline | PRODUCTION (as base) |
| `entities/domain.py` | `DomainRecord` | Domain entity | DomainBuilder, DomainFeaturePipeline | SUPPORTING |
| `entities/motif.py` | `MotifRecord` | Motif occurrence entity | MotifBuilder, MotifFeaturePipeline | SUPPORTING |
| `entities/motif_definition.py` | `MotifDefinition` | General motif-class knowledge | MotifBuilder, MotifDefinitionBuilder, adapter | SUPPORTING |
| `data_sources/base.py` | `BaseDataSource` | download→parse→normalize contract | All 4 data sources | PRODUCTION |
| `data_sources/biogrid.py` | `BioGRIDDataSource` | Ingest raw BioGRID-ALL evidence | build_training_dataset.py | PRODUCTION |
| `data_sources/biogrid_evidence_policy.py` | policy constants/functions | Strong/direct evidence classification | BioGRIDPairAcceptancePipeline | PRODUCTION |
| `data_sources/uniprot.py` | `UniProtDataSource` | Download/normalize UniProt records | ProteinBuilder | PRODUCTION |
| `data_sources/elm.py` | `ELMDataSource` | ELM motif-class regex definitions | MotifDefinitionBuilder (tests) | SUPPORTING |
| `data_sources/elm_interaction_domains.py` | `ELMInteractionDomainsDataSource`, `build_recognized_by_map` | ELM domain↔motif mapping | MotifDefinitionBuilder (tests) | SUPPORTING |
| `mapping/uniprot.py` | `BioGRIDUniProtMapper` | Resolve BioGRID rows to UniProt accessions | BioGRIDPairAcceptancePipeline | PRODUCTION |
| `builders/domain_builder.py` | `DomainBuilder` | Build `DomainRecord`s from UniProt data | DomainFeaturePipeline | SUPPORTING |
| `builders/motif_builder.py` | `MotifBuilder` | Regex-scan protein for motif occurrences | tests only | DORMANT (not in pipeline) |
| `builders/motif_definition_builder.py` | `MotifDefinitionBuilder` | Build validated `MotifDefinition`s | tests only | DORMANT (not in pipeline) |
| `protein_builder/builder.py` | `ProteinBuilder` | Orchestrate canonical protein records | ProteinFeaturePipeline | PRODUCTION |
| `domain_detector/static_detector.py` | `StaticDomainDetector` | Fixed lookup-table domain detector | tests only | TEST UTILITY |
| `matching_domain_motif/domain_motif.py` | `DetectedDomain`, `MotifDefinition`, `DomainMotifRegistry`, `KnownMotifScanner` | Regex-scan partner protein for domain-recognized motifs | ProteinPairAnalyzer, tests | SUPPORTING (unconnected) |
| `matching_domain_motif/motif_definition_adapter.py` | `to_scanner_rule(s)` | Convert general → scanner-rule MotifDefinition | tests | SUPPORTING (unconnected) |
| `matching_domain_motif/pair_analyzer.py` | `ProteinPairAnalyzer`, `PairAnalysisResult` | Bidirectional domain↔motif candidate scan for a pair | tests | SUPPORTING (unconnected) |
| `features/base.py` | `FeatureExtractor` | Extractor ABC | sequence extractors etc. | PRODUCTION (base) |
| `features/registry.py`, `pipeline.py` | `FEATURE_EXTRACTORS`, `FeaturePipeline` | Generic extractor registry (unused) | nobody | DORMANT |
| `features/sequence/*.py` | 5 sequence extractors | Whole-protein numeric features | ProteinFeaturePipeline, RegionFeaturePipeline | PRODUCTION |
| `features/sequence_properties.py` | `SequencePropertyExtractor` | Bundle sequence extractors | nobody in production | DORMANT |
| `features/domains.py` | 4 domain feature extractors | Domain-level features | DomainFeaturePipeline | SUPPORTING (unconnected to training) |
| `features/elm_motifs.py` | Motif geometry/annotation extractors | Motif-occurrence features | MotifFeaturePipeline | SUPPORTING (unconnected) |
| `features/elm_profile.py` | `ELMMotifExtractor` | Alternate protein-level ELM hit summary | own tests only | DORMANT |
| `features/motif_pipeline.py` | `MotifFeaturePipeline` | Assemble per-motif feature table | own tests only | DORMANT |
| `features/motif_properties_of_certain_domain/*` | `WindowInterval`, `MotifSequenceFeatureExtractor`, `MotifContextFeatureExtractor`, `MotifSimilarityFeatureExtractor`, `DomainMotifProfile(Builder)`, `CombinedMotifFeatureExtractor` | Candidate-vs-known-motif feature engineering | own tests only | DORMANT (internally complete) |
| `features/pair/domain_motif.py` | `DomainMotifCompatibilityExtractor`, `DomainMotifMatch` | Label-overlap domain↔motif pairwise features | own tests only | DORMANT |
| `pipelines/biogrid_pair_acceptance.py` | `BioGRIDPairAcceptancePipeline` | Evidence rows → accepted pairs | TrainingDatasetPipeline | PRODUCTION |
| `pipelines/protein_features.py` | `ProteinFeaturePipeline` | Per-protein feature table | TrainingDatasetPipeline | PRODUCTION |
| `pipelines/pair_features.py` | `PairFeaturePipeline`, `MODEL_FEATURE_COLUMNS` | Pair feature table + model-input allowlist | TrainingDatasetPipeline | PRODUCTION |
| `pipelines/training_dataset.py` | `TrainingDatasetPipeline` | Top-level orchestrator | scripts/build_training_dataset.py | PRODUCTION |
| `pipelines/domain_features.py` | `DomainFeaturePipeline` | Per-domain feature table | own tests only | DORMANT |
| `pipelines/region_features.py` | `RegionFeaturePipeline` | Generic per-region sequence features | DomainFeaturePipeline | SUPPORTING (unconnected to training) |
| `knowledge/*.py` | — | Planned GO/KEGG/Reactome/STRING/SynGO integration | nobody | DORMANT (empty) |
