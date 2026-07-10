# ProteinGraph Architecture

## 1. Overview

ProteinGraph is a modular framework for predicting protein-protein interactions (PPIs)
from intrinsic protein properties while keeping external biological knowledge
strictly separated.

The project is designed around three principles:

- modularity
- reproducibility
- extensibility

Each component has a single responsibility and can be developed independently.

---

## 2. Design Philosophy

ProteinGraph separates the complete workflow into independent stages.

```text
External databases
        │
        ▼
Data ingestion
        │
        ▼
Canonical schemas
        │
        ▼
Protein feature extraction
        │
        ▼
Protein feature cache
        │
        ▼
Pair feature generation
        │
        ▼
Machine learning models
```

Each stage communicates through well-defined interfaces.

This architecture makes it possible to replace any component without affecting the
rest of the project.

---

## 3. Core Philosophy

The project distinguishes between intrinsic protein properties and external
biological knowledge.

Intrinsic properties describe the protein itself.

Examples:

- sequence
- domains
- localization
- AlphaFold structure
- surface properties
- embeddings

External biological knowledge describes what is already known about proteins.

Examples:

- STRING
- Reactome
- GO
- KEGG
- SynGO

Keeping these layers separate allows fair benchmarking and reduces information
leakage.

---

## 4. Project Structure

```text
protein_graph/

├── data_sources/
├── schemas/
├── features/
├── knowledge/
├── graph/
├── models/
├── utils/
```

---

## 5. Data Sources

Responsible only for downloading and normalizing public databases.

Examples:

- BioGRID
- UniProt
- AlphaFold
- InterPro
- Pfam
- OrthoDB

A data source should never compute machine-learning features.

Its responsibility ends after producing a normalized DataFrame.

---

## 6. Canonical Schemas

Schemas define the internal representation of every object used by ProteinGraph.

Every external database must be converted into these schemas.

Examples:

- proteins
- interactions
- pairs

No downstream module should depend on database-specific column names.

---

## 7. Feature Extraction

Every feature family is implemented independently.

Feature extractors receive a protein and return only features belonging to
their own domain.

Example:

Protein

↓

LocalizationExtractor

↓

{
    "is_membrane": True,
    "is_secreted": False,
    ...
}

Extractors never depend on each other.

---

## 8. Protein Feature Groups

### Sequence

- protein_id
- sequence
- protein_length
- amino_acid_composition

Source:

UniProt

---

### Embeddings

- sequence_embedding
- structure_embedding

Source:

ESM
ProtT5

---

### Domains

- domains
- domain_count
- domain_architecture

Source:

InterPro
Pfam

---

### Localization

- localization_terms
- is_membrane
- is_secreted
- is_extracellular
- is_cytoplasmic
- is_nuclear
- transmembrane_region_count
- has_transmembrane_region
- signal_peptide
- signal_peptide_length

Source:

UniProt

---

### Secondary Structure

- secondary_structure_fraction_helix
- secondary_structure_fraction_sheet
- secondary_structure_fraction_coil

Source:

AlphaFold

↓

DSSP

---

### Surface

- mean_solvent_accessibility
- exposed_residue_fraction
- surface_charge
- surface_hydrophobicity
- polar_surface_area

Source:

AlphaFold

↓

FreeSASA

---

## 9. Protein Feature Cache

Protein features should be computed only once.

```text
Protein

↓

Feature extractors

↓

Protein feature table

↓

Parquet cache
```

Whenever possible, downstream modules should reuse cached features.

---

## 10. Pair Feature Generation

Pairwise features are generated from cached protein features.

No external databases should be queried during this stage.

### Sequence

- sequence_similarity
- sequence_embedding_cosine_similarity
- structure_embedding_cosine_similarity

### Domains

- domain_overlap
- domain_architecture_compatibility

### Localization

- localization_similarity
- same_localization
- compatible_localization
- both_membrane
- one_secreted_one_membrane
- tm_region_difference
- signal_peptide_match

### Secondary Structure

- secondary_structure_similarity

### Surface

- solvent_accessibility_difference
- surface_charge_difference
- surface_hydrophobicity_difference

---

## 11. Evolution Layer

Optional feature family.

Examples:

- ortholog_count
- paralog_count
- sequence_conservation
- co_evolution_score
- phylogenetic_profile_similarity

---

## 12. Biological Knowledge Layer

Optional feature family.

Examples:

- STRING score
- GO overlap
- Reactome overlap
- KEGG overlap
- SynGO overlap

These features should never be required for the core model.

---

## 13. Model Configurations

ProteinGraph supports multiple experimental configurations.

### Model A

Intrinsic protein features only.

### Model B

Intrinsic features + evolution.

### Model C

Intrinsic features + evolution + biological knowledge.

---

## 14. Implementation Principles

Feature extraction must be modular.

Each feature family should be independently:

- implemented
- tested
- cached
- enabled
- disabled

Missing feature groups must never prevent model training.

---

## 15. Future Modules

Future modules may include:

- graph construction
- graph neural networks
- feature selection
- explainability
- benchmarking
- visualization
- model interpretation

These modules should build upon the existing architecture without requiring changes
to lower layers.
