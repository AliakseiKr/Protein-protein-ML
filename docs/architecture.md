# ProteinGraph Architecture

ProteinGraph separates intrinsic protein features from optional external biological knowledge.

The core model should use features derived from the protein itself: sequence, structure,
domains, localization, surface properties, predicted interfaces, and pairwise physical
compatibility.

External biological knowledge such as STRING, GO, Reactome, KEGG, and SynGO must live in
a separate optional layer to reduce information leakage and allow fair ablation studies.

## Feature layers

### Model A: intrinsic protein features only

Uses sequence, embeddings, domains, localization, structure, surface, interface, and
pairwise compatibility features.

### Model B: intrinsic features + evolution

Adds orthology, paralogy, conservation, phylogenetic profile similarity, and co-evolution.

### Model C: intrinsic features + evolution + biological knowledge

Adds STRING, Reactome, GO, KEGG, SynGO, pathway co-membership, and other external
knowledge-derived features.

## Core protein features

### Sequence features

- protein_id
- sequence
- protein_length
- amino_acid_composition

### Embedding features

- sequence_embedding
- structure_embedding

### Domain features

- domains
- domain_count
- domain_architecture

### Localization features

- localization_terms
- transmembrane_region_count
- has_transmembrane_region
- signal_peptide
- signal_peptide_length

### Structure confidence and disorder-like features

- mean_pLDDT
- low_pLDDT_fraction
- high_pLDDT_fraction
- mean_PAE
- pae_summary

### Secondary structure features

- secondary_structure_fraction_helix
- secondary_structure_fraction_sheet
- secondary_structure_fraction_coil

### Surface features

- mean_solvent_accessibility
- exposed_residue_fraction
- surface_charge
- surface_hydrophobicity
- polar_surface_area

### Compactness and shape features

- radius_of_gyration
- asphericity
- surface_area_to_volume

### Interface features

- predicted_interface_residues
- interface_residue_fraction
- alphafold_interface_score

## Core pair features

### Sequence and embedding similarity

- sequence_similarity
- sequence_embedding_cosine_similarity
- structure_embedding_cosine_similarity

### Domain compatibility

- domain_overlap
- domain_architecture_compatibility

### Localization compatibility

- localization_similarity
- same_localization
- compatible_localization
- both_membrane
- one_secreted_one_membrane
- tm_region_difference
- signal_peptide_match

### Structure compatibility

- plddt_difference
- pae_compatibility
- secondary_structure_similarity

### Surface compatibility

- solvent_accessibility_difference
- surface_charge_difference
- surface_hydrophobicity_difference

### Compactness compatibility

- radius_of_gyration_difference
- compactness_similarity

### Interface compatibility

- interface_residue_compatibility
- alphafold_interface_score

## Evolution layer

- ortholog_count
- has_mouse_ortholog
- has_human_ortholog
- paralog_count
- sequence_conservation
- ortholog_interaction_transfer
- co_evolution_score
- phylogenetic_profile_similarity

## Optional biological knowledge layer

These features are useful but must remain optional because they can introduce leakage from
known biological databases.

- STRING score
- Reactome overlap
- GO overlap
- KEGG overlap
- SynGO overlap
- known pathway co-membership
