## Feature Categories

### Sequence features

Describe the amino acid sequence itself.

- protein_id — Unique protein identifier.
- sequence — Amino acid sequence.
- protein_length — Number of amino acids.
- amino_acid_composition — Frequency of each amino acid in the sequence.

### Embedding features

High-dimensional numerical representations learned by protein language models.

- sequence_embedding — Learned representation of the protein sequence.
- structure_embedding — Learned representation of the protein structure.

These embeddings capture complex biological patterns that are difficult to engineer manually.

### Domain features

Describe the functional building blocks of a protein.

- domains — List of detected domains.
- domain_count — Number of domains.
- domain_architecture — Order and arrangement of domains.

### Localization features

Describe where the protein is expected to function inside the cell.

- localization_terms — Predicted or known cellular locations.
- is_membrane — Whether the protein is membrane-associated.
- is_secreted — Whether the protein is secreted.
- is_extracellular — Whether the protein is extracellular.
- is_cytoplasmic — Whether the protein is cytoplasmic.
- is_nuclear — Whether the protein is nuclear.
- transmembrane_region_count — Number of membrane-spanning regions.
- has_transmembrane_region — Whether the protein has at least one transmembrane region.
- signal_peptide — Presence of a signal peptide.
- signal_peptide_length — Length of the signal peptide.

### Surface features

Describe the physical properties of the protein surface.

- mean_solvent_accessibility — Average exposure to solvent.
- exposed_residue_fraction — Fraction of surface-exposed residues.
- surface_charge — Overall surface electrostatic properties.
- surface_hydrophobicity — Overall surface hydrophobicity.
- polar_surface_area — Accessible polar surface area.

## Pairwise Feature Categories

### Domain compatibility

Evaluate whether domain organization supports interaction.

- domain_overlap — Shared domains.
- domain_architecture_compatibility — Compatibility of domain arrangements.

### Localization compatibility

Determine whether proteins are likely to meet in the same cellular environment.

- localization_similarity — Similarity of localization.
- same_localization — Same compartment.
- compatible_localization — Compatible compartments.
- both_membrane — Both proteins are membrane-associated.
- one_secreted_one_membrane — Secreted protein paired with a membrane protein.
- tm_region_difference — Difference in transmembrane content.
- signal_peptide_match — Similar signal peptide properties.

### Surface compatibility

Compare physicochemical surface properties.

- solvent_accessibility_difference — Difference in solvent exposure.
- surface_charge_difference — Difference in surface charge.
- surface_hydrophobicity_difference — Difference in hydrophobicity.

### Interface compatibility

Estimate whether two proteins have compatible interaction surfaces.

- interface_residue_compatibility — Compatibility of predicted interface regions.
- alphafold_interface_score — Predicted interface compatibility score.

## Evolution Layer

Evolutionary features describe how proteins have changed across species.

- ortholog_count — Number of orthologs.
- has_mouse_ortholog — Presence of a mouse ortholog.
- has_human_ortholog — Presence of a human ortholog.
- paralog_count — Number of paralogs.
- sequence_conservation — Degree of evolutionary conservation.
- ortholog_interaction_transfer — Evidence transferred from ortholog interactions.
- co_evolution_score — Degree of coordinated evolution.
- phylogenetic_profile_similarity — Similarity of evolutionary profiles.

## Biological Knowledge Layer

This layer contains curated biological knowledge rather than intrinsic protein properties.

- STRING score — Evidence-based interaction confidence.
- Reactome overlap — Shared biological pathways.
- GO overlap — Shared Gene Ontology annotations.
- KEGG overlap — Shared KEGG pathways.
- SynGO overlap — Shared synaptic annotations.
- known pathway co-membership — Whether both proteins belong to the same curated pathway.
