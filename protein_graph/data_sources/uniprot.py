"""UniProt data source."""

from __future__ import annotations

from typing import Any

import requests


class UniProtDataSource:
    """Fetch protein records from UniProt."""

    BASE_URL = "https://rest.uniprot.org/uniprotkb"

    def fetch_accession(self, accession: str) -> dict[str, Any]:
        """Fetch one UniProtKB record by accession."""
        url = f"{self.BASE_URL}/{accession}.json"
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        return response.json()

    def normalize_record(self, record: dict[str, Any]) -> dict[str, Any]:
        """Normalize a UniProt record into ProteinGraph protein format."""
        sequence = record.get("sequence", {}).get("value", "")

        genes = record.get("genes", [])
        gene_name = None
        if genes:
            gene_name = genes[0].get("geneName", {}).get("value")

        organism = record.get("organism", {})
        taxon_id = organism.get("taxonId")
        organism_name = organism.get("scientificName")

        return {
            "protein_id": record.get("primaryAccession"),
            "gene_name": gene_name,
            "organism": organism_name,
            "taxon_id": taxon_id,
            "sequence": sequence,
            "protein_length": len(sequence),
        }

    def fetch(self, accessions: list[str]) -> list[dict[str, Any]]:
        """Fetch and normalize multiple UniProt records."""
        proteins = []

        for accession in accessions:
            record = self.fetch_accession(accession)
            proteins.append(self.normalize_record(record))

        return proteins
