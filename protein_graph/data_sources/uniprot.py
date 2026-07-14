"""UniProt data source."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import requests

from protein_graph.data_sources.base import BaseDataSource


class UniProtDataSource(BaseDataSource):
    """Download and normalize UniProtKB protein records."""

    BASE_URL = "https://rest.uniprot.org/uniprotkb"

    def __init__(
        self,
        accessions: list[str],
        raw_dir: str | Path = "data/raw/uniprot",
        timeout: int = 30,
    ) -> None:
        super().__init__(raw_dir=raw_dir)

        cleaned_accessions = {
            accession.strip().upper()
            for accession in accessions
            if accession.strip()
        }

        if not cleaned_accessions:
            raise ValueError("At least one UniProt accession is required.")

        self.accessions = sorted(cleaned_accessions)
        self.timeout = timeout

    def _download(self, force: bool = False) -> None:
        """Download one JSON record for each UniProt accession."""
        for accession in self.accessions:
            destination = self._record_path(accession)

            if destination.exists() and not force:
                continue

            url = f"{self.BASE_URL}/{accession}.json"

            try:
                response = requests.get(
                    url,
                    timeout=self.timeout,
                    headers={"Accept": "application/json"},
                )
                response.raise_for_status()
            except requests.RequestException as error:
                raise RuntimeError(
                    f"Failed to download UniProt accession {accession}."
                ) from error

            try:
                record = response.json()
            except requests.JSONDecodeError as error:
                raise ValueError(
                    f"UniProt returned invalid JSON for {accession}."
                ) from error

            destination.write_text(
                json.dumps(record, indent=2),
                encoding="utf-8",
            )

    def _parse(self) -> list[dict[str, Any]]:
        """Read cached UniProt JSON records."""
        records = []

        for accession in self.accessions:
            path = self._record_path(accession)

            if not path.exists():
                raise FileNotFoundError(
                    f"Cached UniProt record not found: {path}"
                )

            try:
                record = json.loads(path.read_text(encoding="utf-8"))
            except json.JSONDecodeError as error:
                raise ValueError(
                    f"Invalid cached UniProt JSON: {path}"
                ) from error

            records.append(record)

        return records

    def _normalize(
        self,
        records: list[dict[str, Any]],
    ) -> pd.DataFrame:
        """Normalize UniProt records into a protein table."""
        normalized = [
            self._normalize_record(record)
            for record in records
        ]

        return pd.DataFrame(normalized)

    def _normalize_record(
        self,
        record: dict[str, Any],
    ) -> dict[str, Any]:
        """Normalize one UniProtKB JSON record."""
        accession = record.get("primaryAccession")

        if not accession:
            raise ValueError(
                "UniProt record is missing primaryAccession."
            )

        sequence = record.get("sequence", {}).get("value", "")

        if not sequence:
            raise ValueError(
                f"UniProt record {accession} has no sequence."
            )

        organism = record.get("organism", {})
        features = record.get("features", [])

        signal_peptides = self._extract_feature_ranges(
            features,
            feature_types={"Signal"},
        )
        transmembrane_regions = self._extract_feature_ranges(
            features,
            feature_types={"Transmembrane"},
        )
        domains = self._extract_feature_descriptions(
            features,
            feature_types={"Domain"},
        )

        return {
            "protein_id": accession,
            "entry_name": record.get("uniProtkbId"),
            "reviewed": record.get("entryType") == "UniProtKB reviewed (Swiss-Prot)",
            "gene_name": self._extract_gene_name(record),
            "protein_name": self._extract_protein_name(record),
            "organism": organism.get("scientificName"),
            "taxon_id": organism.get("taxonId"),
            "sequence": sequence,
            "protein_length": len(sequence),
            "localization_terms": self._extract_localization_terms(record),
            "domains": domains,
            "domain_count": len(domains),
            "has_signal_peptide": bool(signal_peptides),
            "signal_peptide_ranges": signal_peptides,
            "signal_peptide_length": self._total_feature_length(
                signal_peptides
            ),
            "has_transmembrane_region": bool(transmembrane_regions),
            "transmembrane_region_count": len(transmembrane_regions),
            "transmembrane_region_ranges": transmembrane_regions,
            "interpro_ids": self._extract_cross_references(
                record,
                database="InterPro",
            ),
            "pfam_ids": self._extract_cross_references(
                record,
                database="Pfam",
            ),
            "refseq_ids": self._extract_cross_references(
                record,
                database="RefSeq",
            ),
        }

    def _record_path(self, accession: str) -> Path:
        return self.raw_dir / f"{accession}.json"

    @staticmethod
    def _extract_gene_name(record: dict[str, Any]) -> str | None:
        genes = record.get("genes", [])

        for gene in genes:
            value = gene.get("geneName", {}).get("value")
            if value:
                return value

        return None

    @staticmethod
    def _extract_protein_name(record: dict[str, Any]) -> str | None:
        description = record.get("proteinDescription", {})

        recommended = description.get("recommendedName", {})
        full_name = recommended.get("fullName", {}).get("value")

        if full_name:
            return full_name

        submission_names = description.get("submissionNames", [])

        if submission_names:
            return (
                submission_names[0]
                .get("fullName", {})
                .get("value")
            )

        return None

    @staticmethod
    def _extract_localization_terms(
        record: dict[str, Any],
    ) -> list[str]:
        terms = set()

        for comment in record.get("comments", []):
            if comment.get("commentType") != "SUBCELLULAR LOCATION":
                continue

            for item in comment.get("subcellularLocations", []):
                value = item.get("location", {}).get("value")

                if value:
                    terms.add(value)

        return sorted(terms)

    @staticmethod
    def _extract_feature_descriptions(
        features: list[dict[str, Any]],
        feature_types: set[str],
    ) -> list[str]:
        descriptions = set()

        for feature in features:
            if feature.get("type") not in feature_types:
                continue

            description = feature.get("description")

            if description:
                descriptions.add(description)

        return sorted(descriptions)

    @classmethod
    def _extract_feature_ranges(
        cls,
        features: list[dict[str, Any]],
        feature_types: set[str],
    ) -> list[dict[str, int]]:
        ranges = []

        for feature in features:
            if feature.get("type") not in feature_types:
                continue

            location = feature.get("location", {})
            start = cls._extract_position(location.get("start"))
            end = cls._extract_position(location.get("end"))

            if start is not None and end is not None:
                ranges.append({
                    "start": start,
                    "end": end,
                })

        return ranges

    @staticmethod
    def _extract_position(
        position: dict[str, Any] | None,
    ) -> int | None:
        if not position:
            return None

        value = position.get("value")

        if isinstance(value, int):
            return value

        return None

    @staticmethod
    def _total_feature_length(
        ranges: list[dict[str, int]],
    ) -> int:
        return sum(
            item["end"] - item["start"] + 1
            for item in ranges
        )

    @staticmethod
    def _extract_cross_references(
        record: dict[str, Any],
        database: str,
    ) -> list[str]:
        identifiers = set()

        for cross_reference in record.get(
            "uniProtKBCrossReferences",
            [],
        ):
            if cross_reference.get("database") != database:
                continue

            identifier = cross_reference.get("id")

            if identifier:
                identifiers.add(identifier)

        return sorted(identifiers)

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}("
            f"accessions={len(self.accessions)}, "
            f"raw_dir='{self.raw_dir}')"
        )
