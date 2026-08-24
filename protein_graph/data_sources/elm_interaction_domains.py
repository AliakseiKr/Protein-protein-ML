"""ELM domain-motif interaction knowledge data source.

Structured, authoritative domain <-> motif knowledge, distinct from
``elm.ELMDataSource`` (motif class regex definitions) and from motif
*instances* in specific proteins. This is ELM's own curated "Interaction
Domains" resource: for each ELM motif class, which Pfam domain(s) are known
to recognize it.

Verified schema (2026-08-24), fetched from
``http://elm.eu.org/infos/browse_elm_interactiondomains.tsv``::

    "ELM identifier"    "Interaction Domain Id"    "Interaction Domain Description"    "Interaction Domain Name"
    "CLV_NRD_NRD_1"     "PF00675"                  "Peptidase_M16"                     "Insulinase (Peptidase family M16)"

"ELM identifier" matches the "ELMIdentifier"/motif_id values produced by
``ELMDataSource`` (both come from the same ELM motif-class namespace, e.g.
"CLV_NRD_NRD_1"). "Interaction Domain Id" is a Pfam accession, the same
namespace as ``entities.domain.DomainRecord.pfam_id``.

This is the authoritative source used to populate the previously hardcoded
``recognized_by = ()`` on ELM-derived ``MotifDefinition`` entities - see
``protein_graph.builders.motif_definition_builder.MotifDefinitionBuilder``.
No domain<->motif relationship is ever inferred from motif names or
descriptions; every relationship here comes directly from this curated
ELM export.
"""

from __future__ import annotations

import csv
import re
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any, TextIO

import requests

from protein_graph.data_sources.base import BaseDataSource

#: ELM's "Interaction Domain Id" column is not exclusively Pfam, despite
#: the column's own name suggesting it: real ELM data mixes Pfam ("PF...")
#: with SMART ("SM...") and InterPro ("IPR...") accessions (verified
#: 2026-08-24 against the live export: 343 Pfam, 59 InterPro, 7 SMART out
#: of 409 rows). The namespace is detected from the identifier's own
#: official prefix convention - a structural fact about the ID string
#: itself, not an inference from any name or description text - so a
#: downstream consumer can join against the correct field on DomainRecord
#: (pfam_id / smart_id / interpro_id) instead of everything being
#: mislabeled "Pfam".
_DOMAIN_NAMESPACE_PATTERNS: tuple[tuple[re.Pattern[str], str], ...] = (
    (re.compile(r"^PF\d+$"), "Pfam"),
    (re.compile(r"^SM\d+$"), "SMART"),
    (re.compile(r"^IPR\d+$"), "InterPro"),
    (re.compile(r"^cd\d+$", re.IGNORECASE), "CDD"),
)


def build_recognized_by_map(
    records: Iterable[Mapping[str, Any]],
) -> dict[str, tuple[str, ...]]:
    """Join normalized interaction-domain records into motif_id -> domains.

    Groups the (motif_id, domain_id) pairs produced by
    ``ELMInteractionDomainsDataSource._normalize`` into one recognizing-
    domain tuple per motif, preserving first-seen order and biological
    multiplicity (one motif can legitimately be recognized by several
    domain families; that is represented here as multiple domain_id
    values for one motif_id, not as separate motif definitions).
    """
    by_motif: dict[str, list[str]] = {}

    for record in records:
        motif_id = str(record["motif_id"]).strip()
        domain_id = str(record["domain_id"]).strip()

        if not motif_id or not domain_id:
            continue

        domains = by_motif.setdefault(motif_id, [])

        if domain_id not in domains:
            domains.append(domain_id)

    return {
        motif_id: tuple(domains)
        for motif_id, domains in by_motif.items()
    }


class ELMInteractionDomainsDataSource(BaseDataSource):
    """Load and normalize ELM's motif-class <-> recognizing-domain mapping."""

    RAW_FILENAME = "elm_interaction_domains.tsv"

    _COLUMN_ALIASES: dict[str, tuple[str, ...]] = {
        "motif_id": (
            "ELM identifier",
            "ELMIdentifier",
            "motif_id",
        ),
        "domain_id": (
            "Interaction Domain Id",
            "domain_id",
        ),
        "domain_description": (
            "Interaction Domain Description",
            "domain_description",
        ),
        "domain_name": (
            "Interaction Domain Name",
            "domain_name",
        ),
    }

    _REQUIRED_FIELDS = (
        "motif_id",
        "domain_id",
    )

    def __init__(
        self,
        raw_dir: str | Path = "data/raw/elm",
        *,
        url: str | None = (
            "http://elm.eu.org/infos/browse_elm_interactiondomains.tsv"
        ),
        timeout: int = 30,
    ) -> None:
        super().__init__(raw_dir=raw_dir)

        if url is not None:
            url = url.strip()

            if not url:
                raise ValueError(
                    "ELM interaction-domains source URL cannot be empty."
                )

        if timeout <= 0:
            raise ValueError("timeout must be greater than zero.")

        self.url = url
        self.timeout = timeout

    @property
    def raw_path(self) -> Path:
        """Return the local ELM interaction-domains TSV path."""
        return self.raw_dir / self.RAW_FILENAME

    def _download(self, force: bool = False) -> None:
        """Ensure that the ELM interaction-domains TSV is available locally."""
        if self.raw_path.exists() and not force:
            return

        if self.url is None:
            raise RuntimeError(
                "ELM interaction-domains URL is not configured and no "
                f"cached file is available at {self.raw_path}."
            )

        try:
            response = requests.get(
                self.url,
                timeout=self.timeout,
                headers={
                    "Accept": (
                        "text/tab-separated-values, "
                        "text/plain;q=0.9, */*;q=0.1"
                    ),
                    "User-Agent": "ProteinGraph/0.1",
                },
            )
            response.raise_for_status()
        except requests.RequestException as error:
            raise RuntimeError(
                "Failed to download ELM interaction-domains data "
                f"from {self.url!r}."
            ) from error

        content = response.text

        if not content.strip():
            raise ValueError(
                "ELM returned an empty interaction-domains file."
            )

        self.raw_path.write_text(content, encoding="utf-8")

    def _parse(self) -> list[dict[str, str]]:
        """Parse the locally cached ELM interaction-domains TSV file."""
        if not self.raw_path.exists():
            raise FileNotFoundError(
                f"ELM interaction-domains raw file does not exist: "
                f"{self.raw_path}"
            )

        with self.raw_path.open(
            "r",
            encoding="utf-8-sig",
            newline="",
        ) as handle:
            return self._parse_tsv(handle)

    def _normalize(
        self,
        records: list[dict[str, str]],
    ) -> list[dict[str, Any]]:
        """Convert ELM interaction-domains rows into canonical records."""
        normalized: list[dict[str, Any]] = []
        seen_pairs: set[tuple[str, str]] = set()

        for row_number, record in enumerate(records, start=2):
            motif_id = self._required_value(
                record,
                canonical_name="motif_id",
                row_number=row_number,
            )

            domain_id = self._required_value(
                record,
                canonical_name="domain_id",
                row_number=row_number,
            )

            pair_key = (motif_id, domain_id)

            if pair_key in seen_pairs:
                continue

            seen_pairs.add(pair_key)

            normalized.append(
                {
                    "motif_id": motif_id,
                    "domain_id": domain_id,
                    "domain_namespace": self._detect_domain_namespace(
                        domain_id
                    ),
                    "domain_description": self._optional_value(
                        record,
                        canonical_name="domain_description",
                    ),
                    "domain_name": self._optional_value(
                        record,
                        canonical_name="domain_name",
                    ),
                    "source": "ELM",
                }
            )

        return normalized

    def _parse_tsv(
        self,
        handle: TextIO,
    ) -> list[dict[str, str]]:
        """Parse a TSV stream and ignore empty/comment lines."""
        reader = csv.DictReader(
            self._non_comment_lines(handle),
            delimiter="\t",
        )

        if reader.fieldnames is None:
            raise ValueError(
                "ELM interaction-domains TSV file does not contain a "
                "header."
            )

        fieldnames = {
            field.strip()
            for field in reader.fieldnames
            if field is not None and field.strip()
        }

        if not fieldnames:
            raise ValueError(
                "ELM interaction-domains TSV header is empty."
            )

        self._validate_header(fieldnames)

        records: list[dict[str, str]] = []

        for row in reader:
            cleaned = {
                str(key).strip(): (
                    value.strip()
                    if isinstance(value, str)
                    else ""
                )
                for key, value in row.items()
                if key is not None
            }

            if any(cleaned.values()):
                records.append(cleaned)

        return records

    @staticmethod
    def _detect_domain_namespace(domain_id: str) -> str:
        """Detect a domain identifier's namespace from its own prefix.

        Structural pattern match on the identifier string itself (an
        official convention: Pfam accessions start "PF", SMART "SM",
        InterPro "IPR"), never on the accompanying name/description text.
        Returns "unknown" for an id whose namespace isn't one of the
        conventions seen in real ELM interaction-domains data, rather than
        guessing.
        """
        for pattern, namespace in _DOMAIN_NAMESPACE_PATTERNS:
            if pattern.match(domain_id):
                return namespace

        return "unknown"

    @staticmethod
    def _non_comment_lines(
        handle: TextIO,
    ) -> Iterable[str]:
        """Yield non-empty lines that are not comments."""
        for line in handle:
            if not line.strip():
                continue

            if line.lstrip().startswith("#"):
                continue

            yield line

    def _validate_header(
        self,
        fieldnames: set[str],
    ) -> None:
        """Check that the TSV header contains required fields."""
        missing_fields = [
            canonical_name
            for canonical_name in self._REQUIRED_FIELDS
            if not self._has_alias(fieldnames, canonical_name)
        ]

        if not missing_fields:
            return

        expected = {
            name: self._COLUMN_ALIASES[name]
            for name in missing_fields
        }

        raise ValueError(
            "Unsupported ELM interaction-domains TSV header. Missing "
            f"canonical fields {missing_fields!r}. Accepted aliases: "
            f"{expected!r}. Actual columns: {sorted(fieldnames)!r}."
        )

    def _has_alias(
        self,
        fieldnames: set[str],
        canonical_name: str,
    ) -> bool:
        """Return whether a canonical field has a known column."""
        return any(
            alias in fieldnames
            for alias in self._COLUMN_ALIASES[canonical_name]
        )

    def _required_value(
        self,
        record: dict[str, str],
        *,
        canonical_name: str,
        row_number: int,
    ) -> str:
        """Read a required normalized field from one ELM row."""
        value = self._optional_value(
            record,
            canonical_name=canonical_name,
        )

        if value is None:
            raise ValueError(
                f"ELM interaction-domains row {row_number} has an empty "
                f"{canonical_name!r} value."
            )

        return value

    def _optional_value(
        self,
        record: dict[str, str],
        *,
        canonical_name: str,
    ) -> str | None:
        """Read an optional value through its known aliases."""
        aliases = self._COLUMN_ALIASES[canonical_name]

        for alias in aliases:
            value = record.get(alias)

            if isinstance(value, str) and value.strip():
                return value.strip()

        return None

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}("
            f"raw_dir={str(self.raw_dir)!r}, "
            f"url={self.url!r}, "
            f"timeout={self.timeout}"
            f")"
        )
