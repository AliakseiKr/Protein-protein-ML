"""Build domain-level feature tables."""

from __future__ import annotations

from collections.abc import Iterable
from typing import Any, Protocol

import pandas as pd

from protein_graph.builders import DomainBuilder
from protein_graph.entities import DomainRecord
from protein_graph.features.base import FeatureMap
from protein_graph.features.domains import (
    DomainAnnotationFeatureExtractor,
    DomainFeatureContext,
    DomainFunctionFeatureExtractor,
    DomainGeometryFeatureExtractor,
    DomainStructureFeatureExtractor,
)
from protein_graph.pipelines.region_features import RegionFeaturePipeline


class DomainFeatureExtractor(Protocol):
    """Contract for domain-level feature extractors."""

    def extract(
        self,
        domain: DomainRecord,
        context: DomainFeatureContext,
    ) -> FeatureMap:
        """Extract features from one domain."""
        ...


class DomainFeaturePipeline:
    """Build sequence and annotation features for protein domains."""

    def __init__(
        self,
        domain_builder: DomainBuilder | None = None,
        region_pipeline: RegionFeaturePipeline | None = None,
        extractors: Iterable[DomainFeatureExtractor] | None = None,
    ) -> None:
        self.domain_builder = domain_builder or DomainBuilder()
        self.region_pipeline = (
            region_pipeline or RegionFeaturePipeline()
        )

        if extractors is None:
            extractors = (
                DomainGeometryFeatureExtractor(),
                DomainFunctionFeatureExtractor(),
                DomainAnnotationFeatureExtractor(),
                DomainStructureFeatureExtractor(),
            )

        self.extractors = tuple(extractors)

    def build(
        self,
        protein: dict[str, Any],
    ) -> pd.DataFrame:
        """Build a feature table for all domains in one protein."""
        domains = self.domain_builder.build_many(protein)

        return self.build_from_domains(
            domains=domains,
            protein_length=len(protein["sequence"]),
        )

    def build_from_domains(
        self,
        domains: list[DomainRecord],
        protein_length: int,
    ) -> pd.DataFrame:
        """Build features from already constructed DomainRecord objects."""
        context = DomainFeatureContext(
            protein_length=protein_length,
        )

        rows: list[dict[str, Any]] = []

        for domain_index, domain in enumerate(domains):
            features: dict[str, Any] = self.region_pipeline.build(
                domain
            )

            for extractor in self.extractors:
                extracted_features = extractor.extract(
                    domain=domain,
                    context=context,
                )

                duplicate_names = (
                    features.keys()
                    & extracted_features.keys()
                )

                if duplicate_names:
                    names = ", ".join(
                        sorted(duplicate_names)
                    )

                    raise ValueError(
                        "Duplicate domain feature names: "
                        f"{names}."
                    )

                features.update(extracted_features)

            annotation_features: dict[str, Any] = {
                "domain_index": domain_index,
                "domain_id": domain.domain_id,
                "pfam_id": domain.pfam_id,
                "interpro_id": domain.interpro_id,
                "smart_id": domain.smart_id,
                "cdd_id": domain.cdd_id,
                "family": domain.family,
                "superfamily": domain.superfamily,
                "clan": domain.clan,
                "position_category": domain.position_category(
                    protein_length
                ),
                "relative_position": domain.relative_position(
                    protein_length
                ),
                "functional_roles": list(
                    domain.functional_roles
                ),
                "binding_types": list(
                    domain.binding_types
                ),
                "recognized_motifs": list(
                    domain.recognized_motifs
                ),
                "structural_class": domain.structural_class,
                "structural_architecture": (
                    domain.structural_architecture
                ),
                "structural_fold": domain.structural_fold,
                "known_ligands": list(
                    domain.known_ligands
                ),
                "known_ions": list(
                    domain.known_ions
                ),
                "known_cofactors": list(
                    domain.known_cofactors
                ),
            }

            duplicate_names = (
                features.keys()
                & annotation_features.keys()
            )

            if duplicate_names:
                names = ", ".join(
                    sorted(duplicate_names)
                )

                raise ValueError(
                    "Duplicate domain feature names: "
                    f"{names}."
                )

            features.update(annotation_features)
            rows.append(features)

        return pd.DataFrame(rows)

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}()"
