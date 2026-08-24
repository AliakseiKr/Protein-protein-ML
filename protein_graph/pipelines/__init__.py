"""ProteinGraph processing pipelines."""

from protein_graph.pipelines.biogrid_pair_acceptance import (
    BioGRIDPairAcceptancePipeline,
)
from protein_graph.pipelines.protein_features import ProteinFeaturePipeline
from protein_graph.pipelines.pair_features import PairFeaturePipeline
from protein_graph.pipelines.training_dataset import TrainingDatasetPipeline
from protein_graph.pipelines.domain_features import DomainFeaturePipeline
from protein_graph.pipelines.region_features import RegionFeaturePipeline

__all__ = [
    "RegionFeaturePipeline",
    "DomainFeaturePipeline",
    "ProteinFeaturePipeline",
    "PairFeaturePipeline",
    "BioGRIDPairAcceptancePipeline",
    "TrainingDatasetPipeline",
]
