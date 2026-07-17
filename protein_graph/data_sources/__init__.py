"""External biological data sources."""

from protein_graph.data_sources.base import BaseDataSource
from protein_graph.data_sources.biogrid import BioGRIDDataSource
from protein_graph.data_sources.elm import ELMDataSource
from protein_graph.data_sources.uniprot import UniProtDataSource

__all__ = [
    "BaseDataSource",
    "BioGRIDDataSource",
    "ELMDataSource",
    "UniProtDataSource",
]
