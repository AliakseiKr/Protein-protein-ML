"""Core biological entities."""

from protein_graph.entities.domain import DomainRecord
from protein_graph.entities.motif import MotifRecord
from protein_graph.entities.motif_definition import MotifDefinition
from protein_graph.entities.region import SequenceRegion

__all__ = [
    "DomainRecord",
    "MotifDefinition",
    "MotifRecord",
    "SequenceRegion",
]
