"""S2SignalPathway package.

Refactored package that separates KEGG and Reactome fetching/parsing logic.
"""

from .s2pathway import PathwaySource, SignalingPathway

__all__ = ["PathwaySource", "SignalingPathway"]
