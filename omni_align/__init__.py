"""
omni_align - MUM-anchored sequence alignment.

Uses Maximal Unique Matches as landmarks for optimal global alignment
with reduced search space. Supports NW or WFA for gap alignment.
"""

from .suffix_array import SuffixArray, build_suffix_array, build_lcp_array
from .mum import (
    MUM,
    find_mums,
    can_transition,
    TransitionMatrix,
    build_transition_matrix,
    START,
    END,
)
from .dag import DAG, build_dag
from .segment import (
    Segment,
    AlignmentResult,
    segment,
    get_gap_segments,
    align_segments,
    align_gap,
    align_path,
)

__version__ = "0.1.0"
__all__ = [
    "SuffixArray",
    "build_suffix_array",
    "build_lcp_array",
    "MUM",
    "find_mums",
    "can_transition",
    "TransitionMatrix",
    "build_transition_matrix",
    "START",
    "END",
    "DAG",
    "build_dag",
    "Segment",
    "AlignmentResult",
    "segment",
    "get_gap_segments",
    "align_segments",
    "align_gap",
    "align_path",
]
