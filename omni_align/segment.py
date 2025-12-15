"""Segment extraction and alignment between MUM landmarks."""

from typing import Tuple, Optional, Union, List, Dict, Any
from dataclasses import dataclass

try:
    from .mum import MUM, START, END
except ImportError:
    from mum import MUM, START, END

from nw import needleman_wunsch
from wfa import wfa_align


@dataclass
class Segment:
    """A segment of a sequence between two positions."""
    sequence: str
    start: int
    end: int
    
    @property
    def length(self) -> int:
        return len(self.sequence)


@dataclass
class AlignmentResult:
    """Result of aligning two segments."""
    score: int
    cigar: str
    seg1: Segment
    seg2: Segment
    cells_visited: int = 0


def segment(seq: str, start: int, end: int) -> Segment:
    """Extract a segment from a sequence."""
    start = max(0, start)
    end = min(len(seq), end)
    return Segment(sequence=seq[start:end], start=start, end=end)


def get_gap_segments(seq1: str, seq2: str, mum_from: Optional[MUM], mum_to: Optional[MUM], from_is_start: bool = False, to_is_end: bool = False) -> Tuple[Segment, Segment]:
    """Extract segments between two MUMs from both sequences."""
    n1, n2 = len(seq1), len(seq2)
    
    if from_is_start or mum_from is None:
        start1, start2 = 0, 0
    else:
        start1, start2 = mum_from.end1, mum_from.end2
    
    if to_is_end or mum_to is None:
        end1, end2 = n1, n2
    else:
        end1, end2 = mum_to.pos1, mum_to.pos2
    
    return segment(seq1, start1, end1), segment(seq2, start2, end2)


def align_segments(seg1: Segment, seg2: Segment, use_wfa: bool = False, mismatch: int = 4, gap_open: int = 6, gap_extend: int = 2) -> AlignmentResult:
    """Align two segments using NW or WFA."""
    s1, s2 = seg1.sequence, seg2.sequence
    
    if len(s1) == 0 and len(s2) == 0:
        return AlignmentResult(score=0, cigar="", seg1=seg1, seg2=seg2, cells_visited=0)
    
    if len(s1) == 0:
        score = gap_open + len(s2) * gap_extend
        return AlignmentResult(score=score, cigar=f"{len(s2)}I", seg1=seg1, seg2=seg2, cells_visited=len(s2))
    
    if len(s2) == 0:
        score = gap_open + len(s1) * gap_extend
        return AlignmentResult(score=score, cigar=f"{len(s1)}D", seg1=seg1, seg2=seg2, cells_visited=len(s1))
    
    if use_wfa:
        score, cigar, cells = wfa_align(s1, s2, mismatch=mismatch, gap_open=gap_open, gap_extend=gap_extend, return_cells=True)
    else:
        score, cigar, cells = needleman_wunsch(s1, s2, mismatch=mismatch, gap_open=gap_open, gap_extend=gap_extend, return_cells=True)
    
    return AlignmentResult(score=score, cigar=cigar, seg1=seg1, seg2=seg2, cells_visited=cells)


def align_gap(seq1: str, seq2: str, mum_from: Union[MUM, str, None], mum_to: Union[MUM, str, None], use_wfa: bool = False, mismatch: int = 4, gap_open: int = 6, gap_extend: int = 2) -> AlignmentResult:
    """Align the gap between two MUMs."""
    from_is_start = (mum_from == START or mum_from is None)
    to_is_end = (mum_to == END or mum_to is None)
    
    actual_from = None if from_is_start else mum_from
    actual_to = None if to_is_end else mum_to
    
    seg1, seg2 = get_gap_segments(seq1, seq2, actual_from, actual_to, from_is_start, to_is_end)
    return align_segments(seg1, seg2, use_wfa=use_wfa, mismatch=mismatch, gap_open=gap_open, gap_extend=gap_extend)


def align_path(seq1: str, seq2: str, mum_path: List[MUM], use_wfa: bool = False, mismatch: int = 4, gap_open: int = 6, gap_extend: int = 2) -> Tuple[int, str, List[AlignmentResult], int]:
    """Align all gaps along a path of MUMs. Returns (score, cigar, results, cells)."""
    if not mum_path:
        result = align_gap(seq1, seq2, START, END, use_wfa, mismatch, gap_open, gap_extend)
        return result.score, result.cigar, [result], result.cells_visited
    
    results = []
    total_score = 0
    cigar_parts = []
    
    result = align_gap(seq1, seq2, START, mum_path[0], use_wfa, mismatch, gap_open, gap_extend)
    results.append(result)
    total_score += result.score
    if result.cigar:
        cigar_parts.append(result.cigar)
    cigar_parts.append(f"{mum_path[0].length}M")
    
    for i in range(len(mum_path) - 1):
        result = align_gap(seq1, seq2, mum_path[i], mum_path[i + 1], use_wfa, mismatch, gap_open, gap_extend)
        results.append(result)
        total_score += result.score
        if result.cigar:
            cigar_parts.append(result.cigar)
        cigar_parts.append(f"{mum_path[i + 1].length}M")
    
    result = align_gap(seq1, seq2, mum_path[-1], END, use_wfa, mismatch, gap_open, gap_extend)
    results.append(result)
    total_score += result.score
    if result.cigar:
        cigar_parts.append(result.cigar)
    
    return total_score, "".join(cigar_parts), results, sum(r.cells_visited for r in results)


def compare_algorithms(seq1: str, seq2: str, mum_path: List[MUM]) -> Dict[str, Any]:
    """Compare NW and WFA search space reduction."""
    n, m = len(seq1), len(seq2)
    full_nw_cells = n * m
    
    nw_score, _, _, nw_mum_cells = align_path(seq1, seq2, mum_path, use_wfa=False)
    nw_full = align_gap(seq1, seq2, START, END, use_wfa=False)
    
    wfa_score, _, _, wfa_mum_cells = align_path(seq1, seq2, mum_path, use_wfa=True)
    wfa_full = align_gap(seq1, seq2, START, END, use_wfa=True)
    
    return {
        'seq1_len': n,
        'seq2_len': m,
        'full_nw_cells': full_nw_cells,
        'num_mums': len(mum_path),
        'nw_mum_cells': nw_mum_cells,
        'nw_full_cells': nw_full.cells_visited,
        'nw_mum_reduction_vs_nw': (full_nw_cells - nw_mum_cells) / full_nw_cells * 100,
        'nw_score': nw_score,
        'wfa_mum_cells': wfa_mum_cells,
        'wfa_full_cells': wfa_full.cells_visited,
        'wfa_mum_reduction_vs_nw': (full_nw_cells - wfa_mum_cells) / full_nw_cells * 100,
        'wfa_full_reduction_vs_nw': (full_nw_cells - wfa_full.cells_visited) / full_nw_cells * 100,
        'wfa_score': wfa_score,
        'wfa_vs_nw_mum_reduction': (nw_mum_cells - wfa_mum_cells) / nw_mum_cells * 100 if nw_mum_cells > 0 else 0,
        'wfa_vs_nw_full_reduction': (nw_full.cells_visited - wfa_full.cells_visited) / nw_full.cells_visited * 100 if nw_full.cells_visited > 0 else 0,
        'scores_match': nw_score == wfa_score,
    }
