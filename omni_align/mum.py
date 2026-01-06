"""MUM (Maximal Unique Match) extraction and transition utilities"""

from typing import List, Tuple, Dict
from dataclasses import dataclass

try:
    from .suffix_array import build_suffix_array, build_lcp_array
except ImportError:
    from suffix_array import build_suffix_array, build_lcp_array


START = "START"
END = "END"


@dataclass
class MUM:
    """Maximal Unique Match between two sequences"""
    pos1: int
    pos2: int
    length: int
    
    @property
    def end1(self) -> int:
        return self.pos1 + self.length
    
    @property
    def end2(self) -> int:
        return self.pos2 + self.length
    
    def __repr__(self) -> str:
        return f"MUM(pos1={self.pos1}, pos2={self.pos2}, len={self.length})"
    
    def __hash__(self) -> int:
        return hash((self.pos1, self.pos2, self.length))
    
    def __eq__(self, other) -> bool:
        if not isinstance(other, MUM):
            return False
        return self.pos1 == other.pos1 and self.pos2 == other.pos2 and self.length == other.length


def find_mums(seq1: str, seq2: str, min_length: int = 1) -> List[MUM]:
    """
    Find all Maximal Unique Matches between two sequences
    
    A MUM occurs exactly once in each sequence and cannot be extended
    Uses suffix array with LCP
    """
    if not seq1 or not seq2:
        return []
    
    sep1, sep2 = '$', '#'
    combined = seq1 + sep1 + seq2 + sep2
    n1, n2 = len(seq1), len(seq2)
    n = len(combined)
    
    sa = build_suffix_array(combined)
    lcp = build_lcp_array(combined, sa)
    
    def which_seq(pos: int) -> int:
        if pos < n1:
            return 1
        elif pos == n1:
            return 0
        elif pos < n1 + 1 + n2:
            return 2
        return 0
    
    def pos_in_seq(pos: int) -> int:
        if pos < n1:
            return pos
        elif pos > n1:
            return pos - n1 - 1
        return -1
    
    mums = []
    
    for i in range(1, n):
        pos_i, pos_prev = sa[i], sa[i-1]
        seq_i, seq_prev = which_seq(pos_i), which_seq(pos_prev)
        
        if seq_i == seq_prev or seq_i == 0 or seq_prev == 0:
            continue
        
        match_len = lcp[i]
        if match_len < min_length:
            continue
        
        lcp_before = lcp[i-1] if i > 1 else 0
        lcp_after = lcp[i+1] if i + 1 < n else 0
        
        if lcp_before >= match_len or lcp_after >= match_len:
            continue
        
        if seq_i == 1:
            p1, p2 = pos_in_seq(pos_i), pos_in_seq(pos_prev)
        else:
            p1, p2 = pos_in_seq(pos_prev), pos_in_seq(pos_i)
        
        if p1 + match_len > n1 or p2 + match_len > n2:
            continue
        
        if p1 > 0 and p2 > 0 and seq1[p1 - 1] == seq2[p2 - 1]:
            continue
        
        end1, end2 = p1 + match_len, p2 + match_len
        if end1 < n1 and end2 < n2 and seq1[end1] == seq2[end2]:
            continue
        
        mums.append(MUM(pos1=p1, pos2=p2, length=match_len))
    
    seen = set()
    unique_mums = []
    for mum in mums:
        key = (mum.pos1, mum.pos2, mum.length)
        if key not in seen:
            seen.add(key)
            unique_mums.append(mum)
    
    unique_mums.sort(key=lambda m: (m.pos1, m.pos2))
    return unique_mums


def can_transition(mum_i: MUM, mum_j: MUM) -> bool:
    """Check if transition i -> j is valid (ordered and non-overlapping)"""
    if not (mum_i.pos1 < mum_j.pos1 and mum_i.pos2 < mum_j.pos2):
        return False
    if mum_i.end1 > mum_j.pos1 or mum_i.end2 > mum_j.pos2:
        return False
    return True


class TransitionMatrix:
    """MUM transition graph with START and END nodes"""
    
    def __init__(self, mums: List[MUM], seq1_len: int = 0, seq2_len: int = 0):
        self.mums = mums
        self.n = len(mums)
        self.seq1_len = seq1_len
        self.seq2_len = seq2_len
        self._successors = None
        self._predecessors = None
    
    @property
    def successors(self) -> Dict:
        if self._successors is None:
            adj = {i: [] for i in range(self.n)}
            for i in range(self.n):
                for j in range(self.n):
                    if i != j and can_transition(self.mums[i], self.mums[j]):
                        adj[i].append(j)
                adj[i].append(END)
            adj[START] = list(range(self.n)) + [END]
            adj[END] = []
            self._successors = adj
        return self._successors
    
    @property
    def predecessors(self) -> Dict:
        if self._predecessors is None:
            pred = {j: [] for j in range(self.n)}
            for i in range(self.n):
                for j in range(self.n):
                    if i != j and can_transition(self.mums[i], self.mums[j]):
                        pred[j].append(i)
                pred[i].insert(0, START)
            pred[START] = []
            pred[END] = [START] + list(range(self.n))
            self._predecessors = pred
        return self._predecessors
    
    def can_go(self, i, j) -> bool:
        if i == START:
            return (isinstance(j, int) and 0 <= j < self.n) or j == END
        if i == END:
            return False
        if j == END:
            return isinstance(i, int) and 0 <= i < self.n
        if j == START:
            return False
        return can_transition(self.mums[i], self.mums[j])
    
    def get_edges(self) -> List[Tuple]:
        edges = []
        for i, targets in self.successors.items():
            for j in targets:
                edges.append((i, j))
        return edges
    
    def __len__(self) -> int:
        return self.n
    
    def __repr__(self) -> str:
        return f"TransitionMatrix(n_mums={self.n}, edges={len(self.get_edges())})"


def build_transition_matrix(seq1: str, seq2: str, min_length: int = 1) -> TransitionMatrix:
    """Find MUMs and build transition matrix"""
    mums = find_mums(seq1, seq2, min_length=min_length)
    return TransitionMatrix(mums, len(seq1), len(seq2))
