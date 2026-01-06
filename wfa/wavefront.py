"""Wavefront Alignment (WFA) with affine gap penalties"""

from typing import Tuple, Dict, List


class WavefrontAligner:
    """
    WFA aligner
    
    Offset i = position in seq1
    j = i + k = position in seq2
    """
    
    def __init__(self, pattern: str, text: str, mismatch: int = 4, gap_open: int = 6, gap_extend: int = 2):
        self.P = pattern
        self.T = text
        self.N = len(pattern)
        self.M = len(text)
        self.X = mismatch
        self.O = gap_open
        self.E = gap_extend
        
        self.w_M: List[Dict[int, int]] = []
        self.w_I: List[Dict[int, int]] = []
        self.w_D: List[Dict[int, int]] = []
        
        self.positions_visited = set()
        self.cells_visited = 0
    
    def align(self) -> Tuple[int, str, int]:
        """Returns (score, cigar, cells_visited)"""
        if self.N == 0 and self.M == 0:
            return 0, "", 0
        if self.N == 0:
            return self.O + self.M * self.E, f"{self.M}I", self.M
        if self.M == 0:
            return self.O + self.N * self.E, f"{self.N}D", self.N
        
        self.w_M.append({})
        self.w_I.append({})
        self.w_D.append({})
        
        k = 0
        initial_offset = self._extend(0, k)
        self.w_M[0][k] = initial_offset
        
        if self._finished(initial_offset, k):
            self.cells_visited = len(self.positions_visited)
            return 0, self._backtrack(0, k), self.cells_visited
        
        score = 1
        max_score = (self.N + self.M) * max(self.X, self.O + self.E)
        
        while score <= max_score:
            self.w_M.append({})
            self.w_I.append({})
            self.w_D.append({})
            
            active_k = self._get_active_diagonals(score)
            if not active_k:
                score += 1
                continue
            
            for k in sorted(active_k):
                val_I = self._compute_I(score, k)
                if val_I >= 0:
                    self.w_I[score][k] = val_I
                
                val_D = self._compute_D(score, k)
                if val_D >= 0:
                    self.w_D[score][k] = val_D
                
                val_M = self._compute_M(score, k)
                if val_M >= 0:
                    extended = self._extend(val_M, k)
                    self.w_M[score][k] = extended
                    
                    if self._finished(extended, k):
                        self.cells_visited = len(self.positions_visited)
                        return score, self._backtrack(score, k), self.cells_visited
            
            score += 1
        
        raise Exception(f"WFA failed within score limit {max_score}")
    
    def _get_active_diagonals(self, score: int) -> set:
        active = set()
        
        for s_src in [score - (self.O + self.E), score - self.E]:
            if 0 <= s_src < len(self.w_M):
                for k_prev in self.w_M[s_src]:
                    active.add(k_prev + 1)
                for k_prev in self.w_I[s_src]:
                    active.add(k_prev + 1)
        
        for s_src in [score - (self.O + self.E), score - self.E]:
            if 0 <= s_src < len(self.w_M):
                for k_prev in self.w_M[s_src]:
                    active.add(k_prev - 1)
                for k_prev in self.w_D[s_src]:
                    active.add(k_prev - 1)
        
        s_src = score - self.X
        if 0 <= s_src < len(self.w_M):
            for k_prev in self.w_M[s_src]:
                active.add(k_prev)
        
        return active
    
    def _compute_I(self, score: int, k: int) -> int:
        val = -1
        s_open = score - (self.O + self.E)
        if s_open >= 0:
            prev = self._get(self.w_M, s_open, k - 1)
            if prev >= 0:
                val = max(val, prev)
        s_ext = score - self.E
        if s_ext >= 0:
            prev = self._get(self.w_I, s_ext, k - 1)
            if prev >= 0:
                val = max(val, prev)
        return val
    
    def _compute_D(self, score: int, k: int) -> int:
        val = -1
        s_open = score - (self.O + self.E)
        if s_open >= 0:
            prev = self._get(self.w_M, s_open, k + 1)
            if prev >= 0 and prev + 1 <= self.N:
                val = max(val, prev + 1)
        s_ext = score - self.E
        if s_ext >= 0:
            prev = self._get(self.w_D, s_ext, k + 1)
            if prev >= 0 and prev + 1 <= self.N:
                val = max(val, prev + 1)
        return val
    
    def _compute_M(self, score: int, k: int) -> int:
        val = -1
        s_mis = score - self.X
        if s_mis >= 0:
            prev = self._get(self.w_M, s_mis, k)
            if prev >= 0:
                new_i, new_j = prev + 1, prev + 1 + k
                if new_i <= self.N and new_j <= self.M:
                    val = max(val, new_i)
        if k in self.w_I[score]:
            val = max(val, self.w_I[score][k])
        if k in self.w_D[score]:
            val = max(val, self.w_D[score][k])
        return val
    
    def _extend(self, offset: int, k: int) -> int:
        i, j = offset, offset + k
        while i < self.N and j < self.M:
            self.positions_visited.add((i, j))
            if self.P[i] == self.T[j]:
                i += 1
                j += 1
            else:
                break
        return i
    
    def _finished(self, offset: int, k: int) -> bool:
        return offset >= self.N and (offset + k) >= self.M
    
    def _get(self, w: List[Dict[int, int]], s: int, k: int) -> int:
        if 0 <= s < len(w) and k in w[s]:
            return w[s][k]
        return -1
    
    def _backtrack(self, final_score: int, final_k: int) -> str:
        path = []
        s, k = final_score, final_k
        offset = self.w_M[s][k]
        state = 'M'
        
        while s > 0 or offset > 0:
            if s == 0 and offset == 0:
                break
            
            if state == 'M':
                if s == 0:
                    for _ in range(offset):
                        path.append('M')
                    break
                
                cand_mis = -1
                s_mis = s - self.X
                if s_mis >= 0:
                    prev = self._get(self.w_M, s_mis, k)
                    if prev >= 0:
                        cand_mis = prev + 1
                
                cand_I = self._get(self.w_I, s, k)
                cand_D = self._get(self.w_D, s, k)
                m_init = max(cand_mis, cand_I, cand_D)
                
                for _ in range(offset - m_init):
                    path.append('M')
                offset = m_init
                
                if m_init == cand_mis and cand_mis >= 0:
                    path.append('X')
                    s -= self.X
                    offset -= 1
                elif m_init == cand_I and cand_I >= 0:
                    state = 'I'
                elif m_init == cand_D and cand_D >= 0:
                    state = 'D'
                else:
                    raise Exception(f"Backtrack error M: s={s}, k={k}")
            
            elif state == 'I':
                s_open = s - (self.O + self.E)
                s_ext = s - self.E
                cand_open = self._get(self.w_M, s_open, k - 1) if s_open >= 0 else -1
                cand_ext = self._get(self.w_I, s_ext, k - 1) if s_ext >= 0 else -1
                
                if cand_ext == offset:
                    path.append('I')
                    s, k, state = s_ext, k - 1, 'I'
                elif cand_open == offset:
                    path.append('I')
                    s, k, state = s_open, k - 1, 'M'
                else:
                    raise Exception(f"Backtrack error I: s={s}, k={k}")
            
            elif state == 'D':
                src_offset = offset - 1
                s_open = s - (self.O + self.E)
                s_ext = s - self.E
                cand_open = self._get(self.w_M, s_open, k + 1) if s_open >= 0 else -1
                cand_ext = self._get(self.w_D, s_ext, k + 1) if s_ext >= 0 else -1
                
                if cand_ext == src_offset:
                    path.append('D')
                    s, k, offset, state = s_ext, k + 1, offset - 1, 'D'
                elif cand_open == src_offset:
                    path.append('D')
                    s, k, offset, state = s_open, k + 1, offset - 1, 'M'
                else:
                    raise Exception(f"Backtrack error D: s={s}, k={k}")
        
        return "".join(reversed(path))


def wfa_align(seq1: str, seq2: str, mismatch: int = 4, gap_open: int = 6, gap_extend: int = 2, return_cells: bool = False) -> Tuple[int, str]:
    """
    WFA alignment
    
    Returns
    - (score, cigar) if return_cells == False
    - (score, cigar, cells) if return_cells == True
    """
    aligner = WavefrontAligner(seq1, seq2, mismatch, gap_open, gap_extend)
    score, cigar, cells = aligner.align()
    if return_cells:
        return score, cigar, cells
    return score, cigar

