"""Needleman-Wunsch global alignment with affine gap penalties."""

from typing import Tuple


def needleman_wunsch(seq1: str, seq2: str, match: int = 0, mismatch: int = 4, gap_open: int = 6, gap_extend: int = 2, return_cells: bool = False) -> Tuple[int, str]:
    """Global alignment with affine gaps."""
    n, m = len(seq1), len(seq2)
    
    if n == 0 and m == 0:
        return (0, "", 0) if return_cells else (0, "")
    if n == 0:
        result = (gap_open + m * gap_extend, f"{m}I")
        return (*result, m) if return_cells else result
    if m == 0:
        result = (gap_open + n * gap_extend, f"{n}D")
        return (*result, n) if return_cells else result
    
    cells_visited = n * m
    INF = float('inf')
    
    M = [[INF] * (m + 1) for _ in range(n + 1)]
    I = [[INF] * (m + 1) for _ in range(n + 1)]
    D = [[INF] * (m + 1) for _ in range(n + 1)]
    
    M_ptr = [[None] * (m + 1) for _ in range(n + 1)]
    I_ptr = [[None] * (m + 1) for _ in range(n + 1)]
    D_ptr = [[None] * (m + 1) for _ in range(n + 1)]
    
    M[0][0] = 0
    
    for j in range(1, m + 1):
        I[0][j] = gap_open + j * gap_extend
        I_ptr[0][j] = ('I', 0, j - 1)
    
    for i in range(1, n + 1):
        D[i][0] = gap_open + i * gap_extend
        D_ptr[i][0] = ('D', i - 1, 0)
    
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            cost = match if seq1[i-1] == seq2[j-1] else mismatch
            
            cands_M = []
            if M[i-1][j-1] < INF:
                cands_M.append((M[i-1][j-1] + cost, ('M', i-1, j-1)))
            if I[i-1][j-1] < INF:
                cands_M.append((I[i-1][j-1] + cost, ('I', i-1, j-1)))
            if D[i-1][j-1] < INF:
                cands_M.append((D[i-1][j-1] + cost, ('D', i-1, j-1)))
            if cands_M:
                best = min(cands_M, key=lambda x: x[0])
                M[i][j], M_ptr[i][j] = best[0], best[1]
            
            cands_I = []
            if M[i][j-1] < INF:
                cands_I.append((M[i][j-1] + gap_open + gap_extend, ('M', i, j-1)))
            if I[i][j-1] < INF:
                cands_I.append((I[i][j-1] + gap_extend, ('I', i, j-1)))
            if cands_I:
                best = min(cands_I, key=lambda x: x[0])
                I[i][j], I_ptr[i][j] = best[0], best[1]
            
            cands_D = []
            if M[i-1][j] < INF:
                cands_D.append((M[i-1][j] + gap_open + gap_extend, ('M', i-1, j)))
            if D[i-1][j] < INF:
                cands_D.append((D[i-1][j] + gap_extend, ('D', i-1, j)))
            if cands_D:
                best = min(cands_D, key=lambda x: x[0])
                D[i][j], D_ptr[i][j] = best[0], best[1]
    
    final = [(M[n][m], 'M'), (I[n][m], 'I'), (D[n][m], 'D')]
    best_score, state = min(final, key=lambda x: x[0])
    
    path = []
    i, j = n, m
    
    while i > 0 or j > 0:
        if state == 'M':
            if M_ptr[i][j] is None:
                break
            prev, pi, pj = M_ptr[i][j]
            path.append('M' if seq1[i-1] == seq2[j-1] else 'X')
            state, i, j = prev, pi, pj
        elif state == 'I':
            if I_ptr[i][j] is None:
                break
            prev, pi, pj = I_ptr[i][j]
            path.append('I')
            state, i, j = prev, pi, pj
        elif state == 'D':
            if D_ptr[i][j] is None:
                break
            prev, pi, pj = D_ptr[i][j]
            path.append('D')
            state, i, j = prev, pi, pj
        else:
            break
    
    path.reverse()
    cigar = ''.join(path)
    
    if return_cells:
        return int(best_score), cigar, cells_visited
    return int(best_score), cigar

