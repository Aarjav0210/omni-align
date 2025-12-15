"""Co-linear seed chaining using dynamic programming."""

from dataclasses import dataclass, field
from typing import List
from .seeding import Seed


@dataclass
class Chain:
    """A chain of co-linear seeds."""
    seeds: List[Seed] = field(default_factory=list)
    score: float = 0.0
    
    @property
    def coverage(self) -> int:
        return sum(s.length for s in self.seeds)


def can_chain(seed_i: Seed, seed_j: Seed) -> bool:
    """Check if seed_j can follow seed_i (non-overlapping, co-linear)."""
    return seed_i.end1 <= seed_j.pos1 and seed_i.end2 <= seed_j.pos2


def chain_seeds(seeds: List[Seed], gap_penalty: float = 0.1) -> List[Chain]:
    """Chain seeds using DP. Returns list with best chain."""
    if not seeds:
        return []
    
    sorted_seeds = sorted(seeds, key=lambda s: (s.pos1, s.pos2))
    n = len(sorted_seeds)
    
    dp = [float(s.length) for s in sorted_seeds]
    pred = [-1] * n
    
    for j in range(1, n):
        for i in range(j):
            if can_chain(sorted_seeds[i], sorted_seeds[j]):
                gap1 = sorted_seeds[j].pos1 - sorted_seeds[i].end1
                gap2 = sorted_seeds[j].pos2 - sorted_seeds[i].end2
                new_score = dp[i] + sorted_seeds[j].length - gap_penalty * (gap1 + gap2)
                if new_score > dp[j]:
                    dp[j] = new_score
                    pred[j] = i
    
    best_end = max(range(n), key=lambda i: dp[i])
    
    chain_indices = []
    idx = best_end
    while idx >= 0:
        chain_indices.append(idx)
        idx = pred[idx]
    chain_indices.reverse()
    
    return [Chain(seeds=[sorted_seeds[i] for i in chain_indices], score=dp[best_end])]
