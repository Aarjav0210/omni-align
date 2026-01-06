"""Chain extension using NW or WFA for gap alignment"""

from dataclasses import dataclass
from typing import List, Tuple

from .seeding import Seed, find_kmer_seeds, filter_unique_seeds
from .chaining import Chain, chain_seeds
from nw import needleman_wunsch
from wfa import wfa_align


@dataclass
class ExtensionResult:
    """Result of chain extension"""
    score: int
    cigar: str
    cells_visited: int
    num_seeds: int
    seed_coverage: int


def extend_chain(seq1: str, seq2: str, chain: Chain, use_wfa: bool = False) -> ExtensionResult:
    """Extend chain by aligning gaps between seeds"""
    if not chain.seeds:
        if use_wfa:
            score, cigar, cells = wfa_align(seq1, seq2, return_cells=True)
        else:
            score, cigar, cells = needleman_wunsch(seq1, seq2, return_cells=True)
        return ExtensionResult(score=score, cigar=cigar, cells_visited=cells, num_seeds=0, seed_coverage=0)
    
    total_score, total_cells = 0, 0
    cigar_parts = []
    seeds = chain.seeds
    
    def align_gap(s1: str, s2: str) -> Tuple[int, str, int]:
        if not s1 and not s2:
            return 0, "", 0
        if not s1:
            return 6 + len(s2) * 2, f"{len(s2)}I", len(s2)
        if not s2:
            return 6 + len(s1) * 2, f"{len(s1)}D", len(s1)
        if use_wfa:
            return wfa_align(s1, s2, return_cells=True)
        return needleman_wunsch(s1, s2, return_cells=True)
    
    if seeds[0].pos1 > 0 or seeds[0].pos2 > 0:
        score, cigar, cells = align_gap(seq1[:seeds[0].pos1], seq2[:seeds[0].pos2])
        total_score += score
        total_cells += cells
        if cigar:
            cigar_parts.append(cigar)
    
    cigar_parts.append(f"{seeds[0].length}M")
    
    for i in range(len(seeds) - 1):
        score, cigar, cells = align_gap(seq1[seeds[i].end1:seeds[i+1].pos1], seq2[seeds[i].end2:seeds[i+1].pos2])
        total_score += score
        total_cells += cells
        if cigar:
            cigar_parts.append(cigar)
        cigar_parts.append(f"{seeds[i+1].length}M")
    
    if seeds[-1].end1 < len(seq1) or seeds[-1].end2 < len(seq2):
        score, cigar, cells = align_gap(seq1[seeds[-1].end1:], seq2[seeds[-1].end2:])
        total_score += score
        total_cells += cells
        if cigar:
            cigar_parts.append(cigar)
    
    return ExtensionResult(
        score=total_score,
        cigar="".join(cigar_parts),
        cells_visited=total_cells,
        num_seeds=len(seeds),
        seed_coverage=sum(s.length for s in seeds)
    )


def align_with_seeds(seq1: str, seq2: str, k: int = 4, use_wfa: bool = False, unique_only: bool = True) -> Tuple[ExtensionResult, int, int]:
    """
    Full seed-chain-extend pipeline
    
    Returns (result, seeds_found, seeds_in_chain)
    """
    
    seeds = find_kmer_seeds(seq1, seq2, k=k)
    if unique_only:
        seeds = filter_unique_seeds(seeds, seq1, seq2)
    
    chains = chain_seeds(seeds)
    best_chain = chains[0] if chains else Chain()
    
    result = extend_chain(seq1, seq2, best_chain, use_wfa=use_wfa)
    return result, len(seeds), len(best_chain.seeds)
