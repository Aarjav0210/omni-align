"""K-mer and minimizer seeding strategies"""

from dataclasses import dataclass
from typing import List, Dict, Set, Tuple
from collections import defaultdict


@dataclass
class Seed:
    """A seed match between two sequences"""
    pos1: int
    pos2: int
    length: int
    
    @property
    def end1(self) -> int:
        return self.pos1 + self.length
    
    @property
    def end2(self) -> int:
        return self.pos2 + self.length


def find_kmer_seeds(seq1: str, seq2: str, k: int = 11) -> List[Seed]:
    """Find all k-mer matches between sequences"""
    if k > len(seq1) or k > len(seq2):
        return []
    
    kmer_index: Dict[str, List[int]] = defaultdict(list)
    for i in range(len(seq2) - k + 1):
        kmer_index[seq2[i:i+k]].append(i)
    
    seeds = []
    for i in range(len(seq1) - k + 1):
        kmer = seq1[i:i+k]
        if kmer in kmer_index:
            for j in kmer_index[kmer]:
                seeds.append(Seed(pos1=i, pos2=j, length=k))
    
    return seeds


def find_minimizer_seeds(seq1: str, seq2: str, k: int = 11, w: int = 5) -> List[Seed]:
    """Find minimizer-based seed matches"""
    if k > len(seq1) or k > len(seq2):
        return []
    
    def get_minimizers(seq: str) -> Set[Tuple[str, int]]:
        minimizers = set()
        if len(seq) < k + w - 1:
            for i in range(len(seq) - k + 1):
                minimizers.add((seq[i:i+k], i))
            return minimizers
        
        for window_start in range(len(seq) - k - w + 2):
            min_kmer, min_pos = None, None
            for i in range(window_start, window_start + w):
                if i + k > len(seq):
                    break
                kmer = seq[i:i+k]
                if min_kmer is None or kmer < min_kmer:
                    min_kmer, min_pos = kmer, i
            if min_kmer:
                minimizers.add((min_kmer, min_pos))
        return minimizers
    
    min1 = get_minimizers(seq1)
    min2 = get_minimizers(seq2)
    
    kmer_to_pos2: Dict[str, List[int]] = defaultdict(list)
    for kmer, pos in min2:
        kmer_to_pos2[kmer].append(pos)
    
    seeds = []
    for kmer, pos1 in min1:
        if kmer in kmer_to_pos2:
            for pos2 in kmer_to_pos2[kmer]:
                seeds.append(Seed(pos1=pos1, pos2=pos2, length=k))
    
    return seeds


def filter_unique_seeds(seeds: List[Seed], seq1: str, seq2: str) -> List[Seed]:
    """Keep only seeds that are unique in both sequences"""
    if not seeds:
        return []
    
    k = seeds[0].length
    count1: Dict[str, int] = defaultdict(int)
    count2: Dict[str, int] = defaultdict(int)
    
    for i in range(len(seq1) - k + 1):
        count1[seq1[i:i+k]] += 1
    for i in range(len(seq2) - k + 1):
        count2[seq2[i:i+k]] += 1
    
    return [s for s in seeds if count1[seq1[s.pos1:s.pos1+k]] == 1 and count2[seq2[s.pos2:s.pos2+k]] == 1]
