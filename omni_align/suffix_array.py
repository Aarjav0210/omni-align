"""Suffix Array construction and LCP array utilities."""

from typing import List


def build_suffix_array(text: str) -> List[int]:
    """Build suffix array using prefix doubling. O(n log n) time."""
    n = len(text)
    if n == 0:
        return []
    if n == 1:
        return [0]
    
    sa = list(range(n))
    rank = [ord(c) for c in text]
    tmp = [0] * n
    
    k = 1
    while k < n:
        def sort_key(i):
            return (rank[i], rank[i + k] if i + k < n else -1)
        
        sa.sort(key=sort_key)
        
        tmp[sa[0]] = 0
        for i in range(1, n):
            prev_key = (rank[sa[i-1]], rank[sa[i-1] + k] if sa[i-1] + k < n else -1)
            curr_key = (rank[sa[i]], rank[sa[i] + k] if sa[i] + k < n else -1)
            tmp[sa[i]] = tmp[sa[i-1]] + (1 if curr_key != prev_key else 0)
        
        rank, tmp = tmp, rank
        
        if rank[sa[n-1]] == n - 1:
            break
            
        k *= 2
    
    return sa


def build_lcp_array(text: str, sa: List[int]) -> List[int]:
    """Build LCP array using Kasai's algorithm. O(n) time."""
    n = len(text)
    if n == 0:
        return []
    
    rank = [0] * n
    for i in range(n):
        rank[sa[i]] = i
    
    lcp = [0] * n
    k = 0
    
    for i in range(n):
        if rank[i] == 0:
            k = 0
            continue
            
        j = sa[rank[i] - 1]
        
        while i + k < n and j + k < n and text[i + k] == text[j + k]:
            k += 1
            
        lcp[rank[i]] = k
        
        if k > 0:
            k -= 1
    
    return lcp


class SuffixArray:
    """Suffix Array index for efficient substring operations."""
    
    def __init__(self, text: str, add_sentinel: bool = True):
        self._original_text = text
        
        if add_sentinel and not text.endswith('$'):
            self.text = text + '$'
        else:
            self.text = text
            
        self.sa = build_suffix_array(self.text)
        self.lcp = build_lcp_array(self.text, self.sa)
        
        self._rank = [0] * len(self.text)
        for i, s in enumerate(self.sa):
            self._rank[s] = i
    
    def __len__(self) -> int:
        return len(self.text)
    
    def suffix(self, sa_index: int) -> str:
        return self.text[self.sa[sa_index]:]
    
    def rank(self, text_pos: int) -> int:
        return self._rank[text_pos]
    
    def find_all(self, pattern: str) -> List[int]:
        """Find all occurrences of pattern. O(m log n) time."""
        if not pattern:
            return []
        
        n = len(self.text)
        m = len(pattern)
        
        left = 0
        right = n
        while left < right:
            mid = (left + right) // 2
            suffix = self.text[self.sa[mid]:self.sa[mid] + m]
            if suffix < pattern:
                left = mid + 1
            else:
                right = mid
        start = left
        
        right = n
        while left < right:
            mid = (left + right) // 2
            suffix = self.text[self.sa[mid]:self.sa[mid] + m]
            if suffix <= pattern:
                left = mid + 1
            else:
                right = mid
        end = left
        
        return [self.sa[i] for i in range(start, end)]
    
    def count(self, pattern: str) -> int:
        return len(self.find_all(pattern))
