"""Directed Acyclic Graph for MUM-based alignment."""

from typing import List, Dict, Tuple, Optional, Any
from dataclasses import dataclass

try:
    from .mum import MUM, TransitionMatrix, START, END, find_mums
    from .segment import align_gap, AlignmentResult
except ImportError:
    from mum import MUM, TransitionMatrix, START, END, find_mums
    from segment import align_gap, AlignmentResult


@dataclass
class Node:
    """Node in the MUM DAG."""
    id: Any
    mum: Optional[MUM] = None
    
    @property
    def is_start(self) -> bool:
        return self.id == START
    
    @property
    def is_end(self) -> bool:
        return self.id == END
    
    @property
    def is_mum(self) -> bool:
        return isinstance(self.id, int)
    
    def __hash__(self) -> int:
        return hash(self.id)
    
    def __eq__(self, other) -> bool:
        if not isinstance(other, Node):
            return False
        return self.id == other.id


@dataclass
class Edge:
    """Edge in the MUM DAG with alignment score."""
    source: Any
    target: Any
    weight: float = 0.0
    alignment: Optional[AlignmentResult] = None
    
    def __hash__(self) -> int:
        return hash((self.source, self.target))
    
    def __eq__(self, other) -> bool:
        if not isinstance(other, Edge):
            return False
        return self.source == other.source and self.target == other.target


class DAG:
    """
    DAG for MUM-based alignment.
    
    Nodes: START, MUM indices, END
    Edge weights: alignment scores between gaps
    Optimal path minimizes total score with MUM count as tiebreaker.
    """
    
    def __init__(self, transition_matrix: TransitionMatrix, seq1: str = None, seq2: str = None):
        self.tm = transition_matrix
        self.mums = transition_matrix.mums
        self.seq1 = seq1
        self.seq2 = seq2
        
        self._nodes: Dict[Any, Node] = {}
        self._nodes[START] = Node(id=START)
        self._nodes[END] = Node(id=END)
        for i, mum in enumerate(self.mums):
            self._nodes[i] = Node(id=i, mum=mum)
        
        self._edges: Dict[Tuple[Any, Any], Edge] = {}
        self._adjacency: Dict[Any, List[Any]] = {node_id: [] for node_id in self._nodes}
        self._reverse_adjacency: Dict[Any, List[Any]] = {node_id: [] for node_id in self._nodes}
        
        for src, targets in self.tm.successors.items():
            for tgt in targets:
                edge = Edge(source=src, target=tgt)
                self._edges[(src, tgt)] = edge
                self._adjacency[src].append(tgt)
                self._reverse_adjacency[tgt].append(src)
        
        self._weights_computed = False
    
    @property
    def num_nodes(self) -> int:
        return len(self._nodes)
    
    @property
    def num_edges(self) -> int:
        return len(self._edges)
    
    @property
    def num_mums(self) -> int:
        return len(self.mums)
    
    def get_edge(self, source: Any, target: Any) -> Optional[Edge]:
        return self._edges.get((source, target))
    
    def get_edge_weight(self, source: Any, target: Any) -> float:
        edge = self._edges.get((source, target))
        return edge.weight if edge else float('inf')
    
    def in_degree(self, node_id: Any) -> int:
        return len(self._reverse_adjacency[node_id])
    
    def compute_edge_weights(
        self,
        seq1: str = None,
        seq2: str = None,
        use_wfa: bool = False,
        mismatch: int = 4,
        gap_open: int = 6,
        gap_extend: int = 2,
        verbose: bool = False
    ) -> None:
        """Compute alignment scores for all edges."""
        seq1 = seq1 or self.seq1
        seq2 = seq2 or self.seq2
        
        if seq1 is None or seq2 is None:
            raise ValueError("Sequences required for edge weights")
        
        self.seq1 = seq1
        self.seq2 = seq2
        
        for idx, ((src, tgt), edge) in enumerate(self._edges.items()):
            if verbose and idx % 10 == 0:
                print(f"  Computing edge {idx+1}/{len(self._edges)}...")
            
            mum_from = self.mums[src] if isinstance(src, int) else src
            mum_to = self.mums[tgt] if isinstance(tgt, int) else tgt
            
            result = align_gap(
                seq1, seq2, mum_from, mum_to,
                use_wfa=use_wfa,
                mismatch=mismatch,
                gap_open=gap_open,
                gap_extend=gap_extend
            )
            
            edge.weight = result.score
            edge.alignment = result
        
        self._weights_computed = True
        if verbose:
            print(f"  Computed weights for {len(self._edges)} edges")
    
    def topological_order(self) -> List[Any]:
        """Return nodes in topological order using Kahn's algorithm."""
        in_deg = {node_id: self.in_degree(node_id) for node_id in self._nodes}
        queue = [n for n, d in in_deg.items() if d == 0]
        result = []
        
        while queue:
            node_id = queue.pop(0)
            result.append(node_id)
            for succ in self._adjacency[node_id]:
                in_deg[succ] -= 1
                if in_deg[succ] == 0:
                    queue.append(succ)
        
        return result
    
    def optimal_path(self) -> Tuple[List[Any], int, int]:
        """
        Find optimal path: minimize score, tiebreaker maximize MUMs.
        Returns (path, score, num_mums).
        """
        if not self._weights_computed:
            raise ValueError("Call compute_edge_weights() first")
        
        topo = self.topological_order()
        INF = float('inf')
        dist = {n: (INF, 0, None) for n in self._nodes}
        dist[START] = (0, 0, None)
        
        for node_id in topo:
            if dist[node_id][0] == INF:
                continue
            
            score, neg_mums, _ = dist[node_id]
            
            for succ in self._adjacency[node_id]:
                new_score = score + self.get_edge_weight(node_id, succ)
                mum_bonus = -1 if isinstance(succ, int) else 0
                new_neg_mums = neg_mums + mum_bonus
                
                if (new_score, new_neg_mums) < (dist[succ][0], dist[succ][1]):
                    dist[succ] = (new_score, new_neg_mums, node_id)
        
        path = []
        current = END
        while current is not None:
            path.append(current)
            current = dist[current][2]
        path.reverse()
        
        return path, int(dist[END][0]), -dist[END][1]
    
    def longest_path(self) -> Tuple[List[Any], int]:
        """Find path with most MUMs."""
        topo = self.topological_order()
        dist = {n: (0, None) for n in self._nodes}
        
        for node_id in topo:
            for succ in self._adjacency[node_id]:
                w = 1 if isinstance(succ, int) else 0
                if dist[node_id][0] + w > dist[succ][0]:
                    dist[succ] = (dist[node_id][0] + w, node_id)
        
        path = []
        current = END
        while current is not None:
            path.append(current)
            current = dist[current][1]
        path.reverse()
        
        return path, dist[END][0]
    
    def path_score(self, path: List[Any]) -> int:
        if not self._weights_computed:
            raise ValueError("Call compute_edge_weights() first")
        return int(sum(self.get_edge_weight(path[i], path[i+1]) for i in range(len(path)-1)))
    
    def get_mum_path(self, path: List[Any]) -> List[MUM]:
        return [self.mums[n] for n in path if isinstance(n, int)]
    
    def to_dot(
        self,
        highlight_path: Optional[List[Any]] = None,
        show_weights: bool = False,
        use_positions: bool = False,
        scale: float = 0.5
    ) -> str:
        """Generate DOT format for Graphviz visualization."""
        lines = ["digraph MUM_DAG {"]
        
        if use_positions:
            lines.append("  graph [splines=true];")
            lines.append("  node [pin=true];")
        else:
            lines.append("  rankdir=LR;")
        
        highlight_edges = set()
        if highlight_path:
            for i in range(len(highlight_path) - 1):
                highlight_edges.add((highlight_path[i], highlight_path[i + 1]))
        
        if use_positions and self.mums:
            max_p1 = max(m.pos1 for m in self.mums)
            max_p2 = max(m.pos2 for m in self.mums)
            end_x, end_y = (max_p1 + 10) * scale, (max_p2 + 10) * scale
            lines.append(f'  START [shape=circle, style=filled, fillcolor=green, pos="{-2*scale},{-2*scale}!"];')
            lines.append(f'  END [shape=circle, style=filled, fillcolor=red, pos="{end_x},{end_y}!"];')
        else:
            lines.append('  START [shape=circle, style=filled, fillcolor=green];')
            lines.append('  END [shape=circle, style=filled, fillcolor=red];')
        
        for i, mum in enumerate(self.mums):
            label = f"[{i}]\\n({mum.pos1},{mum.pos2})\\nlen={mum.length}"
            if use_positions:
                x, y = mum.pos1 * scale, mum.pos2 * scale
                lines.append(f'  {i} [shape=box, label="{label}", pos="{x},{y}!"];')
            else:
                lines.append(f'  {i} [shape=box, label="{label}"];')
        
        for (src, tgt), edge in self._edges.items():
            src_str = "START" if src == START else str(src)
            tgt_str = "END" if tgt == END else str(tgt)
            
            attrs = []
            if (src, tgt) in highlight_edges:
                attrs.extend(["color=red", "penwidth=2"])
            if show_weights and self._weights_computed:
                attrs.append(f'label="{int(edge.weight)}"')
            
            attr_str = f" [{', '.join(attrs)}]" if attrs else ""
            lines.append(f'  {src_str} -> {tgt_str}{attr_str};')
        
        lines.append("}")
        return "\n".join(lines)
    
    def __repr__(self) -> str:
        return f"DAG(nodes={self.num_nodes}, edges={self.num_edges}, mums={self.num_mums})"


def build_dag(
    seq1: str,
    seq2: str,
    min_length: int = 1,
    compute_weights: bool = False,
    use_wfa: bool = False,
    mismatch: int = 4,
    gap_open: int = 6,
    gap_extend: int = 2,
    verbose: bool = False
) -> DAG:
    """Build DAG from sequences with optional edge weight computation."""
    if verbose:
        print(f"Finding MUMs (min_length={min_length})...")
    
    mums = find_mums(seq1, seq2, min_length=min_length)
    
    if verbose:
        print(f"  Found {len(mums)} MUMs")
        print("Building transition matrix...")
    
    tm = TransitionMatrix(mums, len(seq1), len(seq2))
    dag = DAG(tm, seq1=seq1, seq2=seq2)
    
    if verbose:
        print(f"  DAG has {dag.num_edges} edges")
    
    if compute_weights:
        if verbose:
            print("Computing edge weights...")
        dag.compute_edge_weights(
            use_wfa=use_wfa,
            mismatch=mismatch,
            gap_open=gap_open,
            gap_extend=gap_extend,
            verbose=verbose
        )
    
    return dag
