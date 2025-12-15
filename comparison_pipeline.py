#!/usr/bin/env python3
"""Comparison pipeline for alignment methods."""

import pandas as pd
from omni_align import find_mums, build_dag, align_path
from nw import needleman_wunsch
from wfa import wfa_align
from seed_chain_extend import align_with_seeds

pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)

SEQ1 = 'TGATGAATCCGAGATTTTCTACGCCGGTTGAATTCGGAAG'
SEQ2 = 'CTCCGCATGATACCGTACGCAGCCGACATTTCCAGCGTTG'


def run_pipeline(seq1: str, seq2: str, min_mum_length: int = 4, sce_k: int = 4):
    """Run all alignment methods and return comparison results."""
    full_cells = len(seq1) * len(seq2)
    
    nw_score, _ = needleman_wunsch(seq1, seq2)
    wfa_score, _, wfa_cells = wfa_align(seq1, seq2, return_cells=True)
    
    dag = build_dag(seq1, seq2, min_length=min_mum_length, compute_weights=True)
    path, omni_score, num_mums = dag.optimal_path()
    mum_path = dag.get_mum_path(path)
    
    _, _, _, omni_nw_cells = align_path(seq1, seq2, mum_path, use_wfa=False)
    _, _, _, omni_wfa_cells = align_path(seq1, seq2, mum_path, use_wfa=True)
    
    sce_nw_result, seeds_found, seeds_in_chain = align_with_seeds(seq1, seq2, k=sce_k, use_wfa=False)
    sce_wfa_result, _, _ = align_with_seeds(seq1, seq2, k=sce_k, use_wfa=True)
    
    results = [
        {
            'Method': 'NW',
            'Score': nw_score,
            'Cells': full_cells,
            'Reduction': 0.0,
            'Optimal': 'Y'
        },
        {
            'Method': 'WFA',
            'Score': wfa_score,
            'Cells': wfa_cells,
            'Reduction': round((1 - wfa_cells / full_cells) * 100, 1),
            'Optimal': 'Y'
        },
        {
            'Method': 'omniNW',
            'Score': omni_score,
            'Cells': omni_nw_cells,
            'Reduction': round((1 - omni_nw_cells / full_cells) * 100, 1),
            'Optimal': 'Y' if omni_score == nw_score else 'N'
        },
        {
            'Method': 'omniWFA',
            'Score': omni_score,
            'Cells': omni_wfa_cells,
            'Reduction': round((1 - omni_wfa_cells / full_cells) * 100, 1),
            'Optimal': 'Y' if omni_score == nw_score else 'N'
        },
        {
            'Method': 'SCE-NW',
            'Score': sce_nw_result.score,
            'Cells': sce_nw_result.cells_visited,
            'Reduction': round((1 - sce_nw_result.cells_visited / full_cells) * 100, 1),
            'Optimal': 'Y' if sce_nw_result.score == nw_score else 'N'
        },
        {
            'Method': 'SCE-WFA',
            'Score': sce_wfa_result.score,
            'Cells': sce_wfa_result.cells_visited,
            'Reduction': round((1 - sce_wfa_result.cells_visited / full_cells) * 100, 1),
            'Optimal': 'Y' if sce_wfa_result.score == nw_score else 'N'
        }
    ]
    
    return pd.DataFrame(results)


def get_mum_details(seq1: str, seq2: str, min_length: int = 4):
    """Get details of MUMs found between sequences."""
    mums = find_mums(seq1, seq2, min_length=min_length)
    data = []
    for i, m in enumerate(mums):
        data.append({
            'ID': i,
            'Sequence': seq1[m.pos1:m.end1],
            'Length': m.length,
            'Pos1': m.pos1,
            'Pos2': m.pos2
        })
    return pd.DataFrame(data)


if __name__ == '__main__':
    print(f'Seq1: {SEQ1}')
    print(f'Seq2: {SEQ2}')
    print(f'Size: {len(SEQ1)} x {len(SEQ2)} = {len(SEQ1)*len(SEQ2)} cells')
    print()
    
    print('MUMs:')
    print(get_mum_details(SEQ1, SEQ2).to_string(index=False))
    print()
    
    print('Results:')
    print(run_pipeline(SEQ1, SEQ2).to_string(index=False))

