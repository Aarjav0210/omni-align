Sequence alignment using Maximal Unique Matches (MUMs) as landmarks to reduce search space.
Compares omniAlign (MUM-based), Needleman-Wunsch, Wavefront Alignment, and Seed-Chain-Extend.
omniAlign achieves optimal alignment scores with 93%+ reduction in cells visited.

1. Clone the repository
2. Create virtual environment: python3 -m venv venv
3. Activate virtual environment: source venv/bin/activate
4. Install dependencies: pip install -r requirements.txt
5. Run comparison: python comparison_pipeline.py