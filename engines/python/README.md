# biodp (Python engines)

Dynamic programming engines for the Bioinformatics DP Engine & Visualizer. Algorithm write-ups live in [`docs/algorithms/`](../../docs/algorithms/).

## Install and test

```bash
cd engines/python
pip install -e ".[dev]"
pytest            # 186 tests: golden fixtures, exhaustive oracles, Biopython cross-checks
mypy              # strict mode over src/
```

## Usage

```python
from biodp.align import needleman_wunsch
from biodp.scoring import GapPenalty, SubstitutionMatrix
from biodp.trace import dumps, to_dptrace

blosum62 = SubstitutionMatrix.from_ncbi_file("../../fixtures/matrices/BLOSUM62.txt")
result = needleman_wunsch("HEAGAWGHEE", "PAWHEAE", matrix=blosum62, gap=GapPenalty(-11, -1))

result.score             # optimal score
result.optimal_count     # exact number of co-optimal alignments (arbitrary precision)
print(result.best)       # canonical alignment, gaps left-aligned
result.best.cigar()      # extended CIGAR, seq_a as reference
result.alignments.take(10)      # first ten, generated lazily
result.alignments.sample()      # one drawn uniformly at random
trace_json = dumps(to_dptrace(result))   # DPTrace for the visualiser
```

For sequences too long for the full table (linear gaps only):

```python
from biodp.align import global_score, hirschberg

alignment = hirschberg(long_a, long_b, matrix=blosum62, gap=GapPenalty.linear(-4))   # O(n + m) memory
score = global_score(long_a, long_b, matrix=blosum62, gap=GapPenalty.linear(-4))     # score only
```

RNA secondary structure (Nussinov, maximum base pairs):

```python
from biodp.rna import FoldingRules, nussinov

result = nussinov("GGGAAAUCC", rules=FoldingRules(min_loop=3, wobble=True))
result.score                 # 3 base pairs
result.optimal_count         # 1
result.best().dot_bracket    # '(((...)))'
```

## Regenerating golden fixtures

```bash
python ../../scripts/generate_nw_golden.py
python ../../scripts/generate_nussinov_golden.py
```

Expected values come from exhaustive enumeration, never from the engine, and are cross-checked against Biopython.
