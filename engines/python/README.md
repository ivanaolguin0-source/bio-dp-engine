# biodp (Python engines)

Dynamic programming engines for the Bioinformatics DP Engine & Visualizer. Algorithm write-ups live in [`docs/algorithms/`](../../docs/algorithms/).

## Install and test

```bash
cd engines/python
pip install -e ".[dev]"
pytest            # 116 tests: golden fixtures, exhaustive oracle, Biopython cross-checks
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

## Regenerating golden fixtures

```bash
python ../../scripts/generate_nw_golden.py
```

Expected values come from exhaustive enumeration, never from the engine, and are cross-checked against Biopython.
