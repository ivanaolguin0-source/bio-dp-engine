# DP Visualizer (web)

React + TypeScript front end that replays DPTrace files. It never computes DP: the Python (and later Java) engines export traces, and this app animates them.

## Run it

Requires Node.js 20.19+ or 22.12+ (check with `node --version`).

```bash
cd web
npm install        # first time only
npm run dev        # then open the URL it prints, usually http://localhost:5173
```

## Test it

```bash
npm test           # 57 unit and contract tests (Vitest)
npm run typecheck  # TypeScript strict mode
npm run build      # production bundle in dist/
```

The contract tests validate every file in `public/traces/` against `contracts/schemas/dptrace.schema.json`, the same schema the Python engine's tests use.

## Regenerate the traces

After changing the engine, re-export from the repo root:

```bash
python scripts/export_traces.py
```

## What it shows

| Element | Meaning |
|---|---|
| Cell colour | Score: blue above zero, red below, grey at zero, pale for unreachable |
| Number | The cell's score (toggle "Scores") |
| Small arrows | Which neighbour(s) produced the maximum; several arrows mean a tie |
| Dashed cells during fill | The three cells the current cell's recurrence reads |
| Dashed cells after fill | Cells on at least one optimal alignment |
| Orange path | The traceback for the selected alignment |
| Matrix menu | Switch between the textbook view H and the raw Gotoh layers M, X, Y |

For RNA folding (Nussinov) traces the same page changes shape:

| Element | Meaning |
|---|---|
| Triangular table | Cell (i, j) is the interval from base i to base j; it fills one diagonal at a time, shortest intervals first |
| Dashed L during fill | Every cell the current interval reads: along its row to the left and down the column to its left |
| Arrows | Left: base j unpaired. Down-left: i pairs with j. Fork mark: j pairs with a base inside, splitting the interval |
| Orange outlines | Intervals the traceback visits; it branches, so there is no single path |
| Arc diagram | Predicted pairs above the line (blue G-C, orange A-U, green G-U); the reference structure, when the trace has one, below it |
| Dot-bracket | The same structures as text, so colour is never the only cue |

Keyboard: Space play/pause, Left/Right step, Home/End jump. A trace can be linked directly with its id in the URL hash, e.g. `#demo-homopolymer-ties`.
