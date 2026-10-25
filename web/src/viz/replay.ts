// Pure replay logic: turns a trace plus a step counter into "what is visible".
// No React and no canvas here, so every rule is unit-testable.
import { isFolding, type Algorithm, type Cell, type DPTrace, type TraceAlignment, type TraceStructure } from "../trace/types";

export type Phase = "fill" | "traceback" | "done";

export interface Timeline {
  fillOrder: Cell[];
  fillSteps: number;
  /** Cells the traceback visits, in the order it visits them. */
  walk: Cell[];
  tracebackSteps: number;
  total: number;
}

/** Cell order the engine used to fill the table, expanded to an explicit list. */
export function fillOrder(trace: DPTrace): Cell[] {
  const [rows, cols] = trace.shape;
  const kind = trace.fill_order.kind;
  if (kind === "explicit") return trace.fill_order.cells ?? [];
  const cells: Cell[] = [];
  if (kind === "anti-diagonal") {
    for (let d = 0; d <= rows + cols - 2; d++) {
      for (let i = Math.max(0, d - cols + 1); i <= Math.min(d, rows - 1); i++) cells.push([i, d - i]);
    }
    return cells;
  }
  if (kind === "interval-length") {
    // Upper triangle, one diagonal at a time: every interval of length d
    // before any interval of length d + 1.
    for (let d = 0; d < cols; d++) for (let i = 0; i + d < rows; i++) cells.push([i, i + d]);
    return cells;
  }
  for (let i = 0; i < rows; i++) for (let j = 0; j < cols; j++) cells.push([i, j]);
  return cells;
}

/** The traceback walk for the selected alignment or structure. */
export function walkFor(trace: DPTrace, index: number): Cell[] {
  if (isFolding(trace)) return trace.result.structures[index]?.visited ?? [];
  const path = trace.result.alignments?.[index]?.path ?? [];
  return [...path].reverse(); // alignment traceback runs from the corner back to (0, 0)
}

export function buildTimeline(trace: DPTrace, walk: Cell[]): Timeline {
  const order = fillOrder(trace);
  return { fillOrder: order, fillSteps: order.length, walk, tracebackSteps: walk.length, total: order.length + walk.length };
}

export interface Frame {
  phase: Phase;
  /** revealed[i][j] is true once the fill has reached that cell. */
  revealed: boolean[][];
  /** The cell being filled right now, if any. */
  cursor: Cell | null;
  /** Traceback cells visited so far, in visiting order. */
  traced: Cell[];
  /** Alignment only: columns [firstColumn, end) have been reconstructed. */
  firstColumn: number;
}

export function frameAt(trace: DPTrace, timeline: Timeline, step: number): Frame {
  const [rows, cols] = trace.shape;
  const s = Math.max(0, Math.min(Math.floor(step), timeline.total));
  const filled = Math.min(s, timeline.fillSteps);
  const revealed = Array.from({ length: rows }, () => new Array<boolean>(cols).fill(false));
  for (let k = 0; k < filled; k++) {
    const [i, j] = timeline.fillOrder[k];
    revealed[i][j] = true;
  }
  const cursor = filled > 0 && filled < timeline.fillSteps ? timeline.fillOrder[filled - 1] : null;

  const tracedCount = Math.max(0, s - timeline.fillSteps);
  const traced = timeline.walk.slice(0, tracedCount);
  const columns = timeline.walk.length - 1;
  // After visiting t cells walking back from the end, t - 1 columns are known.
  const firstColumn = tracedCount === 0 ? columns : Math.max(0, columns - (tracedCount - 1));

  let phase: Phase = "fill";
  if (s >= timeline.total && timeline.total > 0) phase = "done";
  else if (s >= timeline.fillSteps) phase = "traceback";
  return { phase, revealed, cursor, traced, firstColumn };
}

export type ColumnKind = "match" | "mismatch" | "gap";

export function columnKind(top: string, bottom: string): ColumnKind {
  if (top === "-" || bottom === "-") return "gap";
  return top === bottom ? "match" : "mismatch";
}

/**
 * Neighbour offsets (di, dj) that a cell's move bits point to.
 *
 * Alignment: 1 diagonal (up-left), 2 up, 4 left.
 * Nussinov:  1 "i pairs with j" points to the inner interval (i+1, j-1), down-left;
 *            4 "j unpaired" points to (i, j-1), left. Bit 2 splits the interval in
 *            two and has no single neighbour; see hasSplit().
 */
export function moveOffsets(bits: number, algorithm: Algorithm = "needleman-wunsch"): Array<[number, number]> {
  const out: Array<[number, number]> = [];
  if (algorithm === "nussinov") {
    if (bits & 1) out.push([1, -1]);
    if (bits & 4) out.push([0, -1]);
    return out;
  }
  if (bits & 1) out.push([-1, -1]);
  if (bits & 2) out.push([-1, 0]);
  if (bits & 4) out.push([0, -1]);
  return out;
}

export const hasSplit = (bits: number, algorithm: Algorithm) => algorithm === "nussinov" && (bits & 2) !== 0;

/**
 * Cells a Nussinov cell (i, j) reads: (i, j-1) for "j unpaired", plus, for
 * "j pairs with k", the row segment N(i, k-1) and the column segment N(k+1, j-1).
 * Together they form an L along row i and column j-1.
 */
export function foldingInputs(i: number, j: number, minLoop: number): Cell[] {
  const cells: Cell[] = [];
  if (j - 1 >= i) cells.push([i, j - 1]);
  for (let k = i; k <= j - minLoop - 1; k++) {
    if (k - 1 >= i) cells.push([i, k - 1]);
    if (j - 1 >= k + 1) cells.push([k + 1, j - 1]);
  }
  const seen = new Set<string>();
  return cells.filter(([a, b]) => (seen.has(`${a},${b}`) ? false : (seen.add(`${a},${b}`), true)));
}

/**
 * Pairs of a structure that the traceback has decided so far. Each visited
 * interval has a distinct right end j, and the pair (k, j) is decided at the
 * interval ending in j, so a pair is known once that cell has been visited.
 */
export function decidedPairs(structure: TraceStructure, traced: Cell[]): Cell[] {
  const rightEnds = new Set(traced.map(([, j]) => j));
  return structure.pairs.filter(([, j]) => rightEnds.has(j));
}

export type { TraceAlignment };
