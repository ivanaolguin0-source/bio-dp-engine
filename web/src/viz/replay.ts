// Pure replay logic: turns a trace plus a step counter into "what is visible".
// No React and no canvas here, so every rule is unit-testable.
import type { Cell, DPTrace, TraceAlignment } from "../trace/types";

export type Phase = "fill" | "traceback" | "done";

export interface Timeline {
  fillOrder: Cell[];
  fillSteps: number;
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
  // row-major (and the fallback for kinds this view does not draw yet)
  for (let i = 0; i < rows; i++) for (let j = 0; j < cols; j++) cells.push([i, j]);
  return cells;
}

export function buildTimeline(trace: DPTrace, alignment: TraceAlignment | undefined): Timeline {
  const order = fillOrder(trace);
  const tracebackSteps = alignment ? alignment.path.length : 0;
  return { fillOrder: order, fillSteps: order.length, tracebackSteps, total: order.length + tracebackSteps };
}

export interface Frame {
  phase: Phase;
  /** revealed[i][j] is true once the fill has reached that cell. */
  revealed: boolean[][];
  /** The cell being filled right now, if any. */
  cursor: Cell | null;
  /** Traceback cells visited so far, from (n, m) backwards. */
  traced: Cell[];
  /** Alignment columns revealed so far: indices [firstColumn, columns). */
  firstColumn: number;
}

export function frameAt(trace: DPTrace, timeline: Timeline, alignment: TraceAlignment | undefined, step: number): Frame {
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
  const path = alignment?.path ?? [];
  const traced = path.slice(path.length - tracedCount).reverse();
  const columns = path.length - 1;
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

/** Which of the three predecessor moves a cell's arrow bits encode, as (di, dj) offsets. */
export function moveOffsets(bits: number): Array<[number, number]> {
  const out: Array<[number, number]> = [];
  if (bits & 1) out.push([-1, -1]);
  if (bits & 2) out.push([-1, 0]);
  if (bits & 4) out.push([0, -1]);
  return out;
}
