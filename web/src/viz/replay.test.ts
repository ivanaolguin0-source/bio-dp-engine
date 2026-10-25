import { readFileSync } from "node:fs";
import { join } from "node:path";
import { describe, expect, it } from "vitest";
import { isFolding, type DPTrace, type FoldingTrace } from "../trace/types";
import { parseTrace } from "../trace/validate";
import { comparePairs, pairType, toDotBracket } from "./ArcDiagram";
import { buildTimeline, columnKind, decidedPairs, fillOrder, foldingInputs, frameAt, hasSplit, moveOffsets, walkFor } from "./replay";

const load = (name: string): DPTrace => parseTrace(JSON.parse(readFileSync(join(__dirname, `../../public/traces/${name}.json`), "utf8")));

// ---------------------------------------------------------------- alignment replay

const trace = load("nw-005-textbook-gattaca");
const alignment = trace.result.alignments![0];
const timeline = buildTimeline(trace, walkFor(trace, 0));
const [rows, cols] = trace.shape;

describe("fill order", () => {
  it("row-major visits every cell once, dependencies first", () => {
    const order = fillOrder(trace);
    expect(order).toHaveLength(rows * cols);
    const seen = new Set<string>();
    for (const [i, j] of order) {
      for (const [di, dj] of [[-1, -1], [-1, 0], [0, -1]]) {
        if (i + di >= 0 && j + dj >= 0) expect(seen.has(`${i + di},${j + dj}`)).toBe(true);
      }
      seen.add(`${i},${j}`);
    }
  });
  it("anti-diagonal order also respects dependencies", () => {
    const order = fillOrder({ ...trace, fill_order: { kind: "anti-diagonal" } });
    expect(order).toHaveLength(rows * cols);
    const pos = new Map(order.map(([i, j], k) => [`${i},${j}`, k]));
    for (const [i, j] of order) {
      if (i > 0 && j > 0) expect(pos.get(`${i - 1},${j - 1}`)!).toBeLessThan(pos.get(`${i},${j}`)!);
    }
  });
});

describe("alignment frames", () => {
  it("starts empty", () => {
    const f = frameAt(trace, timeline, 0);
    expect(f.phase).toBe("fill");
    expect(f.revealed.flat().some(Boolean)).toBe(false);
  });
  it("reveals exactly k cells after k steps", () => {
    const f = frameAt(trace, timeline, 10);
    expect(f.revealed.flat().filter(Boolean)).toHaveLength(10);
    expect(f.cursor).toEqual(timeline.fillOrder[9]);
  });
  it("the step right after the fill is traceback with nothing traced yet", () => {
    // Regression: the UI once read traced[traced.length - 1] here and crashed.
    const f = frameAt(trace, timeline, timeline.fillSteps);
    expect(f.phase).toBe("traceback");
    expect(f.traced).toEqual([]);
    expect(f.revealed.flat().every(Boolean)).toBe(true);
    expect(f.firstColumn).toBe(alignment.aligned_a.length);
  });
  it("traceback starts at the bottom-right corner and walks back", () => {
    const f = frameAt(trace, timeline, timeline.fillSteps + 2);
    expect(f.phase).toBe("traceback");
    expect(f.traced[0]).toEqual([rows - 1, cols - 1]);
    expect(f.traced).toHaveLength(2);
    expect(f.firstColumn).toBe(alignment.aligned_a.length - 1);
  });
  it("ends with the full path and every column visible", () => {
    const f = frameAt(trace, timeline, timeline.total);
    expect(f.phase).toBe("done");
    expect(f.traced[f.traced.length - 1]).toEqual([0, 0]);
    expect(f.firstColumn).toBe(0);
  });
  it("clamps out-of-range steps", () => {
    expect(frameAt(trace, timeline, -5).phase).toBe("fill");
    expect(frameAt(trace, timeline, 1e9).phase).toBe("done");
  });
});

// ---------------------------------------------------------------- folding replay

const fold = load("nus-009-many-ties") as FoldingTrace;
const n = fold.inputs.sequence.length;

describe("folding traces", () => {
  it("are recognised as folding", () => {
    expect(isFolding(fold)).toBe(true);
    expect(isFolding(trace)).toBe(false);
  });

  it("interval-length order covers the upper triangle, shortest intervals first", () => {
    const order = fillOrder(fold);
    expect(order).toHaveLength((n * (n + 1)) / 2);
    expect(order.every(([i, j]) => i <= j)).toBe(true);
    const lengths = order.map(([i, j]) => j - i);
    expect(lengths).toEqual([...lengths].sort((a, b) => a - b));
  });

  it("every cell a folding cell reads is filled before it", () => {
    const order = fillOrder(fold);
    const pos = new Map(order.map(([i, j], k) => [`${i},${j}`, k]));
    for (const [i, j] of order) {
      for (const [a, b] of foldingInputs(i, j, fold.inputs.folding.min_loop)) {
        expect(pos.get(`${a},${b}`)!).toBeLessThan(pos.get(`${i},${j}`)!);
      }
    }
  });

  it("the traceback walk starts at the whole sequence and stays in the triangle", () => {
    const walk = walkFor(fold, 0);
    expect(walk[0]).toEqual([0, n - 1]);
    expect(walk.every(([i, j]) => i <= j)).toBe(true);
  });

  it("each traceback step uncovers pairs, and all of them by the end", () => {
    const structure = fold.result.structures[0];
    const walk = walkFor(fold, 0);
    const tl = buildTimeline(fold, walk);
    let previous = 0;
    for (let k = 0; k <= walk.length; k++) {
      const shown = decidedPairs(structure, frameAt(fold, tl, tl.fillSteps + k).traced).length;
      expect(shown).toBeGreaterThanOrEqual(previous);
      previous = shown;
    }
    expect(previous).toBe(structure.pairs.length);
  });

  it("Nussinov move bits point inward and left", () => {
    expect(moveOffsets(1, "nussinov")).toEqual([[1, -1]]);
    expect(moveOffsets(4, "nussinov")).toEqual([[0, -1]]);
    expect(moveOffsets(2, "nussinov")).toEqual([]);
    expect(hasSplit(2, "nussinov")).toBe(true);
    expect(hasSplit(2, "needleman-wunsch")).toBe(false);
  });
});

describe("structure helpers", () => {
  it("dot-bracket round trip", () => {
    for (const s of fold.result.structures) expect(toDotBracket(n, s.pairs)).toBe(s.dot_bracket);
  });
  it("pair types ignore orientation", () => {
    expect(pairType("G", "C")).toBe("GC");
    expect(pairType("C", "G")).toBe("GC");
    expect(pairType("U", "A")).toBe("AU");
    expect(pairType("U", "G")).toBe("GU");
  });
  it("tRNA comparison matches the engine's own figures", () => {
    const trna = load("demo-trna-phe") as FoldingTrace;
    const ref = trna.inputs.reference!;
    const c = comparePairs(trna.result.structures[0].pairs, ref.pairs);
    expect(ref.pairs).toHaveLength(21);
    expect(c.shared).toBe(10);
    expect(c.sensitivity).toBeCloseTo((ref as unknown as { canonical_sensitivity: number }).canonical_sensitivity);
  });
});

describe("helpers", () => {
  it("classifies columns", () => {
    expect(columnKind("A", "A")).toBe("match");
    expect(columnKind("A", "C")).toBe("mismatch");
    expect(columnKind("-", "C")).toBe("gap");
  });
  it("decodes alignment move bits", () => {
    expect(moveOffsets(7)).toEqual([[-1, -1], [-1, 0], [0, -1]]);
    expect(moveOffsets(0)).toEqual([]);
  });
});
