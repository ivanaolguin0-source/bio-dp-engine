import { readFileSync } from "node:fs";
import { join } from "node:path";
import { describe, expect, it } from "vitest";
import { parseTrace } from "../trace/validate";
import { buildTimeline, columnKind, fillOrder, frameAt, moveOffsets } from "./replay";

const trace = parseTrace(JSON.parse(readFileSync(join(__dirname, "../../public/traces/nw-005-textbook-gattaca.json"), "utf8")));
const alignment = trace.result.alignments[0];
const timeline = buildTimeline(trace, alignment);
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

describe("frames", () => {
  it("starts empty", () => {
    const f = frameAt(trace, timeline, alignment, 0);
    expect(f.phase).toBe("fill");
    expect(f.revealed.flat().some(Boolean)).toBe(false);
  });
  it("reveals exactly k cells after k steps", () => {
    const f = frameAt(trace, timeline, alignment, 10);
    expect(f.revealed.flat().filter(Boolean)).toHaveLength(10);
    expect(f.cursor).toEqual(timeline.fillOrder[9]);
  });
  it("the step right after the fill is traceback with nothing traced yet", () => {
    // Regression: the UI once read traced[traced.length - 1] here and crashed.
    const f = frameAt(trace, timeline, alignment, timeline.fillSteps);
    expect(f.phase).toBe("traceback");
    expect(f.traced).toEqual([]);
    expect(f.revealed.flat().every(Boolean)).toBe(true);
    expect(f.firstColumn).toBe(alignment.aligned_a.length);
  });
  it("traceback starts at the bottom-right corner and walks back", () => {
    const f = frameAt(trace, timeline, alignment, timeline.fillSteps + 2);
    expect(f.phase).toBe("traceback");
    expect(f.traced[0]).toEqual([rows - 1, cols - 1]);
    expect(f.traced).toHaveLength(2);
    expect(f.firstColumn).toBe(alignment.aligned_a.length - 1);
  });
  it("ends with the full path and every column visible", () => {
    const f = frameAt(trace, timeline, alignment, timeline.total);
    expect(f.phase).toBe("done");
    expect(f.traced[f.traced.length - 1]).toEqual([0, 0]);
    expect(f.firstColumn).toBe(0);
  });
  it("clamps out-of-range steps", () => {
    expect(frameAt(trace, timeline, alignment, -5).phase).toBe("fill");
    expect(frameAt(trace, timeline, alignment, 1e9).phase).toBe("done");
  });
});

describe("helpers", () => {
  it("classifies columns", () => {
    expect(columnKind("A", "A")).toBe("match");
    expect(columnKind("A", "C")).toBe("mismatch");
    expect(columnKind("-", "C")).toBe("gap");
  });
  it("decodes move bits", () => {
    expect(moveOffsets(7)).toEqual([[-1, -1], [-1, 0], [0, -1]]);
    expect(moveOffsets(0)).toEqual([]);
  });
});
