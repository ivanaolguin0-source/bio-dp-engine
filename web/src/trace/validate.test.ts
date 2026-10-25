// Contract test: every trace the Python engine exported must satisfy the
// shared JSON Schema plus the cross-field checks in parseTrace.
import { readFileSync, readdirSync } from "node:fs";
import { join } from "node:path";
import { describe, expect, it } from "vitest";
import { parseTrace, TraceValidationError } from "./validate";

const dir = join(__dirname, "../../public/traces");
const files = readdirSync(dir).filter((f) => f.endsWith(".json") && f !== "index.json");
const load = (f: string) => JSON.parse(readFileSync(join(dir, f), "utf8"));

describe("bundled traces", () => {
  it("exist", () => expect(files.length).toBeGreaterThanOrEqual(8));
  it.each(files)("%s satisfies the DPTrace contract", (f) => {
    expect(() => parseTrace(load(f))).not.toThrow();
  });
  it("manifest lists every trace file", () => {
    const manifest = JSON.parse(readFileSync(join(dir, "index.json"), "utf8"));
    expect(manifest.traces.map((t: { file: string }) => t.file).sort()).toEqual([...files].sort());
  });
});

describe("rejection", () => {
  const good = () => load("nw-005-textbook-gattaca.json");
  it("rejects a missing required field and names it", () => {
    const t = good();
    delete t.projection;
    expect(() => parseTrace(t)).toThrow(/missing required field 'projection'/);
  });
  it("names unexpected fields, one line each, without duplicates", () => {
    try {
      parseTrace({ name: "bio-dp-visualizer", scripts: {}, version: "0.1.0" });
      throw new Error("should have thrown");
    } catch (e) {
      expect(e).toBeInstanceOf(TraceValidationError);
      const problems = (e as TraceValidationError).problems;
      expect(problems).toContain("(top level): unexpected field 'scripts'");
      expect(problems).toContain("(top level): missing required field 'schema_version'");
      expect(new Set(problems).size).toBe(problems.length);
    }
  });
  it("rejects a numeric optimal_count (must be a string for > 2^53 values)", () => {
    const t = good();
    t.result.optimal_count = 3;
    expect(() => parseTrace(t)).toThrow(/optimal_count/);
  });
  it("rejects a grid whose size disagrees with shape", () => {
    const t = good();
    t.projection.moves.pop();
    expect(() => parseTrace(t)).toThrow(/projection.moves/);
  });
  it("rejects a folding trace without its structures", () => {
    const t = load("nus-003-simple-hairpin.json");
    delete t.result.structures;
    expect(() => parseTrace(t)).toThrow(/result.structures/);
  });
  it("rejects a folding trace whose shape disagrees with the sequence", () => {
    const t = load("nus-003-simple-hairpin.json");
    t.inputs.sequence = t.inputs.sequence + "A";
    expect(() => parseTrace(t)).toThrow(/shape must be/);
  });
  it("rejects a path that does not reach the corner", () => {
    const t = good();
    t.result.alignments[0].path[t.result.alignments[0].path.length - 1] = [0, 0];
    expect(() => parseTrace(t)).toThrow(/path must run/);
  });
});
