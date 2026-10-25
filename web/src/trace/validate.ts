// Runtime validation against the shared contract. The schema file is imported
// straight from contracts/, so the frontend cannot drift from the engines.
import Ajv2020, { type ErrorObject } from "ajv/dist/2020";
import schema from "../../../contracts/schemas/dptrace.schema.json";
import type { DPTrace } from "./types";

const ajv = new Ajv2020({ allErrors: true, strict: false });
const validateSchema = ajv.compile<DPTrace>(schema);

export class TraceValidationError extends Error {
  constructor(public readonly problems: string[]) {
    super(`Not a valid DPTrace:\n${problems.join("\n")}`);
    this.name = "TraceValidationError";
  }
}

/** Schema check plus the cross-field invariants JSON Schema cannot express. */
export function parseTrace(data: unknown): DPTrace {
  if (!validateSchema(data)) {
    throw new TraceValidationError(describeSchemaErrors(validateSchema.errors ?? []));
  }
  const trace = data;
  const [rows, cols] = trace.shape;
  const problems: string[] = [];
  const checkGrid = (label: string, grid: unknown[][]) => {
    if (grid.length !== rows || grid.some((row) => row.length !== cols)) {
      problems.push(`${label} is not ${rows} x ${cols}`);
    }
  };
  checkGrid("projection.scores", trace.projection.scores);
  checkGrid("projection.moves", trace.projection.moves);
  checkGrid("projection.on_optimal_path", trace.projection.on_optimal_path);
  trace.layers?.forEach((layer) => {
    checkGrid(`layer ${layer.name} scores`, layer.scores);
    checkGrid(`layer ${layer.name} predecessors`, layer.predecessors);
  });
  if (trace.inputs.seq_a.length + 1 !== rows || trace.inputs.seq_b.length + 1 !== cols) {
    problems.push("shape does not match sequence lengths");
  }
  trace.result.alignments.forEach((a, k) => {
    if (a.aligned_a.length !== a.aligned_b.length) problems.push(`alignment ${k}: rows differ in length`);
    if (a.path.length !== a.aligned_a.length + 1) problems.push(`alignment ${k}: path length mismatch`);
    const end = a.path[a.path.length - 1];
    if (a.path[0][0] !== 0 || a.path[0][1] !== 0 || end[0] !== rows - 1 || end[1] !== cols - 1) {
      problems.push(`alignment ${k}: path must run from (0,0) to (${rows - 1},${cols - 1})`);
    }
  });
  if (problems.length) throw new TraceValidationError(problems);
  return trace;
}

/** One readable line per problem, naming the field involved, without duplicates. */
export function describeSchemaErrors(errors: ErrorObject[]): string[] {
  const lines = errors.map((e) => {
    const where = e.instancePath || "(top level)";
    const params = e.params as Record<string, unknown>;
    if (e.keyword === "required") return `${where}: missing required field '${String(params.missingProperty)}'`;
    if (e.keyword === "additionalProperties") return `${where}: unexpected field '${String(params.additionalProperty)}'`;
    if (e.keyword === "enum") return `${where}: must be one of ${JSON.stringify(params.allowedValues)}`;
    return `${where}: ${e.message ?? "is invalid"}`;
  });
  return [...new Set(lines)];
}
