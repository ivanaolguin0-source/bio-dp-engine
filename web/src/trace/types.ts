// TypeScript view of contracts/schemas/dptrace.schema.json.
// Hand-written for now; replace with generated types once the OpenAPI
// contract lands. Runtime shape is enforced separately by validate.ts.

export type Grid<T> = T[][];
export type Cell = [number, number];

export interface TraceAlignment {
  aligned_a: string;
  aligned_b: string;
  cigar?: string;
  path: Cell[];
}

export interface TraceLayer {
  name: string;
  scores: Grid<number | null>;
  predecessors: Grid<number>;
}

export interface DPTrace {
  schema_version: string;
  algorithm: "needleman-wunsch" | "smith-waterman" | "hirschberg" | "nussinov" | "viterbi";
  engine: { language: "python" | "java"; name: string; version: string };
  inputs: {
    seq_a: string;
    seq_b: string;
    scoring: {
      substitution: { name: string; alphabet: string; scores: Grid<number> };
      gap: { open: number; extend: number };
    };
  };
  shape: [number, number];
  fill_order: { kind: "row-major" | "anti-diagonal" | "interval-length" | "explicit"; cells?: Cell[] };
  projection: {
    scores: Grid<number | null>;
    moves: Grid<number>;
    on_optimal_path: Grid<boolean>;
  };
  layers?: TraceLayer[];
  result: {
    score: number;
    optimal_count: string;
    alignments: TraceAlignment[];
    alignments_truncated: boolean;
  };
}

export interface TraceManifestEntry {
  id: string;
  title: string;
  description: string;
  file: string;
  score: number;
  optimal_count: string;
}

// Move bits in projection.moves (see biodp/trace/dptrace.py).
export const MOVE_DIAGONAL = 1;
export const MOVE_UP = 2;
export const MOVE_LEFT = 4;
