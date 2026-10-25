// TypeScript view of contracts/schemas/dptrace.schema.json (version 1.1).
// Hand-written for now; replace with generated types once the OpenAPI
// contract lands. Runtime shape is enforced separately by validate.ts.

export type Grid<T> = T[][];
export type Cell = [number, number];
export type Algorithm = "needleman-wunsch" | "smith-waterman" | "hirschberg" | "nussinov" | "viterbi";

export interface TraceAlignment {
  aligned_a: string;
  aligned_b: string;
  cigar?: string;
  path: Cell[];
}

export interface TraceStructure {
  dot_bracket: string;
  pairs: Cell[];
  /** Table cells the traceback visits, in order. */
  visited: Cell[];
}

export interface TraceLayer {
  name: string;
  scores: Grid<number | null>;
  predecessors: Grid<number>;
}

export interface ReferenceStructure {
  dot_bracket: string;
  source: string;
  pairs: Cell[];
}

export interface DPTrace {
  schema_version: string;
  algorithm: Algorithm;
  engine: { language: "python" | "java"; name: string; version: string };
  inputs: {
    // Alignment
    seq_a?: string;
    seq_b?: string;
    scoring?: {
      substitution: { name: string; alphabet: string; scores: Grid<number> };
      gap: { open: number; extend: number };
    };
    // RNA folding
    sequence?: string;
    folding?: { min_loop: number; wobble: boolean };
    reference?: ReferenceStructure;
  };
  shape: [number, number];
  fill_order: { kind: "row-major" | "anti-diagonal" | "interval-length" | "explicit"; cells?: Cell[] };
  projection: {
    scores: Grid<number | null>;
    moves: Grid<number>;
    on_optimal_path: Grid<boolean>;
    move_labels?: Record<string, string>;
  };
  layers?: TraceLayer[];
  result: {
    score: number;
    optimal_count: string;
    alignments?: TraceAlignment[];
    alignments_truncated?: boolean;
    structures?: TraceStructure[];
    structures_truncated?: boolean;
  };
}

/** An alignment trace after validation: alignment fields are guaranteed present. */
export type AlignmentTrace = DPTrace & {
  inputs: Required<Pick<DPTrace["inputs"], "seq_a" | "seq_b" | "scoring">> & DPTrace["inputs"];
  result: DPTrace["result"] & { alignments: TraceAlignment[]; alignments_truncated: boolean };
};

/** A folding trace after validation: folding fields are guaranteed present. */
export type FoldingTrace = DPTrace & {
  inputs: Required<Pick<DPTrace["inputs"], "sequence" | "folding">> & DPTrace["inputs"];
  result: DPTrace["result"] & { structures: TraceStructure[]; structures_truncated: boolean };
};

export const isFolding = (t: DPTrace): t is FoldingTrace => t.algorithm === "nussinov";
export const isAlignment = (t: DPTrace): t is AlignmentTrace => !isFolding(t);

export interface TraceManifestEntry {
  id: string;
  algorithm?: Algorithm;
  title: string;
  description: string;
  file: string;
  score: number;
  optimal_count: string;
}

// Move bits in projection.moves for alignment traces (see biodp/trace/dptrace.py).
export const MOVE_DIAGONAL = 1;
export const MOVE_UP = 2;
export const MOVE_LEFT = 4;

// Move bits for Nussinov traces.
export const FOLD_I_PAIRS_J = 1;
export const FOLD_J_PAIRS_INNER = 2;
export const FOLD_J_UNPAIRED = 4;
