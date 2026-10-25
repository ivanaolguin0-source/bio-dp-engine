// The alignment as three monospace rows, revealed column by column as the
// traceback walks from the end of the sequences back to the start.
import type { TraceAlignment } from "../trace/types";
import { columnKind } from "./replay";

interface Props {
  alignment: TraceAlignment | undefined;
  firstColumn: number;
  phase: "fill" | "traceback" | "done";
}

const MIDLINE = { match: "|", mismatch: ".", gap: " " } as const;

export function AlignmentView({ alignment, firstColumn, phase }: Props) {
  if (!alignment) return <p className="muted">No alignment in this trace.</p>;
  if (alignment.aligned_a.length === 0) return <p className="muted">Both sequences are empty: the optimal alignment has no columns.</p>;
  if (phase === "fill") return <p className="muted">The alignment is reconstructed during traceback, after the table is complete.</p>;

  const columns = Array.from(alignment.aligned_a, (top, k) => {
    const bottom = alignment.aligned_b[k];
    return { top, bottom, kind: columnKind(top, bottom), visible: k >= firstColumn, current: phase === "traceback" && k === firstColumn };
  });

  return (
    <div className="alignment" role="table" aria-label="Optimal alignment">
      {(["top", "mid", "bottom"] as const).map((row) => (
        <div className="alignment-row" role="row" key={row}>
          {columns.map((c, k) => {
            const text = row === "top" ? c.top : row === "bottom" ? c.bottom : MIDLINE[c.kind];
            const cls = ["col", row === "mid" ? "mid" : `kind-${c.kind}`, c.visible ? "" : "hidden", c.current ? "current" : ""].join(" ");
            return (
              <span key={k} className={cls} role="cell">
                {c.visible ? text : " "}
              </span>
            );
          })}
        </div>
      ))}
    </div>
  );
}
