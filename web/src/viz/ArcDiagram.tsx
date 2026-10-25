// Arc diagram of an RNA secondary structure.
//
// Bases run along a horizontal line. Each predicted pair is an arc above the
// line, coloured by pair type; when the trace carries a reference structure,
// its pairs are arcs below the line in neutral ink, so prediction and truth
// can be compared directly. Non-crossing arcs are exactly what "no
// pseudoknots" means, which this view makes visible at a glance.
import { useEffect, useRef, useState } from "react";
import type { Cell } from "../trace/types";

export type PairType = "GC" | "AU" | "GU";

export function pairType(x: string, y: string): PairType {
  const key = [x, y].sort().join("");
  if (key === "CG") return "GC";
  if (key === "AU") return "AU";
  return "GU";
}

export interface ArcColors {
  GC: string;
  AU: string;
  GU: string;
  reference: string;
  text: string;
  muted: string;
  accent: string;
}

// First three categorical slots of the reference palette (blue, orange, aqua);
// they pass the colour-vision checks as a set of three in both themes.
export const ARC_LIGHT = { GC: "#2a78d6", AU: "#eb6834", GU: "#1baf7a" };
export const ARC_DARK = { GC: "#3987e5", AU: "#d95926", GU: "#199e70" };

interface Props {
  sequence: string;
  pairs: Cell[];
  reference?: Cell[];
  /** Index of the base the traceback is deciding now, if any. */
  focus?: number | null;
  colors: ArcColors;
}

export function ArcDiagram({ sequence, pairs, reference, focus, colors }: Props) {
  const wrapRef = useRef<HTMLDivElement>(null);
  const [width, setWidth] = useState(720);
  useEffect(() => {
    const el = wrapRef.current;
    if (!el) return;
    const observer = new ResizeObserver(([entry]) => setWidth(entry.contentRect.width));
    observer.observe(el);
    return () => observer.disconnect();
  }, []);

  const n = sequence.length;
  const pad = 14;
  const step = Math.max(6, Math.min(24, (width - 2 * pad) / Math.max(n, 1)));
  const x = (i: number) => pad + step * (i + 0.5);
  const span = (p: Cell) => x(p[1]) - x(p[0]);
  // One shared vertical scale, so nested arcs keep their proportions and never
  // touch: the widest arc (above or below) is MAX_ARC_HEIGHT tall.
  const MAX_ARC_HEIGHT = 140;
  const widest = Math.max(1, ...[...pairs, ...(reference ?? [])].map((p) => span(p) / 2));
  const vScale = Math.min(1, MAX_ARC_HEIGHT / widest);
  const height_of = (p: Cell) => (span(p) / 2) * vScale;
  const maxUp = Math.max(18, ...pairs.map(height_of));
  const maxDown = reference?.length ? Math.max(18, ...reference.map(height_of)) : 0;
  const baseline = maxUp + 10;
  const letterBand = step >= 9 ? 22 : 12;
  const height = baseline + letterBand + (reference?.length ? maxDown + 10 : 4) + 14;
  const svgWidth = pad * 2 + step * n;

  const arc = (p: Cell, below: boolean) => {
    const x1 = x(p[0]);
    const x2 = x(p[1]);
    const rx = (x2 - x1) / 2;
    const ry = rx * vScale;
    const y = below ? baseline + letterBand : baseline;
    // sweep 1 bows the arc upward on screen; sweep 0 bows it downward.
    return `M ${x1} ${y} A ${rx} ${ry} 0 0 ${below ? 0 : 1} ${x2} ${y}`;
  };

  return (
    <div ref={wrapRef} className="arc-wrap">
      <svg
        width={svgWidth}
        height={height}
        role="img"
        aria-label={`Arc diagram: ${pairs.length} predicted base pairs${reference ? `, ${reference.length} reference pairs below the line` : ""}.`}
      >
        {pairs.map((p) => {
          const type = pairType(sequence[p[0]], sequence[p[1]]);
          return (
            <path key={`p${p[0]}-${p[1]}`} d={arc(p, false)} fill="none" stroke={colors[type]} strokeWidth={2} strokeLinecap="round">
              <title>{`Predicted pair ${p[0] + 1}-${p[1] + 1}: ${sequence[p[0]]}-${sequence[p[1]]}`}</title>
            </path>
          );
        })}
        {reference?.map((p) => (
          <path key={`r${p[0]}-${p[1]}`} d={arc(p, true)} fill="none" stroke={colors.reference} strokeWidth={1.5}>
            <title>{`Reference pair ${p[0] + 1}-${p[1] + 1}: ${sequence[p[0]]}-${sequence[p[1]]}`}</title>
          </path>
        ))}
        <line x1={pad} x2={svgWidth - pad} y1={baseline + letterBand / 2} y2={baseline + letterBand / 2} stroke={colors.muted} strokeWidth={0.5} opacity={step >= 9 ? 0 : 1} />
        {Array.from(sequence, (base, i) => (
          <g key={i}>
            {focus === i && (
              <rect x={x(i) - step / 2} y={baseline + 1} width={step} height={letterBand - 2} rx={3} fill="none" stroke={colors.accent} strokeWidth={2} />
            )}
            {step >= 9 ? (
              <text x={x(i)} y={baseline + letterBand / 2 + 1} textAnchor="middle" dominantBaseline="middle" fontSize={Math.min(14, step * 0.75)} fontFamily="ui-monospace, Menlo, monospace" fill={colors.text}>
                {base}
              </text>
            ) : null}
          </g>
        ))}
        {Array.from({ length: Math.floor(n / 10) }, (_, k) => (k + 1) * 10).map((pos) => (
          <text key={pos} x={x(pos - 1)} y={height - 3} textAnchor="middle" fontSize={10} fill={colors.muted}>
            {pos}
          </text>
        ))}
      </svg>
    </div>
  );
}

/** Dot-bracket from a list of pairs. */
export function toDotBracket(n: number, pairs: Cell[]): string {
  const chars = Array.from({ length: n }, () => ".");
  for (const [i, j] of pairs) {
    chars[i] = "(";
    chars[j] = ")";
  }
  return chars.join("");
}

/** Shared pairs, sensitivity (recall) and PPV (precision) of a prediction against a reference. */
export function comparePairs(predicted: Cell[], reference: Cell[]) {
  const ref = new Set(reference.map(([i, j]) => `${i},${j}`));
  const shared = predicted.filter(([i, j]) => ref.has(`${i},${j}`)).length;
  return {
    shared,
    sensitivity: reference.length ? shared / reference.length : NaN,
    ppv: predicted.length ? shared / predicted.length : NaN,
  };
}
