// Canvas heatmap of a DP table, drawn from a DPTrace and a replay Frame.
import { useEffect, useRef, useState } from "react";
import type { Cell, DPTrace } from "../trace/types";
import { inkFor, scoreColor, scoreExtent, type Palette } from "./color";
import { moveOffsets, type Frame } from "./replay";

export interface HeatmapOptions {
  /** -1 for the projected single matrix, otherwise an index into trace.layers. */
  layer: number;
  showArrows: boolean;
  showNumbers: boolean;
  showOptimal: boolean;
}

interface Props {
  trace: DPTrace;
  frame: Frame;
  palette: Palette;
  ink: { primary: string; secondary: string; muted: string; accent: string };
  options: HeatmapOptions;
  onHover: (cell: Cell | null, clientX: number, clientY: number) => void;
}

// Fixed direction for each raw Gotoh layer: M diagonal, X up, Y left.
const LAYER_MOVE = [1, 2, 4];

export function layerGrid(trace: DPTrace, layer: number): Array<Array<number | null>> {
  return layer < 0 ? trace.projection.scores : trace.layers![layer].scores;
}

export function Heatmap({ trace, frame, palette, ink, options, onHover }: Props) {
  const wrapRef = useRef<HTMLDivElement>(null);
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const [width, setWidth] = useState(640);
  const [rows, cols] = trace.shape;

  useEffect(() => {
    const el = wrapRef.current;
    if (!el) return;
    const observer = new ResizeObserver(([entry]) => setWidth(entry.contentRect.width));
    observer.observe(el);
    return () => observer.disconnect();
  }, []);

  const cell = Math.max(8, Math.min(52, Math.floor((width - 4) / (cols + 1))));
  const header = cell;
  const cssW = header + cols * cell;
  const cssH = header + rows * cell;

  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    const dpr = window.devicePixelRatio || 1;
    canvas.width = Math.round(cssW * dpr);
    canvas.height = Math.round(cssH * dpr);
    canvas.style.width = `${cssW}px`;
    canvas.style.height = `${cssH}px`;
    const ctx = canvas.getContext("2d")!;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, cssW, cssH);

    const grid = layerGrid(trace, options.layer);
    const extent = scoreExtent(grid);
    const gap = cell >= 20 ? 2 : 1;
    const x0 = (j: number) => header + j * cell;
    const y0 = (i: number) => header + i * cell;
    const cx = (j: number) => x0(j) + cell / 2;
    const cy = (i: number) => y0(i) + cell / 2;
    const font = (px: number, weight = 500) =>
      `${weight} ${px}px ui-monospace, "SF Mono", Menlo, Consolas, monospace`;

    // Residue headers. Row/column 0 is the empty prefix.
    ctx.textAlign = "center";
    ctx.textBaseline = "middle";
    ctx.font = font(Math.max(9, Math.min(15, cell * 0.5)), 600);
    ctx.fillStyle = ink.secondary;
    const seqA = trace.inputs.seq_a;
    const seqB = trace.inputs.seq_b;
    for (let j = 0; j < cols; j++) ctx.fillText(j === 0 ? "∅" : seqB[j - 1], cx(j), header / 2);
    for (let i = 0; i < rows; i++) ctx.fillText(i === 0 ? "∅" : seqA[i - 1], header / 2, cy(i));

    // Cells (fills only; numbers are drawn after the traceback line).
    ctx.font = font(Math.max(8, Math.min(13, cell * 0.3)));
    const fills: Array<Array<string | null>> = Array.from({ length: rows }, () => new Array<string | null>(cols).fill(null));
    for (let i = 0; i < rows; i++) {
      for (let j = 0; j < cols; j++) {
        const x = x0(j) + gap / 2;
        const y = y0(i) + gap / 2;
        const size = cell - gap;
        if (!frame.revealed[i][j]) {
          ctx.strokeStyle = palette.empty;
          ctx.lineWidth = 1;
          ctx.strokeRect(x + 0.5, y + 0.5, size - 1, size - 1);
          continue;
        }
        const fill = scoreColor(grid[i][j], extent, palette);
        fills[i][j] = fill;
        ctx.fillStyle = fill;
        ctx.beginPath();
        ctx.roundRect(x, y, size, size, Math.min(4, size / 4));
        ctx.fill();
      }
    }

    // The traceback line goes under the numbers so it never hides a score.
    const numbersShown = options.showNumbers && cell >= 24;
    if (frame.traced.length > 1) {
      ctx.strokeStyle = ink.accent;
      ctx.lineWidth = Math.max(2, cell * 0.07);
      ctx.lineCap = "round";
      ctx.lineJoin = "round";
      ctx.beginPath();
      frame.traced.forEach(([i, j], k) => (k === 0 ? ctx.moveTo(cx(j), cy(i)) : ctx.lineTo(cx(j), cy(i))));
      ctx.stroke();
    }

    // Scores, each with a halo in its own cell colour so the line passes behind the digits.
    if (numbersShown) {
      ctx.lineJoin = "round";
      ctx.lineWidth = Math.max(3, cell * 0.12);
      for (let i = 0; i < rows; i++) {
        for (let j = 0; j < cols; j++) {
          const fill = fills[i][j];
          if (!fill) continue;
          const value = grid[i][j];
          const label = value === null ? "\u00b7" : String(value);
          ctx.strokeStyle = fill;
          ctx.strokeText(label, cx(j), cy(i));
          ctx.fillStyle = inkFor(fill, palette);
          ctx.fillText(label, cx(j), cy(i));
        }
      }
    }

    // Arrows toward each optimal predecessor.
    if (options.showArrows && cell >= 16) {
      ctx.strokeStyle = ink.muted;
      ctx.fillStyle = ink.muted;
      ctx.lineWidth = 1.5;
      for (let i = 0; i < rows; i++) {
        for (let j = 0; j < cols; j++) {
          if (!frame.revealed[i][j]) continue;
          const bits =
            options.layer < 0
              ? trace.projection.moves[i][j]
              : trace.layers![options.layer].predecessors[i][j]
                ? LAYER_MOVE[options.layer]
                : 0;
          for (const [di, dj] of moveOffsets(bits)) drawArrow(ctx, cx(j), cy(i), di, dj, cell);
        }
      }
    }

    // Cells lying on at least one optimal alignment.
    if (options.showOptimal && frame.phase !== "fill") {
      ctx.strokeStyle = ink.secondary;
      ctx.lineWidth = 1.5;
      ctx.setLineDash([3, 2]);
      for (let i = 0; i < rows; i++) {
        for (let j = 0; j < cols; j++) {
          if (trace.projection.on_optimal_path[i][j]) ctx.strokeRect(x0(j) + 2.5, y0(i) + 2.5, cell - 5, cell - 5);
        }
      }
      ctx.setLineDash([]);
    }

    // The cell being filled, and the three cells its recurrence reads.
    if (frame.cursor) {
      const [i, j] = frame.cursor;
      ctx.strokeStyle = ink.secondary;
      ctx.lineWidth = 1.5;
      ctx.setLineDash([2, 2]);
      for (const [di, dj] of [[-1, -1], [-1, 0], [0, -1]]) {
        if (i + di >= 0 && j + dj >= 0) ctx.strokeRect(x0(j + dj) + 1.5, y0(i + di) + 1.5, cell - 3, cell - 3);
      }
      ctx.setLineDash([]);
      ctx.strokeStyle = ink.accent;
      ctx.lineWidth = 3;
      ctx.strokeRect(x0(j) + 1.5, y0(i) + 1.5, cell - 3, cell - 3);
    }

    // Traceback markers. With numbers visible, an outline around each visited
    // cell (thicker for the current one) keeps the score readable; without
    // numbers, small dots on the line are clearer.
    for (const [k, [i, j]] of frame.traced.entries()) {
      const head = k === frame.traced.length - 1;
      ctx.strokeStyle = ink.accent;
      if (numbersShown) {
        const inset = head ? 2 : 3;
        ctx.lineWidth = head ? 3.5 : 2;
        ctx.beginPath();
        ctx.roundRect(x0(j) + inset, y0(i) + inset, cell - 2 * inset, cell - 2 * inset, 4);
        ctx.stroke();
      } else {
        ctx.beginPath();
        ctx.arc(cx(j), cy(i), head ? Math.max(4, cell * 0.22) : Math.max(2.5, cell * 0.12), 0, Math.PI * 2);
        ctx.fillStyle = head ? ink.accent : palette.surface;
        ctx.fill();
        ctx.lineWidth = 2;
        ctx.stroke();
      }
    }
  }, [trace, frame, palette, ink, options, cell, header, cssW, cssH, rows, cols]);

  const locate = (e: React.MouseEvent<HTMLCanvasElement>): Cell | null => {
    const rect = e.currentTarget.getBoundingClientRect();
    const j = Math.floor((e.clientX - rect.left - header) / cell);
    const i = Math.floor((e.clientY - rect.top - header) / cell);
    return i >= 0 && j >= 0 && i < rows && j < cols ? [i, j] : null;
  };

  return (
    <div ref={wrapRef} className="heatmap-wrap">
      <canvas
        ref={canvasRef}
        role="img"
        aria-label={`Dynamic programming table, ${rows} rows by ${cols} columns. Hover a cell for its values.`}
        onMouseMove={(e) => onHover(locate(e), e.clientX, e.clientY)}
        onMouseLeave={(e) => onHover(null, e.clientX, e.clientY)}
      />
    </div>
  );
}

function drawArrow(ctx: CanvasRenderingContext2D, x: number, y: number, di: number, dj: number, cell: number) {
  // A short arrow in the band between the centred score and the cell edge,
  // pointing at the predecessor: left edge, top edge, or top-left corner.
  const reach = cell * 0.46;
  const start = cell * 0.3; // clear of the centred number
  const sx = x + dj * start;
  const sy = y + di * start;
  const ex = x + dj * reach;
  const ey = y + di * reach;
  ctx.beginPath();
  ctx.moveTo(sx, sy);
  ctx.lineTo(ex, ey);
  ctx.stroke();
  const head = Math.max(2.5, cell * 0.08);
  const angle = Math.atan2(di, dj);
  ctx.beginPath();
  ctx.moveTo(ex, ey);
  ctx.lineTo(ex - head * Math.cos(angle - 0.55), ey - head * Math.sin(angle - 0.55));
  ctx.lineTo(ex - head * Math.cos(angle + 0.55), ey - head * Math.sin(angle + 0.55));
  ctx.closePath();
  ctx.fill();
}
