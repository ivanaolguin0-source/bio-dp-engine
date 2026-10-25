import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import type { Cell, DPTrace, TraceManifestEntry } from "./trace/types";
import { parseTrace, TraceValidationError } from "./trace/validate";
import { DARK, LIGHT, scoreColor, scoreExtent } from "./viz/color";
import { Heatmap, layerGrid, type HeatmapOptions } from "./viz/Heatmap";
import { AlignmentView } from "./viz/AlignmentView";
import { buildTimeline, frameAt, moveOffsets } from "./viz/replay";

const SPEEDS = [2, 8, 30, 120, 600];
const DEFAULT_TRACE = "nw-005-textbook-gattaca";
const INK = {
  light: { primary: "#0b0b0b", secondary: "#52514e", muted: "rgba(11,11,11,0.62)", accent: "#c26a00" },
  dark: { primary: "#f0efec", secondary: "#c3c2b7", muted: "rgba(240,239,236,0.66)", accent: "#f2a33a" },
};

type Theme = "light" | "dark";

function useSystemTheme(): Theme {
  const query = window.matchMedia?.("(prefers-color-scheme: dark)");
  const [theme, setTheme] = useState<Theme>(query?.matches ? "dark" : "light");
  useEffect(() => {
    if (!query) return;
    const update = () => setTheme(query.matches ? "dark" : "light");
    query.addEventListener("change", update);
    return () => query.removeEventListener("change", update);
  }, [query]);
  return theme;
}

export default function App() {
  const [manifest, setManifest] = useState<TraceManifestEntry[]>([]);
  const [selected, setSelected] = useState<string>(() => decodeURIComponent(location.hash.slice(1)) || DEFAULT_TRACE);
  const [trace, setTrace] = useState<DPTrace | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [alignmentIndex, setAlignmentIndex] = useState(0);
  const [step, setStep] = useState(0);
  const [playing, setPlaying] = useState(false);
  const [speed, setSpeed] = useState(SPEEDS[1]);
  const [options, setOptions] = useState<HeatmapOptions>({ layer: -1, showArrows: true, showNumbers: true, showOptimal: true });
  const [hover, setHover] = useState<{ cell: Cell; x: number; y: number } | null>(null);
  const systemTheme = useSystemTheme();
  const [themeOverride, setThemeOverride] = useState<Theme | null>(null);
  const theme = themeOverride ?? systemTheme;
  const palette = theme === "dark" ? DARK : LIGHT;
  const ink = INK[theme];

  useEffect(() => {
    document.documentElement.dataset.theme = theme;
  }, [theme]);

  // Manifest of bundled traces.
  useEffect(() => {
    fetch(`${import.meta.env.BASE_URL}traces/index.json`)
      .then((r) => r.json())
      .then((data: { traces: TraceManifestEntry[] }) => setManifest(data.traces))
      .catch((e) => setError(`Could not load traces/index.json: ${e}`));
  }, []);

  const showTrace = useCallback((data: unknown) => {
    try {
      setTrace(parseTrace(data));
      setError(null);
    } catch (e) {
      setTrace(null);
      setError(e instanceof TraceValidationError ? e.message : String(e));
    }
    setAlignmentIndex(0);
    setStep(0);
    setPlaying(false);
  }, []);

  useEffect(() => {
    const entry = manifest.find((m) => m.id === selected);
    if (!entry) return;
    location.hash = entry.id;
    fetch(`${import.meta.env.BASE_URL}traces/${entry.file}`)
      .then((r) => r.json())
      .then(showTrace)
      .catch((e) => setError(`Could not load ${entry.file}: ${e}`));
  }, [manifest, selected, showTrace]);

  const alignment = trace?.result.alignments[alignmentIndex];
  const timeline = useMemo(() => (trace ? buildTimeline(trace, alignment) : null), [trace, alignment]);
  const frame = useMemo(
    () => (trace && timeline ? frameAt(trace, timeline, alignment, step) : null),
    [trace, timeline, alignment, step],
  );

  // Animation loop: advance `speed` steps per second while playing.
  const last = useRef<number | null>(null);
  useEffect(() => {
    if (!playing || !timeline) return;
    let raf = 0;
    const tick = (now: number) => {
      const dt = last.current === null ? 0 : (now - last.current) / 1000;
      last.current = now;
      setStep((s) => {
        const next = Math.min(timeline.total, s + dt * speed);
        if (next >= timeline.total) setPlaying(false);
        return next;
      });
      raf = requestAnimationFrame(tick);
    };
    raf = requestAnimationFrame(tick);
    return () => {
      cancelAnimationFrame(raf);
      last.current = null;
    };
  }, [playing, speed, timeline]);

  // Keyboard: space play/pause, arrows step, End jumps to the end.
  useEffect(() => {
    const onKey = (e: KeyboardEvent) => {
      if (!timeline || (e.target as HTMLElement).tagName === "SELECT" || (e.target as HTMLElement).tagName === "INPUT") return;
      if (e.key === " ") {
        e.preventDefault();
        if (Math.floor(step) >= timeline.total) setStep(0);
        setPlaying((p) => !p);
      } else if (e.key === "ArrowRight") setStep((s) => Math.min(timeline.total, Math.floor(s) + 1));
      else if (e.key === "ArrowLeft") setStep((s) => Math.max(0, Math.ceil(s) - 1));
      else if (e.key === "End") setStep(timeline.total);
      else if (e.key === "Home") setStep(0);
    };
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [timeline, step]);

  const onFile = async (file: File | undefined) => {
    if (!file) return;
    try {
      showTrace(JSON.parse(await file.text()));
    } catch (e) {
      setError(`Could not read ${file.name}: ${e}`);
    }
  };

  const entry = manifest.find((m) => m.id === selected);
  const whole = Math.floor(step);

  return (
    <div className="app">
      <header className="top">
        <div>
          <h1>DP Alignment Visualizer</h1>
          <p className="muted">Replays DPTrace files produced by the engines. The browser never computes DP itself.</p>
        </div>
        <button className="ghost" onClick={() => setThemeOverride(theme === "dark" ? "light" : "dark")} aria-label="Toggle colour theme">
          {theme === "dark" ? "Light theme" : "Dark theme"}
        </button>
      </header>

      <section className="toolbar" aria-label="Trace selection">
        <label>
          Trace{" "}
          <select value={selected} onChange={(e) => setSelected(e.target.value)}>
            {manifest.map((m) => (
              <option key={m.id} value={m.id}>
                {m.title}
              </option>
            ))}
          </select>
        </label>
        <label className="file">
          Load DPTrace file
          <input type="file" accept="application/json,.json" onChange={(e) => onFile(e.target.files?.[0])} />
        </label>
      </section>

      {error && (
        <pre className="error" role="alert">
          {error}
        </pre>
      )}

      {trace && timeline && frame && (
        <>
          <section className="summary">
            {entry && <p>{entry.description}</p>}
            <dl>
              <dt>seq_a (rows)</dt>
              <dd className="mono">{trace.inputs.seq_a || "(empty)"}</dd>
              <dt>seq_b (columns)</dt>
              <dd className="mono">{trace.inputs.seq_b || "(empty)"}</dd>
              <dt>Scoring</dt>
              <dd>
                {trace.inputs.scoring.substitution.name}; gap open {trace.inputs.scoring.gap.open}, extend{" "}
                {trace.inputs.scoring.gap.extend}
              </dd>
              <dt>Optimal score</dt>
              <dd className="mono">{trace.result.score}</dd>
              <dt>Co-optimal alignments</dt>
              <dd className="mono">
                {trace.result.optimal_count}
                {trace.result.alignments_truncated && ` (first ${trace.result.alignments.length} included in this trace)`}
              </dd>
            </dl>
          </section>

          <section className="controls" aria-label="Playback">
            <button onClick={() => { if (whole >= timeline.total) setStep(0); setPlaying(!playing); }}>
              {playing ? "Pause" : whole >= timeline.total ? "Replay" : "Play"}
            </button>
            <button onClick={() => setStep(Math.max(0, Math.ceil(step) - 1))} aria-label="Step back">&#9664;</button>
            <button onClick={() => setStep(Math.min(timeline.total, whole + 1))} aria-label="Step forward">&#9654;</button>
            <button className="ghost" onClick={() => setStep(timeline.fillSteps)}>Skip to traceback</button>
            <button className="ghost" onClick={() => { setStep(0); setPlaying(false); }}>Reset</button>
            <label>
              Speed{" "}
              <select value={speed} onChange={(e) => setSpeed(Number(e.target.value))}>
                {SPEEDS.map((s) => (
                  <option key={s} value={s}>
                    {s} steps/s
                  </option>
                ))}
              </select>
            </label>
            <input
              className="scrubber"
              type="range"
              min={0}
              max={timeline.total}
              value={whole}
              onChange={(e) => { setPlaying(false); setStep(Number(e.target.value)); }}
              aria-label="Replay position"
            />
            <span className="mono muted">
              {whole}/{timeline.total}
            </span>
          </section>

          <p className="status" aria-live="polite">
            {statusLine(trace, frame.phase, frame.cursor, frame.traced, timeline.tracebackSteps)}
          </p>

          <section className="view-options" aria-label="View options">
            <label>
              Matrix{" "}
              <select value={options.layer} onChange={(e) => setOptions({ ...options, layer: Number(e.target.value) })}>
                <option value={-1}>H = max(M, X, Y)</option>
                {trace.layers?.map((l, k) => (
                  <option key={l.name} value={k}>
                    {l.name} ({["last column pairs residues", "gap in seq_b", "gap in seq_a"][k] ?? "layer"})
                  </option>
                ))}
              </select>
            </label>
            {(["showNumbers", "showArrows", "showOptimal"] as const).map((key) => (
              <label key={key} className="check">
                <input type="checkbox" checked={options[key]} onChange={(e) => setOptions({ ...options, [key]: e.target.checked })} />
                {{ showNumbers: "Scores", showArrows: "Arrows", showOptimal: "Optimal cells" }[key]}
              </label>
            ))}
            {trace.result.alignments.length > 1 && (
              <label>
                Alignment{" "}
                <select value={alignmentIndex} onChange={(e) => setAlignmentIndex(Number(e.target.value))}>
                  {trace.result.alignments.map((_, k) => (
                    <option key={k} value={k}>
                      {k + 1} of {trace.result.optimal_count}
                      {k === 0 ? " (canonical)" : ""}
                    </option>
                  ))}
                </select>
              </label>
            )}
          </section>

          <section className="stage">
            <Heatmap trace={trace} frame={frame} palette={palette} ink={ink} options={options}
              onHover={(cell, x, y) => setHover(cell ? { cell, x, y } : null)} />
            <Legend trace={trace} layer={options.layer} palette={palette} />
          </section>

          <section className="alignment-panel" aria-label="Alignment">
            <h2>Alignment</h2>
            <AlignmentView alignment={alignment} firstColumn={frame.firstColumn} phase={frame.phase} />
            <p className="muted small">
              <span className="swatch kind-match" /> identical &nbsp; <span className="swatch kind-mismatch" /> substitution &nbsp;
              <span className="swatch kind-gap" /> gap. The midline repeats this as | . and blank, so colour is never the only cue.
            </p>
          </section>

          {hover && frame.revealed[hover.cell[0]][hover.cell[1]] && (
            <CellTooltip trace={trace} cell={hover.cell} x={hover.x} y={hover.y} />
          )}
        </>
      )}

      <footer className="muted small">
        Keyboard: Space play/pause, arrow keys step, Home/End jump. Colour: blue above zero, red below, grey at zero.
      </footer>
    </div>
  );
}

function statusLine(trace: DPTrace, phase: string, cursor: Cell | null, traced: Cell[], total: number): string {
  if (phase === "fill") {
    if (!cursor) return "Press Play (or Space) to fill the table cell by cell.";
    const [i, j] = cursor;
    if (i === 0 && j === 0) return "Cell (0, 0): aligning two empty prefixes scores 0.";
    if (i === 0 || j === 0) return `Cell (${i}, ${j}): boundary. A prefix against nothing is one leading gap.`;
    const a = trace.inputs.seq_a[i - 1];
    const b = trace.inputs.seq_b[j - 1];
    return `Cell (${i}, ${j}): best of diagonal + s(${a}, ${b}), up + gap, left + gap. The three dashed cells are its inputs.`;
  }
  if (phase === "traceback") {
    const head = traced[traced.length - 1];
    if (!head) {
      const [rows, cols] = trace.shape;
      return `Table complete. Traceback starts at the bottom-right cell (${rows - 1}, ${cols - 1}), where the optimal score sits.`;
    }
    return `Traceback ${traced.length} of ${total}: at (${head[0]}, ${head[1]}), following an arrow back toward (0, 0).`;
  }
  return `Done. Score ${trace.result.score}; ${trace.result.optimal_count} optimal alignment(s). Dashed cells lie on at least one of them.`;
}

function CellTooltip({ trace, cell: [i, j], x, y }: { trace: DPTrace; cell: Cell; x: number; y: number }) {
  const { seq_a, seq_b, scoring } = trace.inputs;
  const sub = scoring.substitution;
  const a = i > 0 ? seq_a[i - 1] : null;
  const b = j > 0 ? seq_b[j - 1] : null;
  const s = a && b ? sub.scores[sub.alphabet.indexOf(a)]?.[sub.alphabet.indexOf(b)] : undefined;
  const moves = moveOffsets(trace.projection.moves[i][j]).map(([di, dj]) => (di && dj ? "diagonal" : di ? "up" : "left"));
  const fmt = (v: number | null | undefined) => (v === null || v === undefined ? "unreachable" : String(v));
  const left = Math.min(x + 14, window.innerWidth - 280);
  return (
    <div className="tooltip" style={{ left, top: y + 14 }} role="tooltip">
      <div className="mono strong">
        ({i}, {j}) {a ?? "∅"} / {b ?? "∅"}
      </div>
      <div>Prefixes: "{seq_a.slice(0, i)}" vs "{seq_b.slice(0, j)}"</div>
      <div>H = {fmt(trace.projection.scores[i][j])}</div>
      {trace.layers && <div className="mono">{trace.layers.map((l) => `${l.name}=${fmt(l.scores[i][j])}`).join("  ")}</div>}
      {s !== undefined && <div>s({a}, {b}) = {s}</div>}
      <div>Best from: {moves.length ? moves.join(", ") : "none (origin)"}</div>
      <div>{trace.projection.on_optimal_path[i][j] ? "On an optimal path" : "Not on any optimal path"}</div>
    </div>
  );
}

function Legend({ trace, layer, palette }: { trace: DPTrace; layer: number; palette: typeof LIGHT }) {
  const extent = scoreExtent(layerGrid(trace, layer));
  const stops = Array.from({ length: 11 }, (_, k) => -extent + (2 * extent * k) / 10);
  return (
    <div className="legend" aria-label="Colour legend">
      <div className="ramp">
        {stops.map((v) => (
          <span key={v} style={{ background: scoreColor(v, extent, palette) }} />
        ))}
      </div>
      <div className="ramp-labels mono small">
        <span>-{extent}</span>
        <span>0</span>
        <span>+{extent}</span>
      </div>
      <div className="small muted legend-keys">
        <span><i className="key key-empty" style={{ background: palette.empty }} /> unreachable</span>
        <span><i className="key key-optimal" /> on an optimal path</span>
        <span><i className="key key-trace" /> traceback</span>
      </div>
    </div>
  );
}
