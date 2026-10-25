import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import {
  isFolding,
  type AlignmentTrace,
  type Cell,
  type DPTrace,
  type FoldingTrace,
  type TraceManifestEntry,
} from "./trace/types";
import { parseTrace, TraceValidationError } from "./trace/validate";
import { DARK, LIGHT, scoreColor, scoreExtent, type Palette } from "./viz/color";
import { Heatmap, layerGrid, type HeatmapOptions } from "./viz/Heatmap";
import { AlignmentView } from "./viz/AlignmentView";
import { ARC_DARK, ARC_LIGHT, ArcDiagram, comparePairs, toDotBracket } from "./viz/ArcDiagram";
import { buildTimeline, decidedPairs, frameAt, hasSplit, moveOffsets, walkFor, type Frame } from "./viz/replay";

const SPEEDS = [2, 8, 30, 120, 600];
const DEFAULT_TRACE = "nw-005-textbook-gattaca";
const INK = {
  light: { primary: "#0b0b0b", secondary: "#52514e", muted: "rgba(11,11,11,0.62)", accent: "#c26a00" },
  dark: { primary: "#f0efec", secondary: "#c3c2b7", muted: "rgba(240,239,236,0.66)", accent: "#f2a33a" },
};
const GROUPS: Array<{ label: string; test: (m: TraceManifestEntry) => boolean }> = [
  { label: "Sequence alignment", test: (m) => (m.algorithm ?? "needleman-wunsch") !== "nussinov" },
  { label: "RNA folding (Nussinov)", test: (m) => m.algorithm === "nussinov" },
];

type Theme = "light" | "dark";
type Ink = (typeof INK)["light"];

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
  const [resultIndex, setResultIndex] = useState(0);
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

  // Follow the URL hash after load too (Back/Forward, or a pasted trace link).
  useEffect(() => {
    const onHash = () => {
      const id = decodeURIComponent(location.hash.slice(1));
      if (id) setSelected(id);
    };
    window.addEventListener("hashchange", onHash);
    return () => window.removeEventListener("hashchange", onHash);
  }, []);

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
    setResultIndex(0);
    setStep(0);
    setPlaying(false);
    setOptions((o) => ({ ...o, layer: -1 }));
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

  const walk = useMemo(() => (trace ? walkFor(trace, resultIndex) : []), [trace, resultIndex]);
  const timeline = useMemo(() => (trace ? buildTimeline(trace, walk) : null), [trace, walk]);
  const frame = useMemo(() => (trace && timeline ? frameAt(trace, timeline, step) : null), [trace, timeline, step]);

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

  // Keyboard: Space play/pause, arrows step, Home/End jump.
  useEffect(() => {
    const onKey = (e: KeyboardEvent) => {
      const tag = (e.target as HTMLElement).tagName;
      if (!timeline || tag === "SELECT" || tag === "INPUT") return;
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
  const folding = trace ? isFolding(trace) : false;
  const results = trace ? (isFolding(trace) ? trace.result.structures : trace.result.alignments ?? []) : [];
  const resultNoun = folding ? "Structure" : "Alignment";

  return (
    <div className="app">
      <header className="top">
        <div>
          <h1>DP Visualizer</h1>
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
            {GROUPS.map((g) => (
              <optgroup key={g.label} label={g.label}>
                {manifest.filter(g.test).map((m) => (
                  <option key={m.id} value={m.id}>
                    {m.title}
                  </option>
                ))}
              </optgroup>
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
            {isFolding(trace) ? <FoldingSummary trace={trace} /> : <AlignmentSummary trace={trace as AlignmentTrace} />}
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
            {isFolding(trace)
              ? foldingStatus(trace, frame, timeline.tracebackSteps)
              : alignmentStatus(trace as AlignmentTrace, frame, timeline.tracebackSteps)}
          </p>

          <section className="view-options" aria-label="View options">
            {!folding && (
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
            )}
            {(["showNumbers", "showArrows", "showOptimal"] as const).map((key) => (
              <label key={key} className="check">
                <input type="checkbox" checked={options[key]} onChange={(e) => setOptions({ ...options, [key]: e.target.checked })} />
                {{ showNumbers: "Scores", showArrows: "Arrows", showOptimal: "Optimal cells" }[key]}
              </label>
            ))}
            {results.length > 1 && (
              <label>
                {resultNoun}{" "}
                <select value={resultIndex} onChange={(e) => { setResultIndex(Number(e.target.value)); }}>
                  {results.map((_, k) => (
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

          {isFolding(trace) ? (
            <StructurePanel trace={trace} index={resultIndex} frame={frame} theme={theme} ink={ink} />
          ) : (
            <section className="alignment-panel" aria-label="Alignment">
              <h2>Alignment</h2>
              <AlignmentView alignment={(trace as AlignmentTrace).result.alignments[resultIndex]} firstColumn={frame.firstColumn} phase={frame.phase} />
              <p className="muted small">
                <span className="swatch kind-match" /> identical &nbsp; <span className="swatch kind-mismatch" /> substitution &nbsp;
                <span className="swatch kind-gap" /> gap. The midline repeats this as | . and blank, so colour is never the only cue.
              </p>
            </section>
          )}

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

// ------------------------------------------------------------------ summaries

function AlignmentSummary({ trace }: { trace: AlignmentTrace }) {
  const { seq_a, seq_b, scoring } = trace.inputs;
  return (
    <dl>
      <dt>seq_a (rows)</dt>
      <dd className="mono">{seq_a || "(empty)"}</dd>
      <dt>seq_b (columns)</dt>
      <dd className="mono">{seq_b || "(empty)"}</dd>
      <dt>Scoring</dt>
      <dd>
        {scoring.substitution.name}; gap open {scoring.gap.open}, extend {scoring.gap.extend}
      </dd>
      <dt>Optimal score</dt>
      <dd className="mono">{trace.result.score}</dd>
      <dt>Co-optimal alignments</dt>
      <dd className="mono">
        {trace.result.optimal_count}
        {trace.result.alignments_truncated && ` (first ${trace.result.alignments.length} included in this trace)`}
      </dd>
    </dl>
  );
}

function FoldingSummary({ trace }: { trace: FoldingTrace }) {
  const { sequence, folding, reference } = trace.inputs;
  return (
    <dl>
      <dt>Sequence</dt>
      <dd className="mono">{sequence}</dd>
      <dt>Length</dt>
      <dd className="mono">{sequence.length} nt</dd>
      <dt>Rules</dt>
      <dd>
        hairpin loops of at least {folding.min_loop} bases; G-U wobble {folding.wobble ? "allowed" : "not allowed"}
      </dd>
      <dt>Maximum base pairs</dt>
      <dd className="mono">{trace.result.score}</dd>
      <dt>Optimal structures</dt>
      <dd className="mono">
        {trace.result.optimal_count}
        {trace.result.structures_truncated && ` (first ${trace.result.structures.length} included in this trace)`}
      </dd>
      {reference && (
        <>
          <dt>Reference</dt>
          <dd>
            {reference.pairs.length} pairs; {reference.source}
          </dd>
        </>
      )}
    </dl>
  );
}

// ------------------------------------------------------------------ status lines

function alignmentStatus(trace: AlignmentTrace, frame: Frame, total: number): string {
  const { phase, cursor, traced } = frame;
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

function foldingStatus(trace: FoldingTrace, frame: Frame, total: number): string {
  const { phase, cursor, traced } = frame;
  const seq = trace.inputs.sequence;
  const h = trace.inputs.folding.min_loop;
  if (phase === "fill") {
    if (!cursor) return "Press Play (or Space) to fill the table one diagonal at a time, shortest intervals first.";
    const [i, j] = cursor;
    const len = j - i + 1;
    if (j - i - 1 < h) {
      return `Interval ${i + 1}-${j + 1} (${len} nt): too short to close a hairpin loop of at least ${h} bases, so 0 pairs.`;
    }
    return `Interval ${i + 1}-${j + 1} (${len} nt): best of "base ${j + 1} (${seq[j]}) unpaired" and "base ${j + 1} pairs with some k". The dashed L is every cell it reads.`;
  }
  if (phase === "traceback") {
    const head = traced[traced.length - 1];
    if (!head) return `Table complete. Traceback starts at the top-right cell, the whole sequence (1-${seq.length}).`;
    return `Traceback ${traced.length} of ${total}: deciding base ${head[1] + 1} within interval ${head[0] + 1}-${head[1] + 1}.`;
  }
  return `Done. ${trace.result.score} base pairs; ${trace.result.optimal_count} optimal structure(s). Dashed cells lie on at least one of them.`;
}

// ------------------------------------------------------------------ structure panel

function StructurePanel({ trace, index, frame, theme, ink }: { trace: FoldingTrace; index: number; frame: Frame; theme: Theme; ink: Ink }) {
  const structure = trace.result.structures[index];
  if (!structure) return null;
  const seq = trace.inputs.sequence;
  const reference = trace.inputs.reference;
  const shown = frame.phase === "fill" ? [] : frame.phase === "done" ? structure.pairs : decidedPairs(structure, frame.traced);
  const head = frame.phase === "traceback" ? frame.traced[frame.traced.length - 1] : undefined;
  const colors = { ...(theme === "dark" ? ARC_DARK : ARC_LIGHT), reference: ink.secondary, text: ink.primary, muted: ink.muted, accent: ink.accent };
  const counts = { GC: 0, AU: 0, GU: 0 };
  for (const [i, j] of structure.pairs) {
    const key = [seq[i], seq[j]].sort().join("");
    counts[key === "CG" ? "GC" : key === "AU" ? "AU" : "GU"] += 1;
  }
  const comparison = reference ? comparePairs(structure.pairs, reference.pairs) : null;
  const pct = (v: number) => `${Math.round(v * 100)}%`;

  return (
    <section className="alignment-panel" aria-label="Structure">
      <h2>Structure</h2>
      {frame.phase === "fill" ? (
        <p className="muted">Base pairs are recovered during traceback, after the table is complete.</p>
      ) : (
        <>
          <ArcDiagram sequence={seq} pairs={shown} reference={reference?.pairs} focus={head ? head[1] : null} colors={colors} />
          <div className="dot-bracket mono" aria-label="Dot-bracket notation">
            <div><span className="db-label">sequence</span>{seq}</div>
            <div><span className="db-label">{frame.phase === "done" ? "predicted" : "so far"}</span>{toDotBracket(seq.length, shown)}</div>
            {reference && <div><span className="db-label">reference</span>{reference.dot_bracket}</div>}
          </div>
          <p className="small legend-keys">
            <span><i className="key" style={{ background: colors.GC }} /> G-C ({counts.GC})</span>
            <span><i className="key" style={{ background: colors.AU }} /> A-U ({counts.AU})</span>
            <span><i className="key" style={{ background: colors.GU }} /> G-U wobble ({counts.GU})</span>
            {reference && <span><i className="key" style={{ background: colors.reference }} /> reference, below the line</span>}
          </p>
          {comparison && frame.phase === "done" && (
            <p>
              This structure shares <strong>{comparison.shared}</strong> of the reference's {reference!.pairs.length} pairs
              (sensitivity {pct(comparison.sensitivity)}), and {comparison.shared} of its own {structure.pairs.length} pairs
              are in the reference (PPV {pct(comparison.ppv)}). Maximising the pair count over-pairs real RNA, which is why
              practical predictors minimise free energy instead (Zuker's algorithm, as in ViennaRNA).
            </p>
          )}
        </>
      )}
    </section>
  );
}

// ------------------------------------------------------------------ tooltip and legend

function CellTooltip({ trace, cell, x, y }: { trace: DPTrace; cell: Cell; x: number; y: number }) {
  const left = Math.min(x + 14, window.innerWidth - 300);
  return (
    <div className="tooltip" style={{ left, top: y + 14 }} role="tooltip">
      {isFolding(trace) ? <FoldingTip trace={trace} cell={cell} /> : <AlignmentTip trace={trace as AlignmentTrace} cell={cell} />}
    </div>
  );
}

function AlignmentTip({ trace, cell: [i, j] }: { trace: AlignmentTrace; cell: Cell }) {
  const { seq_a, seq_b, scoring } = trace.inputs;
  const sub = scoring.substitution;
  const a = i > 0 ? seq_a[i - 1] : null;
  const b = j > 0 ? seq_b[j - 1] : null;
  const s = a && b ? sub.scores[sub.alphabet.indexOf(a)]?.[sub.alphabet.indexOf(b)] : undefined;
  const moves = moveOffsets(trace.projection.moves[i][j]).map(([di, dj]) => (di && dj ? "diagonal" : di ? "up" : "left"));
  const fmt = (v: number | null | undefined) => (v === null || v === undefined ? "unreachable" : String(v));
  return (
    <>
      <div className="mono strong">
        ({i}, {j}) {a ?? "∅"} / {b ?? "∅"}
      </div>
      <div>Prefixes: "{seq_a.slice(0, i)}" vs "{seq_b.slice(0, j)}"</div>
      <div>H = {fmt(trace.projection.scores[i][j])}</div>
      {trace.layers && <div className="mono">{trace.layers.map((l) => `${l.name}=${fmt(l.scores[i][j])}`).join("  ")}</div>}
      {s !== undefined && <div>s({a}, {b}) = {s}</div>}
      <div>Best from: {moves.length ? moves.join(", ") : "none (origin)"}</div>
      <div>{trace.projection.on_optimal_path[i][j] ? "On an optimal path" : "Not on any optimal path"}</div>
    </>
  );
}

function FoldingTip({ trace, cell: [i, j] }: { trace: FoldingTrace; cell: Cell }) {
  const seq = trace.inputs.sequence;
  if (j < i) return <div className="muted">Below the diagonal: not an interval.</div>;
  const bits = trace.projection.moves[i][j];
  const labels = trace.projection.move_labels ?? {};
  const decisions = [1, 4, 2].filter((b) => bits & b).map((b) => labels[String(b)] ?? `move ${b}`);
  const interval = seq.slice(i, j + 1);
  return (
    <>
      <div className="mono strong">
        bases {i + 1}-{j + 1} ({j - i + 1} nt)
      </div>
      <div className="mono">{interval.length > 36 ? `${interval.slice(0, 16)}...${interval.slice(-16)}` : interval}</div>
      <div>Maximum pairs in this interval: {trace.projection.scores[i][j]}</div>
      <div>Best by: {decisions.length ? decisions.join("; ") : "nothing to decide (no pair possible)"}</div>
      {hasSplit(bits, trace.algorithm) && <div className="muted small">The fork mark means the interval splits in two.</div>}
      <div>{trace.projection.on_optimal_path[i][j] ? "On an optimal structure's traceback" : "Not used by any optimal structure"}</div>
    </>
  );
}

function Legend({ trace, layer, palette }: { trace: DPTrace; layer: number; palette: Palette }) {
  const extent = scoreExtent(layerGrid(trace, layer));
  const folding = isFolding(trace);
  // Folding scores are pair counts, never negative: show only the 0..max half.
  const stops = folding
    ? Array.from({ length: 6 }, (_, k) => (extent * k) / 5)
    : Array.from({ length: 11 }, (_, k) => -extent + (2 * extent * k) / 10);
  return (
    <div className="legend" aria-label="Colour legend">
      <div className="ramp">
        {stops.map((v) => (
          <span key={v} style={{ background: scoreColor(v, extent, palette) }} />
        ))}
      </div>
      <div className="ramp-labels mono small">
        {folding ? (
          <>
            <span>0 pairs</span>
            <span>{extent} pairs</span>
          </>
        ) : (
          <>
            <span>-{extent}</span>
            <span>0</span>
            <span>+{extent}</span>
          </>
        )}
      </div>
      <div className="small muted legend-keys">
        {!folding && <span><i className="key key-empty" style={{ background: palette.empty }} /> unreachable</span>}
        <span><i className="key key-optimal" /> {folding ? "used by an optimal structure" : "on an optimal path"}</span>
        <span><i className="key key-trace" /> traceback</span>
      </div>
      {folding && (
        <div className="small muted legend-keys">
          <span>&larr; base j unpaired</span>
          <span>&#8601; i pairs with j</span>
          <span>fork: j pairs inside, splitting the interval</span>
        </div>
      )}
    </div>
  );
}
