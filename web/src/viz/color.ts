// Diverging colour scale for alignment scores.
//
// Scores have a meaningful sign (positive = net similarity, negative = net
// cost), so the scale is diverging: red arm below zero, blue arm above, a
// neutral grey at exactly zero. Interpolation runs in OKLab so equal score
// steps look like equal colour steps. Poles and midpoints come from the
// reference palette in the dataviz guidance (blue/red diverging pair).

export type Rgb = [number, number, number];

export interface Palette {
  surface: string;
  mid: string;
  positive: string;
  negative: string;
  empty: string;
  inkDark: string;
  inkLight: string;
}

export const LIGHT: Palette = {
  surface: "#fcfcfb",
  mid: "#f0efec",
  positive: "#184f95",
  negative: "#b3302f",
  empty: "#e4e3df",
  inkDark: "#000000",
  inkLight: "#ffffff",
};

export const DARK: Palette = {
  surface: "#1a1a19",
  mid: "#383835",
  positive: "#6da7ec",
  negative: "#e66767",
  empty: "#262624",
  inkDark: "#000000",
  inkLight: "#ffffff",
};

export function hexToRgb(hex: string): Rgb {
  const h = hex.replace("#", "");
  return [0, 2, 4].map((k) => parseInt(h.slice(k, k + 2), 16) / 255) as Rgb;
}

export function rgbToHex([r, g, b]: Rgb): string {
  const c = (v: number) => Math.round(Math.min(1, Math.max(0, v)) * 255).toString(16).padStart(2, "0");
  return `#${c(r)}${c(g)}${c(b)}`;
}

const toLinear = (v: number) => (v <= 0.04045 ? v / 12.92 : ((v + 0.055) / 1.055) ** 2.4);
const toGamma = (v: number) => (v <= 0.0031308 ? 12.92 * v : 1.055 * v ** (1 / 2.4) - 0.055);

export function rgbToOklab([r, g, b]: Rgb): Rgb {
  const [lr, lg, lb] = [toLinear(r), toLinear(g), toLinear(b)];
  const l = Math.cbrt(0.4122214708 * lr + 0.5363325363 * lg + 0.0514459929 * lb);
  const m = Math.cbrt(0.2119034982 * lr + 0.6806995451 * lg + 0.1073969566 * lb);
  const s = Math.cbrt(0.0883024619 * lr + 0.2817188376 * lg + 0.6299787005 * lb);
  return [
    0.2104542553 * l + 0.793617785 * m - 0.0040720468 * s,
    1.9779984951 * l - 2.428592205 * m + 0.4505937099 * s,
    0.0259040371 * l + 0.7827717662 * m - 0.808675766 * s,
  ];
}

export function oklabToRgb([L, a, b]: Rgb): Rgb {
  const l = (L + 0.3963377774 * a + 0.2158037573 * b) ** 3;
  const m = (L - 0.1055613458 * a - 0.0638541728 * b) ** 3;
  const s = (L - 0.0894841775 * a - 1.291485548 * b) ** 3;
  return [
    toGamma(4.0767416621 * l - 3.3077115913 * m + 0.2309699292 * s),
    toGamma(-1.2684380046 * l + 2.6097574011 * m - 0.3413193965 * s),
    toGamma(-0.0041960863 * l - 0.7034186147 * m + 1.707614701 * s),
  ];
}

function mix(fromHex: string, toHex: string, t: number): string {
  const a = rgbToOklab(hexToRgb(fromHex));
  const b = rgbToOklab(hexToRgb(toHex));
  return rgbToHex(oklabToRgb([0, 1, 2].map((k) => a[k] + (b[k] - a[k]) * t) as Rgb));
}

/** Largest |score| among reachable cells; the scale is symmetric around zero. */
export function scoreExtent(grid: Array<Array<number | null>>): number {
  let max = 0;
  for (const row of grid) for (const v of row) if (v !== null) max = Math.max(max, Math.abs(v));
  return max || 1;
}

export function scoreColor(value: number | null, extent: number, palette: Palette): string {
  if (value === null) return palette.empty;
  const t = Math.min(1, Math.abs(value) / extent);
  if (value === 0) return palette.mid;
  return mix(palette.mid, value > 0 ? palette.positive : palette.negative, t);
}

/** WCAG relative luminance, used to pick legible text over a fill.
 *
 * Cell numbers use pure black or pure white ink. For any fill, the better of
 * the two is at least 4.58:1 (the worst case is a fill of luminance ~0.18),
 * so every point on the scale clears the WCAG AA 4.5:1 text threshold.
 * Slightly tinted inks would drop that floor to about 4.45:1.
 */
export function luminance(hex: string): number {
  const [r, g, b] = hexToRgb(hex).map(toLinear);
  return 0.2126 * r + 0.7152 * g + 0.0722 * b;
}

export function contrastRatio(a: string, b: string): number {
  const [hi, lo] = [luminance(a), luminance(b)].sort((x, y) => y - x);
  return (hi + 0.05) / (lo + 0.05);
}

export function inkFor(fill: string, palette: Palette): string {
  return contrastRatio(fill, palette.inkDark) >= contrastRatio(fill, palette.inkLight)
    ? palette.inkDark
    : palette.inkLight;
}
