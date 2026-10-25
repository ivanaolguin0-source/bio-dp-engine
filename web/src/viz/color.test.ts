import { describe, expect, it } from "vitest";
import { contrastRatio, DARK, hexToRgb, inkFor, LIGHT, oklabToRgb, rgbToHex, rgbToOklab, scoreColor, scoreExtent } from "./color";

describe("OKLab round trip", () => {
  it.each(["#184f95", "#b3302f", "#f0efec", "#000000", "#ffffff"])("%s", (hex) => {
    expect(rgbToHex(oklabToRgb(rgbToOklab(hexToRgb(hex))))).toBe(hex);
  });
});

describe("diverging scale", () => {
  it("is neutral at zero and hits the poles at the extent", () => {
    expect(scoreColor(0, 5, LIGHT)).toBe(LIGHT.mid);
    expect(scoreColor(5, 5, LIGHT)).toBe(LIGHT.positive);
    expect(scoreColor(-5, 5, LIGHT)).toBe(LIGHT.negative);
    expect(scoreColor(null, 5, LIGHT)).toBe(LIGHT.empty);
  });
  it("uses a symmetric extent and ignores unreachable cells", () => {
    expect(scoreExtent([[0, -7, null], [3, 2, 1]])).toBe(7);
    expect(scoreExtent([[0]])).toBe(1);
  });
  it.each([["light", LIGHT], ["dark", DARK]] as const)("cell numbers stay legible in %s mode", (_, palette) => {
    for (let v = -20; v <= 20; v++) {
      const fill = scoreColor(v, 20, palette);
      expect(contrastRatio(fill, inkFor(fill, palette))).toBeGreaterThanOrEqual(4.5);
    }
  });
});
