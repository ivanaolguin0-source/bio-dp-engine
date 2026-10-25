import { defineConfig } from "vitest/config";
import react from "@vitejs/plugin-react";

export default defineConfig({
  plugins: [react()],
  // The DPTrace schema is imported from ../contracts, outside this folder.
  server: { fs: { allow: [".."] } },
  test: { include: ["src/**/*.test.ts"] },
});
