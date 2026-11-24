# plspower news

Version history. The current version and its validation status are in the [README](README.md).

## v0.8.2 (November 24, 2026): hosted-server limits

- **Publication tier is run offline when the app is hosted.** On a public server (Posit Connect Cloud,
  shinyapps.io) the tier menu offers Quick and Standard; a note explains how to run the
  Publication tier with `plspower_batch.R` on your own computer, with the exact command.
- **Standard runs on the hosted app are limited to about 3 hours** of estimated runtime. Longer runs
  show a message with the two-step tip (Quick on a wide range, then Standard on a narrow range) or the batch runner.
- Hosted mode is detected automatically (non-interactive R) and can be set with the environment
  variable `PLSPOWER_HOSTED=true/false`. Run locally in RStudio, every tier is available as before.
- **No change to the estimation engine, the tiers or any power result;** the smoke test and the
  seminr cross-validation still apply.
- Maintainer email updated; the test workflow renamed to `tests.yml`; version history moved here.

## v0.8.1 (export fix and validation wording)

- **JSON export fix:** when Find my N does not reach the target power, only
  the fast scan runs (200 replications, no bootstrap). The exported settings
  now record those values, with the requested tier kept separately under
  `requested`. Earlier versions recorded the tier's values instead.
- Validation documentation now states exactly what the seminr check
  establishes (see Validation below).
- No change to the estimation engine or to any power result.

## v0.8.0 (exact values and c limit)

- **Type exact values:** a new tick box under "Standardized paths" swaps the
  sliders for number boxes with 0.01 steps (for example a = 0.13, c = 0.12),
  with the same ranges: paths 0 to 0.70, c 0 to 0.50, loading 0.50 to 0.95.
  Typed values outside the range are rejected with a message.
- **Limit for c shown live:** when the paths into Y are large, the model note
  says how high the direct effect c can go (for example "With these paths,
  c can be at most 0.37"), including when a setting is blocked.
- Per-construct loadings now take 0.01 steps.
- The engine is unchanged, so earlier results and the smoke test still apply.

## v0.7.1 (phone display fixes)

- Formulas split into shorter lines and scrollable, so nothing is cut off
  on a phone; three stray "\;" symbols removed from the text.
- Result box, methods text and tables wrap or scroll on small screens.
- Power curve zooms to the sample sizes that were analysed, with a title
  and legend that fit a phone screen.
- Diagnostics show whole numbers for N; the model note reads
  "True indirect effect: a x b = 0.165".
- **Smoke-test reference updated:** now computed with seminr 2.6.0 on the
  same seed-2026 datasets, so the test checks the current engine (the app
  matches it exactly). The old reference came from the v0.4.1 engine,
  before the path-weighting fix.

## v0.7.0

- **Title bar** at the top of the app: *plspower — Sample-size and power
  planner for PLS-SEM*, version and author.
- **Progress bar** during a run: replications done, percent and time left
  (Find my N shows step 1 "scanning" and step 2 with the bar).
- **"How long will it take?"** guide in the app: Quick / Standard /
  Publication times for *your* model on *this* server, for a needed N of
  about 150, 300 and 600.
- **Publication tier 5× faster**: 5,000 Monte Carlo replications × 2,000
  bootstrap resamples (was 10,000; 2,000 resamples is ample for a
  significance decision at α = .05 or .01).
- **All formulas shown as rendered math** in the Methods & Validation tab
  (data generation, PLS-SEM estimation, tests, power, minimum N, G*Power
  and inverse square-root comparisons).

## v0.6.1 (mobile and usability fixes)

- **Live status while a run is going**: elapsed time, the estimate, and a
  reminder to keep the page open. Errors are now shown instead of failing
  silently. The smoke test shows a "running" message too.
- **Download buttons appear only after a run** (tapping them before gave an
  error), with a short explanation of the CSV and the JSON run record.
- **Correct CPU count on cloud hosting**: uses `parallelly::availableCores()`,
  which respects container limits. `detectCores()` could report the host's
  cores, so runs took much longer than the estimate.
- **Quick tier is 5× faster**: 200 MC × 100 bootstrap resamples (was 500);
  still a preview — use Standard or Publication for reporting.
- **Phones can reconnect** to the same session after a short interruption
  (switching apps, screen lock).

**Using plspower on a phone:** keep the page open and the screen on until
the results appear; use Quick tier to explore and run Standard or
Publication tier on a computer.

## v0.6.0

- **1–3 mediators**, **parallel** (X → each M → Y) or **serial**
  (X → M1 → M2 → M3 → Y), always with the direct path X → Y.
- **Different loading and number of indicators per construct** (or one
  setting for all).
- **Effects**: direct, total indirect, specific indirect via each mediator,
  and serial indirect; tested with the delta-method z test (= Sobel for a
  single a × b; verified against the empirical SD) and the percentile bootstrap.
- Model check in the app: impossible path combinations (explained variance
  ≥ 100%) are blocked with a clear message; R² of Y and the true indirect
  effects are shown.
- G*Power-equivalent comparison generalised to any of these models.
- **Fix:** minimum-N interpolation returned the lower grid point when power
  reached exactly 100% at the next point (present since v0.4.1).
- Cross-check against seminr extended to six model types (simple, parallel
  with 2–3 mediators, serial with 2–3 mediators; equal and unequal measurement).
- Simple mediation gives bit-for-bit the same results as v0.5.0.

## v0.5.0

- **About 10× faster engine**, numerically identical to v0.4.2 (max
  difference ~1e-15). Grid points run in parallel; the default Quick tier
  now takes minutes instead of hours.
- **Find my N**: enter a minimum and maximum N (e.g. 50–1,000). The app
  scans the range quickly (log-spaced, 200 replications, no bootstrap), then
  runs the full tier only on a fine grid where power reaches your target —
  large N ranges in minutes instead of hours. "Custom N grid" is still
  available.
- **Indirect-effect power** (a × b): Sobel test and percentile bootstrap,
  with minimum N and its own power curve. Choose "Target effect" in the app.
- **Cross-validated against seminr** (`validation/cross_validate_seminr.R`):
  identical estimates on the same data and identical power curves.
- **Methods text export**: a ready-to-paste sample-size justification
  paragraph with references (download as TXT).
- **How to cite** card, citation in the footer, plot caption, JSON export
  and methods text; `CITATION.cff` so GitHub shows "Cite this repository".
- Reproducible parallel runs: each N uses its own L'Ecuyer-CMRG stream
  derived from the seed, so results are identical with 1 or many cores.

### v0.4.2

- Path weighting now includes the X → Y link (v0.4.1 ignored it; its
  direct-effect estimates differed from seminr by up to 0.09).
- Honest runtime estimate shown before running; validated N-grid input;
  background runs; package versions in the JSON export.
