# plspower — A Monte Carlo Simulation Planner for PLS-SEM Power Analysis

**Status: v0.6.0 — models with 1–3 mediators (parallel or serial) and
construct-specific measurement; engine cross-validated against `seminr` 2.6.0
for six model types (all paths agree within 2e-7; power identical — see
[Validation](#validation-protocol-read-before-reporting-anything)).
Built-in smoke test PASSED on Posit Connect Cloud (Nov 16, 2026).** Quick-preview results are never
publication-grade by design.

**Created by Fritz Jeffry Briones.** If you use plspower, please
[cite it](#how-to-cite).

## 🚀 Live Application

**Try the app:** [Launch plspower on Posit Connect Cloud](https://01a1dfae-8fd9-812a-b8f2-f8ffdfb1f80c.share.connect.posit.cloud/)


## What it does

`plspower` is a **prospective** power-analysis planner for PLS-SEM studies.
You specify a mediation model — X → 1 to 3 mediators → Y (parallel or
serial) plus the direct effect X → Y — with measurement quality (loadings
and indicator counts, per construct if you like) and an inference rule; the
app runs a seed-controlled Monte Carlo simulation in the background and
reports:

- power for the **direct effect**, the **total indirect effect**, each
  **specific indirect effect** (parallel) or the **serial indirect effect**,
- the power curve over a sample-size grid (normal-theory / Sobel and
  percentile-bootstrap PLS-SEM inference),
- the minimum N for your target power (80% or 90%) — use **Find my N** to
  search a wide range (e.g. 50–1,000) automatically, or set your own N grid,
- and a **comparison** against two common planning shortcuts:
  the G*Power-equivalent regression calculation (exact noncentral-F) and
  the Kock–Hadaya (2018) inverse square-root rule (direct effect),
- a **ready-to-paste methods paragraph** (sample-size justification with
  references) you can download for your manuscript.

Differences between methods are reported in neutral wording as a
**"planning-method discrepancy"** — this tool does not claim that any
method is wrong. Results apply **only** to the model and settings you
specify; the 10-times rule is not a power criterion.

**Note on the inverse square-root rule:** it is computed **two-tailed**,
consistent with this app's two-sided CI inference. Kock & Hadaya's (2018)
published tables are **one-tailed**, so the same inputs give different
numbers (e.g., b = 0.20, α = 0.05, power = 0.80 → 155 one-tailed vs
197 two-tailed).

## What's different from existing tools

Existing PLS-SEM planning aids each implement a single method: the
inverse square-root rule (Kock & Hadaya, 2018; the PLS-SEM-power R
package/Shiny app), G*Power regression recipes (Cohen/Faul), or scripted
Monte Carlo simulation (matrixpls; Aguirre-Urreta & Rönkkö, 2015).
`plspower` is not a new estimator. Its contribution is a transparent
**comparison and calibration layer**: the same study, planned three ways,
with measurement quality (loadings, indicator counts) as an explicit
design input.

Why the methods disagree: regression shortcuts and the inverse
square-root rule implicitly assume error-free composites. With imperfect
indicators, measurement-error attenuation reduces the effective effect
size, so simulation-based planning typically requires a larger N than
the shortcuts suggest — and the size of the gap depends on your
measurement quality, which the shortcuts cannot see. `plspower` makes
that cost visible.

## Quick start

```r
# Requires R with shiny >= 1.8.1 and bslib >= 0.7.0 (for ExtendedTask):
install.packages(c("shiny", "bslib", "ggplot2", "future",
                   "promises", "jsonlite"))

shiny::runApp("app.R")
```

Or deploy `app.R` as a single-file Shiny app to shinyapps.io
(`rsconnect::deployApp(appDir = ".", appFiles = "app.R")`).

## What's new in v0.8.0 (exact values and c limit)

- **Type exact values:** a new tick box under "Standardized paths" swaps the
  sliders for number boxes with 0.01 steps (for example a = 0.13, c = 0.12),
  with the same ranges: paths 0 to 0.70, c 0 to 0.50, loading 0.50 to 0.95.
  Typed values outside the range are rejected with a message.
- **Limit for c shown live:** when the paths into Y are large, the model note
  says how high the direct effect c can go (for example "With these paths,
  c can be at most 0.37"), including when a setting is blocked.
- Per-construct loadings now take 0.01 steps.
- The engine is unchanged, so earlier results and the smoke test still apply.

## What's new in v0.7.1 (phone display fixes)

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

## What's new in v0.7.0

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

### How long will it take?

Typical **Find my N** run (scan 50–1,000), simple mediation with 3
indicators per construct, measured on a 2-core server. Your model and server
change these numbers — the app shows the estimate for your settings before
you run.

| Tier | N about 150 | N about 300 | N about 600 |
|---|---|---|---|
| Quick (preview) | ~2 min | ~3 min | ~6 min |
| Standard (report) | ~2.5 h | ~5 h | ~10 h |
| Publication | ~12 h | ~1 day | ~2 days |

Explore with Quick on any device; run Standard or Publication on a computer
and keep the page open until it finishes.

## What's new in v0.6.1 (mobile and usability fixes)

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

## What's new in v0.6.0

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

## What's new in v0.5.0

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

## Model scope (v1, deliberately narrow)

- X → 1–3 mediators (parallel or serial) → Y, with direct effect X → Y
- Reflective (Mode A) constructs: one setting for all (2–8 indicators,
  loading 0.50–0.95) or per construct (1–12 indicators, loading 0.30–0.97)
- Direct, total indirect, specific indirect or serial indirect power; two-sided inference
- Normal-theory and percentile-bootstrap tests (bootstrap power is
  computed over valid bootstrap decisions only)
- Single built-in Mode A engine (cross-validated against seminr) — no engine selection in-app
- No user-data fitting in v1 — prospective power only

## Validation protocol (read before reporting anything)

1. **Built-in smoke test** (Methods & Validation tab): reproduces a
   seed-2026 reference fixture (power computed with seminr 2.6.0 on the
   same datasets, from v0.7.1) with explicit acceptance criteria
   (power: within 2.5 combined Monte Carlo standard errors;
   f²: 1e-10 against a full-precision reference; inverse-sqrt: exact
   integer). At 400 Monte Carlo replications this verifies
   *implementation agreement*, not publication-grade precision.
   **Status: PASSED on Posit Connect Cloud, v0.7.0, Nov 17, 2026** (re-run after deploying v0.7.1).
2. **Independent cross-validation against seminr**
   (`Rscript validation/cross_validate_seminr.R`): the built-in engine and
   `seminr::estimate_pls()` (Mode A, path weighting) fit the same simulated
   datasets. Acceptance rules: every structural path within 1e-4 and power
   for every effect within ±0.03 at every N.
   **Status: PASSED (v0.6.0 vs seminr 2.6.0, Nov 16, 2026)** — 6 models
   (simple; parallel with 2 and 3 mediators; serial with 2 and 3 mediators;
   equal and unequal loadings/indicator counts) × N = 100, 200, 300 × 250
   datasets: largest path difference 1.3e-7, largest power difference 0.000
   across 51 effect-by-condition cells. The v0.4.1 engine fails the simple
   check (difference 0.09). Results: `validation/results/`.
3. **Publication settings**: the Publication tier (5,000 MC replications ×
   2,000 bootstrap resamples); for a final reported N, repeat with at least
   3 seeds. Quick-tier results are never citable.

Exported JSON contains seed, all settings, package versions, and results,
so every run is reproducible.

## Known limitations

- The built-in engine agrees with seminr for all supported models; it has
  not been compared with cSEM or SmartPLS.
- One X and one Y; up to three mediators; no moderators, formative or
  higher-order constructs, or multi-group analysis.
- No custom path diagrams: models are built from the mediator count and
  structure options.
- Percentile bootstrap only (BCa planned).
- Results apply only to the specified model, settings, and inference rule.

## Roadmap

- [x] First hosted run + smoke test PASS (Posit Cloud, Nov 15, 2026)
- [x] `seminr` cross-validation via separate validation script
      (±0.03 agreement rule) — PASSED, v0.5.0
- [x] Indirect-effect (a × b) power — v0.5.0
- [x] Construct-specific loadings and indicator counts — v0.6.0
- [x] 2–3 mediators, parallel or serial — v0.6.0
- [ ] Moderation (interaction terms)
- [ ] BCa bootstrap
- [ ] Publication-tier fixture (5,000 MC × 2,000 boot, ≥3 seeds)
- [ ] Condition sweep: loadings, indicator counts, path sizes,
      inference rules (power, Type I error, bias, N-at-80%)
- [ ] testthat suite in this repo; modular package structure; renv.lock

## Related work

Kock & Hadaya (2018) inverse square-root method; `matrixpls` Monte Carlo
simulation; the **PLS-SEM-power** R package and Shiny app (inverse
square-root approach); `semPower` for covariance-based SEM; G*Power.
This tool's intended contribution is a transparent **comparison and
calibration layer** across planning methods, not a new estimator.

## License

MIT License. See LICENSE.

## How to cite

If you use plspower in a thesis, article or report, please cite:

> Briones, F. J. (2026). *plspower: A Monte Carlo simulation planner for
> PLS-SEM power analysis* (Version 0.8.0) [Computer software]. GitHub.
> https://github.com/fritzbriones/plspower

In-text: (Briones, 2026). GitHub's **"Cite this repository"** button
(from `CITATION.cff`) gives APA and BibTeX.

```bibtex
@software{briones2026plspower,
  author  = {Briones, Fritz Jeffry},
  title   = {plspower: A Monte Carlo simulation planner for PLS-SEM power analysis},
  year    = {2026},
  version = {0.8.0},
  url     = {https://github.com/fritzbriones/plspower}
}
```
