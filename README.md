# plspower — A Monte Carlo Simulation Planner for PLS-SEM Power Analysis

**Status: v0.6.0 — models with 1–3 mediators (parallel or serial) and
construct-specific measurement; engine cross-validated against `seminr` 2.6.0
for six model types (all paths agree within 2e-7; power identical — see
[Validation](#validation-protocol-read-before-reporting-anything)).
Built-in smoke test PASSED (Nov 16, 2026).** Quick-preview results are never
publication-grade by design.

**Created by Fritz Jeffry Briones.** If you use plspower, please
[cite it](#how-to-cite).

## 🚀 Live Application

**Try the app:** [Launch plspower on Posit Connect Cloud](https://01a1dfae-8fd9-812a-b8f2-f8ffdfb1f80c.share.connect.posit.cloud/)


## What it does

`plspower` is a **prospective** power-analysis planner for PLS-SEM studies.
You specify a mediation model — X → 1 to 3 mediators → Y (parallel or
serial) plus the direct effect X → Y — with loadings and indicator counts
per construct,
measurement quality, and inference rule; the app runs a seed-controlled
Monte Carlo simulation in the background and reports:

- power for the **direct effect**, the **total indirect effect**, each
  **specific indirect effect** (parallel) or the **serial indirect effect**,
- the power curve over a sample-size grid (normal-theory / Sobel and
  percentile-bootstrap PLS-SEM inference),
- the minimum N for your target power (80% or 90%),
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
- Reflective (Mode A) constructs, 1–12 indicators each, loadings 0.30–0.97
- Direct, total indirect, specific indirect or serial indirect power; two-sided inference
- Normal-theory and percentile-bootstrap tests (bootstrap power is
  computed over valid bootstrap decisions only)
- Single built-in Mode A engine (cross-validated against seminr) — no engine selection in-app
- No user-data fitting in v1 — prospective power only

## Validation protocol (read before reporting anything)

1. **Built-in smoke test** (Methods & Validation tab): reproduces a
   seed-2026 reference fixture with explicit acceptance criteria
   (power: within 2.5 combined Monte Carlo standard errors;
   f²: 1e-10 against a full-precision reference; inverse-sqrt: exact
   integer). At 400 Monte Carlo replications this verifies
   *implementation agreement*, not publication-grade precision.
   **Status: PASSED on Posit Cloud (R 4.6.2), Nov 15, 2026.**
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
3. **Publication settings**: 5,000–10,000 MC replications × 10,000
   bootstrap resamples, at least 3 seeds. Only results at this tier are
   citable.

Exported JSON contains seed, all settings, package versions, and results,
so every run is reproducible.

## Known limitations

- The built-in engine agrees with seminr for this model; it has not been
  compared with cSEM or SmartPLS.
- One X and one Y; up to three mediators; no moderators, formative or
  higher-order constructs, or multi-group analysis.
- One model family only (v1); no formative constructs, multi-group
  analysis, higher-order constructs, or custom path diagrams.
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
- [ ] Publication-tier fixture (5,000–10,000 MC × 10,000 boot, ≥3 seeds)
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
> PLS-SEM power analysis* (Version 0.6.0) [Computer software]. GitHub.
> https://github.com/fritzbriones/plspower

In-text: (Briones, 2026). GitHub's **"Cite this repository"** button
(from `CITATION.cff`) gives APA and BibTeX.

```bibtex
@software{briones2026plspower,
  author  = {Briones, Fritz Jeffry},
  title   = {plspower: A Monte Carlo simulation planner for PLS-SEM power analysis},
  year    = {2026},
  version = {0.6.0},
  url     = {https://github.com/fritzbriones/plspower}
}
```
