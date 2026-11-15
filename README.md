# plspower — A Monte Carlo Simulation Planner for PLS-SEM Power Analysis

**Status: v0.3 — code-complete, engine validation pending.** This is a
research-planning tool under active development. The simulation engine has
not yet been executed on a hosted R runtime; run the built-in smoke test
(see below) before using any output. Quick-preview results are never
publication-grade by design.

## What it does

`plspower` is a **prospective** power-analysis planner for PLS-SEM studies.
You specify a mediation model (X → M → Y plus direct effect X → Y),
measurement quality, and inference rule; the app runs a seed-controlled
Monte Carlo simulation in the background and reports:

- the power curve over a sample-size grid (normal-theory and
  percentile-bootstrap PLS-SEM inference),
- the minimum N for your target power (80% or 90%),
- and a **comparison** against two common planning shortcuts:
  the G*Power-equivalent regression calculation (exact noncentral-F) and
  the Kock–Hadaya (2018) inverse square-root rule.

Differences between methods are reported in neutral wording as a
**"planning-method discrepancy"** — this tool does not claim that any
method is wrong. Results apply **only** to the model and settings you
specify; the 10-times rule is not a power criterion.

## Quick start

```r
# Requires R with:
install.packages(c("shiny", "bslib", "ggplot2", "future",
                   "promises", "jsonlite"))   # seminr is optional

shiny::runApp("app.R")
```

Or deploy `app.R` as a single-file Shiny app to shinyapps.io
(`rsconnect::deployApp(appDir = ".", appFiles = "app.R")`).

## Model scope (v1, deliberately narrow)

- One mediation template: X → M → Y with direct effect X → Y
- Three reflective (Mode A) constructs, 3–6 indicators each
- Direct-effect power analysis; two-sided 95% CI inference
- Normal-theory and percentile-bootstrap tests
- No user-data fitting in v1 — prospective power only

## Validation protocol (read before reporting anything)

1. **Built-in smoke test** (Methods & Validation tab): reproduces a
   seed-2026 reference fixture with explicit acceptance criteria
   (power: ±5% relative; f²: 1e-6; inverse-sqrt: exact integer).
   At 400 Monte Carlo replications this verifies *implementation
   agreement*, not publication-grade precision.
2. **Dual-engine cross-check** (when `seminr` is installed): the built-in
   Mode A reference engine and `seminr::estimate_pls()` should produce
   power curves within ±0.03 absolute at every N. Engines differ by
   design (RNG paths, scaling, stopping rules), so exact agreement is
   not expected.
3. **Publication settings**: 5,000–10,000 MC replications × 10,000
   bootstrap resamples, at least 3 seeds. Only results at this tier are
   citable, and only as preliminary until independently validated.

Exported JSON contains seed, all settings, package versions, and results,
so every run is reproducible.

## Known limitations

- The built-in PLS-SEM engine is a hand-rolled Wold/Lohmöller Mode A
  implementation with path weighting; it has **not yet been cross-validated
  against `seminr`/`cSEM` at production Monte Carlo settings**.
- One model family only (v1); no formative constructs, multi-group
  analysis, higher-order constructs, or custom path diagrams.
- Percentile bootstrap only (BCa planned).
- Results apply only to the specified model, settings, and inference rule.

## Roadmap

- [ ] First hosted run + smoke test PASS (gate for all further claims)
- [ ] `seminr` cross-validation (±0.03 agreement rule)
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

## Citation (temporary)

Briones, F. J. (2026). *plspower: A Monte Carlo simulation planner for
PLS-SEM power analysis* (Version 0.3) [Computer software]. GitHub.
https://github.com/fritzbriones/plspower