# plspower — A Monte Carlo Simulation Planner for PLS-SEM Power Analysis

**Status: v0.4.1 — smoke test PASSED on a hosted R runtime (Posit Cloud,
Nov 15, 2026).** This is a research-planning tool under active development.
The built-in engine is the single supported estimator; independent
cross-validation against established PLS-SEM software is planned as a
separate validation exercise. Quick-preview results are never
publication-grade by design.

## 🚀 Live Application

**Try the app:** [Launch plspower on Posit Connect Cloud](https://01a1dfae-8fd9-812a-b8f2-f8ffdfb1f80c.share.connect.posit.cloud/)

**Source code:** [GitHub repository](https://github.com/fritzbriones/plspower)

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
# Requires R with:
install.packages(c("shiny", "bslib", "ggplot2", "future",
                   "promises", "jsonlite"))

shiny::runApp("app.R")
```

Or deploy `app.R` as a single-file Shiny app to shinyapps.io
(`rsconnect::deployApp(appDir = ".", appFiles = "app.R")`).

## Model scope (v1, deliberately narrow)

- One mediation template: X → M → Y with direct effect X → Y
- Three reflective (Mode A) constructs, 3–6 indicators each
- Direct-effect power analysis; two-sided 95% CI inference
- Normal-theory and percentile-bootstrap tests (bootstrap power is
  computed over valid bootstrap decisions only)
- Single built-in Mode A engine — no engine selection in-app
- No user-data fitting in v1 — prospective power only

## Validation protocol (read before reporting anything)

1. **Built-in smoke test** (Methods & Validation tab): reproduces a
   seed-2026 reference fixture with explicit acceptance criteria
   (power: within 2.5 combined Monte Carlo standard errors;
   f²: 1e-10 against a full-precision reference; inverse-sqrt: exact
   integer). At 400 Monte Carlo replications this verifies
   *implementation agreement*, not publication-grade precision.
   **Status: PASSED on Posit Cloud (R 4.6.2), Nov 15, 2026.**
2. **Independent cross-validation** (planned, separate script): the
   built-in engine vs `seminr`/`cSEM` power curves within ±0.03 absolute
   at every N, under a prespecified protocol. Neither package is an
   in-app estimation backend.
3. **Publication settings**: 5,000–10,000 MC replications × 10,000
   bootstrap resamples, at least 3 seeds. Only results at this tier are
   citable, and only as preliminary until independently validated.

Exported JSON contains seed, all settings, package versions, and results,
so every run is reproducible.

## Known limitations

- The built-in PLS-SEM engine is a hand-rolled Wold/Lohmöller Mode A
  implementation with path weighting; it has **not yet been independently
  cross-validated against `seminr`/`cSEM`** under a prespecified protocol.
- One model family only (v1); no formative constructs, multi-group
  analysis, higher-order constructs, or custom path diagrams.
- Percentile bootstrap only (BCa planned).
- Results apply only to the specified model, settings, and inference rule.

## Roadmap

- [x] First hosted run + smoke test PASS (Posit Cloud, Nov 15, 2026)
- [ ] `seminr`/`cSEM` cross-validation via separate validation script
      (±0.03 agreement rule)
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
PLS-SEM power analysis* (Version 0.4.1) [Computer software]. GitHub.
https://github.com/fritzbriones/plspower
