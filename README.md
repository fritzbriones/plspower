# plspower — A Monte Carlo Simulation Planner for PLS-SEM Power Analysis

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24291129.svg)](https://doi.org/10.5281/zenodo.24291129)

**Status: v0.8.2 — models with 1–3 mediators (parallel or serial) and
construct-specific measurement; engine cross-validated against `seminr` 2.6.0
for six model types (all paths agree within 2e-7; power identical — see
[Validation](#validation-protocol-read-before-reporting-anything)).
Built-in smoke test PASSED on v0.8.2 (Nov 24, 2026); seminr cross-validation PASSED on the same engine (v0.8.0, Nov 17, 2026).** Quick-preview results are never
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

## Online app limits (to protect the free server)

The public app runs on a free Posit Connect Cloud plan, where every running hour counts against the
owner's monthly compute credits, and a long run is lost if the browser tab closes. So the online app:

- offers the **Quick** and **Standard** tiers;
- allows Standard runs of up to about **3 hours** of estimated runtime;
- runs the **Publication** tier only on your own computer, with the same engine and identical results:

```bash
# download app.R and plspower_batch.R into one folder, then:
Rscript plspower_batch.R --tier publication --mode find --n-min 50 --n-max 1000
```

Running the app locally (`shiny::runApp("app.R")` in RStudio) has no limits. To force a mode, set the
environment variable `PLSPOWER_HOSTED=true` or `false`.

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

Version history: see [NEWS.md](NEWS.md).

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
   **Status: PASSED, v0.8.2, Nov 24, 2026** (local R run: all 7 checks pass;
   v0.8.0 also passed on Nov 17, 2026, minimum N 254 vs reference 254). Last confirmed on Posit Connect
   Cloud with v0.7.1 (Nov 17, 2026); the estimation engine is unchanged
   in v0.8.0.
2. **Independent cross-validation against seminr**
   (`Rscript validation/cross_validate_seminr.R`): the built-in engine and
   `seminr::estimate_pls()` (Mode A, path weighting) fit the same simulated
   datasets. Acceptance rules: every structural path within 1e-4 and power
   for every effect within ±0.03 at every N.
   **Status: PASSED (v0.8.0 vs seminr 2.6.0, Nov 17, 2026)** — 6 models
   (simple; parallel with 2 and 3 mediators; serial with 2 and 3 mediators;
   equal and unequal loadings/indicator counts) × N = 100, 200, 300 × 250
   datasets: largest path difference 1.3e-7, largest power difference 0.000
   across 51 effect-by-condition cells. The v0.4.1 engine fails the simple
   check (difference 0.09). Results: `validation/results/`.
   **What this establishes:** the PLS estimation (weights, scores and path
   coefficients) matches seminr. The significance tests (t test, delta-method
   z) are applied identically to both sets of estimates, so the matching
   power confirms the estimates, not the tests themselves; percentile-bootstrap
   power is not part of this comparison. Independent checks of the inference
   (Type I error, coverage, bootstrap) are planned for v0.9.
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

The approved v0.9.0 plan (tests, independent validation, then features) is in
[ROADMAP.md](ROADMAP.md).

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

## References

Aguirre-Urreta, M. I., & Rönkkö, M. (2015). Sample size determination and
statistical power analysis in PLS using R: An annotated tutorial.
*Communications of the Association for Information Systems, 36*, 33–51.
https://doi.org/10.17705/1CAIS.03603

Ansani, A., & Rinallo, E. (2025). *PLS-SEM-power: A Shiny app and R package
for computing required sample size and minimum detectable effect size in
PLS-SEMs* [Preprint]. arXiv. https://arxiv.org/abs/2511.14546

Cohen, J. (1988). *Statistical power analysis for the behavioral sciences*
(2nd ed.). Lawrence Erlbaum Associates.

Efron, B., & Tibshirani, R. J. (1993). *An introduction to the bootstrap*.
Chapman & Hall.

Faul, F., Erdfelder, E., Buchner, A., & Lang, A.-G. (2009). Statistical power
analyses using G*Power 3.1: Tests for correlation and regression analyses.
*Behavior Research Methods, 41*(4), 1149–1160.
https://doi.org/10.3758/BRM.41.4.1149

Faul, F., Erdfelder, E., Lang, A.-G., & Buchner, A. (2007). G*Power 3: A
flexible statistical power analysis program for the social, behavioral, and
biomedical sciences. *Behavior Research Methods, 39*(2), 175–191.
https://doi.org/10.3758/BF03193146

Hair, J. F., Hult, G. T. M., Ringle, C. M., & Sarstedt, M. (2022). *A primer
on partial least squares structural equation modeling (PLS-SEM)* (3rd ed.).
Sage.

Kock, N., & Hadaya, P. (2018). Minimum sample size estimation in PLS-SEM: The
inverse square root and gamma-exponential methods. *Information Systems
Journal, 28*(1), 227–261. https://doi.org/10.1111/isj.12131

Lohmöller, J.-B. (1989). *Latent variable path modeling with partial least
squares*. Physica-Verlag.

Muthén, L. K., & Muthén, B. O. (2002). How to use a Monte Carlo study to
decide on sample size and determine power. *Structural Equation Modeling,
9*(4), 599–620. https://doi.org/10.1207/S15328007SEM0904_8

Preacher, K. J., & Hayes, A. F. (2008). Asymptotic and resampling strategies
for assessing and comparing indirect effects in multiple mediator models.
*Behavior Research Methods, 40*(3), 879–891.
https://doi.org/10.3758/BRM.40.3.879

Ray, S., Danks, N. P., & Calero Valdez, A. (2026). *seminr: Building and
estimating structural equation models* (Version 2.6.0) [R package].
https://CRAN.R-project.org/package=seminr

Sobel, M. E. (1982). Asymptotic confidence intervals for indirect effects in
structural equation models. *Sociological Methodology, 13*, 290–312.
https://doi.org/10.2307/270723

## License

MIT License. See LICENSE.

## How to cite

If you use plspower in a thesis, article or report, please cite:

> Briones, F. J. (2026). *plspower: A Monte Carlo simulation planner for
> PLS-SEM power analysis* (Version 0.8.2) [Computer software]. Zenodo.
> https://doi.org/10.5281/zenodo.24291129

In-text: (Briones, 2026). GitHub's **"Cite this repository"** button
(from `CITATION.cff`) gives APA and BibTeX.

```bibtex
@software{briones2026plspower,
  author  = {Briones, Fritz Jeffry},
  title   = {plspower: A Monte Carlo simulation planner for PLS-SEM power analysis},
  year    = {2026},
  version = {0.8.2},
  publisher = {Zenodo},
  doi     = {10.5281/zenodo.24291129},
  url     = {https://doi.org/10.5281/zenodo.24291129}
}
```
