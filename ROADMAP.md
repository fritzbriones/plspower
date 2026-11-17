# plspower v0.9.0 development plan

Status: **approved, not started** (finalized November 17, 2026).
Current release: v0.8.1 (live app). Citable version: v0.8.0,
DOI [10.5281/zenodo.24291129](https://doi.org/10.5281/zenodo.24291129).

Guiding rule: **build it, test it, validate it independently, then expand it.**
Features are added only after validation shows the current engine behaves
correctly. If validation finds a problem (for example in the delta method or
the minimum-N search), that is fixed before any new feature.

Each stage is delivered as small pull requests (about 4–6 in total), each
reviewed and merged by the author.

---

## Stage 1: v0.9.0-alpha (code safety)

Purpose: make sure the code does what it is intended to do, and keeps doing so.
This checks the code, not the statistical method (that is Stage 2).

Tests (`tests/testthat/`):

| Area | What is checked |
|---|---|
| Data generation | Simulated latent and indicator correlations match the model (large N, fixed tolerance) |
| PLS estimation | Path estimates match seminr on fixed datasets (tolerance 1e-4) |
| Admissibility | Impossible models are blocked with a clear message (for example R² of Y ≥ 1) |
| Input limits | Out-of-range or empty inputs are rejected; negative paths are rejected until Stage 3 adds them |
| Boundary cases | Zero effects, loadings 0.30 and 0.97, 1 and 12 indicators |
| Reproducibility | Same inputs and seed give identical results |
| Minimum N | Interpolation is correct, including power of exactly 0 or 1 |
| Exports | JSON records the settings that actually ran (including the Find-my-N scan-only case) |
| Error handling | Invalid inputs produce messages, not crashes (shiny::testServer) |

Infrastructure:

- GitHub Actions workflow runs all tests on every pull request.
- `renv.lock` records exact package versions.
- The existing smoke test and seminr cross-validation stay as they are.

Done when: all tests pass in GitHub Actions on a clean machine.

---

## Stage 2: v0.9.0-beta (independent statistical validation)

Purpose: show that the statistical procedures behave correctly under known
data-generating conditions. This is the main evidence for the software paper.

Measures:

| Measure | Definition | Target |
|---|---|---|
| Type I error | Rejection rate when the true effect is 0 | ≈ .05 (10,000 replications, MCSE ≈ .002) |
| Bias | Mean and relative bias of path and indirect estimates | Small and decreasing with N |
| RMSE | Root mean squared error of estimates | Decreasing with N |
| Coverage | Share of 95% percentile-bootstrap intervals containing the true value | ≈ .95 |
| Bootstrap power | Agreement with an independent implementation (seminr bootstrap) under identical data, N, resamples and α | Agreement within Monte Carlo error; differences documented |
| Delta method vs bootstrap | Standard errors, rejection rates and power for direct, specific indirect, serial indirect and total indirect effects | Size of the cross-equation approximation error reported, especially for serial models |
| Minimum-N accuracy | Fresh, independent simulation at the recommended N and just below it | Power ≈ target at N, below target just under N |
| Convergence | Share of replications that estimate and bootstrap successfully | Reported for every condition |

Every result is reported with its Monte Carlo standard error or confidence
interval, so real discrepancies can be told apart from simulation noise.

Validation matrix (fractional factorial design, to keep computing time
reasonable):

| Dimension | Conditions |
|---|---|
| Loadings | 0.50, 0.60, 0.70, 0.80, 0.90, 0.95 |
| Indicators per construct | 1, 2, 3, 4, 6, 8 |
| Effect sizes | 0 (Type I error), small, medium, large |
| Model types | Simple mediation; parallel with 2 and 3 mediators; serial with 2 and 3 mediators |
| Sample sizes | Several values across 50–1,000 |
| Inference | Delta method / t test, and percentile bootstrap |

Stored in `validation/`: design file, seeds, package versions, raw summaries,
and a summary table. Simulations run in stages because of computing time.

Done when: all measures are reported for the full design, and any problem found
is fixed or documented before Stage 3.

---

## Stage 3: v0.9.0-rc (features)

Only after Stage 2 passes. Each feature gets its own tests and validation runs.

1. **Correlated mediators.** New input: the correlation between mediator
   *disturbances* after accounting for X (not the raw mediator correlation,
   which already includes their shared dependence on X).
2. **Negative paths and wider ranges.** Negative coefficients allowed. Instead
   of a fixed cap such as 0.90, the permitted range is documented and every
   model is checked for an admissible correlation matrix.
3. **Power for every path.**
   - a dropdown to choose the target path or effect;
   - a table of power for all structural paths;
   - a flag on any path below the target power;
   - direct and indirect effects reported separately.

Done when: features pass their tests and validation, and the methods text
and Methods tab describe them.

---

## Stage 4: v0.9.0 (public release)

- Update README, changelog, Methods tab and generated methods text.
- Publish the validation results and summary tables.
- Tag v0.9.0 and create a GitHub release; Zenodo issues a new version DOI.
  The v0.8.0 DOI remains the citation for that version.
- Update CITATION.cff and the in-app citation with the new DOI.
- Republish on Posit Connect Cloud.
- Use the Stage 2 results for a software paper (*SoftwareX* or *Journal of
  Open Source Software*).

---

## Later (only if users ask)

Moderation, formative or higher-order X, BCa intervals, Likert and skewed data,
missing-data adjustment, nested data, sensitivity analysis (smallest detectable
effect at a given N), "start from my SmartPLS results", PDF report, saved models.
