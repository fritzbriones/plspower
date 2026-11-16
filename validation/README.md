# Validation

`cross_validate_seminr.R` checks plspower's built-in PLS-SEM engine against the
established R package [seminr](https://cran.r-project.org/package=seminr).

```r
install.packages("seminr")
# from the repository root:
#   Rscript validation/cross_validate_seminr.R          # standard (~1 min)
#   Rscript validation/cross_validate_seminr.R --quick  # fast check
```

Both engines fit the same simulated datasets (Mode A, path weighting) for simple,
parallel (2-3 mediators) and serial (2-3 mediators) models with equal and unequal measurement.
Pass rules: path estimates within 1e-4; direct- and indirect-effect power
within ±0.03 at every N. Exit status 0 = PASS.

Latest result (v0.6.0 vs seminr 2.6.0, 2026-11-16; 6 model types): **PASS** — see
`results/summary.txt`, `results/estimate_agreement.csv` and
`results/power_agreement.csv`.
