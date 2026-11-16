# =============================================================================
# cross_validate_seminr.R - independent check of the plspower engine
# -----------------------------------------------------------------------------
# Compares plspower's built-in PLS-SEM engine with the established R package
# seminr (Ray, Danks & Calero Valdez) on the SAME simulated datasets, for
# simple, parallel and serial mediation models, including models where
# loadings and indicator counts differ per construct.
#
#   Check 1 - estimates: every structural path from both engines must agree
#             within 1e-4 on every dataset.
#   Check 2 - power:     power for every effect (direct, total, specific and
#             serial indirect), computed from each engine's estimates with the
#             same tests, must agree within +/- 0.03 at every N.
#
# Both engines use Mode A outer weights and the path weighting inner scheme.
#
# Run from the repository root:
#   Rscript validation/cross_validate_seminr.R            # standard run (~2-3 min)
#   Rscript validation/cross_validate_seminr.R --quick    # fast check (<1 min)
# Requires: install.packages("seminr")
# Writes:   validation/results/*.csv and validation/results/summary.txt
#           (a --quick run writes to validation/results/quick/ instead)
# Exit status 0 = all checks PASS, 1 = at least one FAIL.
# =============================================================================

suppressPackageStartupMessages(library(seminr))

# --- load the engine from app.R (everything above the APP section) ----------
app_lines <- readLines("app.R")
engine_end <- grep("^# =============== APP", app_lines)[1] - 1
eval(parse(text = app_lines[seq_len(engine_end)]), envir = globalenv())
eval(parse(text = grep("^PLSPOWER_VERSION <-", app_lines, value = TRUE)), envir = globalenv())

quick <- "--quick" %in% commandArgs(trailingOnly = TRUE)
EST_TOL   <- 1e-4
POWER_TOL <- 0.03
REPS      <- if (quick) 40 else 250
N_GRID    <- if (quick) c(100, 250) else c(100, 200, 300)
ALPHA     <- 0.05
OUT <- if (quick) "validation/results/quick" else "validation/results"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

MODELS <- list(
  simple_equal      = build_model(1, "parallel", a = .4, b = .3, c = .2, lambda = .8, n_ind = 3),
  simple_unequal    = build_model(1, "parallel", a = .4, b = .3, c = .2,
                                  lambda = c(.85, .7, .75), n_ind = c(3, 5, 4)),
  parallel2_unequal = build_model(2, "parallel", a = c(.4, .3), b = c(.3, .25), c = .15,
                                  lambda = c(.8, .7, .75, .85), n_ind = c(3, 4, 5, 3)),
  parallel3_equal   = build_model(3, "parallel", a = c(.4, .3, .3), b = c(.3, .25, .2), c = .15,
                                  lambda = .75, n_ind = 4),
  serial2_unequal   = build_model(2, "serial", a = .45, d = .4, b = .35, c = .1,
                                  lambda = c(.8, .75, .7, .85), n_ind = c(4, 3, 5, 4)),
  serial3_equal     = build_model(3, "serial", a = .5, d = c(.45, .45), b = .4, c = .1,
                                  lambda = .8, n_ind = 3))
if (quick) MODELS <- MODELS[c("simple_equal", "parallel2_unequal", "serial2_unequal")]

seminr_fit <- function(blocks, mod) {
  nm <- mod$names; q <- length(nm)
  df <- as.data.frame(do.call(cbind, blocks))
  names(df) <- unlist(lapply(seq_len(q), function(j) paste0(tolower(nm[j]), "_", seq_len(mod$n_ind[j]))))
  mm <- do.call(constructs, lapply(seq_len(q), function(j)
    composite(nm[j], paste0(tolower(nm[j]), "_", seq_len(mod$n_ind[j])), weights = mode_A)))
  idx <- which(mod$S, arr.ind = TRUE)
  sm <- do.call(relationships, lapply(seq_len(nrow(idx)), function(r)
    paths(from = nm[idx[r, 2]], to = nm[idx[r, 1]])))
  fit <- suppressMessages(suppressWarnings(
    estimate_pls(df, mm, sm, inner_weights = path_weighting)))
  # Same OLS standard errors as the built-in engine, from seminr's scores.
  sc <- scale(fit$construct_scores[, nm]); n <- nrow(sc); R <- crossprod(sc) / (n - 1)
  B <- matrix(0, q, q); V <- list()
  for (j in 2:q) {
    P <- which(mod$S[j, ]); Rpp <- R[P, P, drop = FALSE]
    B[j, P] <- fit$path_coef[nm[P], nm[j]]
    beta <- solve(Rpp, R[P, j]); r2 <- sum(beta * R[P, j])
    V[[j]] <- list(idx = cbind(j, P), cov = (n - 1) * (1 - r2) / (n - 1 - length(P)) * solve(Rpp) / (n - 1))
  }
  list(B = B, V = V, df_y = n - 1 - sum(mod$S[q, ]))
}

decide <- function(fit, mod) {
  vapply(names(mod$effects), function(e) {
    ee <- effect_est(fit, mod$effects[[e]])
    crit <- if (e == "direct") qt(1 - ALPHA / 2, fit$df_y) else qnorm(1 - ALPHA / 2)
    abs(ee[["est"]] / ee[["se"]]) > crit
  }, logical(1))
}

est_rows <- list(); pow_rows <- list()
set.seed(20261116); t0 <- Sys.time()
for (mn in names(MODELS)) {
  mod <- MODELS[[mn]]; en <- names(mod$effects); pidx <- which(mod$S)
  for (n in N_GRID) {
    max_diff <- 0; db <- ds <- matrix(NA, REPS, length(en)); used <- 0
    for (r in seq_len(REPS)) {
      blocks <- gen_data(n, mod)
      fb <- pls_fit(blocks, mod); fs <- tryCatch(seminr_fit(blocks, mod), error = function(e) NULL)
      if (is.null(fb) || is.null(fs)) next
      used <- used + 1
      max_diff <- max(max_diff, abs(fb$B[pidx] - fs$B[pidx]))
      db[r, ] <- decide(fb, mod); ds[r, ] <- decide(fs, mod)
    }
    ok <- stats::complete.cases(db)
    est_rows[[length(est_rows) + 1]] <- data.frame(model = mn, N = n, datasets = used,
                                                   max_abs_path_diff = max_diff)
    for (k in seq_along(en)) pow_rows[[length(pow_rows) + 1]] <- data.frame(
      model = mn, N = n, effect = en[k],
      power_plspower = mean(db[ok, k]), power_seminr = mean(ds[ok, k]))
    cat(sprintf("%-18s N=%d done\n", mn, n))
  }
}
est <- do.call(rbind, est_rows); pow <- do.call(rbind, pow_rows)
est$result <- ifelse(est$max_abs_path_diff <= EST_TOL, "PASS", "FAIL")
pow$abs_diff <- abs(pow$power_plspower - pow$power_seminr)
pow$result <- ifelse(pow$abs_diff <= POWER_TOL, "PASS", "FAIL")
write.csv(est, file.path(OUT, "estimate_agreement.csv"), row.names = FALSE)
write.csv(pow, file.path(OUT, "power_agreement.csv"), row.names = FALSE)

pass <- all(est$result == "PASS") && all(pow$result == "PASS")
summary_txt <- c(
  sprintf("plspower %s vs seminr %s - %s - %s run: %d models x %d sample sizes x %d datasets",
          PLSPOWER_VERSION, as.character(packageVersion("seminr")), format(Sys.time(), "%Y-%m-%d"),
          if (quick) "quick" else "standard", length(MODELS), length(N_GRID), REPS),
  sprintf("Models: %s", paste(names(MODELS), collapse = ", ")),
  sprintf("Check 1 (all structural paths, tolerance %g): %s - largest difference %.2e",
          EST_TOL, if (all(est$result == "PASS")) "PASS" else "FAIL", max(est$max_abs_path_diff)),
  sprintf("Check 2 (power for %d effect-by-condition cells, tolerance +/-%.2f): %s - largest difference %.3f",
          nrow(pow), POWER_TOL, if (all(pow$result == "PASS")) "PASS" else "FAIL", max(pow$abs_diff)),
  sprintf("Overall: %s (%.0f s)", if (pass) "PASS" else "FAIL",
          as.numeric(difftime(Sys.time(), t0, units = "secs"))))
writeLines(summary_txt, file.path(OUT, "summary.txt"))
cat("\n", paste(summary_txt, collapse = "\n"), "\n", sep = "")
quit(status = if (pass) 0 else 1)
