# =============================================================================
# cross_validate_seminr.R - independent check of the plspower engine
# -----------------------------------------------------------------------------
# Compares plspower's built-in PLS-SEM engine with the established R package
# seminr (Ray, Danks & Calero Valdez) on the SAME simulated datasets.
#
#   Check 1 - estimates: for every dataset, the paths a (X -> M), b (M -> Y)
#             and c (X -> Y) from both engines must agree within 1e-4.
#   Check 2 - power:     power curves for the direct and indirect effect,
#             computed from each engine's estimates with the same tests,
#             must agree within +/- 0.03 at every N (the README's rule).
#
# Both engines use Mode A outer weights and the path weighting inner scheme.
#
# Run from the repository root:
#   Rscript validation/cross_validate_seminr.R            # standard run (~5-10 min)
#   Rscript validation/cross_validate_seminr.R --quick    # fast check (~1 min)
# Requires: install.packages("seminr")
# Writes:   validation/results/*.csv and validation/results/summary.txt
# Exit status 0 = all checks PASS, 1 = at least one FAIL.
# =============================================================================

suppressPackageStartupMessages(library(seminr))

# --- load the engine from app.R (everything above the APP section) ----------
app_lines <- readLines("app.R")
engine_end <- grep("^# =============== APP", app_lines)[1] - 1
eval(parse(text = app_lines[seq_len(engine_end)]), envir = globalenv())
eval(parse(text = grep("^PLSPOWER_VERSION <-", app_lines, value = TRUE)), envir = globalenv())

quick <- "--quick" %in% commandArgs(trailingOnly = TRUE)
EST_TOL   <- 1e-4    # estimates: max absolute difference
POWER_TOL <- 0.03    # power: max absolute difference at any N
REPS      <- if (quick) 60 else 300
N_GRID    <- if (quick) c(100, 250) else c(100, 150, 200, 300)
CONDITIONS <- expand.grid(lambda = if (quick) 0.8 else c(0.7, 0.85),
                          n_ind  = if (quick) 3 else c(3, 5))
PATHS <- list(a = 0.40, b = 0.30, c = 0.20)
ALPHA <- 0.05
dir.create("validation/results", showWarnings = FALSE, recursive = TRUE)

seminr_fit <- function(blocks) {
  k <- ncol(blocks[[1]])
  df <- as.data.frame(do.call(cbind, blocks))
  names(df) <- c(paste0("x", 1:k), paste0("m", 1:k), paste0("y", 1:k))
  mm <- constructs(composite("X", multi_items("x", 1:k), weights = mode_A),
                   composite("M", multi_items("m", 1:k), weights = mode_A),
                   composite("Y", multi_items("y", 1:k), weights = mode_A))
  sm <- relationships(paths(from = "X", to = c("M", "Y")),
                      paths(from = "M", to = "Y"))
  fit <- suppressMessages(suppressWarnings(
    estimate_pls(df, mm, sm, inner_weights = path_weighting)))
  pc <- fit$path_coef
  # Same OLS standard errors as the built-in engine, from seminr's scores.
  sc <- scale(fit$construct_scores[, c("X", "M", "Y")])
  n <- nrow(sc); R <- crossprod(sc) / (n - 1)
  beta <- solve(R[1:2, 1:2], R[1:2, 3]); r2 <- sum(beta * R[1:2, 3])
  se <- sqrt(diag((n - 1) * (1 - r2) / (n - 3) * solve(R[1:2, 1:2]) / (n - 1)))
  list(a = pc["X", "M"], b_m = pc["M", "Y"], b_direct = pc["X", "Y"],
       se = se[1], se_a = sqrt((1 - R[1, 2]^2) / (n - 2)), se_b = se[2])
}

decide <- function(f, n) {
  ab <- f$a * f$b_m
  c(direct   = abs(f$b_direct / f$se) > qt(1 - ALPHA / 2, n - 3),
    indirect = abs(ab / sqrt(f$b_m^2 * f$se_a^2 + f$a^2 * f$se_b^2)) > qnorm(1 - ALPHA / 2))
}

est_rows <- list(); pow_rows <- list()
set.seed(20261116)
t0 <- Sys.time()
for (ci in seq_len(nrow(CONDITIONS))) {
  cond <- CONDITIONS[ci, ]
  p <- c(PATHS, list(lambda = cond$lambda, n_ind = cond$n_ind))
  for (n in N_GRID) {
    dec_b <- dec_s <- matrix(NA, REPS, 2); diffs <- matrix(NA, REPS, 3)
    for (r in seq_len(REPS)) {
      blocks <- gen_data(n, p)
      fb <- pls_fit(blocks); fs <- tryCatch(seminr_fit(blocks), error = function(e) NULL)
      if (is.null(fb) || is.null(fs)) next
      diffs[r, ] <- c(fb$a - fs$a, fb$b_m - fs$b_m, fb$b_direct - fs$b_direct)
      dec_b[r, ] <- decide(fb, n); dec_s[r, ] <- decide(fs, n)
    }
    ok <- stats::complete.cases(diffs)
    est_rows[[length(est_rows) + 1]] <- data.frame(
      lambda = cond$lambda, n_ind = cond$n_ind, N = n, datasets = sum(ok),
      max_abs_diff_a = max(abs(diffs[ok, 1])), max_abs_diff_b = max(abs(diffs[ok, 2])),
      max_abs_diff_c = max(abs(diffs[ok, 3])))
    pow_rows[[length(pow_rows) + 1]] <- data.frame(
      lambda = cond$lambda, n_ind = cond$n_ind, N = n,
      power_direct_plspower = mean(dec_b[ok, 1]), power_direct_seminr = mean(dec_s[ok, 1]),
      power_indirect_plspower = mean(dec_b[ok, 2]), power_indirect_seminr = mean(dec_s[ok, 2]))
    cat(sprintf("lambda=%.2f k=%d N=%d done\n", cond$lambda, cond$n_ind, n))
  }
}
est <- do.call(rbind, est_rows); pow <- do.call(rbind, pow_rows)
est$result <- ifelse(pmax(est$max_abs_diff_a, est$max_abs_diff_b, est$max_abs_diff_c) <= EST_TOL, "PASS", "FAIL")
pow$max_abs_diff <- pmax(abs(pow$power_direct_plspower - pow$power_direct_seminr),
                         abs(pow$power_indirect_plspower - pow$power_indirect_seminr))
pow$result <- ifelse(pow$max_abs_diff <= POWER_TOL, "PASS", "FAIL")
write.csv(est, "validation/results/estimate_agreement.csv", row.names = FALSE)
write.csv(pow, "validation/results/power_agreement.csv", row.names = FALSE)

pass <- all(est$result == "PASS") && all(pow$result == "PASS")
summary_txt <- c(
  sprintf("plspower %s vs seminr %s - %s - %s run, %d datasets per cell",
          if (exists("PLSPOWER_VERSION")) PLSPOWER_VERSION else "(engine)",
          as.character(packageVersion("seminr")), format(Sys.time(), "%Y-%m-%d"),
          if (quick) "quick" else "standard", REPS),
  sprintf("Check 1 (estimates, tolerance %g): %s - largest difference %.2e",
          EST_TOL, if (all(est$result == "PASS")) "PASS" else "FAIL",
          max(est$max_abs_diff_a, est$max_abs_diff_b, est$max_abs_diff_c)),
  sprintf("Check 2 (power, tolerance +/-%.2f): %s - largest difference %.3f",
          POWER_TOL, if (all(pow$result == "PASS")) "PASS" else "FAIL", max(pow$max_abs_diff)),
  sprintf("Overall: %s (%.0f s)", if (pass) "PASS" else "FAIL",
          as.numeric(difftime(Sys.time(), t0, units = "secs"))))
writeLines(summary_txt, "validation/results/summary.txt")
cat("\n", paste(summary_txt, collapse = "\n"), "\n", sep = "")
quit(status = if (pass) 0 else 1)
