# plspower - single-file deployable v0.5.0 (concatenated from modular project)
# Author: Fritz Jeffry Briones - https://github.com/fritzbriones/plspower
# Cite: Briones, F. J. (2026). plspower (Version 0.5.0) [Computer software].
# Deploy: paste into shinyapps.io as app.R, or run with shiny::runApp()
# Requires: shiny, bslib, ggplot2, future, promises, jsonlite

# =============== ENGINE ===============

# run_monte_carlo.R ------------------------------------------------------------
# Core simulation engine for PLS-SEM power analysis (mediation template).
# DGP: X ~ N(0,1); M = aX + eM; Y = bM + cX + eY (all latents unit variance);
#      reflective indicators x_ij = lambda*L + sqrt(1-lambda^2)*e_ij.
# Estimator: PLS-SEM Mode A, path-weighting inner scheme, fixed-point
#      iteration (1e-7, max 100). Paths: OLS of Y composite on X, M composites.
# Inference (direct effect X -> Y): (a) normal-theory t; (b) percentile
#      bootstrap CI excluding zero. Power = share of "significant" reps.

library(stats)

gen_data <- function(n, p) {
  e_m_sd  <- sqrt(1 - p$a^2)
  var_y   <- (p$c + p$a * p$b)^2 + p$b^2 * e_m_sd^2
  e_y_sd  <- sqrt(max(1e-9, 1 - var_y))
  lam_err <- sqrt(1 - p$lambda^2)
  X <- rnorm(n)
  M <- p$a * X + e_m_sd * rnorm(n)
  Y <- p$b * M + p$c * X + e_y_sd * rnorm(n)
  lats <- list(X, M, Y)
  lapply(lats, function(L)
    replicate(p$n_ind, p$lambda * L + lam_err * rnorm(n)))
}

# ---- PLS-SEM Mode A, path weighting (v0.5.0: matrix implementation) --------
# Numerically identical to v0.4.2 (max difference ~1e-15) and ~10x faster:
# indicators are standardised once, composite scores and inner weights are
# computed with matrix products, and the structural paths come from the
# composite correlation matrix (same estimates and SEs as lm()).
std_cols <- function(B) {
  n <- nrow(B)
  B <- B - rep(colMeans(B), each = n)
  s <- sqrt(colSums(B * B) / (n - 1))
  if (any(!is.finite(s) | s < 1e-12)) return(NULL)
  B / rep(s, each = n)
}

pls_fit <- function(blocks, tol = 1e-7, max_it = 100) {
  n <- nrow(blocks[[1]])
  blocks <- lapply(blocks, std_cols)
  if (any(vapply(blocks, is.null, logical(1)))) return(NULL)
  w <- lapply(blocks, function(B) rep(1, ncol(B)))
  nw <- function(v) v / sqrt(sum(v * v))
  S <- matrix(0, n, 3); conv <- FALSE
  for (it in seq_len(max_it)) {
    for (i in 1:3) {
      s <- blocks[[i]] %*% w[[i]]
      sd_s <- sqrt(sum(s * s) / (n - 1))   # indicators centred -> mean(s) = 0
      if (!is.finite(sd_s) || sd_s < 1e-12) return(NULL)
      S[, i] <- s / sd_s
    }
    R <- crossprod(S) / (n - 1)
    # Path weighting (Lohmoeller) for X -> M, X -> Y, M -> Y: successors get
    # correlation weights, predecessors get regression weights.
    b_y <- tryCatch(solve(R[1:2, 1:2], R[1:2, 3]), error = function(e) NULL)
    if (is.null(b_y)) return(NULL)
    E <- matrix(c(0,       R[1, 2], R[1, 3],    # X: successors M, Y
                  R[1, 2], 0,       R[2, 3],    # M: predecessor X, successor Y
                  b_y[1],  b_y[2],  0), 3, 3)   # Y: predecessors X, M
    Z <- S %*% E
    w_new <- lapply(1:3, function(i) drop(crossprod(blocks[[i]], Z[, i])) / (n - 1))
    delta <- max(abs(unlist(Map(function(a, b) nw(a) - nw(b), w_new, w))))
    w <- w_new
    if (delta < tol) { conv <- TRUE; break }
  }
  if (!conv) return(NULL)
  Rxx <- R[1:2, 1:2]; rxy <- R[1:2, 3]
  beta <- solve(Rxx, rxy)                       # Y ~ X + M
  r2 <- sum(beta * rxy)
  s2 <- (n - 1) * (1 - r2) / (n - 3)
  se_beta <- sqrt(diag(s2 * solve(Rxx) / (n - 1)))
  a <- R[1, 2]                                  # M ~ X
  list(b_direct = beta[1], se = se_beta[1], r2 = r2,
       a = a, se_a = sqrt((1 - a^2) / (n - 2)),
       b_m = beta[2], se_b = se_beta[2], converged = TRUE)
}

# One Monte Carlo replication. Tests (two-sided, level alpha):
#   direct effect c : normal-theory t test; percentile bootstrap CI
#   indirect a x b  : Sobel z test;         percentile bootstrap CI
one_rep <- function(n, p, boot_reps, alpha) {
  blocks <- gen_data(n, p)
  fit <- pls_fit(blocks)
  na <- list(nt = NA, boot = NA, b = NA, ind_sobel = NA, ind_boot = NA, ab = NA)
  if (is.null(fit)) return(na)
  nt <- abs(fit$b_direct / fit$se) > qt(1 - alpha / 2, n - 3)
  ab <- fit$a * fit$b_m
  sobel_se <- sqrt(fit$b_m^2 * fit$se_a^2 + fit$a^2 * fit$se_b^2)
  ind_sobel <- abs(ab / sobel_se) > qnorm(1 - alpha / 2)
  boot <- ind_boot <- NA
  if (boot_reps > 0) {
    bs <- bab <- numeric(boot_reps); ok <- 0
    for (bb in seq_len(boot_reps)) {
      rows <- sample.int(n, n, replace = TRUE)
      fb <- pls_fit(lapply(blocks, function(B) B[rows, , drop = FALSE]))
      if (!is.null(fb)) { ok <- ok + 1; bs[ok] <- fb$b_direct; bab[ok] <- fb$a * fb$b_m }
    }
    if (ok >= boot_reps / 2) {
      pr <- c(alpha / 2, 1 - alpha / 2)
      q  <- quantile(bs[seq_len(ok)], pr, names = FALSE)
      qi <- quantile(bab[seq_len(ok)], pr, names = FALSE)
      boot     <- !(q[1]  <= 0 && q[2]  >= 0)
      ind_boot <- !(qi[1] <= 0 && qi[2] >= 0)
    }
  }
  list(nt = nt, boot = boot, b = fit$b_direct,
       ind_sobel = ind_sobel, ind_boot = ind_boot, ab = ab)
}

# Power at one sample size.
power_at_n <- function(n, p, mc_reps, boot_reps, alpha) {
  r <- lapply(seq_len(mc_reps), function(i) one_rep(n, p, boot_reps, alpha))
  ok <- !vapply(r, function(x) is.na(x$b), logical(1))
  g <- function(k) { v <- unlist(lapply(r[ok], `[[`, k)); v <- v[!is.na(v)]
                     if (length(v)) mean(v) else NA_real_ }
  data.frame(N = n,
             power_nt = g("nt"), power_boot = g("boot"),
             power_ind_sobel = g("ind_sobel"), power_ind_boot = g("ind_boot"),
             bias = if (any(ok)) mean(unlist(lapply(r[ok], `[[`, "b"))) - p$c else NA_real_,
             bias_ind = if (any(ok)) mean(unlist(lapply(r[ok], `[[`, "ab"))) - p$a * p$b else NA_real_,
             nonconv = sum(!ok), reps = sum(ok))
}

# Number of CPU cores used for a run (forking is not available on Windows).
plspower_cores <- function() {
  if (.Platform$OS.type == "windows") return(1L)
  max(1L, min(4L, parallel::detectCores(logical = FALSE), na.rm = TRUE))
}

# Power curve over the N grid. v0.5.0: each grid point gets its own
# L'Ecuyer-CMRG random-number stream derived from `seed`, so results are
# identical whether the grid points run one after another or in parallel.
power_curve <- function(p, n_grid, mc_reps, boot_reps, alpha, seed,
                        cores = plspower_cores()) {
  old_kind <- RNGkind()[1]
  on.exit(RNGkind(old_kind), add = TRUE)
  RNGkind("L'Ecuyer-CMRG"); set.seed(seed)
  streams <- vector("list", length(n_grid)); s <- .Random.seed
  for (i in seq_along(n_grid)) { streams[[i]] <- s; s <- parallel::nextRNGStream(s) }
  job <- function(i) {
    assign(".Random.seed", streams[[i]], envir = globalenv())
    power_at_n(n_grid[i], p, mc_reps, boot_reps, alpha)
  }
  idx <- seq_along(n_grid)
  res <- if (cores > 1 && length(idx) > 1) {
    parallel::mclapply(idx, job, mc.cores = cores, mc.preschedule = FALSE)
  } else lapply(idx, job)
  bad <- vapply(res, inherits, logical(1), "try-error")
  if (any(bad)) stop("Simulation failed: ", conditionMessage(attr(res[[which(bad)[1]]], "condition")))
  do.call(rbind, res)
}

min_n_for_power <- function(n_grid, pw, target) {
  if (all(is.na(pw))) return(NA)
  lg <- function(x) log(x / (1 - x))
  for (i in seq_along(pw)) {
    if (!is.na(pw[i]) && pw[i] >= target) {
      if (i == 1) return(n_grid[i])
      p0 <- pw[i - 1]; n0 <- n_grid[i - 1]; p1 <- pw[i]; n1 <- n_grid[i]
      if (is.na(p0) || p1 == p0) return(n_grid[i])
      t <- (lg(target) - lg(p0)) / (lg(p1) - lg(p0))
      return(round(n0 + t * (n1 - n0)))
    }
  }
  NA
}
# run_gpower_equivalent.R --------------------------------------------------------
# Exact power for the F test of the R-squared INCREASE contributed by the
# target predictor (direct effect X -> Y) in the observed-score regression
# of Y on (X, M). Answers the regression-shortcut planning question only --
# NOT a claim about PLS-SEM power. Reported as a comparison method.

gpower_required_n <- function(a, b, c, alpha = 0.05, target = 0.80,
                              n_max = 5000) {
  r_xm <- a
  r_xy <- c + a * b
  r_my <- b + a * c
  b_x  <- (r_xy - r_xm * r_my) / (1 - r_xm^2)
  b_m  <- (r_my - r_xm * r_xy) / (1 - r_xm^2)
  r2_full <- b_x * r_xy + b_m * r_my
  inc <- b_x^2 * (1 - r_xm^2)
  f2 <- inc / (1 - r2_full)
  for (n in 4:n_max) {
    df2 <- n - 3
    crit <- qf(1 - alpha, 1, df2)
    power <- 1 - pf(crit, 1, df2, f2 * n)
    if (power >= target) return(list(N = n, f2 = f2, r2_full = r2_full))
  }
  list(N = NA, f2 = f2, r2_full = r2_full)
}
# run_inverse_sqrt.R --------------------------------------------------------------
# Kock & Hadaya (2018) inverse square-root rule. Comparison method only.

inverse_sqrt_n <- function(b_min, alpha = 0.05, target = 0.80) {
  z <- qnorm(1 - alpha / 2) + qnorm(target)
  ceiling(z^2 / b_min^2)
}
# Independent cross-validation against seminr/cSEM is planned as a separate
# validation script; neither package is an in-app estimation backend.
# self_test.R --------------------------------------------------------------------
# SMOKE TEST (not publication-grade validation).
# Verifies that the deployed R engine reproduces the reference fixture
# (Python twin implementation, same DGP, same algorithm class) within
# EXPLICIT tolerance rules. 400 MC reps means MC noise alone can move power
# by ~+/-0.025; this test detects implementation/environment failures, not
# precise minimum-N estimates.
#
# Acceptance criteria (fixed a priori):
#   - power (each N, normal theory): difference <= 2.5 combined Monte Carlo SE
#   - f2 (G*Power-equivalent, deterministic): full-precision reference,
#     absolute tolerance 1e-10
#   - inverse-sqrt N (deterministic): exact integer match
#   - min-N is reported INFORMATIONALLY (interpolation artifact), no pass/fail
# PASS requires ALL tolerance checks to pass.
#
# Future independent cross-validation should compare this implementation
# with established PLS-SEM software using a separately documented protocol.

FIXTURE <- list(
  a = 0.40, b = 0.30, c = 0.20, lambda = 0.80, n_ind = 3,
  seed = 2026, reps = 400, grid = c(100, 150, 200, 250, 300),
  ref_power = c(0.432, 0.585, 0.690, 0.795, 0.902),
  ref_min_n = 252, ref_f2 = 0.04087591240875911, ref_inv_sqrt = 197
)

run_self_test <- function() {
  p <- list(a = FIXTURE$a, b = FIXTURE$b, c = FIXTURE$c,
            lambda = FIXTURE$lambda, n_ind = FIXTURE$n_ind)
  curve <- power_curve(p, FIXTURE$grid, FIXTURE$reps,
                       boot_reps = 0, alpha = 0.05, FIXTURE$seed)
  pw <- curve$power_nt
  gp <- gpower_required_n(p$a, p$b, p$c, 0.05, 0.80)
  isq <- inverse_sqrt_n(p$c, 0.05, 0.80)

  rows <- list()
  add_row <- function(check, ref, app, tol, ok) rows[[length(rows) + 1]] <<-
    data.frame(check = check, reference = ref, app = app,
               abs_diff = abs(app - ref), tolerance = tol,
               result = ifelse(ok, "PASS", "FAIL"))

  for (i in seq_along(FIXTURE$grid)) {
    ref_pw <- FIXTURE$ref_power[i]
    app_pw <- pw[i]
    # Both the reference and app power estimates use finite Monte Carlo
    # samples. Compare them using 2.5 times the combined binomial MC standard
    # error, rather than an arbitrary relative percentage tolerance.
    tol <- 2.5 * sqrt(
      ref_pw * (1 - ref_pw) / FIXTURE$reps +
      app_pw * (1 - app_pw) / FIXTURE$reps
    )
    add_row(sprintf("Power NT, N=%d", FIXTURE$grid[i]),
            ref_pw, app_pw, tol,
            is.finite(app_pw) && abs(app_pw - ref_pw) <= tol)
  }
  # Full-precision deterministic reference; small floating-point tolerance.
  f2_tol <- 1e-10
  f2_ok <- abs(gp$f2 - FIXTURE$ref_f2) <= f2_tol
  add_row("f2 (G*Power-equivalent)", FIXTURE$ref_f2, gp$f2,
          f2_tol, f2_ok)
  isq_ok <- isq == FIXTURE$ref_inv_sqrt
  add_row("Inverse-sqrt N", FIXTURE$ref_inv_sqrt, isq, 0, isq_ok)

  table <- do.call(rbind, rows)
  list(
    curve = curve, table = table,
    pass = all(table$result == "PASS"),
    min_n = min_n_for_power(FIXTURE$grid, pw, 0.80),
    min_n_ref = FIXTURE$ref_min_n,
    gpower_n = gp$N, f2 = gp$f2, inv_sqrt = isq,
    label = "Smoke test: verifies implementation agreement (400 MC reps),
             not publication-grade precision."
  )
}

# =============== APP ===============

# app.R ---------------------------------------------------------------------------
# plspower: a simulation PLANNER for PLS-SEM power analysis (v0.5.0)
# Prospective power only. v1 scope: one mediation template, reflective Mode A.
# Comparison methods reported with neutral wording ("planning-method
# discrepancy"). In-app self-test is a SMOKE TEST (implementation agreement),
# not publication-grade validation.

REQ_PKGS <- c("shiny", "bslib", "ggplot2", "future", "promises", "jsonlite")
missing_pkgs <- REQ_PKGS[!vapply(REQ_PKGS, requireNamespace, logical(1),
                                 quietly = TRUE)]
if (length(missing_pkgs))
  stop("plspower requires: ", paste(missing_pkgs, collapse = ", "),
       ". Install the missing packages with install.packages().")

library(shiny)
library(bslib)
library(ggplot2)
library(future)
library(promises)
library(jsonlite)

# (engine inlined above)

PLSPOWER_VERSION <- "0.5.0"
APP_VERSION <- paste("plspower", PLSPOWER_VERSION)
AUTHOR <- "Fritz Jeffry Briones"
REPO_URL <- "https://github.com/fritzbriones/plspower"
CITATION_APA <- sprintf(paste0(
  "Briones, F. J. (2026). plspower: A Monte Carlo simulation planner for ",
  "PLS-SEM power analysis (Version %s) [Computer software]. %s"),
  PLSPOWER_VERSION, REPO_URL)
CITATION_INTEXT <- "(Briones, 2026)"
CITATION_BIBTEX <- sprintf(paste0(
  "@software{briones2026plspower,\n",
  "  author  = {Briones, Fritz Jeffry},\n",
  "  title   = {plspower: A Monte Carlo simulation planner for PLS-SEM power analysis},\n",
  "  year    = {2026},\n  version = {%s},\n  url     = {%s}\n}"),
  PLSPOWER_VERSION, REPO_URL)


parse_grid <- function(txt) {
  rng <- suppressWarnings(as.numeric(trimws(strsplit(txt, ",")[[1]])))
  if (length(rng) != 3 || anyNA(rng) || rng[1] < 20 || rng[2] < rng[1] || rng[3] <= 0)
    stop("N grid must be 'from,to,by' with from >= 20, to >= from and by > 0")
  seq(rng[1], rng[2], by = rng[3])
}

check_range <- function(lo, hi) {
  if (!is.numeric(lo) || !is.numeric(hi) || is.na(lo) || is.na(hi) ||
      lo < 20 || hi <= lo || hi > 5000)
    stop("Find my N: minimum N must be at least 20 and below the maximum N (maximum 5,000)")
  invisible(TRUE)
}

fmt_duration <- function(secs) {
  if (secs < 90) return(sprintf("%d seconds", max(1, round(secs))))
  if (secs < 5400) return(sprintf("%d minutes", round(secs / 60)))
  if (secs < 172800) return(sprintf("%.1f hours", secs / 3600))
  sprintf("%.1f days", secs / 86400)
}

# Time PLS fits at N = 200 once at startup (a short bootstrap replication, so
# resampling overhead is included), for the runtime estimate. ~0.1 s.
FIT_MS_N200 <- local({
  p0 <- list(a = .4, b = .3, c = .2, lambda = .8, n_ind = 3)
  set.seed(1); one_rep(200, p0, 5, .05)                    # warm-up
  t <- system.time(one_rep(200, p0, 60, .05))[["elapsed"]]
  max(0.05, 1000 * t / 61)
})
N_CORES <- plspower_cores()

# Run simulations in a background R process so the app stays responsive.
# (Set after the timing above so starting the workers does not slow it.)
future::plan(future::multisession, workers = 2)

TIERS <- list(
  quick        = list(mc = 200,  boot = 500,   label = "Quick preview (200 MC x 500 boot)"),
  standard     = list(mc = 1000, boot = 2000,  label = "Standard (1,000 MC x 2,000 boot)"),
  publication  = list(mc = 5000, boot = 10000, label = "Publication (5,000 MC x 10,000 boot)")
)

EFFECT_LABEL <- c(direct = "direct effect c (X -> Y)",
                  indirect = "indirect effect a x b (X -> M -> Y)")

# ---- Find my N: two-stage search ------------------------------------------
# Stage 1 scans a wide N range on a log-spaced grid with fast normal-theory
# tests (Sobel for a x b) and 200 replications. Stage 2 runs the full tier
# (bootstrap included) on a fine grid around the N where stage-1 power first
# reached the target.
COARSE_POINTS <- 8
COARSE_MC <- 200
FINE_POINTS <- 5
coarse_grid <- function(lo, hi)
  unique(round(exp(seq(log(lo), log(hi), length.out = COARSE_POINTS))))
fine_grid <- function(lower, upper)
  unique(round(seq(max(20, lower), upper * 1.1, length.out = FINE_POINTS)))

run_power_simulation <- function(settings) {
  p <- list(a = settings$a, b = settings$b, c = settings$c,
            lambda = settings$lambda, n_ind = settings$n_ind)
  if (!identical(settings$engine, "built-in")) {
    stop("Only the built-in engine is currently integrated and supported.")
  }
  t0 <- Sys.time()
  coarse <- NULL; search_note <- NULL
  if (identical(settings$mode, "find")) {
    cg <- coarse_grid(settings$n_min, settings$n_max)
    coarse <- power_curve(p, cg, COARSE_MC, 0, settings$alpha, settings$seed,
                          cores = settings$cores)
    pw <- if (identical(settings$effect, "indirect")) coarse$power_ind_sobel else coarse$power_nt
    i <- which(pw >= settings$target)[1]
    if (is.na(i)) {
      settings$n_grid <- cg
      search_note <- sprintf(paste0("Target power was not reached anywhere up to N = %d ",
        "(highest stage-1 power %.2f). Increase the maximum N or check the effect sizes."),
        settings$n_max, max(pw, na.rm = TRUE))
      curve <- coarse
    } else {
      lower <- if (i == 1) max(20, round(cg[1] / 2)) else cg[i - 1]
      settings$n_grid <- fine_grid(lower, cg[i])
      if (i == 1) search_note <- sprintf(paste0("Power already reached the target at the ",
        "minimum N (%d); the fine grid starts lower. Lower the minimum N to scan further down."), cg[1])
      curve <- power_curve(p, settings$n_grid, settings$mc_reps, settings$boot_reps,
                           settings$alpha, settings$seed + 1L, cores = settings$cores)
    }
  } else {
    curve <- power_curve(p, settings$n_grid, settings$mc_reps,
                         settings$boot_reps, settings$alpha, settings$seed,
                         cores = settings$cores)
  }
  gp <- gpower_required_n(p$a, p$b, p$c, settings$alpha, settings$target)
  isq <- inverse_sqrt_n(p$c, settings$alpha, settings$target)
  g <- settings$n_grid; tg <- settings$target
  list(
    settings = settings, curve = curve,
    min_n_nt = min_n_for_power(g, curve$power_nt, tg),
    min_n_boot = min_n_for_power(g, curve$power_boot, tg),
    min_n_ind_sobel = min_n_for_power(g, curve$power_ind_sobel, tg),
    min_n_ind_boot = min_n_for_power(g, curve$power_ind_boot, tg),
    gpower = gp, inverse_sqrt = isq,
    runtime = settings$runtime,
    elapsed_sec = as.numeric(difftime(Sys.time(), t0, units = "secs")),
    version = APP_VERSION, engine = settings$engine,
    citation = CITATION_APA,
    coarse = coarse, search_note = search_note
  )
}

fmt_n <- function(x) if (is.na(x)) "not reached in grid" else format(x)
num <- function(x, d = 2) formatC(x, format = "f", digits = d)

# Ready-to-paste sample-size justification for a manuscript.
methods_text <- function(r) {
  s <- r$settings; ind <- identical(s$effect, "indirect")
  scan_only <- !is.null(r$coarse) && identical(r$curve, r$coarse)
  if (scan_only) { s$boot_reps <- 0; s$mc_reps <- COARSE_MC }   # only the fast scan ran
  boot <- s$boot_reps > 0
  min_n <- if (ind) { if (boot) r$min_n_ind_boot else r$min_n_ind_sobel }
           else      { if (boot) r$min_n_boot     else r$min_n_nt }
  pw <- if (ind) { if (boot) r$curve$power_ind_boot else r$curve$power_ind_sobel }
        else      { if (boot) r$curve$power_boot     else r$curve$power_nt }
  if (all(is.na(pw)))   # e.g. Find my N stopped after the scan (no bootstrap run)
    pw <- if (ind) r$curve$power_ind_sobel else r$curve$power_nt
  test <- if (boot) sprintf(paste0("the two-sided %d%% percentile bootstrap confidence ",
                                   "interval (%s resamples) excluded zero"),
                            round(100 * (1 - s$alpha)), format(s$boot_reps, big.mark = ","))
          else if (ind) "the Sobel test was significant (two-sided)"
          else "the two-sided t test was significant"
  grid <- s$n_grid
  out <- paste0(
    "An a priori Monte Carlo power analysis was conducted with plspower (Version ",
    PLSPOWER_VERSION, "; Briones, 2026). Data were generated from a simple mediation ",
    "model (X -> M -> Y with a direct path X -> Y) with standardized path coefficients ",
    "a = ", num(s$a), ", b = ", num(s$b), " and c = ", num(s$c), ". Each construct ",
    "was measured reflectively by ", s$n_ind, " indicators with loadings of ",
    num(s$lambda), ". Every simulated sample was estimated with PLS-SEM (Mode A ",
    "outer weights, path weighting scheme). Statistical power for the ",
    EFFECT_LABEL[[if (ind) "indirect" else "direct"]],
    " was defined as the proportion of ", format(s$mc_reps, big.mark = ","),
    " Monte Carlo replications in which ", test, " at alpha = ",
    sub("^0", "", num(s$alpha)), ". Sample sizes from N = ", min(grid), " to N = ",
    max(grid), " in steps of ", if (length(grid) > 1) diff(grid)[1] else 0,
    " were evaluated (random seed ", s$seed, "). ")
  if (!is.null(r$coarse) && identical(r$curve, r$coarse)) out <- paste0(out, sprintf(
    "These sample sizes come from a log-spaced scan of N = %d to %d (%d replications, %s); no bootstrap stage was run because the target was not reached. ",
    s$n_min, s$n_max, COARSE_MC, if (ind) "Sobel test" else "normal-theory t test"))
  else if (!is.null(r$coarse)) out <- paste0(out, sprintf(paste0(
    "These sample sizes were selected by a two-stage search: a log-spaced scan ",
    "of N = %d to %d (%d replications, %s) located the region where power reached ",
    "the target, and the full analysis was then run on a fine grid in that region. "),
    s$n_min, s$n_max, COARSE_MC,
    if (ind) "Sobel test" else "normal-theory t test"))
  out <- paste0(out, if (!is.na(min_n))
    sprintf("A sample of approximately N = %d was required to reach %d%% power (interpolated from the simulated power curve). ",
            min_n, round(100 * s$target))
    else sprintf("The target power of %d%% was not reached within the evaluated range (maximum simulated power = %s). ",
                 round(100 * s$target), num(max(pw, na.rm = TRUE))))
  if (!ind) out <- paste0(out, sprintf(paste0(
    "For comparison, a regression-based calculation equivalent to G*Power ",
    "(f-squared = %s) indicated N = %s, and the inverse square-root method ",
    "(Kock & Hadaya, 2018; two-tailed) indicated N = %s; these shortcuts assume ",
    "error-free composites."), num(r$gpower$f2, 3), fmt_n(r$gpower$N), r$inverse_sqrt))
  if (ind && !boot) out <- paste0(out, "The Sobel test (Sobel, 1982) was used for the indirect effect.")
  if (identical(s$tier, "quick"))
    out <- paste0(out, "\n\n[NOTE: Quick-preview tier - re-run with the Standard or Publication tier before reporting.]")
  refs <- c(paste0("Briones, F. J. (2026). plspower: A Monte Carlo simulation planner for ",
                   "PLS-SEM power analysis (Version ", PLSPOWER_VERSION, ") [Computer software]. ", REPO_URL),
            if (!ind) paste0("Kock, N., & Hadaya, P. (2018). Minimum sample size estimation in PLS-SEM: ",
                             "The inverse square root and gamma-exponential methods. Information ",
                             "Systems Journal, 28(1), 227-261."),
            if (ind && !boot) paste0("Sobel, M. E. (1982). Asymptotic confidence intervals for indirect ",
                                     "effects in structural equation models. Sociological Methodology, 13, 290-312."))
  paste0(out, "\n\nReferences\n", paste(refs, collapse = "\n"))
}

cite_card <- function() card(
  card_header("How to cite plspower"),
  card_body(
    tags$p("If you use plspower in a thesis, article or report, please cite:"),
    tags$blockquote(style = "font-size:0.95em;", CITATION_APA),
    tags$p(tags$small("In-text: ", CITATION_INTEXT)),
    tags$details(tags$summary("BibTeX"), tags$pre(CITATION_BIBTEX)),
    tags$p(tags$small("Created by ", tags$strong(AUTHOR), " - ",
                      tags$a(href = REPO_URL, target = "_blank", REPO_URL)))
  ))

ui <- bslib::page_fillable(
  title = "plspower - PLS-SEM Monte Carlo power planner",
  theme = NULL,
  bslib::layout_columns(
    col_widths = c(4, 8),
    card(
      card_header("1. Specify the model"),
      card_body(
        selectInput("template", "Model template",
                    choices = c("Simple mediation (X -> M -> Y + X -> Y)")),
        radioButtons("effect", "Target effect",
                     choices = c("Direct effect c (X -> Y)" = "direct",
                                 "Indirect effect a x b (X -> M -> Y)" = "indirect")),
        sliderInput("a", "Path a (X -> M)", 0.05, 0.60, 0.40, 0.05),
        sliderInput("b", "Path b (M -> Y)", 0.05, 0.60, 0.30, 0.05),
        sliderInput("c", "Direct effect c (X -> Y)", 0.05, 0.40, 0.20, 0.05),
        sliderInput("lambda", "Indicator loading (all)", 0.50, 0.95, 0.80, 0.05),
        sliderInput("n_ind", "Indicators per construct", 3, 6, 3, 1),
        tags$p(strong("Estimation engine: "), "Built-in PLS Mode A, path weighting"),
        tags$p(class = "text-muted", "Cross-check against seminr: see validation/ in the GitHub repository.")
      ),
      card_header("2. Analysis settings"),
      card_body(
        selectInput("alpha", "Significance level", c(0.05, 0.01, 0.10)),
        selectInput("target", "Target power", c(0.80, 0.90)),
        selectInput("tier", "Compute tier", choices = names(TIERS), selected = "quick"),
        textOutput("tier_note"),
        uiOutput("tier_banner"),
        numericInput("seed", "Random seed", 2026),
        radioButtons("mode", "Sample sizes",
                     choices = c("Find my N (automatic search)" = "find",
                                 "Custom N grid" = "grid")),
        conditionalPanel("input.mode == 'find'",
          tags$div(style = "display:flex;gap:8px;",
            numericInput("n_min", "Minimum N", 50, min = 20, step = 10),
            numericInput("n_max", "Maximum N", 1000, min = 30, step = 50)),
          tags$p(class = "text-muted", style = "font-size:0.85em;",
                 "Scans this range quickly, then runs the full analysis only where power reaches your target.")),
        conditionalPanel("input.mode == 'grid'",
          textInput("n_grid_text", "N grid (from,to,by)", "100,300,50")),
        textOutput("runtime_note")
      )
    ),
    navset_card_tab(
      nav_panel("3. Results",
        card(card_header("Run"),
             card_body(actionButton("run_analysis", "Run power analysis",
                                    class = "btn-primary"),
                       textOutput("elapsed_note"))),
        uiOutput("required_n_card"),
        plotOutput("power_curve"),
        card(card_header("Simulation diagnostics"),
             card_body(tableOutput("diag_table"))),
        card(card_header("4. Methods text for your manuscript"),
             card_body(
               tags$p(class = "text-muted",
                      "A ready-to-paste sample-size justification. Check it against your study before use."),
               verbatimTextOutput("methods_text"),
               downloadButton("dl_methods", "Methods text (TXT)")
             )),
        card(card_header("5. Export"),
             card_body(
               downloadButton("dl_csv", "Power curve (CSV)"),
               downloadButton("dl_json", "Full run (JSON - seed, settings, package versions, citation)")
             )),
        cite_card()
      ),
      nav_panel("Methods & Validation",
        card(card_body(
          tags$h4("Data-generating model"),
          tags$code("M = a*X + eM;  Y = b*M + c*X + eY (all latents standardized)"),
          tags$p("Indicators: x_ij = lambda*L + sqrt(1-lambda^2)*e_ij; Gaussian errors."),
          tags$h4("Estimator"),
          tags$p("Built-in PLS Mode A outer weights, path weighting inner scheme
                 (including the direct path X -> Y), fixed-point iteration (1e-7),
                 structural paths by OLS on standardized composites."),
          tags$h4("Definition of power"),
          tags$p("Direct effect c: proportion of replications in which the two-sided
                 t test (normal theory) or the percentile bootstrap CI excludes zero."),
          tags$p("Indirect effect a x b: proportion of replications in which the Sobel
                 test (normal theory) or the percentile bootstrap CI of a x b excludes
                 zero. The bootstrap is generally preferred for indirect effects."),
          tags$h4("Comparison methods (direct effect only)"),
          tags$p("G*Power-equivalent (exact noncentral F, observed-score regression)
                 and the Kock & Hadaya (2018) inverse square-root rule. Both apply to
                 a single path, so they are not shown for the indirect effect.
                 Differences are reported as 'planning-method discrepancy'."),
          tags$p("Note: the inverse square-root rule is computed two-tailed,
                 consistent with this app's two-sided inference; Kock &
                 Hadaya's published tables are one-tailed (e.g., 155 vs 197
                 at b = 0.20, alpha = 0.05, power = 0.80)."),
          tags$h4("Positioning"),
          tags$p("This tool is not a new estimator. Its contribution is a
                 transparent comparison and calibration layer across
                 planning methods for PLS-SEM: the same study, planned
                 several ways, with measurement quality as an explicit
                 design input. Regression shortcuts and the inverse
                 square-root rule implicitly assume error-free composites;
                 with imperfect indicators the simulation typically requires
                 a larger N."),
          tags$h4("Speed and reproducibility"),
          tags$p(sprintf("Each N in the grid runs on its own random-number stream
                 (L'Ecuyer-CMRG, derived from the seed), so results are identical
                 whether grid points run in sequence or in parallel. This server
                 uses %d core(s).", N_CORES)),
          tags$h4("Built-in smoke test"),
          tags$p("Reproduces the seed-2026 reference fixture with explicit tolerance
                 rules (power: 2.5 combined Monte Carlo standard errors;
                 f2: 1e-10; inverse-sqrt: exact). SMOKE TEST: verifies
                 implementation agreement, not publication-grade precision
                 (400 MC reps)."),
          actionButton("self_test", "Run engine smoke test"),
          tableOutput("self_test_table"),
          textOutput("self_test_verdict"),
          tags$h4("Known limitations (v1)"),
          tags$ul(
            tags$li("One mediation template; no formative constructs, MGA,
                    custom diagrams, or user-data fitting (prospective power only)."),
            tags$li("Equal loadings and indicator counts for all constructs."),
            tags$li("Results apply ONLY to the specified model and settings.")
          )
        )),
        cite_card()
      )
    )
  ),
  tags$footer(style = "margin-top:12px;color:#888;font-size:0.85em;",
              APP_VERSION, " - created by ", AUTHOR, " - ",
              tags$a(href = REPO_URL, target = "_blank", "GitHub"),
              " - please cite: ", CITATION_INTEXT)
)

server <- function(input, output, session) {
  settings <- eventReactive(input$run_analysis, {
    grid <- tryCatch(
      if (identical(input$mode, "find")) { check_range(input$n_min, input$n_max)
                                           coarse_grid(input$n_min, input$n_max) }
      else parse_grid(input$n_grid_text),
      error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL })
    req(grid)
    tier <- TIERS[[input$tier]]
    list(a = input$a, b = input$b, c = input$c, lambda = input$lambda,
         n_ind = input$n_ind, effect = input$effect,
         mode = input$mode, n_min = input$n_min, n_max = input$n_max,
         alpha = as.numeric(input$alpha), target = as.numeric(input$target),
         n_grid = grid, mc_reps = tier$mc, boot_reps = tier$boot,
         seed = input$seed, runtime = Sys.time(), cores = N_CORES,
         tier = input$tier, engine = "built-in",
         plspower_version = PLSPOWER_VERSION,
         versions = c(list(R = R.version.string),
                      lapply(setNames(REQ_PKGS, REQ_PKGS),
                             function(p) as.character(utils::packageVersion(p)))))
  })

  power_task <- ExtendedTask$new(function(s) future_promise({
    run_power_simulation(s)
  }, seed = TRUE)) |> bind_task_button("run_analysis")
  observeEvent(input$run_analysis, power_task$invoke(settings()))

  self_task <- ExtendedTask$new(function(...) future_promise({
    run_self_test()
  }, seed = TRUE)) |> bind_task_button("self_test")
  observeEvent(input$self_test, self_task$invoke(list()))

  res <- power_task$result

  output$tier_note <- renderText(paste("Tier:", TIERS[[input$tier]]$label))
  output$tier_banner <- renderUI({
    if (input$tier == "quick")
      tags$p(style = "color:#c0392b;font-weight:600;",
             "Quick tier: DO NOT USE FOR REPORTING. Use Standard or
              Publication tier for any reported estimate.")
  })

  # Runtime estimate shown BEFORE the run: N grid x MC reps x (1 + bootstrap)
  # PLS fits, fit time measured on this server, grid points split over cores.
  output$runtime_note <- renderText({
    tier <- TIERS[[input$tier]]
    cost <- function(grid, mc, boot) sum(grid / 200) * mc * (1 + boot) *
      FIT_MS_N200 / 1000 * (input$n_ind / 3) / min(N_CORES, length(grid))
    if (identical(input$mode, "find")) {
      ok <- tryCatch({ check_range(input$n_min, input$n_max); TRUE }, error = function(e) FALSE)
      if (!ok) return("Enter a minimum N of at least 20 and a maximum N above it (up to 5,000).")
      cg <- coarse_grid(input$n_min, input$n_max)
      s1 <- cost(cg, COARSE_MC, 0)
      lo <- cost(fine_grid(cg[1], cg[2]), tier$mc, tier$boot)
      hi <- cost(fine_grid(cg[length(cg) - 1], cg[length(cg)]), tier$mc, tier$boot)
      return(sprintf("Estimated runtime: ~%s to ~%s, depending on where N lands (%d core(s)).",
                     fmt_duration(s1 + lo), fmt_duration(s1 + hi), N_CORES))
    }
    grid <- tryCatch(parse_grid(input$n_grid_text), error = function(e) NULL)
    if (is.null(grid)) return("N grid: enter from,to,by (e.g. 100,300,50)")
    sprintf("Estimated runtime: ~%s (%s PLS fits, %d core(s))",
            fmt_duration(cost(grid, tier$mc, tier$boot)),
            format(round(length(grid) * tier$mc * (1 + tier$boot)), big.mark = ","),
            N_CORES)
  })

  output$elapsed_note <- renderText({
    r <- res(); req(r); sprintf("Finished in %s.", fmt_duration(r$elapsed_sec))
  })

  output$required_n_card <- renderUI({
    r <- res(); req(r)
    hdr <- sprintf("Estimated minimum N for %d%% power - %s",
                   round(100 * r$settings$target), EFFECT_LABEL[[r$settings$effect]])
    body <- if (identical(r$settings$effect, "indirect")) sprintf(
"Built-in bootstrap (a x b):  %s
Sobel test (a x b):          %s

G*Power and the inverse square-root rule apply to single paths,
so they are not shown for the indirect effect.",
      fmt_n(r$min_n_ind_boot), fmt_n(r$min_n_ind_sobel))
    else {
      gp <- r$gpower$N
      dif <- if (!is.na(r$min_n_boot) && !is.na(gp)) sprintf("%+d participants", r$min_n_boot - gp) else "n/a"
      sprintf(
"Built-in bootstrap:        %s
PLS-SEM normal theory:     %s
G*Power regression:        %s
Inverse square-root rule:  %d

Planning-method discrepancy (bootstrap vs G*Power): %s",
        fmt_n(r$min_n_boot), fmt_n(r$min_n_nt), fmt_n(gp), r$inverse_sqrt, dif)
    }
    extra <- list()
    if (!is.null(r$coarse)) extra <- c(extra, list(tags$p(class = "text-muted", sprintf(
      "Find my N: scanned N = %s; full analysis at N = %s.",
      paste(r$coarse$N, collapse = ", "), paste(r$settings$n_grid, collapse = ", ")))))
    if (!is.null(r$search_note)) extra <- c(extra, list(tags$p(
      style = "color:#c0392b;font-weight:600;", r$search_note)))
    card(card_header(hdr), card_body(tags$pre(body), extra))
  })

  output$power_curve <- renderPlot({
    r <- res(); req(r)
    ind <- identical(r$settings$effect, "indirect")
    d <- r$curve; tgt <- r$settings$target
    d$p_nt   <- if (ind) d$power_ind_sobel else d$power_nt
    d$p_boot <- if (ind) d$power_ind_boot  else d$power_boot
    se <- sqrt(d$p_nt * (1 - d$p_nt) / pmax(d$reps, 1))
    ggplot(d, aes(N, p_nt)) +
      geom_ribbon(aes(ymin = pmax(0, p_nt - 1.96 * se), ymax = pmin(1, p_nt + 1.96 * se)),
                  fill = "#4a69bd", alpha = .15) +
      geom_line(color = "#4a69bd") + geom_point(color = "#4a69bd") +
      geom_line(aes(y = p_boot), color = "#e58e26", linetype = "dashed", na.rm = TRUE) +
      geom_hline(yintercept = tgt, color = "#c0392b", linetype = "dashed") +
      { if (!is.null(r$coarse) && !identical(r$curve, r$coarse))
          geom_point(data = transform(r$coarse,
                       p_nt = if (ind) power_ind_sobel else power_nt),
                     aes(N, p_nt), colour = "grey55", shape = 1, size = 2.5) } +
      scale_y_continuous(limits = c(0, 1)) +
      labs(y = "Power", x = "Sample size (N)",
           title = paste("Monte Carlo power curve -", EFFECT_LABEL[[r$settings$effect]]),
           subtitle = if (ind) "Blue: Sobel test; orange dashed: percentile bootstrap of a x b"
                      else "Blue: normal-theory t test; orange dashed: percentile bootstrap",
           caption = paste0(if (!is.null(r$coarse)) "Grey circles: Find-my-N scan (200 replications, no bootstrap).  " else "",
                            APP_VERSION, " - cite: Briones (2026)")) +
      theme_minimal()
  })

  output$diag_table <- renderTable({ r <- res(); req(r); r$curve }, digits = 3)

  output$methods_text <- renderText({ r <- res(); req(r); methods_text(r) })
  output$dl_methods <- downloadHandler(
    filename = "plspower_methods_text.txt",
    content = function(file) writeLines(methods_text(res()), file, useBytes = TRUE))
  output$dl_csv <- downloadHandler(
    filename = "plspower_curve.csv",
    content = function(file) write.csv(res()$curve, file, row.names = FALSE))
  output$dl_json <- downloadHandler(
    filename = "plspower_run.json",
    content = function(file) jsonlite::write_json(res(), file,
                                                  auto_unbox = TRUE, force = TRUE))

  st <- reactive(self_task$result())
  output$self_test_table <- renderTable({
    s <- st(); req(s); s$table
  }, digits = 4)
  output$self_test_verdict <- renderText({
    s <- st(); req(s)
    sprintf("%s | min N (informational): this engine %s, reference %d | %s",
            if (s$pass) "SMOKE TEST PASSED" else "SMOKE TEST FAILED",
            format(s$min_n), s$min_n_ref, s$label)
  })
}

shinyApp(ui, server)
