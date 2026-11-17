# plspower - single-file deployable v0.8.1 (concatenated from modular project)
# Author: Fritz Jeffry Briones - https://github.com/fritzbriones/plspower
# Cite: Briones, F. J. (2026). plspower (Version 0.8.1) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.24291129
# Deploy: paste into shinyapps.io as app.R, or run with shiny::runApp()
# Requires: shiny, bslib, ggplot2, future, promises, jsonlite

# =============== ENGINE ===============

# run_monte_carlo.R ------------------------------------------------------------
# Core simulation engine for PLS-SEM power analysis (mediation models).
# Models (v0.6.0): X -> 1-3 mediators -> Y, parallel or serial, plus X -> Y.
# DGP: standardized latents in causal order; reflective indicators
#      x_ij = lambda_j*L_j + sqrt(1-lambda_j^2)*e_ij (loadings and indicator
#      counts may differ per construct).
# Estimator: PLS-SEM Mode A, path-weighting inner scheme, fixed-point
#      iteration (1e-7, max 100); paths by OLS on standardized composites.
# Inference per effect: normal theory (t test for the direct effect,
#      delta-method/Sobel z for indirect effects) and percentile bootstrap.
#      Power = share of replications with a significant result.

library(stats)

# ---- Model specification (v0.6.0) -------------------------------------------
# Constructs are ordered X, M1..Mk, Y (k = 1..3 mediators).
#   parallel: X -> Mj (a_j), Mj -> Y (b_j), X -> Y (c)
#   serial:   X -> M1 (a), Mj -> Mj+1 (d), Mk -> Y (b), X -> Y (c)
# B[j, i] is the standardized path i -> j. Loadings and indicator counts may
# differ per construct (scalars are recycled).
build_model <- function(k = 1, structure = c("parallel", "serial"),
                        a, b, c, d = 0.3, lambda = 0.8, n_ind = 3) {
  structure <- match.arg(structure)
  if (k == 1) structure <- "parallel"
  q <- k + 2; X <- 1; Y <- q; M <- 2:(k + 1)
  nm <- c("X", paste0("M", seq_len(k)), "Y")
  B <- matrix(0, q, q, dimnames = list(nm, nm))
  if (structure == "parallel") {
    a <- rep_len(a, k); b <- rep_len(b, k)
    B[cbind(M, X)] <- a; B[cbind(Y, M)] <- b
  } else {
    B[M[1], X] <- a[1]
    if (k > 1) for (j in 2:k) B[M[j], M[j - 1]] <- rep_len(d, k - 1)[j - 1]
    B[Y, M[k]] <- b[1]
  }
  B[Y, X] <- c
  # Paths present in the estimated model (structure), whatever their values.
  S <- matrix(FALSE, q, q, dimnames = dimnames(B))
  if (structure == "parallel") { S[cbind(M, X)] <- TRUE; S[cbind(Y, M)] <- TRUE }
  else { S[M[1], X] <- TRUE; if (k > 1) S[cbind(M[-1], M[-k])] <- TRUE; S[Y, M[k]] <- TRUE }
  S[Y, X] <- TRUE
  # Effects: each is a list of chains (construct index paths from X to Y).
  eff <- list(direct = list(c(X, Y)))
  if (structure == "parallel") {
    if (k == 1) eff$indirect <- list(c(X, M, Y))
    else {
      eff$total <- lapply(M, function(m) c(X, m, Y))
      for (j in seq_len(k)) eff[[paste0("via_M", j)]] <- list(c(X, M[j], Y))
    }
  } else eff$indirect <- list(c(X, M, Y))
  mod <- list(k = k, structure = structure, names = nm, B = B, S = S,
              lambda = rep_len(lambda, q), n_ind = as.integer(rep_len(n_ind, q)),
              effects = eff)
  mod$Sigma <- implied_cor(mod)
  mod$true_effects <- vapply(eff, effect_value, numeric(1), B = B)
  mod
}

# Implied latent correlation matrix; fails if a construct's paths explain
# more than 99% of its variance (impossible standardized parameter values).
implied_cor <- function(mod) {
  B <- mod$B; q <- nrow(B); Sg <- diag(q); psi <- numeric(q); psi[1] <- 1
  for (j in 2:q) {
    for (i in 1:(j - 1)) Sg[j, i] <- Sg[i, j] <- sum(B[j, 1:(j - 1)] * Sg[1:(j - 1), i])
    expl <- drop(t(B[j, 1:(j - 1)]) %*% Sg[1:(j - 1), 1:(j - 1)] %*% B[j, 1:(j - 1)])
    psi[j] <- 1 - expl
    if (psi[j] < 0.01) stop(sprintf(paste0("The paths into %s explain %.0f%% of its variance; ",
      "reduce the path values (standardized paths must leave residual variance)."),
      rownames(B)[j], 100 * expl))
  }
  attr(Sg, "psi") <- psi
  Sg
}

effect_value <- function(chains, B)
  sum(vapply(chains, function(ch) prod(B[cbind(ch[-1], ch[-length(ch)])]), numeric(1)))

# Backward compatibility: v0.5.0 parameter lists (a, b, c, lambda, n_ind).
as_model <- function(p) if (!is.null(p$B)) p else
  build_model(1, "parallel", p$a, p$b, p$c, lambda = p$lambda, n_ind = p$n_ind)

# DGP: standardized latents built in causal order, reflective indicators
# x = lambda * L + sqrt(1 - lambda^2) * e. For simple mediation this draws
# random numbers in the same order as v0.5.0.
gen_data <- function(n, p) {
  mod <- as_model(p); B <- mod$B; q <- nrow(B); psi <- attr(mod$Sigma, "psi")
  L <- matrix(0, n, q); L[, 1] <- rnorm(n)
  for (j in 2:q) {
    pr <- which(B[j, ] != 0)
    L[, j] <- (if (length(pr)) L[, pr, drop = FALSE] %*% B[j, pr] else 0) + sqrt(psi[j]) * rnorm(n)
  }
  lapply(seq_len(q), function(j)
    replicate(mod$n_ind[j], mod$lambda[j] * L[, j] + sqrt(1 - mod$lambda[j]^2) * rnorm(n)))
}

std_cols <- function(B) {
  n <- nrow(B)
  B <- B - rep(colMeans(B), each = n)
  s <- sqrt(colSums(B * B) / (n - 1))
  if (any(!is.finite(s) | s < 1e-12)) return(NULL)
  B / rep(s, each = n)
}

# PLS-SEM: Mode A outer weights, path weighting inner scheme (Lohmoeller):
# successors weighted by correlations, predecessors by regression weights.
# Structural paths by OLS on standardized composites (same as lm()).
pls_fit <- function(blocks, p = NULL, tol = 1e-7, max_it = 100) {
  mod <- if (is.null(p)) as_model(list(a = .1, b = .1, c = .1, lambda = .8, n_ind = 1)) else as_model(p)
  S_str <- mod$S; q <- nrow(S_str); n <- nrow(blocks[[1]])
  preds <- lapply(1:q, function(j) which(S_str[j, ]))
  succs <- lapply(1:q, function(j) which(S_str[, j]))
  blocks <- lapply(blocks, std_cols)
  if (any(vapply(blocks, is.null, logical(1)))) return(NULL)
  w <- lapply(blocks, function(B) rep(1, ncol(B)))
  nw <- function(v) v / sqrt(sum(v * v))
  Sc <- matrix(0, n, q); conv <- FALSE
  for (it in seq_len(max_it)) {
    for (i in 1:q) {
      s <- blocks[[i]] %*% w[[i]]
      sd_s <- sqrt(sum(s * s) / (n - 1))
      if (!is.finite(sd_s) || sd_s < 1e-12) return(NULL)
      Sc[, i] <- s / sd_s
    }
    R <- crossprod(Sc) / (n - 1)
    E <- matrix(0, q, q)
    for (j in 1:q) {
      if (length(succs[[j]])) E[succs[[j]], j] <- R[succs[[j]], j]
      P <- preds[[j]]
      if (length(P)) {
        bj <- tryCatch(solve(R[P, P, drop = FALSE], R[P, j]), error = function(e) NULL)
        if (is.null(bj)) return(NULL)
        E[P, j] <- bj
      }
    }
    Z <- Sc %*% E
    w_new <- lapply(1:q, function(i) drop(crossprod(blocks[[i]], Z[, i])) / (n - 1))
    delta <- max(abs(unlist(Map(function(a, b) nw(a) - nw(b), w_new, w))))
    w <- w_new
    if (delta < tol) { conv <- TRUE; break }
  }
  if (!conv) return(NULL)
  Bhat <- matrix(0, q, q); V <- list(); r2 <- numeric(q)
  for (j in 2:q) {
    P <- preds[[j]]; Rpp <- R[P, P, drop = FALSE]
    beta <- solve(Rpp, R[P, j]); r2[j] <- sum(beta * R[P, j])
    s2 <- (n - 1) * (1 - r2[j]) / (n - 1 - length(P))
    Bhat[j, P] <- beta
    V[[j]] <- list(idx = cbind(j, P), cov = s2 * solve(Rpp) / (n - 1))
  }
  list(B = Bhat, V = V, r2 = r2, n = n, df_y = n - 1 - length(preds[[q]]),
       b_direct = Bhat[q, 1], se = sqrt(V[[q]]$cov[match(1, preds[[q]]), match(1, preds[[q]])]),
       converged = TRUE)
}

# Effect estimate and delta-method SE (exact Sobel SE for a single a x b).
effect_est <- function(fit, chains) {
  B <- fit$B; est <- effect_value(chains, B)
  g <- matrix(0, nrow(B), ncol(B))
  for (ch in chains) {
    e <- cbind(ch[-1], ch[-length(ch)]); v <- B[e]
    for (t in seq_len(nrow(e))) g[e[t, 1], e[t, 2]] <- g[e[t, 1], e[t, 2]] + prod(v[-t])
  }
  var <- 0
  for (blk in fit$V) if (!is.null(blk)) { gv <- g[blk$idx]; var <- var + drop(t(gv) %*% blk$cov %*% gv) }
  c(est = est, se = sqrt(var))
}

# One Monte Carlo replication: for every effect, a normal-theory decision
# (t test for the direct effect, delta-method/Sobel z for indirect effects)
# and a percentile-bootstrap decision.
one_rep <- function(n, p, boot_reps, alpha) {
  mod <- as_model(p); effs <- mod$effects; ne <- length(effs)
  blocks <- gen_data(n, mod)
  fit <- pls_fit(blocks, mod)
  if (is.null(fit)) return(NULL)
  ee <- vapply(effs, function(ch) effect_est(fit, ch), numeric(2))
  crit <- ifelse(names(effs) == "direct", qt(1 - alpha / 2, fit$df_y), qnorm(1 - alpha / 2))
  nt <- abs(ee["est", ] / ee["se", ]) > crit
  boot <- rep(NA, ne)
  if (boot_reps > 0) {
    bs <- matrix(NA_real_, boot_reps, ne); ok <- 0
    for (bb in seq_len(boot_reps)) {
      rows <- sample.int(n, n, replace = TRUE)
      fb <- pls_fit(lapply(blocks, function(B) B[rows, , drop = FALSE]), mod)
      if (!is.null(fb)) { ok <- ok + 1; bs[ok, ] <- vapply(effs, effect_value, numeric(1), B = fb$B) }
    }
    if (ok >= boot_reps / 2) for (e in seq_len(ne)) {
      qq <- quantile(bs[seq_len(ok), e], c(alpha / 2, 1 - alpha / 2), names = FALSE)
      boot[e] <- !(qq[1] <= 0 && qq[2] >= 0)
    }
  }
  list(est = ee["est", ], nt = nt, boot = boot)
}

# Power at one N: for each effect e, columns <e>_nt, <e>_boot, <e>_bias.
power_at_n <- function(n, p, mc_reps, boot_reps, alpha) {
  mod <- as_model(p); en <- names(mod$effects)
  # v0.7.0: when the app sets option plspower.progress_file, every finished
  # replication appends one byte to that file; the app shows a progress bar
  # from the file size. (No effect on results or random numbers.)
  pf <- getOption("plspower.progress_file")
  r <- lapply(seq_len(mc_reps), function(i) {
    out <- one_rep(n, mod, boot_reps, alpha)
    if (!is.null(pf)) try(cat(".", file = pf, append = TRUE), silent = TRUE)
    out
  })
  r <- r[!vapply(r, is.null, logical(1))]
  out <- data.frame(N = n)
  for (k in seq_along(en)) {
    nt <- vapply(r, function(x) x$nt[k], logical(1))
    bt <- vapply(r, function(x) x$boot[k], logical(1))
    es <- vapply(r, function(x) x$est[k], numeric(1))
    out[[paste0(en[k], "_nt")]]   <- if (length(nt)) mean(nt) else NA_real_
    out[[paste0(en[k], "_boot")]] <- if (any(!is.na(bt))) mean(bt, na.rm = TRUE) else NA_real_
    out[[paste0(en[k], "_bias")]] <- if (length(es)) mean(es) - mod$true_effects[[k]] else NA_real_
  }
  out$nonconv <- mc_reps - length(r); out$reps <- length(r)
  out
}

# Number of CPU cores used for a run (forking is not available on Windows).
# v0.6.1: parallelly::availableCores() respects container CPU limits (cgroups);
# detectCores() reports the host machine's cores, which on cloud hosting can
# be far more than the app may use - runs were then much slower than estimated.
plspower_cores <- function() {
  if (.Platform$OS.type == "windows") return(1L)
  n <- if (requireNamespace("parallelly", quietly = TRUE))
         parallelly::availableCores(omit = 0) else parallel::detectCores(logical = FALSE)
  max(1L, min(4L, as.integer(n), na.rm = TRUE))
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
  # Power estimates of exactly 0 or 1 have infinite logits; clamp them so the
  # interpolation stays between the two grid points (v0.6.0 fix: before, a
  # power of 1.00 at the next grid point returned the lower N).
  lg <- function(x) { x <- min(max(x, 0.001), 0.999); log(x / (1 - x)) }
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
# General version for any v0.6.0 model: F test of the R-squared increase from
# X in the population regression of Y on its predictors (X and the mediators
# that point to Y). Identical to gpower_required_n() for simple mediation.
gpower_required_n_model <- function(mod, alpha = 0.05, target = 0.80, n_max = 5000) {
  q <- nrow(mod$S); P <- which(mod$S[q, ]); Sg <- mod$Sigma
  r2 <- function(P) drop(t(Sg[P, q]) %*% solve(Sg[P, P, drop = FALSE], Sg[P, q]))
  r2_full <- r2(P); P0 <- setdiff(P, 1)
  r2_red <- if (length(P0)) r2(P0) else 0
  f2 <- (r2_full - r2_red) / (1 - r2_full); u <- length(P)
  for (n in (u + 2):n_max) {
    df2 <- n - u - 1
    if (1 - pf(qf(1 - alpha, 1, df2), 1, df2, f2 * n) >= target)
      return(list(N = n, f2 = f2, r2_full = r2_full))
  }
  list(N = NA, f2 = f2, r2_full = r2_full)
}

# run_inverse_sqrt.R --------------------------------------------------------------
# Kock & Hadaya (2018) inverse square-root rule. Comparison method only.

inverse_sqrt_n <- function(b_min, alpha = 0.05, target = 0.80) {
  z <- qnorm(1 - alpha / 2) + qnorm(target)
  ceiling(z^2 / b_min^2)
}
# Independent cross-validation against seminr: validation/cross_validate_seminr.R
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
  # v0.7.1: reference power computed with seminr 2.6.0 on the same seed-2026
  # datasets (estimates agree with this engine within 1e-7). The earlier
  # reference (0.432, 0.585, 0.690, 0.795, 0.902; min N 252) came from a Python
  # twin of the v0.4.1 engine, before the path-weighting fix.
  ref_power = c(0.4875, 0.6275, 0.745, 0.7875, 0.8975),
  ref_min_n = 254, ref_f2 = 0.04087591240875911, ref_inv_sqrt = 197
)

run_self_test <- function() {
  p <- list(a = FIXTURE$a, b = FIXTURE$b, c = FIXTURE$c,
            lambda = FIXTURE$lambda, n_ind = FIXTURE$n_ind)
  curve <- power_curve(p, FIXTURE$grid, FIXTURE$reps,
                       boot_reps = 0, alpha = 0.05, FIXTURE$seed)
  pw <- curve$direct_nt
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
    label = "Smoke test: verifies implementation agreement (400 MC reps), not publication-grade precision."
  )
}

# =============== APP ===============

# app.R ---------------------------------------------------------------------------
# plspower: a simulation PLANNER for PLS-SEM power analysis (v0.6.0)
# Prospective power only. Models: X -> 1-3 mediators (parallel or serial) -> Y
# plus X -> Y; reflective Mode A constructs with per-construct loadings and
# indicator counts. Comparison methods reported with neutral wording
# ("planning-method discrepancy"). In-app self-test is a SMOKE TEST.

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

if (!exists("%||%", baseenv())) `%||%` <- function(x, y) if (is.null(x)) y else x

PLSPOWER_VERSION <- "0.8.1"
APP_VERSION <- paste("plspower", PLSPOWER_VERSION)
AUTHOR <- "Fritz Jeffry Briones"
REPO_URL <- "https://github.com/fritzbriones/plspower"
DOI <- "10.5281/zenodo.24291129"
DOI_URL <- paste0("https://doi.org/", DOI)
CITATION_APA <- sprintf(paste0(
  "Briones, F. J. (2026). plspower: A Monte Carlo simulation planner for ",
  "PLS-SEM power analysis (Version %s) [Computer software]. Zenodo. %s"),
  PLSPOWER_VERSION, DOI_URL)
CITATION_INTEXT <- "(Briones, 2026)"
CITATION_BIBTEX <- sprintf(paste0(
  "@software{briones2026plspower,\n",
  "  author  = {Briones, Fritz Jeffry},\n",
  "  title   = {plspower: A Monte Carlo simulation planner for PLS-SEM power analysis},\n",
  "  year    = {2026},\n  version = {%s},\n  publisher = {Zenodo},\n",
  "  doi     = {%s},\n  url     = {%s}\n}"),
  PLSPOWER_VERSION, DOI, DOI_URL)

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

# Time PLS fits for the simple model at N = 200 (9 indicators) once at
# startup, including bootstrap overhead, for the runtime estimate. ~0.1 s.
FIT_MS_N200 <- local({
  p0 <- list(a = .4, b = .3, c = .2, lambda = .8, n_ind = 3)
  set.seed(1); one_rep(200, p0, 5, .05)                    # warm-up
  t <- system.time(one_rep(200, p0, 60, .05))[["elapsed"]]
  max(0.05, 1000 * t / 61)
})
N_CORES <- plspower_cores()

# Run simulations in a background R process so the app stays responsive.
future::plan(future::multisession, workers = 2)

TIERS <- list(
  quick        = list(mc = 200,  boot = 100,   label = "Quick preview (200 MC x 100 boot)"),
  standard     = list(mc = 1000, boot = 2000,  label = "Standard (1,000 MC x 2,000 boot)"),
  publication  = list(mc = 5000, boot = 2000,  label = "Publication (5,000 MC x 2,000 boot)")
)

# ---- Effects and labels -----------------------------------------------------
chain_text <- function(mod, ch) paste(mod$names[ch], collapse = " -> ")
effect_labels <- function(mod) {
  e <- names(mod$effects)
  setNames(vapply(e, function(id) switch(id,
    direct   = "Direct effect c (X -> Y)",
    indirect = if (mod$k == 1) "Indirect effect a x b (X -> M1 -> Y)"
               else sprintf("Serial indirect effect (%s)", chain_text(mod, mod$effects$indirect[[1]])),
    total    = sprintf("Total indirect effect (sum over %d mediators)", mod$k),
    sprintf("Specific indirect effect via %s (%s)", sub("via_", "", id),
            chain_text(mod, mod$effects[[id]][[1]]))), character(1)), e)
}
is_indirect <- function(id) !identical(id, "direct")

# Model from app inputs (also used for validation messages).
model_from <- function(s) {
  q <- s$k + 2
  lam <- if (isTRUE(s$same_meas)) s$lambda else unlist(s$lambda_each)[seq_len(q)]
  ind <- if (isTRUE(s$same_meas)) s$n_ind  else unlist(s$n_ind_each)[seq_len(q)]
  if (anyNA(lam) || any(lam < 0.3 | lam > 0.97)) stop("Loadings must be between 0.30 and 0.97.")
  if (anyNA(ind) || any(ind < 1 | ind > 12)) stop("Indicator counts must be between 1 and 12.")
  pa <- c(s$a, s$b, s$d)
  if (anyNA(pa) || any(pa < 0 | pa > 0.70)) stop("Paths a, b and d must be between 0 and 0.70.")
  if (is.null(s$c) || is.na(s$c) || s$c < 0 || s$c > 0.50) stop("Direct effect c must be between 0 and 0.50.")
  if (identical(s$structure, "serial") && s$k > 1)
    build_model(s$k, "serial", a = s$a[1], b = s$b[1], c = s$c, d = s$d[seq_len(s$k - 1)],
                lambda = lam, n_ind = ind)
  else build_model(s$k, "parallel", a = s$a[seq_len(s$k)], b = s$b[seq_len(s$k)], c = s$c,
                   lambda = lam, n_ind = ind)
}

# Largest direct effect c (0.01 steps, up to 0.50) that keeps the model valid
# with the other settings fixed; NA if no c works (other paths too large).
max_c <- function(s) {
  for (cc in seq(0.50, 0, by = -0.01))
    if (!inherits(tryCatch(model_from(modifyList(s, list(c = cc))), error = function(e) e), "error"))
      return(cc)
  NA_real_
}

# ---- Find my N: two-stage search ------------------------------------------
# Stage 1 scans a wide N range on a log-spaced grid with fast normal-theory
# tests and 200 replications. Stage 2 runs the full tier (bootstrap
# included) on a fine grid around the N where stage-1 power first reached
# the target.
COARSE_POINTS <- 8
COARSE_MC <- 200
FINE_POINTS <- 5
coarse_grid <- function(lo, hi)
  unique(round(exp(seq(log(lo), log(hi), length.out = COARSE_POINTS))))
fine_grid <- function(lower, upper)
  unique(round(seq(max(20, lower), upper * 1.1, length.out = FINE_POINTS)))

col_nt   <- function(e) paste0(e, "_nt")
col_boot <- function(e) paste0(e, "_boot")

run_power_simulation <- function(settings) {
  mod <- settings$model
  if (!identical(settings$engine, "built-in"))
    stop("Only the built-in engine is currently integrated and supported.")
  e <- settings$effect
  t0 <- Sys.time()
  on.exit(options(plspower.progress_file = NULL), add = TRUE)
  coarse <- NULL; search_note <- NULL
  if (identical(settings$mode, "find")) {
    cg <- coarse_grid(settings$n_min, settings$n_max)
    coarse <- power_curve(mod, cg, COARSE_MC, 0, settings$alpha, settings$seed,
                          cores = settings$cores)
    pw <- coarse[[col_nt(e)]]
    i <- which(pw >= settings$target)[1]
    if (is.na(i)) {
      settings$n_grid <- cg
      # Only the fast scan ran: record the settings actually used, keep the
      # requested tier separately so exports (JSON, methods) are accurate.
      settings$requested <- list(mc_reps = settings$mc_reps, boot_reps = settings$boot_reps)
      settings$mc_reps <- COARSE_MC; settings$boot_reps <- 0L
      search_note <- sprintf(paste0("Target power was not reached anywhere up to N = %d ",
        "(highest stage-1 power %.2f). Increase the maximum N or check the effect sizes."),
        settings$n_max, max(pw, na.rm = TRUE))
      curve <- coarse
    } else {
      lower <- if (i == 1) max(20, round(cg[1] / 2)) else cg[i - 1]
      settings$n_grid <- fine_grid(lower, cg[i])
      if (i == 1) search_note <- sprintf(paste0("Power already reached the target at the ",
        "minimum N (%d); the fine grid starts lower. Lower the minimum N to scan further down."), cg[1])
      if (!is.null(settings$progress_file)) options(plspower.progress_file = settings$progress_file)
      curve <- power_curve(mod, settings$n_grid, settings$mc_reps, settings$boot_reps,
                           settings$alpha, settings$seed + 1L, cores = settings$cores)
    }
  } else {
    if (!is.null(settings$progress_file)) options(plspower.progress_file = settings$progress_file)
    curve <- power_curve(mod, settings$n_grid, settings$mc_reps,
                         settings$boot_reps, settings$alpha, settings$seed,
                         cores = settings$cores)
  }
  options(plspower.progress_file = NULL)
  g <- settings$n_grid; tg <- settings$target
  gp <- gpower_required_n_model(mod, settings$alpha, tg)
  list(
    settings = settings, curve = curve,
    min_n_nt   = min_n_for_power(g, curve[[col_nt(e)]], tg),
    min_n_boot = min_n_for_power(g, curve[[col_boot(e)]], tg),
    gpower = gp,
    inverse_sqrt = if (mod$B["Y", "X"] > 0) inverse_sqrt_n(mod$B["Y", "X"], settings$alpha, tg) else NA,
    runtime = settings$runtime,
    elapsed_sec = as.numeric(difftime(Sys.time(), t0, units = "secs")),
    version = APP_VERSION, engine = settings$engine,
    citation = CITATION_APA,
    coarse = coarse, search_note = search_note
  )
}

fmt_n <- function(x) if (is.na(x)) "not reached in grid" else format(x)
num <- function(x, d = 2) formatC(x, format = "f", digits = d)

model_sentence <- function(mod) {
  B <- mod$B; idx <- which(mod$S, arr.ind = TRUE)
  q <- nrow(B); last <- idx[, 1] == q & idx[, 2] == 1
  idx <- idx[order(last, idx[, 1], idx[, 2]), , drop = FALSE]
  paths <- paste(sprintf("%s -> %s = %s", mod$names[idx[, 2]], mod$names[idx[, 1]],
                         num(B[idx])), collapse = ", ")
  kind <- if (mod$k == 1) "a simple mediation model (X -> M1 -> Y with a direct path X -> Y)"
          else sprintf("a %s mediation model with %d mediators and a direct path X -> Y",
                       mod$structure, mod$k)
  same <- length(unique(mod$lambda)) == 1 && length(unique(mod$n_ind)) == 1
  meas <- if (same) sprintf("Each construct was measured reflectively by %d indicators with loadings of %s.",
                            mod$n_ind[1], num(mod$lambda[1]))
          else paste0("Constructs were measured reflectively (",
                      paste(sprintf("%s: %d indicators, loadings %s", mod$names, mod$n_ind,
                                    num(mod$lambda)), collapse = "; "), ").")
  sprintf("Data were generated from %s with standardized path coefficients %s. %s",
          kind, paths, meas)
}

# Ready-to-paste sample-size justification for a manuscript.
methods_text <- function(r) {
  s <- r$settings; mod <- s$model; e <- s$effect; ind <- is_indirect(e)
  scan_only <- !is.null(r$coarse) && identical(r$curve, r$coarse)
  if (scan_only) { s$boot_reps <- 0; s$mc_reps <- COARSE_MC }   # only the fast scan ran
  boot <- s$boot_reps > 0
  min_n <- if (boot) r$min_n_boot else r$min_n_nt
  pw <- r$curve[[if (boot) col_boot(e) else col_nt(e)]]
  if (all(is.na(pw))) pw <- r$curve[[col_nt(e)]]
  test <- if (boot) sprintf(paste0("the two-sided %d%% percentile bootstrap confidence ",
                                   "interval (%s resamples) excluded zero"),
                            round(100 * (1 - s$alpha)), format(s$boot_reps, big.mark = ","))
          else if (ind) "the two-sided delta-method (Sobel) z test was significant"
          else "the two-sided t test was significant"
  grid <- s$n_grid
  out <- paste0(
    "An a priori Monte Carlo power analysis was conducted with plspower (Version ",
    PLSPOWER_VERSION, "; Briones, 2026). ", model_sentence(mod), " Every simulated ",
    "sample was estimated with PLS-SEM (Mode A outer weights, path weighting scheme). ",
    "Statistical power for the ", tolower(substr(effect_labels(mod)[[e]], 1, 1)),
    substring(effect_labels(mod)[[e]], 2),
    " was defined as the proportion of ", format(s$mc_reps, big.mark = ","),
    " Monte Carlo replications in which ", test, " at alpha = ",
    sub("^0", "", num(s$alpha)), ". Sample sizes from N = ", min(grid), " to N = ",
    max(grid), if (length(grid) > 1) paste0(" in steps of ", diff(grid)[1]),
    " were evaluated (random seed ", s$seed, "). ")
  if (scan_only) out <- paste0(out, sprintf(
    "These sample sizes come from a log-spaced scan of N = %d to %d (%d replications, %s); no bootstrap stage was run because the target was not reached. ",
    s$n_min, s$n_max, COARSE_MC, if (ind) "delta-method z test" else "t test"))
  else if (!is.null(r$coarse)) out <- paste0(out, sprintf(paste0(
    "These sample sizes were selected by a two-stage search: a log-spaced scan ",
    "of N = %d to %d (%d replications, %s) located the region where power reached ",
    "the target, and the full analysis was then run on a fine grid in that region. "),
    s$n_min, s$n_max, COARSE_MC, if (ind) "delta-method z test" else "t test"))
  out <- paste0(out, if (!is.na(min_n))
    sprintf("A sample of approximately N = %d was required to reach %d%% power (interpolated from the simulated power curve). ",
            min_n, round(100 * s$target))
    else sprintf("The target power of %d%% was not reached within the evaluated range (maximum simulated power = %s). ",
                 round(100 * s$target), num(max(pw, na.rm = TRUE))))
  if (!ind) out <- paste0(out, sprintf(paste0(
    "For comparison, a regression-based calculation equivalent to G*Power ",
    "(f-squared = %s) indicated N = %s, and the inverse square-root method ",
    "(Kock & Hadaya, 2018; two-tailed) indicated N = %s; these shortcuts assume ",
    "error-free composites."), num(r$gpower$f2, 3), fmt_n(r$gpower$N), fmt_n(r$inverse_sqrt)))
  if (identical(s$tier, "quick"))
    out <- paste0(out, "\n\n[NOTE: Quick-preview tier - re-run with the Standard or Publication tier before reporting.]")
  refs <- c(CITATION_APA,
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
    tags$p(tags$small("DOI: ", tags$a(href = DOI_URL, target = "_blank", DOI))),
    tags$p(tags$small("Created by ", tags$strong(AUTHOR), " - ",
                      tags$a(href = REPO_URL, target = "_blank", REPO_URL)))
  ))

app_header <- tags$div(
  style = paste("display:flex;flex-wrap:wrap;align-items:baseline;gap:4px 14px;",
                "padding:12px 16px;margin-bottom:10px;border-radius:10px;",
                "background:#1f3a5f;color:#fff;"),
  tags$span(style = "font-size:1.6em;font-weight:700;letter-spacing:.5px;", "plspower"),
  tags$span(style = "font-size:1em;opacity:.95;", "Sample-size and power planner for PLS-SEM"),
  tags$span(style = "margin-left:auto;font-size:.85em;opacity:.8;",
            paste0("v", PLSPOWER_VERSION, " \u00b7 by ", AUTHOR)))

ui <- bslib::page_fillable(
  title = "plspower - PLS-SEM power planner",
  theme = NULL,
  fillable_mobile = FALSE,
  tags$style(HTML(paste(
    ".MathJax_Display, mjx-container[display='true'] { overflow-x: auto; overflow-y: hidden; max-width: 100%; }",
    "pre, .shiny-text-output { white-space: pre-wrap !important; word-break: break-word; }",
    ".table-scroll { overflow-x: auto; max-width: 100%; }"))),
  app_header,
  bslib::layout_columns(
    col_widths = c(4, 8),
    card(
      card_header("1. Specify the model"),
      card_body(
        selectInput("k", "Number of mediators", c("1" = 1, "2" = 2, "3" = 3)),
        conditionalPanel("input.k != '1'",
          radioButtons("structure", "Mediator structure",
                       c("Parallel (X -> each M -> Y)" = "parallel",
                         "Serial (X -> M1 -> M2 ... -> Y)" = "serial"))),
        tags$h6("Standardized paths"),
        checkboxInput("exact", "Type exact values (0.01 steps) instead of sliders", FALSE),
        uiOutput("path_inputs"),
        tags$h6("Measurement"),
        checkboxInput("same_meas", "Same loading and number of indicators for all constructs", TRUE),
        conditionalPanel("input.same_meas", uiOutput("same_meas_inputs")),
        conditionalPanel("!input.same_meas", uiOutput("meas_inputs")),
        uiOutput("model_note"),
        uiOutput("effect_select"),
        tags$p(class = "text-muted", style = "font-size:0.85em;",
               "Engine: built-in PLS Mode A, path weighting - cross-checked against seminr (see validation/ on GitHub).")
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
        textOutput("runtime_note"),
        tags$details(style = "margin-top:8px;",
          tags$summary(tags$strong("How long will it take?")),
          tags$p(class = "text-muted", style = "font-size:0.85em;margin:6px 0;",
                 "Typical Find my N run (scan 50-1,000) for the model above, measured on this server.",
                 "The needed N is not known in advance, so three cases are shown."),
          tableOutput("time_guide"),
          tags$p(class = "text-muted", style = "font-size:0.85em;",
                 "Tip: explore with Quick (any device), then run Standard or Publication on a computer",
                 "and keep the page open until it finishes."))
      )
    ),
    navset_card_tab(
      nav_panel("3. Results",
        card(card_header("Run"),
             card_body(actionButton("run_analysis", "Run power analysis",
                                    class = "btn-primary"),
                       uiOutput("run_status"),
                       textOutput("elapsed_note"))),
        uiOutput("required_n_card"),
        plotOutput("power_curve"),
        card(card_header("Simulation diagnostics"),
             card_body(tags$div(class = "table-scroll", tableOutput("diag_table")))),
        card(card_header("4. Methods text for your manuscript"),
             card_body(
               tags$p(class = "text-muted",
                      "A ready-to-paste sample-size justification. Check it against your study before use."),
               verbatimTextOutput("methods_text"),
               uiOutput("dl_methods_ui")
             )),
        card(card_header("5. Export"),
             card_body(uiOutput("export_ui"))),
        cite_card()
      ),
      nav_panel("Methods & Validation",
        card(card_body(
          withMathJax(),
          tags$p(class = "text-muted", "Formulas used by this app, in the order they are applied."),
          tags$h4("1. Model (data generation)"),
          tags$p("Structural model, all latent variables standardized (variance 1), in causal order X, M1..Mk, Y:"),
          helpText("$$\\eta_j = \\sum_{i \\to j} \\beta_{ji}\\,\\eta_i + \\zeta_j$$"),
          helpText("$$\\operatorname{Var}(\\zeta_j) = 1 - \\boldsymbol\\beta_j^{\\top}\\,\\Sigma_{PP}\\,\\boldsymbol\\beta_j$$"),
          tags$p("Simple mediation: \\(M = aX + \\zeta_M\\), \\(Y = cX + bM + \\zeta_Y\\).
                  Parallel: \\(X \\to M_j\\) (\\(a_j\\)), \\(M_j \\to Y\\) (\\(b_j\\)).
                  Serial: \\(X \\to M_1\\) (\\(a\\)), \\(M_j \\to M_{j+1}\\) (\\(d_j\\)), \\(M_k \\to Y\\) (\\(b\\)).
                  A run is blocked if any \\(\\operatorname{Var}(\\zeta_j) < 0.01\\)."),
          tags$p("Reflective indicators (loadings \\(\\lambda_j\\), Gaussian errors):"),
          helpText("$$x_{ij} = \\lambda_j\\,\\eta_j + \\sqrt{1-\\lambda_j^{2}}\\;\\varepsilon_{ij}$$"),
          tags$p("with independent standard normal errors \\(\\varepsilon_{ij}\\)."),
          tags$h4("2. PLS-SEM estimation"),
          tags$p("Indicators are standardized. Composite scores from outer weights \\(\\mathbf w_j\\) (start: all 1):"),
          helpText("$$\\hat\\eta_j = \\frac{X_j \\mathbf w_j}{\\operatorname{sd}(X_j \\mathbf w_j)}$$"),
          helpText("$$\\mathbf z_j = \\sum_i e_{ij}\\,\\hat\\eta_i$$"),
          helpText("$$\\mathbf w_j \\leftarrow \\operatorname{Cov}(X_j, \\mathbf z_j)\\quad(\\text{Mode A})$$"),
          tags$p("Path weighting scheme: \\(e_{ij} = \\operatorname{cor}(\\hat\\eta_i, \\hat\\eta_j)\\) if \\(\\eta_i\\) is a
                  successor of \\(\\eta_j\\); \\(e_{ij}\\) = coefficient of \\(\\hat\\eta_i\\) in the regression of
                  \\(\\hat\\eta_j\\) on its predecessors. Repeat until the weights change by less than \\(10^{-7}\\)
                  (maximum 100 iterations)."),
          tags$p("Structural paths: OLS on the composite correlation matrix \\(R\\), with \\(P\\) the \\(p\\) predictors of \\(\\eta_j\\):"),
          helpText("$$\\hat{\\boldsymbol\\beta}_j = R_{PP}^{-1}\\,\\mathbf r_{Pj}$$"),
          helpText("$$R^2_j = \\hat{\\boldsymbol\\beta}_j^{\\top}\\mathbf r_{Pj}$$"),
          helpText("$$\\widehat{\\operatorname{Var}}(\\hat{\\boldsymbol\\beta}_j) = \\frac{1-R^2_j}{n-1-p}\\,\\,R_{PP}^{-1}$$"),
          tags$h4("3. Significance tests in each simulated sample"),
          tags$p("Direct effect \\(c\\):"),
          helpText("$$t = \\frac{\\hat c}{\\operatorname{SE}(\\hat c)}$$"),
          tags$p("Significant if \\(|t| > t_{1-\\alpha/2,\\,n-1-p}\\)."),
          tags$p("Indirect effect = sum over chains of products of paths (e.g. \\(ab\\), \\(a\\,d\\,b\\), \\(\\sum_j a_j b_j\\)).
                  Delta-method standard error, with \\(\\mathbf g\\) the gradient of the effect and \\(V\\) the path covariance
                  (equal to the Sobel test for a single \\(ab\\)):"),
          helpText("$$\\operatorname{SE}(\\hat a\\hat b) = \\sqrt{\\hat b^{2}\\,\\operatorname{SE}(\\hat a)^{2} + \\hat a^{2}\\,\\operatorname{SE}(\\hat b)^{2}}$$"),
          helpText("$$\\operatorname{SE} = \\sqrt{\\mathbf g^{\\top} V \\mathbf g}$$"),
          tags$p("Significant if \\(|\\text{effect}/\\operatorname{SE}| > z_{1-\\alpha/2}\\)."),
          tags$p("Percentile bootstrap: refit the model on \\(B\\) resamples (rows drawn with replacement);
                  significant if the interval \\([q_{\\alpha/2},\\, q_{1-\\alpha/2}]\\) of the \\(B\\) estimates excludes 0."),
          tags$h4("4. Power and minimum N"),
          helpText("$$\\widehat{\\text{Power}}(N) = \\frac{1}{R}\\sum_{r=1}^{R} \\mathbb 1[\\text{sig}_r]$$"),
          helpText("$$\\operatorname{SE}_{MC} = \\sqrt{\\widehat{\\text{Power}}\\,(1-\\widehat{\\text{Power}})/R}$$"),
          tags$p("Minimum N: first grid point \\(N_1\\) with power \\(\\ge\\) target \\(\\pi^*\\), interpolated on the logit scale
                  from the previous point \\(N_0\\) (powers clamped to [0.001, 0.999]):"),
          helpText("$$N^* = N_0 + t\\,(N_1 - N_0)$$"),
          helpText("$$t = \\frac{\\operatorname{logit}\\pi^* - \\operatorname{logit}\\pi_0}{\\operatorname{logit}\\pi_1 - \\operatorname{logit}\\pi_0}$$"),
          tags$p("Find my N: step 1 scans 8 log-spaced N from the minimum to the maximum (200 replications, no bootstrap);
                  step 2 runs the chosen tier on 5 N values between the last scan point below the target and 1.1 times
                  the first point above it."),
          tags$h4("5. Comparison methods (direct effect only)"),
          tags$p("G*Power-equivalent: exact noncentral F for the increase in \\(R^2\\) from X in the population regression of
                  Y on its \\(u\\) predictors:"),
          helpText("$$f^2 = \\frac{R^2_{\\text{full}} - R^2_{\\text{without }X}}{1 - R^2_{\\text{full}}}$$"),
          helpText("$$F \\sim F(1,\\; N-u-1,\\; \\lambda = f^2 N)$$"),
          tags$p("Inverse square-root rule (Kock & Hadaya, 2018), computed two-tailed:"),
          helpText("$$N_{\\min} = \\left\\lceil \\left(\\frac{z_{1-\\alpha/2} + z_{\\pi^*}}{|\\beta_{\\min}|}\\right)^{2} \\right\\rceil$$"),
          tags$p("Kock & Hadaya's published tables are one-tailed (e.g., 155 vs 197 at \\(\\beta = 0.20\\), \\(\\alpha = 0.05\\),
                  power 0.80). Both shortcuts apply to a single path and assume error-free composites, so they are not
                  shown for indirect effects."),
          tags$h4("6. Speed and reproducibility"),
          tags$p(sprintf("Each N in the grid runs on its own random-number stream (L'Ecuyer-CMRG, derived from the seed),
                 so results are identical whether grid points run in sequence or in parallel. This server uses %d core(s).
                 The engine is cross-checked against seminr (validation/ on GitHub): every path agrees within 2e-7.", N_CORES)),
          tags$h4("Built-in smoke test"),
          tags$p("Reproduces the seed-2026 reference fixture for simple mediation (reference computed with seminr) with
                 explicit tolerance rules (power: 2.5 combined Monte Carlo standard
                 errors; f2: 1e-10; inverse-sqrt: exact). SMOKE TEST: verifies
                 implementation agreement, not publication-grade precision."),
          actionButton("self_test", "Run engine smoke test"),
          uiOutput("self_status"),
          tags$div(class = "table-scroll", tableOutput("self_test_table")),
          textOutput("self_test_verdict"),
          tags$h4("Known limitations"),
          tags$ul(
            tags$li("One X and one Y; up to three mediators; no moderators, formative
                    or higher-order constructs, MGA, or user-data fitting."),
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
  # Let phones reconnect to the same session after a short disconnect
  # (switching apps, screen lock), so a running analysis is not lost.
  session$allowReconnect(TRUE)

  k <- reactive(as.integer(input$k %||% 1))
  structure <- reactive(if (k() > 1) input$structure %||% "parallel" else "parallel")
  val <- function(id, default) { v <- isolate(input[[id]]); if (is.null(v)) default else v }

  # A slider (0.05 steps) or, with "Type exact values", a number box (0.01
  # steps) with the same range. Previous values are kept on redraw.
  num_in <- function(id, label, lo, hi, def, step = 0.05) {
    v <- val(id, def)
    if (isTRUE(input$exact))
      numericInput(id, sprintf("%s [%s to %s]", label, lo, hi), v, min = lo, max = hi, step = 0.01)
    else sliderInput(id, label, lo, hi, v, step)
  }

  # Path inputs for the chosen model, then the direct effect c.
  output$path_inputs <- renderUI({
    kk <- k(); st <- structure()
    sl <- function(id, label, def) num_in(id, label, 0.00, 0.70, def)
    cc <- num_in("c", "Direct effect c (X -> Y)", 0.00, 0.50, 0.20)
    if (kk == 1) return(tagList(sl("a1", "a (X -> M1)", 0.40), sl("b1", "b (M1 -> Y)", 0.30), cc))
    if (st == "parallel") return(tagList(
      lapply(seq_len(kk), function(j) sl(paste0("a", j), sprintf("a%d (X -> M%d)", j, j), 0.40 - 0.1 * (j - 1))),
      lapply(seq_len(kk), function(j) sl(paste0("b", j), sprintf("b%d (M%d -> Y)", j, j), 0.30)), cc))
    tagList(sl("a1", "a (X -> M1)", 0.40),
            lapply(seq_len(kk - 1), function(j) sl(paste0("d", j), sprintf("d%d (M%d -> M%d)", j, j, j + 1), 0.35)),
            sl("b1", sprintf("b (M%d -> Y)", kk), 0.30), cc)
  })

  output$same_meas_inputs <- renderUI(tagList(
    num_in("lambda", "Indicator loading", 0.50, 0.95, 0.80),
    sliderInput("n_ind", "Indicators per construct", 2, 8, val("n_ind", 3), 1)))

  # Per-construct loadings and indicator counts.
  output$meas_inputs <- renderUI({
    nm <- c("X", paste0("M", seq_len(k())), "Y")
    tagList(lapply(nm, function(v) tags$div(style = "display:flex;gap:8px;",
      numericInput(paste0("lam_", v), paste(v, "loading"), val(paste0("lam_", v), 0.80),
                   min = 0.3, max = 0.97, step = 0.01),
      numericInput(paste0("ind_", v), paste(v, "indicators"), val(paste0("ind_", v), 3),
                   min = 1, max = 12, step = 1))))
  })

  model_inputs <- reactive({
    kk <- k(); nm <- c("X", paste0("M", seq_len(kk)), "Y")
    getv <- function(ids, def) vapply(ids, function(id) { v <- input[[id]]; if (is.null(v)) def else as.numeric(v) }, numeric(1))
    list(k = kk, structure = structure(),
         a = getv(paste0("a", 1:3), 0.4), b = getv(paste0("b", 1:3), 0.3),
         d = getv(paste0("d", 1:2), 0.35), c = as.numeric(input$c %||% 0.2),
         same_meas = isTRUE(input$same_meas), lambda = input$lambda %||% 0.8, n_ind = input$n_ind %||% 3,
         lambda_each = getv(paste0("lam_", nm), 0.8), n_ind_each = getv(paste0("ind_", nm), 3))
  })
  model <- reactive(tryCatch(model_from(model_inputs()), error = function(e) e))

  output$model_note <- renderUI({
    m <- model(); s <- model_inputs()
    # Hint on the largest allowed c: shown when the model is fine or blocked
    # only by Y's variance, not for input-range errors.
    bad <- inherits(m, "error")
    mc <- if (!bad || grepl("variance", conditionMessage(m))) max_c(s) else NULL
    c_hint <- if (is.null(mc) || (!is.na(mc) && mc >= 0.50)) NULL
      else if (is.na(mc)) " Even with c = 0 this is not possible: lower a, b or d."
      else sprintf(" With these paths, c can be at most %.2f.", mc)
    if (bad) return(tags$p(style = "color:#c0392b;font-weight:600;", paste0(conditionMessage(m), c_hint)))
    q <- nrow(m$S); P <- which(m$S[q, ])
    r2y <- drop(t(m$Sigma[P, q]) %*% solve(m$Sigma[P, P, drop = FALSE], m$Sigma[P, q]))
    eff <- m$true_effects[names(m$true_effects) != "direct"]
    nice <- function(id) if (id == "indirect") (if (m$k == 1) "a x b" else "serial") else sub("via_", "via ", id)
    names(eff) <- vapply(names(eff), nice, character(1))
    tags$p(class = "text-muted", style = "font-size:0.85em;",
           sprintf("Model OK. R-squared of Y = %.2f. True indirect effect%s: %s.", r2y,
                   if (length(eff) > 1) "s" else "",
                   paste(sprintf("%s = %.3f", names(eff), eff), collapse = ", ")),
           c_hint)
  })

  output$effect_select <- renderUI({
    m <- model(); if (inherits(m, "error")) return(NULL)
    labs <- effect_labels(m)
    sel <- isolate(input$effect); if (is.null(sel) || !sel %in% names(labs)) sel <- names(labs)[1]
    radioButtons("effect", "Target effect", setNames(names(labs), labs), selected = sel)
  })

  settings <- eventReactive(input$run_analysis, {
    m <- model()
    if (inherits(m, "error")) { showNotification(conditionMessage(m), type = "error"); req(FALSE) }
    grid <- tryCatch(
      if (identical(input$mode, "find")) { check_range(input$n_min, input$n_max)
                                           coarse_grid(input$n_min, input$n_max) }
      else parse_grid(input$n_grid_text),
      error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL })
    req(grid)
    eff <- input$effect; if (is.null(eff) || !eff %in% names(m$effects)) eff <- "direct"
    tier <- TIERS[[input$tier]]
    list(model = m, model_inputs = model_inputs(), effect = eff,
         mode = input$mode, n_min = input$n_min, n_max = input$n_max,
         alpha = as.numeric(input$alpha), target = as.numeric(input$target),
         n_grid = grid, mc_reps = tier$mc, boot_reps = tier$boot,
         seed = input$seed, runtime = Sys.time(), cores = N_CORES,
         progress_file = progress_file,
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

  # Live status while a run is going (results arrive only at the end).
  run_started <- reactiveVal(NULL)
  observeEvent(input$run_analysis, run_started(Sys.time()), priority = 10)
  progress_file <- tempfile("plspower_progress_", tmpdir = tempdir())
  stage2_started <- reactiveVal(NULL)
  observeEvent(input$run_analysis, { unlink(progress_file); stage2_started(NULL) }, priority = 20)
  session$onSessionEnded(function() unlink(progress_file))

  progress_bar <- function(frac, label) tags$div(style = "margin:8px 0 4px;",
    tags$div(style = "height:14px;background:#dde5f2;border-radius:7px;overflow:hidden;",
      tags$div(style = sprintf("height:100%%;width:%.1f%%;background:#4a69bd;transition:width .8s;", 100 * frac))),
    tags$small(label))
  output$run_status <- renderUI({
    st <- power_task$status()
    if (st == "running") {
      invalidateLater(1000)
      el <- as.numeric(difftime(Sys.time(), run_started() %||% Sys.time(), units = "secs"))
      s <- isolate(tryCatch(settings(), error = function(e) NULL))
      done <- if (file.exists(progress_file)) file.size(progress_file) else 0
      total <- if (!is.null(s)) (if (identical(s$mode, "find")) FINE_POINTS else length(s$n_grid)) * s$mc_reps else NA
      bar <- if (is.na(total)) NULL
        else if (done == 0 && identical(s$mode, "find"))
          progress_bar(0, "Step 1 of 2: scanning the N range (fast)...")
        else {
          if (is.null(isolate(stage2_started()))) stage2_started(Sys.time())
          frac <- min(1, done / total)
          t2 <- as.numeric(difftime(Sys.time(), isolate(stage2_started()), units = "secs"))
          eta <- if (frac > 0.03) sprintf(" - about %s left", fmt_duration(t2 * (1 - frac) / frac)) else ""
          progress_bar(frac, sprintf("%s%s replications done (%d%%)%s",
            if (identical(s$mode, "find")) "Step 2 of 2: " else "",
            format(done, big.mark = ","), round(100 * frac), eta))
        }
      tags$div(style = "margin-top:8px;padding:8px;border-radius:6px;background:#eef3fb;",
        tags$strong(sprintf("Running... %d:%02d elapsed.", floor(el / 60), round(el %% 60))),
        bar, isolate(output_runtime_text()),
        tags$br(), tags$small("Keep this page open and your screen on. On phones, switching apps or",
                              "locking the screen can pause the page; results appear here when done."))
    } else if (st == "error") {
      msg <- tryCatch({ power_task$result(); "" }, error = function(e) conditionMessage(e))
      tags$p(style = "color:#c0392b;font-weight:600;margin-top:8px;", "The run stopped with an error: ", msg)
    }
  })
  output$self_status <- renderUI({
    st <- self_task$status()
    if (st == "running") tags$p(style = "margin-top:6px;", tags$strong("Smoke test running... usually 15-60 seconds. Keep this page open."))
    else if (st == "error") tags$p(style = "color:#c0392b;", "Smoke test error: ",
      tryCatch({ self_task$result(); "" }, error = function(e) conditionMessage(e)))
  })

  # Downloads appear only once there is a result (tapping them earlier gave an error).
  output$dl_methods_ui <- renderUI({
    if (power_task$status() != "success") return(tags$p(class = "text-muted", "Run the analysis to get the methods text."))
    downloadButton("dl_methods", "Download methods text (TXT)")
  })
  output$export_ui <- renderUI({
    if (power_task$status() != "success") return(tags$p(class = "text-muted", "Run the analysis to enable downloads."))
    tagList(
      downloadButton("dl_csv", "Power curve (CSV)"),
      tags$p(class = "text-muted", style = "font-size:0.85em;margin-top:6px;",
             "CSV: the power values for each N, for Excel or your own charts."),
      downloadButton("dl_json", "Full run record (JSON)"),
      tags$p(class = "text-muted", style = "font-size:0.85em;margin-top:6px;",
             "JSON: an exact record of the model, settings, random seed, software versions and results,",
             "so anyone can reproduce your analysis. Keep it with your data or add it to supplementary materials."))
  })

  output$tier_note <- renderText(paste("Tier:", TIERS[[input$tier]]$label))
  output$tier_banner <- renderUI({
    if (input$tier == "quick")
      tags$p(style = "color:#c0392b;font-weight:600;",
             "Quick tier: DO NOT USE FOR REPORTING. Use Standard or
              Publication tier for any reported estimate.")
  })

  # Runtime estimate shown BEFORE the run: N grid x MC reps x (1 + bootstrap)
  # PLS fits, fit time measured on this server and scaled by model size.
  output_runtime_text <- reactive({
    m <- model(); if (inherits(m, "error")) return("")
    tier <- TIERS[[input$tier]]
    size <- sum(m$n_ind) / 9 * (nrow(m$S) / 3)
    cost <- function(grid, mc, boot) sum(grid / 200) * mc * (1 + boot) *
      FIT_MS_N200 / 1000 * size / min(N_CORES, length(grid))
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
  output$runtime_note <- renderText(output_runtime_text())

  # "How long will it take?" - Find my N (scan 50-1,000) for the current model,
  # with the needed N around 150, 300 or 600, for each tier.
  output$time_guide <- renderTable({
    m <- model(); if (inherits(m, "error")) return(NULL)
    size <- sum(m$n_ind) / 9 * (nrow(m$S) / 3)
    cost <- function(grid, mc, boot) sum(grid / 200) * mc * (1 + boot) *
      FIT_MS_N200 / 1000 * size / min(N_CORES, length(grid))
    scan <- cost(coarse_grid(50, 1000), COARSE_MC, 0)
    cases <- c(150, 300, 600)
    out <- data.frame(Tier = c("Quick (preview)", "Standard (report)", "Publication"))
    for (nn in cases) out[[sprintf("N about %d", nn)]] <- vapply(TIERS, function(t)
      fmt_duration(scan + cost(fine_grid(round(nn * 0.8), nn), t$mc, t$boot)), character(1))
    out
  }, striped = TRUE, spacing = "xs")

  output$elapsed_note <- renderText({
    r <- res(); req(r); sprintf("Finished in %s.", fmt_duration(r$elapsed_sec))
  })

  output$required_n_card <- renderUI({
    r <- res(); req(r)
    m <- r$settings$model; e <- r$settings$effect
    hdr <- sprintf("Estimated minimum N for %d%% power - %s",
                   round(100 * r$settings$target), effect_labels(m)[[e]])
    body <- if (is_indirect(e)) sprintf(
"Built-in bootstrap:          %s
Delta-method (Sobel) z test: %s

G*Power and the inverse square-root rule apply to single paths,
so they are not shown for indirect effects.",
      fmt_n(r$min_n_boot), fmt_n(r$min_n_nt))
    else {
      gp <- r$gpower$N
      dif <- if (!is.na(r$min_n_boot) && !is.na(gp)) sprintf("%+d participants", r$min_n_boot - gp) else "n/a"
      sprintf(
"Built-in bootstrap:        %s
PLS-SEM normal theory:     %s
G*Power regression:        %s
Inverse square-root rule:  %s

Planning-method discrepancy (bootstrap vs G*Power): %s",
        fmt_n(r$min_n_boot), fmt_n(r$min_n_nt), fmt_n(gp), fmt_n(r$inverse_sqrt), dif)
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
    e <- r$settings$effect; m <- r$settings$model; ind <- is_indirect(e)
    d <- r$curve; tgt <- r$settings$target
    d$p_nt <- d[[col_nt(e)]]; d$p_boot <- d[[col_boot(e)]]
    se <- sqrt(d$p_nt * (1 - d$p_nt) / pmax(d$reps, 1))
    ggplot(d, aes(N, p_nt)) +
      geom_ribbon(aes(ymin = pmax(0, p_nt - 1.96 * se), ymax = pmin(1, p_nt + 1.96 * se)),
                  fill = "#4a69bd", alpha = .15) +
      geom_line(color = "#4a69bd") + geom_point(color = "#4a69bd") +
      geom_line(aes(y = p_boot), color = "#e58e26", linetype = "dashed", na.rm = TRUE) +
      geom_hline(yintercept = tgt, color = "#c0392b", linetype = "dashed") +
      { if (!is.null(r$coarse) && !identical(r$curve, r$coarse))
          geom_point(data = data.frame(N = r$coarse$N, p_nt = r$coarse[[col_nt(e)]]),
                     aes(N, p_nt), colour = "grey55", shape = 1, size = 2.5) } +
      scale_y_continuous(limits = c(0, 1)) +
      { if (!is.null(r$coarse) && !identical(r$curve, r$coarse)) {
          cn <- r$coarse$N; lo <- max(c(cn[cn < min(d$N)], min(d$N)))
          hi <- min(c(cn[cn > max(d$N)], max(d$N)))
          coord_cartesian(xlim = c(lo, hi)) } } +
      labs(y = "Power", x = "Sample size (N)",
           title = paste(strwrap(paste("Power curve -", effect_labels(m)[[e]]), 42), collapse = "\n"),
           subtitle = if (ind) "Blue: delta-method (Sobel) z test\nOrange dashed: percentile bootstrap"
                      else "Blue: normal-theory t test\nOrange dashed: percentile bootstrap",
           caption = paste0(if (!is.null(r$coarse)) "Grey circles: Find-my-N scan (200 reps, no bootstrap)\n" else "",
                            APP_VERSION, " - cite: Briones (2026)")) +
      theme_minimal()
  })

  output$diag_table <- renderTable({
    r <- res(); req(r); d <- r$curve
    d <- d[, c("N", grep("_(nt|boot|bias)$", names(d), value = TRUE), "nonconv", "reps")]
    for (k in c("N", "nonconv", "reps")) d[[k]] <- as.integer(round(d[[k]]))
    d
  }, digits = 3)

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
