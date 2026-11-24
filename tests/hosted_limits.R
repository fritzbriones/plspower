# Checks the hosted-server limits added in v0.8.2 (run from the repo root):
#   Rscript tests/hosted_limits.R
# Hosted mode: Publication tier not offered and refused by the server;
# Quick allowed; Standard allowed only up to HOSTED_MAX_HOURS of estimated runtime.
# Local mode: all tiers offered. Exits with status 1 on any failure.
fails <- 0
check <- function(ok, what) {
  cat(if (isTRUE(ok)) "PASS" else "FAIL", "-", what, "\n")
  if (!isTRUE(ok)) fails <<- fails + 1
}

run_mode <- function(hosted) {
  env <- new.env()
  Sys.setenv(PLSPOWER_HOSTED = if (hosted) "true" else "false")
  suppressMessages(sys.source("app.R", envir = env))
  html <- as.character(env$ui)
  check(identical(env$HOSTED, hosted), sprintf("HOSTED flag is %s", hosted))
  check(identical(names(env$TIERS), c("quick", "standard", "publication")), "TIERS unchanged (batch runner relies on it)")
  check(grepl('<option value="publication"', html, fixed = TRUE) != hosted,
        sprintf("Publication tier %s in the tier menu", if (hosted) "hidden" else "shown"))
  check(grepl("plspower_batch.R --tier publication", html, fixed = TRUE) == hosted,
        sprintf("batch-runner note %s", if (hosted) "shown" else "not shown"))
  if (!hosted) return(invisible())

  try_run <- function(tier, mode, grid = "100,300,50") {
    out <- NULL
    shiny::testServer(env$server, {
      session$setInputs(k = 1, structure = "parallel", a1 = 0.4, b1 = 0.3, c = 0.2,
                        same_meas = TRUE, lambda = 0.8, n_ind = 3, alpha = "0.05",
                        target = "0.8", tier = tier, seed = 2026, mode = mode,
                        n_min = 50, n_max = 1000, n_grid_text = grid, effect = "direct")
      up <- runtime_upper()
      session$setInputs(run_analysis = 1)
      allowed <- tryCatch({ settings(); TRUE }, shiny.silent.error = function(e) FALSE)
      out <<- list(allowed = allowed, hours = up / 3600)
    })
    out
  }
  q <- try_run("quick", "find")
  check(q$allowed, "Quick tier run allowed")
  p <- try_run("publication", "grid")
  check(!p$allowed, "Publication tier run refused")
  for (cfg in list(list("grid", "100,300,50"), list("find", "100,300,50"), list("grid", "100,2000,100"))) {
    s <- try_run("standard", cfg[[1]], cfg[[2]])
    check(s$allowed == (s$hours <= env$HOSTED_MAX_HOURS),
          sprintf("Standard %s %s: est. %.1f h -> %s", cfg[[1]], cfg[[2]], s$hours,
                  if (s$allowed) "allowed" else "refused"))
  }
}

run_mode(TRUE)
run_mode(FALSE)
if (fails > 0) { cat(fails, "check(s) failed\n"); quit(status = 1) }
cat("ALL HOSTED-LIMIT CHECKS PASSED\n")
