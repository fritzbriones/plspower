#!/usr/bin/env Rscript
# plspower_batch.R - headless plspower runner (v0.8.2; engine unchanged since v0.8.0).
#
# Runs the same Monte Carlo power engine as the Shiny app, without any
# hosting limits: use it for Standard or Publication tier runs on your
# own laptop/cluster. Sources app.R for the verified engine
# (cross-validated against seminr, see validation/), then writes the
# same exports the app produces.
#
# Usage:
#   Rscript plspower_batch.R [--k 1] [--structure parallel|serial]
#        [--a 0.40] [--b 0.30] [--d 0.35] [--c 0.20]
#        [--lambda 0.80] [--n-ind 3]
#        [--effect direct|indirect|total|via_M1,...]
#        [--mode find|grid] [--n-min 50] [--n-max 1000] [--grid 100,300,50]
#        [--alpha 0.05] [--target 0.80]
#        [--tier quick|standard|publication|custom] [--mc 200] [--boot 100]
#        [--seed 2026] [--cores -1] [--out plspower_run]
#
#   Vector inputs take comma-separated values, e.g. --a 0.40,0.30,0.20
#   for k = 3 parallel. --tier custom uses --mc and --boot directly.
#   Default effect: direct. Default settings equal the app defaults.
#
# Outputs (in the working directory, unless --out gives another prefix):
#   <out>_curve.csv   power curve (same columns as the app's CSV export)
#   <out>_methods.txt ready-to-paste sample-size justification
#   <out>_run.json    full reproducible run record (same as the app's JSON)

# ---- locate app.R next to this script --------------------------------------
.argv0 <- commandArgs(trailingOnly = FALSE)
.f <- grep("^--file=", .argv0, value = TRUE)
.app_dir <- if (length(.f)) dirname(normalizePath(sub("^--file=", "", .f[1]))) else getwd()
if (!file.exists(file.path(.app_dir, "app.R"))) {
  cat("Error: app.R not found next to plspower_batch.R.\n",
      "Keep the two files in the same directory.\n", sep = "")
  quit(status = 1)
}
source(file.path(.app_dir, "app.R"), echo = FALSE)

# ---- argument parsing -------------------------------------------------------
.args <- commandArgs(trailingOnly = TRUE)
.arg_val <- function(flag) {
  i <- match(flag, .args)
  if (is.na(i) || i == length(.args)) return(NULL)
  v <- .args[i + 1]
  if (grepl("^--", v)) return(NULL)
  v
}
.get_num <- function(flag, default) { v <- .arg_val(flag); if (is.null(v)) default else suppressWarnings(as.numeric(v)) }
.get_chr <- function(flag, default) { v <- .arg_val(flag); if (is.null(v)) default else v }
.get_vec <- function(flag, default) {
  v <- .arg_val(flag)
  if (is.null(v)) return(default)
  suppressWarnings(as.numeric(strsplit(v, ",", fixed = TRUE)[[1]]))
}

k         <- as.integer(.get_num("--k", 1))
structure <- .get_chr("--structure", "parallel")
a_vec     <- .get_vec("--a", 0.40)
b_vec     <- .get_vec("--b", 0.30)
d_vec     <- .get_vec("--d", 0.35)
c_dir     <- .get_num("--c", 0.20)
lambda    <- .get_num("--lambda", 0.80)
n_ind     <- as.integer(.get_num("--n-ind", 3))
effect_id <- .get_chr("--effect", "direct")
mode      <- .get_chr("--mode", "find")
n_min     <- as.integer(.get_num("--n-min", 50))
n_max     <- as.integer(.get_num("--n-max", 1000))
grid_txt  <- .get_chr("--grid", "100,300,50")
alpha     <- .get_num("--alpha", 0.05)
target    <- .get_num("--target", 0.80)
tier      <- .get_chr("--tier", "quick")
mc_reps   <- as.integer(.get_num("--mc", 200))
boot_reps <- as.integer(.get_num("--boot", 100))
seed      <- as.integer(.get_num("--seed", 2026))
out       <- .get_chr("--out", "plspower_run")

# ---- build the model exactly as the app does --------------------------------
spec <- list(k = k, structure = structure, same_meas = TRUE,
             lambda = lambda, n_ind = n_ind,
             a = a_vec, b = b_vec, d = d_vec, c = c_dir)
mod <- tryCatch(model_from(spec),
                error = function(e) {
                  cat("Model error: ", conditionMessage(e), "\n", sep = "")
                  quit(status = 1)
                })
if (!effect_id %in% names(mod$effects)) {
  cat("Error: unknown --effect '", effect_id, "'. This model supports: ",
      paste(names(mod$effects), collapse = ", "), "\n", sep = "")
  quit(status = 1)
}

# ---- tier bookkeeping (same TIER table as the app) ---------------------------
tiers_note <- switch(tier,
  quick      = "Quick preview",
  standard   = "Standard",
  publication = "Publication",
  custom     = "Custom",
  { cat("Error: unknown --tier '", tier,
        "'. Use quick, standard, publication or custom.\n", sep = ""); quit(status = 1) })
if (tier %in% c("quick", "standard", "publication")) {
  mc_reps  <- as.integer(TIERS[[tier]]$mc)
  boot_reps <- as.integer(TIERS[[tier]]$boot)
}
if (tier == "custom" && (mc_reps < 1 || boot_reps < 0)) {
  cat("Error: custom tier needs --mc >= 1 and --boot >= 0.\n"); quit(status = 1)
}

# assemble settings (same contract run_power_simulation expects) ---------------
.get_grid <- function() {
  if (identical(mode, "find")) { check_range(n_min, n_max); coarse_grid(n_min, n_max) }
  else if (identical(mode, "grid")) parse_grid(grid_txt)
  else { cat("Error: --mode must be 'find' or 'grid'.\n"); quit(status = 1) }
}
n_grid <- tryCatch(.get_grid(),
                   error = function(e) {
                     cat("Sample-size error: ", conditionMessage(e), "\n", sep = "")
                     quit(status = 1)
                   })

settings <- list(
  model = mod, effect = effect_id, mode = mode,
  n_min = n_min, n_max = n_max, n_grid = n_grid,
  alpha = alpha, target = target,
  mc_reps = mc_reps, boot_reps = boot_reps,
  seed = seed, runtime = Sys.time(), cores = plspower_cores(),
  tier = tier, engine = "built-in",
  plspower_version = PLSPOWER_VERSION,
  versions = c(list(R = R.version.string),
              lapply(setNames(REQ_PKGS, REQ_PKGS),
                     function(p) as.character(utils::packageVersion(p))))
)

cat(sprintf("plspower %s | %s | %s tier (%d MC x %d boot) | %s | seed %d\n",
            PLSPOWER_VERSION, model_sentence(mod), tiers_note,
            mc_reps, boot_reps,
            if (identical(mode, "find")) sprintf("Find my N: %d-%d", n_min, n_max)
            else paste(n_grid, collapse = ","),
            seed))
t0 <- Sys.time()

r <- tryCatch(run_power_simulation(settings),
              error = function(e) {
                cat("Simulation error: ", conditionMessage(e), "\n", sep = "")
                quit(status = 1)
              })

# ---- report (mirrors the app's results card) --------------------------------
cat(sprintf("Finished in %s.\n\n", fmt_duration(as.numeric(difftime(Sys.time(), t0, units = "secs")))))
if (!is.null(r$search_note)) cat(r$search_note, "\n\n")
e <- effect_id
cat("Target effect: ", effect_labels(mod)[[e]], "\n", sep = "")
cat("  Normal-theory power (N, power): ",
    paste(sprintf("%d: %s", r$curve$N, num(r$curve[[col_nt(e)]])), collapse = ", "), "\n", sep = "")
if (boot_reps > 0)
  cat("  Bootstrap power     (N, power): ",
      paste(sprintf("%d: %s", r$curve$N, num(r$curve[[col_boot(e)]])), collapse = ", "), "\n", sep = "")
cat("  Minimum N for ", sprintf("%.0f%%", 100 * target), " power: ", fmt_n(r$min_n_nt),
    " (normal-theory)", sep = "")
if (boot_reps > 0) cat(", ", fmt_n(r$min_n_boot), " (bootstrap)", sep = "")
cat("\n")
if (identical(e, "direct")) {
  cat("  Planning comparison (direct effect only):\n")
  cat("    G*Power-equivalent N: ", fmt_n(r$gpower$N),
      " | inverse square-root rule N: ", fmt_n(r$inverse_sqrt), "\n\n", sep = "")
  
} else cat("\n")

cat("---- Methods text (same paragraph the app exports) ----\n")
cat(methods_text(r), "\n", sep = "")
cat("----\n")

# ---- exports (same formats the app exports) ----------------------------------
write.csv(r$curve, paste0(out, "_curve.csv"), row.names = FALSE)
writeLines(methods_text(r), paste0(out, "_methods.txt"), useBytes = TRUE)
jsonlite::write_json(r, paste0(out, "_run.json"), auto_unbox = TRUE, force = TRUE)

cat("Wrote: ", out, "_curve.csv, ", out, "_methods.txt, ", out, "_run.json\n", sep = "")
cat("BATCH RUN COMPLETED\n")
