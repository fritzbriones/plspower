# plspower - single-file deployable MVP 0.4.1 (concatenated from modular project)
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

pls_fit <- function(blocks) {
  n <- nrow(blocks[[1]])
  blocks <- lapply(blocks, function(B)
    structure(scale(B), "scaled:centers" = NULL, "scaled:scale" = NULL))
  k <- length(blocks)
  w <- lapply(blocks, function(B) rep(1, ncol(B)))
  scores <- NULL; conv <- FALSE
  for (it in 1:100) {
    scores <- lapply(seq_len(k), function(i) {
      s <- as.numeric(blocks[[i]] %*% w[[i]])
      sd_s <- sd(s)
      if (!is.finite(sd_s) || sd_s < 1e-12) return(NULL)
      (s - mean(s)) / sd_s
    })
    if (any(vapply(scores, is.null, logical(1)))) return(NULL)
    R <- matrix(0, k, k)
    for (i in seq_len(k)) for (j in seq_len(k))
      R[i, j] <- cor(scores[[i]], scores[[j]])
    prox <- list(
      R[1, 2] * scores[[2]],
      R[1, 2] * scores[[1]] + R[2, 3] * scores[[3]],
      R[2, 3] * scores[[2]]
    )
    w_new <- lapply(seq_len(k), function(i)
      apply(blocks[[i]], 2, function(x) cov(x, prox[[i]])))
    nw <- function(v) v / sqrt(sum(v^2))
    delta <- max(abs(unlist(mapply(function(v1, v2)
      nw(v1) - nw(v2), w_new, w, SIMPLIFY = FALSE))))
    w <- w_new
    if (delta < 1e-7) { conv <- TRUE; break }
  }
  if (!conv) return(NULL)
  sc <- scores
  df <- data.frame(x = sc[[1]], m = sc[[2]], y = sc[[3]])
  fit <- lm(y ~ x + m, data = df)
  list(b_direct = unname(coef(fit)[2]),
       se = unname(sqrt(diag(vcov(fit)))[2]),
       r2 = unname(summary(fit)$r.squared),
       converged = conv)
}

one_rep <- function(n, p, boot_reps, alpha) {
  blocks <- gen_data(n, p)
  fit <- pls_fit(blocks)
  if (is.null(fit)) return(list(nt = NA, boot = NA, b = NA))
  nt <- abs(fit$b_direct / fit$se) > qt(1 - alpha / 2, n - 3)
  boot <- NA
  if (boot_reps > 0) {
    idx <- seq_len(n)
    bs <- numeric(boot_reps)
    ok <- 0
    for (bb in seq_len(boot_reps)) {
      rows <- sample(idx, n, replace = TRUE)
      fb <- pls_fit(lapply(blocks, function(B) B[rows, , drop = FALSE]))
      if (!is.null(fb)) { ok <- ok + 1; bs[ok] <- fb$b_direct }
    }
    if (ok >= boot_reps / 2) {
      q <- quantile(bs[seq_len(ok)], c(alpha / 2, 1 - alpha / 2))
      boot <- !(q[1] <= 0 && q[2] >= 0)
    }
  }
  list(nt = nt, boot = boot, b = fit$b_direct)
}

power_curve <- function(p, n_grid, mc_reps, boot_reps, alpha, seed) {
  set.seed(seed)
  out <- data.frame(N = n_grid, power_nt = NA_real_, power_boot = NA_real_,
                    bias = NA_real_, nonconv = 0, reps = 0L)
  for (i in seq_along(n_grid)) {
    n <- n_grid[i]
    nt <- boot <- 0; bs <- c(); nc <- 0; done <- 0; boot_done <- 0
    for (r in seq_len(mc_reps)) {
      res <- one_rep(n, p, boot_reps, alpha)
      if (is.na(res$b)) { nc <- nc + 1; next }
      done <- done + 1
      nt <- nt + as.integer(res$nt)
      if (!is.na(res$boot)) { boot <- boot + as.integer(res$boot); boot_done <- boot_done + 1L }
      bs <- c(bs, res$b)
    }
    out$power_nt[i]   <- if (done) nt / done else NA
    out$power_boot[i] <- if (boot_done) boot / boot_done else NA
    out$bias[i]       <- if (done) mean(bs) - p$c else NA
    out$nonconv[i]    <- nc
    out$reps[i]       <- done
  }
  out
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
# plspower: a simulation PLANNER for PLS-SEM power analysis (MVP 0.4.1)
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

APP_VERSION <- "plspower MVP 0.4.1"
TIERS <- list(
  quick        = list(mc = 200,  boot = 500,   label = "Quick preview (200 MC x 500 boot)"),
  standard     = list(mc = 1000, boot = 2000,  label = "Standard (1,000 MC x 2,000 boot)"),
  publication  = list(mc = 5000, boot = 10000, label = "Publication (5,000 MC x 10,000 boot)")
)

run_power_simulation <- function(settings) {
  p <- list(a = settings$a, b = settings$b, c = settings$c,
            lambda = settings$lambda, n_ind = settings$n_ind)
  if (!identical(settings$engine, "built-in")) {
    stop("Only the built-in engine is currently integrated and supported.")
  }
  curve <- power_curve(p, settings$n_grid, settings$mc_reps,
                       settings$boot_reps, settings$alpha, settings$seed)
  gp <- gpower_required_n(p$a, p$b, p$c, settings$alpha, settings$target)
  isq <- inverse_sqrt_n(p$c, settings$alpha, settings$target)
  list(
    settings = settings, curve = curve,
    min_n_nt = min_n_for_power(settings$n_grid, curve$power_nt, settings$target),
    min_n_boot = min_n_for_power(settings$n_grid, curve$power_boot, settings$target),
    gpower = gp, inverse_sqrt = isq,
    runtime = settings$runtime, version = APP_VERSION,
    engine = settings$engine
  )
}

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
        sliderInput("a", "Indirect path a (X -> M)", 0.05, 0.60, 0.40, 0.05),
        sliderInput("b", "Indirect path b (M -> Y)", 0.05, 0.60, 0.30, 0.05),
        sliderInput("c", "Direct effect c (X -> Y, TARGET)", 0.05, 0.40, 0.20, 0.05),
        sliderInput("lambda", "Indicator loading (all)", 0.50, 0.95, 0.80, 0.05),
        sliderInput("n_ind", "Indicators per construct", 3, 6, 3, 1),
        tags$p(strong("Estimation engine: "), "Built-in Mode A reference engine"),
        tags$p(class = "text-muted", "Independent comparison with seminr/cSEM is planned as a separate validation exercise, not as an in-app engine option.")
      ),
      card_header("2. Analysis settings"),
      card_body(
        selectInput("alpha", "Significance level", c(0.05, 0.01, 0.10)),
        selectInput("target", "Target power", c(0.80, 0.90)),
        selectInput(
          inputId = "tier",
          label = "Compute tier",
          choices = names(TIERS),
          selected = "quick"
        ),
        textOutput("tier_note"),
        uiOutput("tier_banner"),
        numericInput("seed", "Random seed", 2026),
        textInput("n_grid_text", "N grid (from,to,by)", "100,300,10"),
        textOutput("runtime_note")
      )
    ),
    navset_card_tab(
      nav_panel("3. Results",
        card(card_header("Run"),
             card_body(actionButton("run_analysis", "Run power analysis",
                                    class = "btn-primary"))),
        uiOutput("required_n_card"),
        plotOutput("power_curve"),
        card(card_header("Simulation diagnostics"),
             card_body(tableOutput("diag_table"))),
        card(card_header("4. Export"),
             card_body(
               downloadButton("dl_csv", "Power curve (CSV)"),
               downloadButton("dl_json", "Full run (JSON - includes seed,
                              settings, package versions)")
             ))
      ),
      nav_panel("Methods & Validation",
        card(card_body(
          tags$h4("Data-generating model"),
          tags$code("M = a*X + eM;  Y = b*M + c*X + eY (all latents standardized)"),
          tags$p("Indicators: x_ij = lambda*L + sqrt(1-lambda^2)*e_ij; Gaussian errors."),
          tags$h4("Estimator"),
          tags$p("Active estimator: built-in PLS Mode A, path weighting, fixed-point
                 iteration (1e-7), paths via OLS on composites, and percentile
                 bootstrap. seminr/cSEM comparisons are planned separately."),
          tags$h4("Definition of power"),
          tags$p("Proportion of MC replications whose two-sided 95% CI for
                 the direct effect excludes zero."),
          tags$h4("Comparison methods"),
          tags$p("G*Power-equivalent (exact noncentral F, observed-score
                 regression - the regression-shortcut question) and
                 Kock & Hadaya (2018) inverse square-root rule. Differences
                 are reported as 'planning-method discrepancy'."),
          tags$p("Note: the inverse square-root rule is computed two-tailed,
                 consistent with this app's two-sided inference; Kock &
                 Hadaya's published tables are one-tailed (e.g., 155 vs 197
                 at b = 0.20, alpha = 0.05, power = 0.80)."),
          tags$h4("Positioning"),
          tags$p("This tool is not a new estimator. Its contribution is a
                 transparent comparison and calibration layer across
                 planning methods for PLS-SEM: the same study, planned
                 three ways, with measurement quality as an explicit
                 design input."),
          tags$p("Planning methods disagree largely because regression
                 shortcuts and the inverse square-root rule implicitly
                 assume error-free composites. With imperfect indicators
                 (e.g., loading 0.80 with 3 indicators), measurement-error
                 attenuation reduces the effective effect size, so the
                 simulation typically requires a larger N than the
                 shortcuts suggest. The discrepancy you see here is the
                 cost of ignoring measurement quality at the planning stage."),
          tags$h4("Built-in smoke test"),
          tags$p("Reproduces the seed-2026 fixture with explicit tolerance
                 rules (power: 2.5 combined Monte Carlo standard errors;
                 f2: 1e-10 against a full-precision reference;
                 inverse-sqrt: exact). SMOKE TEST: verifies implementation
                 agreement, not
                 publication-grade precision (400 MC reps)."),
          actionButton("self_test", "Run engine smoke test"),
          tableOutput("self_test_table"),
          textOutput("self_test_verdict"),
          tags$h4("Known limitations (v1)"),
          tags$ul(
            tags$li("One mediation template; no formative constructs, MGA,
                    custom diagrams, or user-data fitting (prospective power only)."),
            tags$li("The built-in engine has not yet been independently cross-validated
                    against seminr/cSEM using a prespecified validation protocol."),
            tags$li("Results apply ONLY to the specified model and settings.")
          )
        ))
      )
    )
  ),
  tags$footer(style = "margin-top:12px;color:#888;font-size:0.85em;", APP_VERSION)
)

server <- function(input, output, session) {
  settings <- eventReactive(input$run_analysis, {
    rng <- as.numeric(strsplit(input$n_grid_text, ",")[[1]])
    tier <- TIERS[[input$tier]]
    grid <- seq(rng[1], rng[2], by = rng[3])
    engine <- "built-in"
    list(a = input$a, b = input$b, c = input$c, lambda = input$lambda,
         n_ind = input$n_ind, alpha = as.numeric(input$alpha),
         target = as.numeric(input$target),
         n_grid = grid, mc_reps = tier$mc, boot_reps = tier$boot,
         seed = input$seed, runtime = Sys.time(),
         tier = input$tier, engine = engine,
         versions = as.list(R.version.string))
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

  output$runtime_note <- renderText({
    s <- settings(); req(s)
    est <- length(s$n_grid) * s$mc_reps *
      (1 + if (s$boot_reps > 0) 0.02 * s$boot_reps else 0) * 2e-4
    sprintf("Estimated runtime: ~%d minute(s) on a laptop CPU",
            max(1, round(est / 60)))
  })

  output$required_n_card <- renderUI({
    r <- res(); req(r)
    gp <- r$gpower$N; isq <- r$inverse_sqrt
    dif <- if (!is.na(r$min_n_boot) && !is.na(gp))
      sprintf("%+d participants", r$min_n_boot - gp) else "n/a"
    card(card_header(sprintf("Estimated minimum N for %d%% power",
                             round(100 * r$settings$target))),
         card_body(tags$pre(sprintf(
"Built-in bootstrap:       %s
PLS-SEM normal theory:     %s
G*Power regression:        %d
Inverse square-root rule:  %d

Planning-method discrepancy (bootstrap vs G*Power): %s",
  format(r$min_n_boot), format(r$min_n_nt), gp, isq, dif))))
  })

  output$power_curve <- renderPlot({
    r <- res(); req(r)
    d <- r$curve; tgt <- r$settings$target
    ggplot(d, aes(N, power_nt)) +
      geom_ribbon(aes(ymin = pmax(0, power_nt - 1.96 * sqrt(power_nt *
        (1 - power_nt) / pmax(reps, 1))),
        ymax = pmin(1, power_nt + 1.96 * sqrt(power_nt *
        (1 - power_nt) / pmax(reps, 1)))), fill = "#4a69bd", alpha = .15) +
      geom_line(color = "#4a69bd") + geom_point(color = "#4a69bd") +
      geom_line(aes(y = power_boot), color = "#e58e26", linetype = "dashed") +
      geom_hline(yintercept = tgt, color = "#c0392b", linetype = "dashed") +
      scale_y_continuous(limits = c(0, 1)) +
      labs(y = "Power", x = "Sample size (N)",
           title = "Monte Carlo power curve (simulation-based)",
           subtitle = "Blue: normal-theory test; orange dashed: built-in percentile bootstrap") +
      theme_minimal()
  })

  output$diag_table <- renderTable({ r <- res(); req(r); r$curve }, digits = 3)

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
