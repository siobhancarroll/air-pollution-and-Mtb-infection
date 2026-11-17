# Synthetic checks only; no participant data or imputation.
source(here::here("outcome_adjustment_helpers.R"), local = TRUE)
stopifnot(identical(names(resolve_hap_analysis_scope("primary")),
  c("solid_fuel", "wood", "kerosene", "solid_primary")),
  identical(resolve_hap_analysis_scope("sensitivity"), c(solid_fuel = "solid_fuel_use")),
  identical(names(resolve_hap_analysis_scope("primary", c("wood", "kerosene"))),
    c("wood", "kerosene")),
  identical(resolve_hap_analysis_scope("sensitivity", c("solid_fuel", "solid_primary")),
    c(solid_fuel = "solid_fuel_use")),
  inherits(tryCatch(resolve_hap_analysis_scope("sensitivity", "wood"),
    error = function(e) e), "error"))
stopifnot(identical(resolve_outcome_covariates(c("age", "sex", "entry_season"),
  c("site", "end_season"), "current"), c("site", "end_season")),
  setequal(resolve_outcome_covariates(c("age", "sex", "entry_season"),
    c("site", "end_season"), "full"), c("age", "sex", "entry_season", "site", "end_season")))
# A cached random input must be identical AND leave the same subsequent RNG
# state as computing it. This protects pairing of future PM draws.
outcome_model_comparison <- TRUE
outcome_adjustment_mode <- "current"
outcome_comparison_cache <- tempfile("outcome-pair-test-")
set.seed(812)
first <- cache_comparison_input("synthetic", list(n = 4), function() stats::rnorm(4))
next_first <- stats::runif(3)
set.seed(812)
second <- cache_comparison_input("synthetic", list(n = 4), function() stop("Cache was not reused"))
next_second <- stats::runif(3)
stopifnot(identical(first, second), identical(next_first, next_second))
unlink(outcome_comparison_cache, recursive = TRUE)
# Full outcome adjustment must flag rank deficiency, not drop a coefficient
# and present the reduced model as the requested full specification.
d <- data.frame(PID = 1:12, y = rep(c(0, 1, 0), 4), x = 1:12)
d$z <- 2 * d$x
design <- survey::svydesign(ids = ~1, weights = ~1, data = d)
# A valid model must pass, ensuring that the strict checks are active without
# rejecting every fit. The requested rank-deficient model must then fail.
valid <- fit_adjustment_outcome(y ~ x, design, stats::quasibinomial(), d)
stopifnot(inherits(valid, "svyglm"),
  attr(valid, "outcome_diagnostics")$status == "ok",
  attr(valid, "outcome_diagnostics")$n_aliased == 0L)
error <- tryCatch(fit_adjustment_outcome(y ~ x + z, design,
  stats::quasibinomial(), d), error = function(e) e)
stopifnot(inherits(error, "error"), !is.null(attr(error, "outcome_diagnostics")),
  attr(error, "outcome_diagnostics")$status == "failed",
  attr(error, "outcome_diagnostics")$n_aliased > 0,
  attr(error, "outcome_diagnostics")$n_requested_coefficients == 3L)
# Exercise the survival initialization and independent numerical checks on
# unequal weights; logistic baseline models above retain their original fit.
d$w <- rep(c(0.7, 1.3, 2), 4)
weighted_design <- survey::svydesign(ids = ~1, weights = ~w, data = d)
stable <- fit_adjustment_outcome(y ~ x, weighted_design,
  stats::binomial(link = "cloglog"), d)
stable_diagnostics <- attr(stable, "outcome_diagnostics")
stopifnot(stable_diagnostics$status == "ok",
  stable_diagnostics$initialization == "weighted_intercept_only",
  stable_diagnostics$max_iterations == 100L,
  stable_diagnostics$original_verification_ok,
  stable_diagnostics$alternate_start_verification_ok,
  stable_diagnostics$tight_tolerance_verification_ok)
message("Outcome-adjustment comparison synthetic checks passed")
