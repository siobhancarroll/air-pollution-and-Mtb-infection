# Shared support for a planned current-versus-full outcome-adjustment comparison.
# The cache contains participant-level inputs and must stay under ignored Data/.
resolve_hap_analysis_scope <- function(conversion_definition, subset = NULL) {
  exposures <- c(solid_fuel = "solid_fuel_use", wood = "wood_use",
    kerosene = "kerosene_use", solid_primary = "solid_fuel_primary")
  stopifnot(conversion_definition %in% c("primary", "sensitivity"))
  if (!is.null(subset)) {
    stopifnot(length(subset) > 0L, all(subset %in% names(exposures)))
    exposures <- exposures[subset]
  }
  if (conversion_definition == "sensitivity") {
    exposures <- exposures[names(exposures) == "solid_fuel"]
    if (!length(exposures)) stop("Stricter conversion analyses retain only solid_fuel (any polluting fuel); alternative HAP definitions are analysed with primary conversion only")
  }
  exposures
}

resolve_outcome_covariates <- function(weight_covariates, current_covariates,
                                        mode = "current", mandatory = character()) {
  stopifnot(mode %in% c("current", "full"))
  unique(c(mandatory, current_covariates, if (mode == "full") weight_covariates else character()))
}

# Cache both the result and RNG state. Replaying an input must consume the same
# random stream as computing it, including if WeightIt uses random starts.
cache_comparison_input <- function(kind, payload, compute) {
  if (!isTRUE(outcome_model_comparison)) return(compute())
  if (!requireNamespace("digest", quietly = TRUE)) stop("The digest package is required for paired input verification")
  dir.create(outcome_comparison_cache, recursive = TRUE, showWarnings = FALSE)
  before <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  key <- digest::digest(list(kind = kind, payload = payload, seed = before,
    WeightIt = as.character(utils::packageVersion("WeightIt"))), algo = "sha256")
  path <- file.path(outcome_comparison_cache, paste0(kind, "-", key, ".rds"))
  if (file.exists(path)) {
    entry <- readRDS(path)
  } else {
    caught <- character()
    value <- withCallingHandlers(compute(), warning = function(w) {
      caught <<- c(caught, conditionMessage(w)); invokeRestart("muffleWarning")
    })
    entry <- list(value = value, seed_after = get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE),
                  warnings = unique(caught))
    tmp <- tempfile(pattern = "paired-input-", tmpdir = outcome_comparison_cache)
    saveRDS(entry, tmp)
    if (!file.rename(tmp, path)) stop("Could not save comparison input cache")
  }
  if (!is.null(entry$seed_after)) assign(".Random.seed", entry$seed_after, envir = .GlobalEnv)
  for (w in entry$warnings) warning(w, call. = FALSE)
  entry$value
}

comparison_weightit <- function(args) {
  if (!isTRUE(outcome_model_comparison)) return(do.call(WeightIt::weightit, args))
  payload <- args
  payload$formula <- paste(deparse(args$formula), collapse = " ")
  payload$data <- as.list(args$data)
  cache_comparison_input("weights", payload, function() {
    w <- do.call(WeightIt::weightit, args)
    list(weights = w$weights)
  })
}

comparison_lognormal_draw <- function(meanlog, sdlog) {
  cache_comparison_input("pm-draw", list(meanlog = meanlog, sdlog = sdlog),
    function() stats::rlnorm(length(meanlog), meanlog = meanlog, sdlog = sdlog))
}

# Fixed initialization for all cloglog survival models, including ordinary runs.
# Weights and the formula are unchanged; only numerical fitting controls differ.
fit_stable_cloglog <- function(formula, design, family, verify = FALSE) {
  mf <- stats::model.frame(formula, data = design$variables, na.action = stats::na.pass)
  X <- stats::model.matrix(attr(mf, "terms"), mf)
  y <- stats::model.response(mf)
  w <- as.numeric(stats::weights(design, type = "sampling"))
  keep <- !is.na(y) & is.finite(w) & w > 0
  p <- stats::weighted.mean(y[keep], w[keep])
  if (!is.finite(p) || p <= 0 || p >= 1) stop("Cloglog outcome requires both events and non-events")
  start <- setNames(rep(0, ncol(X)), colnames(X))
  if (!"(Intercept)" %in% names(start)) stop("Stable cloglog fitting requires an intercept")
  start["(Intercept)"] <- family$linkfun(p)
  fit <- survey::svyglm(formula, design = design, family = family, start = start,
    control = stats::glm.control(maxit = 100))
  if (!verify || !isTRUE(fit$converged)) return(fit)
  checks <- list()
  # Compare with the original algorithm, a different initialization, and a
  # tighter tolerance. No probe replaces the prespecified primary fit.
  for (probe in c("original", "alternate_start", "tight_tolerance")) {
    caught <- character()
    candidate <- tryCatch(withCallingHandlers({
      args <- list(formula = formula, design = design, family = family)
      if (probe != "original") {
        args$start <- start
        if (probe == "alternate_start") args$start["(Intercept)"] <- start["(Intercept)"] + 0.25
        args$control <- stats::glm.control(maxit = 500,
          epsilon = if (probe == "tight_tolerance") 1e-10 else 1e-8)
      }
      do.call(survey::svyglm, args)
    }, warning = function(w) {
      caught <<- c(caught, conditionMessage(w)); invokeRestart("muffleWarning")
    }), error = function(e) e)
    converged <- !inherits(candidate, "error") && isTRUE(candidate$converged)
    delta_beta <- delta_p <- NA_real_
    if (converged) {
      b <- fit$coefficients; b2 <- candidate$coefficients
      if (identical(names(b), names(b2)) && !anyNA(b) && !anyNA(b2))
        delta_beta <- max(abs(b - b2))
      delta_p <- max(abs(stats::fitted(fit) - stats::fitted(candidate)))
    }
    # Original failures are expected and informational. Every converged
    # original fit must agree; both new numerical probes must pass.
    agrees <- converged && is.finite(delta_beta) && is.finite(delta_p) &&
      delta_beta <= 0.01 && delta_p <= 0.001
    checks[[paste0(probe, "_converged")]] <- converged
    checks[[paste0(probe, "_max_abs_coefficient_change")]] <- delta_beta
    checks[[paste0(probe, "_max_abs_fitted_probability_change")]] <- delta_p
    checks[[paste0(probe, "_warnings")]] <- paste(unique(caught), collapse = " | ")
    checks[[paste0(probe, "_error")]] <- if (inherits(candidate, "error")) conditionMessage(candidate) else ""
    checks[[paste0(probe, "_verification_ok")]] <- if (probe == "original" && !converged) TRUE else agrees
  }
  attr(fit, "numerical_verification") <- checks
  fit
}

# Diagnostics are recorded before strict checks, so an unidentified full model
# is visible in the exports rather than silently dropping aliased covariates.
fit_adjustment_outcome <- function(formula, design, family, prediction_data) {
  fit_model <- function() {
    if (identical(family$link, "cloglog"))
      fit_stable_cloglog(formula, design, family, verify = isTRUE(outcome_model_comparison)) else
      survey::svyglm(formula, design = design, family = family)
  }
  if (!isTRUE(outcome_model_comparison))
    return(fit_model())
  warnings <- character()
  result <- tryCatch(withCallingHandlers(
    fit_model(),
    warning = function(w) {warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")}
  ), error = function(e) e)
  hash <- function(x) digest::digest(x, algo = "sha256")
  diag <- tibble::tibble(outcome_adjustment = outcome_adjustment_mode,
    formula = paste(deparse(formula), collapse = " "),
    # Hash column values rather than model-frame terms/environments, which
    # can differ solely because they reference the current or full formula.
    fit_input_hash = hash(list(data = as.list(design$variables),
      weights = as.numeric(stats::weights(design, type = "sampling")))),
    prediction_population_hash = hash(as.list(prediction_data)),
    n_input_rows = nrow(design$variables), n_prediction_rows = nrow(prediction_data),
    n_prediction_people = if ("PID" %in% names(prediction_data)) dplyr::n_distinct(prediction_data$PID) else nrow(prediction_data),
    warnings = paste(unique(warnings), collapse = " | "))
  fail <- function(msg) {
    error <- simpleError(msg)
    attr(error, "outcome_diagnostics") <- dplyr::mutate(diag, status = "failed", error = msg)
    stop(error)
  }
  if (inherits(result, "error")) fail(conditionMessage(result))
  verification <- attr(result, "numerical_verification")
  if (length(verification)) diag <- dplyr::bind_cols(diag, tibble::as_tibble(verification))
  diag$initialization <- if (identical(family$link, "cloglog")) "weighted_intercept_only" else "default"
  diag$max_iterations <- if (identical(family$link, "cloglog")) 100L else 25L
  diag$iterations <- result$iter
  beta <- stats::coef(result); V <- stats::vcov(result)
  # coef.svyglm() can omit aliased terms. Inspect the requested model matrix
  # and the underlying glm coefficient vector, rather than relying on NA
  # values in the returned coefficient subset to detect rank deficiency.
  X <- tryCatch(stats::model.matrix(stats::delete.response(stats::terms(result)),
    data = stats::model.frame(result), contrasts.arg = result$contrasts), error = function(e) NULL)
  p <- stats::fitted(result)
  diag$n_fitted_rows <- stats::nobs(result)
  observed_response <- stats::model.response(stats::model.frame(result))
  if (is.numeric(observed_response) && all(observed_response %in% c(0, 1))) {
    diag$n_events <- sum(observed_response == 1)
    diag$n_nonevents <- sum(observed_response == 0)
  }
  fitted_indices <- match(rownames(stats::model.frame(result)), rownames(design$variables))
  if (!anyNA(fitted_indices)) {
    fitted_weights <- as.numeric(stats::weights(design, type = "sampling"))[fitted_indices]
    diag$fitted_weight_ess <- sum(fitted_weights)^2 / sum(fitted_weights^2)
    diag$max_fitted_weight_share <- max(fitted_weights) / sum(fitted_weights)
  }
  diag$n_coefficients <- length(beta)
  diag$rank <- result$rank
  diag$n_requested_coefficients <- if (is.null(X)) length(result$coefficients) else ncol(X)
  diag$n_aliased <- max(sum(is.na(beta)), sum(is.na(result$coefficients)),
    diag$n_requested_coefficients - result$rank)
  diag$aliased_terms <- paste(unique(c(names(result$coefficients)[is.na(result$coefficients)],
    if (!is.null(X)) setdiff(colnames(X), names(beta)) else character())), collapse = " | ")
  diag$converged <- isTRUE(result$converged)
  diag$max_abs_coefficient <- if (any(is.finite(beta))) max(abs(beta[is.finite(beta)])) else NA_real_
  se <- sqrt(pmax(base::diag(V), 0))
  diag$max_coefficient_se <- if (any(is.finite(se))) max(se[is.finite(se)]) else NA_real_
  diag$fitted_probability_min <- min(p)
  diag$fitted_probability_max <- max(p)
  diag$n_fitted_near_boundary <- sum(p < 1e-6 | p > 1 - 1e-6)
  diag$condition_number <- if (is.null(X)) NA_real_ else tryCatch(kappa(X), error = function(e) Inf)
  # Opt-in exploration only; exported probes never replace this fitted model.
  exploration_dir <- Sys.getenv("FUEL_CONVERGENCE_DIAGNOSTICS", unset = "")
  if (nzchar(exploration_dir)) {
    exploration <- new.env(parent = environment())
    sys.source(here::here("fuel_convergence_diagnostics.R"), envir = exploration)
    exploration$export_fuel_convergence(result, design, formula, family, diag, exploration_dir)
  }
  if (!diag$converged) fail("Outcome model did not converge")
  if (is.null(X)) fail("Could not reconstruct the outcome model matrix to verify rank")
  if (diag$n_aliased > 0L || result$rank < ncol(X) ||
      anyNA(beta) || any(!is.finite(beta)) || any(!is.finite(V)))
    fail("Outcome model has aliased or non-finite coefficients/covariance; no coefficients were silently dropped")
  checks <- verification[grepl("_verification_ok$", names(verification))]
  if (length(checks) && !all(unlist(checks)))
    fail("Cloglog numerical verification failed: inspect initialization and tolerance comparisons before pooling")
  # Missing outcomes can be excluded during fitting; missing predictors cannot
  # justify changing the standardization population between specifications.
  predictors <- stats::delete.response(stats::terms(result))
  frame <- stats::model.frame(predictors, data = design$variables, na.action = stats::na.pass)
  if (anyNA(frame)) fail("Missing outcome-model predictors would change the analysis population")
  attr(result, "outcome_diagnostics") <- dplyr::mutate(diag, status = "ok", error = NA_character_)
  result
}

prediction_range_monitor <- function() {
  values <- list()
  add <- function(p, domain = "contrasts") {
    if (!isTRUE(outcome_model_comparison)) return(invisible(NULL))
    if (any(!is.finite(p))) stop("Non-finite standardized scenario predictions")
    old <- values[[domain]]
    if (is.null(old)) old <- c(min = Inf, max = -Inf, near_boundary = 0, n_predictions = 0)
    values[[domain]] <<- c(min = min(old[["min"]], p), max = max(old[["max"]], p),
      near_boundary = old[["near_boundary"]] + sum(p < 1e-6 | p > 1 - 1e-6),
      n_predictions = old[["n_predictions"]] + length(p))
    invisible(NULL)
  }
  finish <- function(model) {
    row <- attr(model, "outcome_diagnostics")
    if (is.null(row)) return(NULL)
    for (domain in names(values)) for (name in names(values[[domain]]))
      row[[paste(domain, name, sep = "_")]] <- unname(values[[domain]][[name]])
    row
  }
  list(add = add, finish = finish)
}

save_outcome_diagnostics <- function(fits, label, output_dir, nested = FALSE) {
  if (!isTRUE(outcome_model_comparison)) return(invisible(NULL))
  get_row <- function(x, i, k = NA_integer_) {
    row <- attr(x, "outcome_diagnostics")
    if (is.null(row)) {
      msg <- if (inherits(x, "error")) conditionMessage(x) else attr(x, "fit_error")
      row <- tibble::tibble(status = "failed", error = if (is.null(msg)) "Fit diagnostics missing" else msg)
    }
    dplyr::mutate(row, imputation = i, pm_draw = k, .before = 1)
  }
  rows <- dplyr::bind_rows(lapply(seq_along(fits), function(i) {
    if (nested) dplyr::bind_rows(lapply(seq_along(fits[[i]]), function(k) get_row(fits[[i]][[k]], i, k))) else
      get_row(fits[[i]], i)
  }))
  readr::write_csv(rows, file.path(output_dir, paste0(label, "_outcome_model_diagnostics.csv")))
  invisible(rows)
}

record_outcome_failure <- function(label, error, output_dir) {
  path <- file.path(output_dir, "model_fit_failures.csv")
  readr::write_csv(tibble::tibble(run_label = label, error = conditionMessage(error)),
                    path, append = file.exists(path))
}
