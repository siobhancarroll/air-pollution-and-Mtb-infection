# Aggregate diagnostics only. Called before the original strict convergence guard.
export_fuel_convergence <- function(model, design, formula, family, diagnostics, directory) {
  if (!"wood_use" %in% all.vars(formula)) return(invisible(NULL))
  key <- diagnostics$fit_input_hash
  prefix <- file.path(directory, key)
  mf <- stats::model.frame(model)
  y <- stats::model.response(mf)
  idx <- match(rownames(mf), rownames(design$variables))
  if (anyNA(idx)) stop("Cannot align fitted rows for convergence diagnostics")
  w <- as.numeric(stats::weights(design, type = "sampling"))[idx]
  coefficients <- function(fit, probe) {
    b <- fit$coefficients
    v <- stats::vcov(fit)
    se <- rep(NA_real_, length(b)); names(se) <- names(b)
    se[rownames(v)] <- sqrt(pmax(diag(v), 0))
    data.frame(term = names(b), estimate = unname(b), se = unname(se), probe = probe)
  }
  readr::write_csv(coefficients(model, "original"), paste0(prefix, "_coefficients.csv"))
  counts <- lapply(setdiff(names(mf), names(mf)[1]), function(variable) {
    x <- mf[[variable]]
    if (is.numeric(x) && length(unique(x)) > 10L) return(NULL)
    do.call(rbind, lapply(unique(as.character(x)), function(level) {
      take <- as.character(x) == level
      data.frame(variable = variable, level = level, rows = sum(take),
        events = sum(y[take] == 1), nonevents = sum(y[take] == 0),
        weight_sum = sum(w[take]), event_weight = sum(w[take] * y[take]),
        ess = sum(w[take])^2 / sum(w[take]^2))
    }))
  })
  readr::write_csv(dplyr::bind_rows(counts), paste0(prefix, "_category_counts.csv"))
  readr::write_csv(diagnostics, paste0(prefix, "_original.csv"))
  if (isTRUE(model$converged)) return(invisible(NULL))
  # Probe both iteration limit and initialization, keeping family, design,
  # formula and weights fixed. These fits are never used for pooling.
  X <- stats::model.matrix(model)
  start <- rep(0, ncol(X)); names(start) <- colnames(X)
  if ("(Intercept)" %in% names(start))
    start["(Intercept)"] <- family$linkfun(stats::weighted.mean(y, w))
  for (probe in c("default_start_100", "intercept_start_100")) {
    warnings <- character()
    fit <- tryCatch(withCallingHandlers({
      args <- list(formula = formula, design = design, family = family,
        control = stats::glm.control(maxit = 100))
      if (probe == "intercept_start_100") args$start <- start
      do.call(survey::svyglm, args)
    }, warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    }), error = function(e) e)
    record <- data.frame(probe = probe, converged = FALSE, iterations = NA_integer_,
      max_abs_coefficient = NA_real_, near_boundary = NA_integer_,
      warnings = paste(unique(warnings), collapse = " | "), error = "")
    if (inherits(fit, "error")) record$error <- conditionMessage(fit) else {
      record$converged <- isTRUE(fit$converged)
      record$iterations <- fit$iter
      record$max_abs_coefficient <- max(abs(fit$coefficients), na.rm = TRUE)
      p <- stats::fitted(fit)
      record$near_boundary <- sum(p < 1e-6 | p > 1 - 1e-6)
      readr::write_csv(coefficients(fit, probe), paste0(prefix, "_", probe, "_coefficients.csv"))
    }
    readr::write_csv(record, paste0(prefix, "_", probe, "_probe.csv"))
  }
}
