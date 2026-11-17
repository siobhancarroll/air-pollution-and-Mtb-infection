# Helpers shared by propensity_scores.R (baseline models) and
# discrete_time_survival.R (conversion models): data loading, covariate
# collapsing, covariate diagnostics, balance checks and multiple-imputation
# pooling.

library(here)
library(dplyr)
library(forcats)
library(tibble)
library(WeightIt)
library(cobalt)
library(car)

##############################################################################
### Data loading                                                          ###
##############################################################################
# Loads the most recently dated file in Data/ that matches `pattern`.
load_latest <- function(pattern, dir = here("Data")) {
  files <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) stop("No files matching pattern: ", pattern)

  dates <- as.Date(sub(".*_(\\d{4}-\\d{2}-\\d{2})\\.Rdata$", "\\1", files))
  latest <- files[which.max(dates)]
  message("Loading most recent match for '", pattern, "': ", basename(latest))
  load(latest, envir = .GlobalEnv)
}

##############################################################################
### Fuel-specific comparison populations                                  ###
##############################################################################
# These restrictions are applied separately in every completed dataset,
# before fitting weights, censoring models, outcome models or predictions.
# A user of the selected fuel can also use other fuels (fuel stacking).
# The clean-reference trial excludes contradictory completed classifications
# rather than silently deciding which independently imputed indicator wins.
select_hap_comparison <- function(df, exposure_var,
                                  comparator = c("nonusers", "exclusive_clean")) {
  comparator <- match.arg(comparator)
  if (comparator == "nonusers" ||
      !exposure_var %in% c("wood_use", "kerosene_use")) return(df)
  fuel_cols <- c("solid_fuel_use", "wood_use", "kerosene_use", "solid_fuel_primary")
  if (!all(fuel_cols %in% names(df))) stop("Completed fuel indicators are missing")
  if (any(vapply(df[fuel_cols], function(x) anyNA(x) || !all(x %in% c(0, 1)), logical(1))))
    stop("Fuel comparisons require complete numeric 0/1 fuel indicators")
  clean <- df$solid_fuel_use == 0
  inconsistent <- clean & (df$wood_use == 1 | df$kerosene_use == 1 |
                            df$solid_fuel_primary == 1)
  keep <- !inconsistent & (df[[exposure_var]] == 1 | clean)
  df[keep, , drop = FALSE]
}

# Save counts before attempting a fit, including fits which subsequently fail.
# Survival counts include censored rows; baseline counts use observed outcomes.
prepare_hap_comparisons <- function(imputations, exposure_var, comparator,
                                    run_label, output_dir, outcome_var = NULL) {
  counts <- vector("list", length(imputations))
  selected <- lapply(seq_along(imputations), function(i) {
    source_df <- imputations[[i]]
    observed <- if (is.null(outcome_var)) source_df else
      source_df[!is.na(source_df[[outcome_var]]), , drop = FALSE]
    d <- select_hap_comparison(observed, exposure_var, comparator)
    id_n <- function(x) if ("PID" %in% names(x)) n_distinct(x$PID) else nrow(x)
    clean <- observed$solid_fuel_use == 0
    inconsistent <- clean & (observed$wood_use == 1 | observed$kerosene_use == 1 |
                              observed$solid_fuel_primary == 1)
    restricted <- comparator == "exclusive_clean" &&
      exposure_var %in% c("wood_use", "kerosene_use")
    ds <- d
    if (!"interval" %in% names(ds)) ds$interval <- 0L
    if (!"PID" %in% names(ds)) ds$PID <- seq_len(nrow(ds))
    ds$comparison_arm <- ifelse(ds[[exposure_var]] == 1, "exposed",
                                if (restricted) "exclusive_clean" else "nonusers")
    ds$.event <- if (is.null(outcome_var)) ds$event else ds[[outcome_var]]
    cells <- ds %>% group_by(interval, comparison_arm) %>%
      summarise(n_people = n_distinct(PID), n_observed = sum(!is.na(.event)),
                n_events = sum(.event == 1, na.rm = TRUE),
                n_censored = sum(is.na(.event)), .groups = "drop")
    # Explicit zeros make an absent exposure arm visible in diagnostic exports.
    template <- expand.grid(interval = sort(unique(if ("interval" %in% names(observed))
      observed$interval else 0L)), comparison_arm = c("exposed",
        if (restricted) "exclusive_clean" else "nonusers"), stringsAsFactors = FALSE)
    cells <- left_join(as_tibble(template), cells, by = c("interval", "comparison_arm"))
    for (v in c("n_people", "n_observed", "n_events", "n_censored")) cells[[v]][is.na(cells[[v]])] <- 0L
    baseline <- if ("interval" %in% names(d)) d[d$interval == 1, , drop = FALSE] else d
    counts[[i]] <<- cells %>% mutate(
      imputation = i, exposure = exposure_var, comparator = comparator,
      n_input_people = id_n(source_df),
      n_missing_outcome_excluded = id_n(source_df) - id_n(observed),
      n_inconsistent_excluded = if (restricted) id_n(observed[inconsistent, , drop = FALSE]) else 0L,
      n_other_fuel_excluded = if (restricted)
        id_n(observed[!inconsistent & observed[[exposure_var]] == 0 & !clean, , drop = FALSE]) else 0L,
      n_prediction_people = id_n(baseline), .before = 1)
    d
  })
  readr::write_csv(bind_rows(counts), file.path(output_dir, paste0(run_label, "_comparison_counts.csv")))
  selected
}

record_hap_failure <- function(run_label, error, output_dir) {
  path <- file.path(output_dir, "hap_fit_failures.csv")
  row <- tibble(run_label = run_label, error = conditionMessage(error))
  readr::write_csv(row, path, append = file.exists(path))
}

capture_hap_warnings <- function(expr) {
  warnings <- character(0)
  out <- withCallingHandlers(expr, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  attr(out, "fit_warnings") <- unique(warnings)
  out
}

empty_hap_summary <- function(survival = FALSE) {
  keys <- if (survival) tibble(cohort = character(), exposure = character(), horizon = character()) else
    tibble(population = character(), exposure = character(), outcome = character())
  bind_cols(keys, tibble(measure = character(), estimate = double(), se = double(),
                        ci_low = double(), ci_high = double(), df = double()))
}

##############################################################################
### Covariate collapsing                                                  ###
##############################################################################
# Collapses sparse covariate categories and derives BMI. Applied to every
# completed dataset after loading, so both analyses use the same definitions.
collapse_covariates <- function(df) {
  df %>%
    mutate(
      B_Gender_numeric = case_when(
        B_Gender == "Male" ~ 1,
        B_Gender == "Female" ~ 2,
        TRUE ~ NA_real_
      ),
      B_Age_months = B_Age * 12,
      Education_collapsed =
        fct_collapse(B_Education,
                     Low = c("Illiterate",
                             "PrimarySchool0to5years",
                             "NotApplicable"),
                     Medium = c("HighSchool6to10years",
                                "JuniorCollege11to12years"),
                     High = c("College12plusyears",
                              "Graduate",
                              "Professional")
      ),
      Employment_collapsed =
        fct_collapse(B_Employment,
                     Not_in_Labor = c("Housewife",
                                      "Nonworkingunemployed",
                                      "NotApplicable"),
                     Manual_Trade = c("Unskilledmanualworker",
                                      "Skilledmanualworkerfarmer",
                                      "Tradesbusinessshopowner"),
                     Professional = c("Professional",
                                      "Semiprofessionclerk",
                                      "Other")
      ),
      Income_collapsed =
        fct_collapse(B_HH_Income.Categories,
                     Low = c("1000to4999", "5000to9999"),
                     Medium = c("10000to14999", "15000to19999"),
                     High = c("20000to24999", "25000to29999", "morethan30000")
      ),
      # Adult BMI (kg/m^2) from B_Weight and B_Height (cm); B_Height is also kept
      # as a separate covariate.
      BMI = B_Weight / (B_Height / 100)^2,
      # Two-level Living.Area (Urban/Rural), with peri-urban towns grouped as
      # Urban. Levels are matched by pattern rather than exact spelling.
      Living.Area_collapsed = case_when(
        grepl("Rural", Living.Area) ~ "Rural",
        grepl("Urban|PeriUrban", Living.Area) ~ "Urban",
        TRUE ~ NA_character_
      ),
      Living.Area_collapsed = factor(Living.Area_collapsed)
    ) %>%
    # Two-level season (Monsoon vs non-Monsoon): monsoon rainfall lowers ambient
    # PM2.5, while the other three seasons are drier with higher PM2.5. Applied
    # to whichever of season / season_Entry exist in df.
    mutate(across(any_of(c("season", "season_Entry")),
                 ~factor(if_else(.x == "Monsoon", "Monsoon", "Non-Monsoon")),
                 .names = "{.col}_collapsed"))
}

# Replaces one covariate name (`placeholder`) in a list with `new_col`, so a
# base list can be reused by an exposure model that needs a differently timed
# version of the same confounder (e.g. time-varying season_collapsed in place
# of season_Entry_collapsed).
substitute_confounder <- function(cov_list, new_col, placeholder) {
  sub(paste0("^", placeholder, "$"), new_col, cov_list)
}

##############################################################################
### Covariate-list capacity diagnostics                                   ###
##############################################################################
# Reports, for a fixed covariate list, the number of model-matrix parameters
# (a k-level factor counts k-1), the observations per parameter, and whether
# that meets min_obs_per_param. Diagnostic only: the covariate list is never
# modified, so every imputation and draw is fitted with the same model.
#
# min_obs_per_param = 15 follows the usual 10-20 observations-per-parameter
# guidance for logistic regression (e.g. Peduzzi et al. 1996).
covariate_capacity <- function(cov_list, data_frame, min_obs_per_param = 15, label = NULL) {
  n <- nrow(data_frame)
  if (length(cov_list) == 0) {
    n_params <- 0L
  } else {
    form <- as.formula(paste("~", paste(cov_list, collapse = " + ")))
    n_params <- tryCatch(
      ncol(model.matrix(form, data = data_frame)) - 1L, # drop intercept
      error = function(e) NA_integer_
    )
  }
  obs_per_param <- if (is.na(n_params) || n_params == 0) Inf else n / n_params
  meets_threshold <- is.finite(obs_per_param) && obs_per_param >= min_obs_per_param

  if (!is.na(n_params)) {
    message(if (!is.null(label)) paste0("[", label, "] ") else "",
           "n = ", n, ", ", length(cov_list), " covariates -> ", n_params,
           " model-matrix params (", round(obs_per_param, 1), " obs/param)",
           if (!meets_threshold) {
             paste0(" -- BELOW min_obs_per_param = ", min_obs_per_param,
                    "; consider manually shortening this covariate list")
           } else "")
  } else {
    message(if (!is.null(label)) paste0("[", label, "] ") else "",
           "could not build a model matrix for this covariate list on this ",
           "subset (check for a factor level with zero observations here)")
  }

  tibble(label = label, n = n, n_covariates = length(cov_list),
        n_params = n_params, obs_per_param = round(obs_per_param, 2),
        min_obs_per_param = min_obs_per_param, meets_threshold = meets_threshold)
}

# Interactive tool, not used in the pipeline: how many covariates, dropped
# from the end of the list, would be needed to reach min_obs_per_param.
thin_covariates_for_n <- function(cov_list, data_frame, min_obs_per_param = 15,
                                  min_covariates = 3, label = NULL) {
  covs <- cov_list
  n <- nrow(data_frame)
  repeat {
    if (length(covs) == 0) break
    form <- as.formula(paste("~", paste(covs, collapse = " + ")))
    n_params <- ncol(model.matrix(form, data = data_frame)) - 1 # drop intercept
    if (length(covs) <= min_covariates || n / max(n_params, 1) >= min_obs_per_param) break
    dropped <- covs[length(covs)]
    message(if (!is.null(label)) paste0("[", label, "] ") else "",
           "thinning covariates for small n: n = ", n, ", ", n_params,
           " model-matrix params (", round(n / max(n_params, 1), 1),
           " obs/param) - dropping lowest-priority covariate '", dropped, "'")
    covs <- covs[-length(covs)]
  }
  covs
}

##############################################################################
### Multicollinearity among confounders                                   ###
##############################################################################
# Generalized variance inflation factors (GVIF; Fox & Monette 1992) via
# car::vif(), which handle multi-level factors. gvif_adj = GVIF^(1/(2*Df)) is
# on the sqrt(VIF) scale and is compared with `threshold` (2, i.e. VIF = 4
# for a single-df term).
#
# GVIF depends only on the covariates, so it is checked once per
# population/cohort. The lm() response is random noise, used only to give
# car::vif() a model object.
check_multicollinearity <- function(cov_list, data_frame, label = NULL, threshold = 2) {
  if (length(cov_list) < 2) {
    message(if (!is.null(label)) paste0("[", label, "] ") else "",
           "fewer than 2 covariates - nothing to check for multicollinearity.")
    return(invisible(tibble(label = label, term = character(0), gvif = numeric(0),
                            df = integer(0), gvif_adj = numeric(0), flagged = logical(0))))
  }

  fit_data <- data_frame
  fit_data$..dummy_y <- rnorm(nrow(fit_data))
  form <- as.formula(paste("..dummy_y ~", paste(cov_list, collapse = " + ")))

  fit <- tryCatch(lm(form, data = fit_data), error = function(e) e)
  if (inherits(fit, "error")) {
    message(if (!is.null(label)) paste0("[", label, "] ") else "",
           "could not fit a model for this covariate list (", conditionMessage(fit), ")")
    return(invisible(tibble(label = label, term = character(0), gvif = numeric(0),
                            df = integer(0), gvif_adj = numeric(0), flagged = logical(0))))
  }
  # Perfectly collinear terms give NA coefficients (rank-deficient fit), for
  # which VIF is undefined; these are reported and the check is skipped.
  aliased_terms <- names(which(is.na(coef(fit))))
  if (length(aliased_terms) > 0) {
    message(if (!is.null(label)) paste0("[", label, "] ") else "",
           "PERFECTLY COLLINEAR terms (model coefficients undefined, VIF not computable): ",
           paste(aliased_terms, collapse = ", "),
           " - two or more covariates are redundant given the others; drop one by hand.")
  }

  v <- car::vif(fit)
  if (is.matrix(v)) {
    result <- tibble(term = rownames(v), gvif = v[, "GVIF"], df = v[, "Df"],
                     gvif_adj = v[, "GVIF^(1/(2*Df))"])
  } else {
    result <- tibble(term = names(v), gvif = v, df = 1L, gvif_adj = sqrt(v))
  }
  result <- result %>%
    mutate(label = label, flagged = gvif_adj > threshold, .before = 1)

  flagged_terms <- result %>% filter(flagged)
  if (nrow(flagged_terms) > 0) {
    message(if (!is.null(label)) paste0("[", label, "] ") else "",
           nrow(flagged_terms), " of ", nrow(result),
           " covariate(s) above GVIF^(1/(2*Df)) = ", threshold, ": ",
           paste0(flagged_terms$term, " (", round(flagged_terms$gvif_adj, 2), ")", collapse = ", "))
  } else {
    message(if (!is.null(label)) paste0("[", label, "] ") else "",
           "no covariates above GVIF^(1/(2*Df)) = ", threshold, " - no strong multicollinearity flagged.")
  }

  invisible(result)
}

##############################################################################
### Covariate balance diagnostics                                        ###
##############################################################################
# Fits the exposure weights on one completed dataset and saves a balance
# table and love plot. Requires an output folder `fp` in the global
# environment.
#
# `method` is passed to WeightIt::weightit() (stabilize = TRUE when
# method = "glm"). The covariate list is used as given and returned as
# attr(., "cov_list_used").
#
# Also reports the covariate_capacity() diagnostic, the effective sample size
# (ESS) as a fraction of n, the maximum weight and any fitting warnings. A
# fit that errors is caught and reported, and attr(., "fit_failed") is set,
# so that a loop over populations and exposures can continue.
assess_balance <- function(data_frame, cov_list, exposure_var,
                           exposure_type = c("binary", "continuous"),
                           method = "glm", label, min_obs_per_param = 15,
                           fuel_comparator = "nonusers") {
  exposure_type <- match.arg(exposure_type)
  if (exposure_type == "binary")
    data_frame <- select_hap_comparison(data_frame, exposure_var, fuel_comparator)

  capacity <- covariate_capacity(cov_list, data_frame,
                                 min_obs_per_param = min_obs_per_param, label = label)

  bal_formula <- as.formula(paste(exposure_var, "~",
                                  paste(cov_list, collapse = " + ")))

  weightit_args <- list(formula = bal_formula, data = data_frame, method = method)
  if (method == "glm") weightit_args$stabilize <- TRUE

  fit_warnings <- character(0)
  fit_result <- tryCatch(
    withCallingHandlers(
      do.call(weightit, weightit_args),
      warning = function(wrn) {
        fit_warnings <<- c(fit_warnings, conditionMessage(wrn))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )

  if (inherits(fit_result, "error")) {
    message("[", label, "] FIT FAILED with the ", length(cov_list),
           "-covariate list as given (", capacity$obs_per_param, " obs/param): ",
           conditionMessage(fit_result),
           " -- shorten this population/cohort's covariate list by hand and re-run.")
    bal <- structure(list(), class = "assess_balance_failure")
    attr(bal, "cov_list_used") <- cov_list
    attr(bal, "capacity") <- capacity
    attr(bal, "fit_failed") <- TRUE
    attr(bal, "fit_error") <- conditionMessage(fit_result)
    return(invisible(bal))
  }
  w <- fit_result

  ess <- sum(w$weights)^2 / sum(w$weights^2)
  message("[", label, "] n = ", nrow(data_frame),
         ", covariates used = ", length(cov_list),
         ", ESS = ", round(ess, 1), " (", round(100 * ess / nrow(data_frame), 1), "% of n)",
         ", max weight = ", round(max(w$weights), 2))
  if (length(fit_warnings) > 0) {
    message("[", label, "] fitting warnings: ", paste(unique(fit_warnings), collapse = "; "))
  }
  if (ess / nrow(data_frame) < 0.5) {
    warning(label, ": effective sample size is less than 50% of n (",
           round(100 * ess / nrow(data_frame), 1), "%) after weighting - ",
           "weights may be unstable; consider manually shortening cov_list.")
  }

  thresh <- if (exposure_type == "binary") c(m = .1) else c(cor = .1)

  bal <- bal.tab(w, un = TRUE, thresholds = thresh)
  print(bal)

  love_path <- paste0(fp, "/love_plot_", label, ".jpg")
  jpeg(file = love_path, height = 7, width = 9, units = "in", res = 300)
  print(love.plot(w, thresholds = thresh, stars = "raw",
                  title = paste("Covariate Balance:", label)))
  dev.off()

  attr(bal, "cov_list_used") <- cov_list
  attr(bal, "capacity") <- capacity
  attr(bal, "ess") <- ess
  attr(bal, "fit_failed") <- FALSE
  invisible(bal)
}

##############################################################################
### Diagnostics of weights actually used in survival fits                  ###
##############################################################################
# One row per interval, weight stage and covariate. Exposure-stage rows use
# the full observed risk set; combined-stage rows use the uncensored outcome
# rows after truncation. These are distinct diagnostic populations. No model
# is refitted and the diagnostics do not change any weights or estimates.
fitted_survival_balance <- function(df, cov_list, exposure_var, weight_col, stage) {
  bind_rows(lapply(sort(unique(df$interval)), function(iv) {
    d <- df[df$interval == iv, , drop = FALSE]
    w <- d[[weight_col]]
    if (any(!is.finite(w)) || any(w < 0) || sum(w) <= 0)
      stop("Invalid fitted weights in interval ", iv, " at stage ", stage)
    ess <- function(x) if (length(x) && sum(x^2) > 0) sum(x)^2 / sum(x^2) else NA_real_
    a <- d[[exposure_var]]
    binary <- all(a %in% c(0, 1))
    meta <- tibble(interval = iv, stage = stage, n = nrow(d),
                   n_people = n_distinct(d$PID), ess = ess(w),
                   ess_fraction = ess(w) / nrow(d), max_weight = max(w),
                   max_weight_share = max(w) / sum(w),
                   ess_unexposed = if (binary) ess(w[a == 0]) else NA_real_,
                   ess_exposed = if (binary) ess(w[a == 1]) else NA_real_)
    # A constant exposure cannot support a within-interval balance comparison.
    if (length(unique(a)) < 2L)
      return(bind_cols(meta, tibble(covariate = NA_character_,
                                   balance_metric = NA_character_,
                                   imbalance = NA_real_,
                                   diagnostic_error = "Exposure constant in interval")))
    covs <- cov_list[vapply(cov_list, function(v)
      length(unique(stats::na.omit(d[[v]]))) > 1L, logical(1))]
    if (!length(covs))
      return(bind_cols(meta, tibble(covariate = NA_character_,
                                   balance_metric = NA_character_,
                                   imbalance = NA_real_, diagnostic_error = NA_character_)))
    bal <- tryCatch(cobalt::bal.tab(
      reformulate(covs, response = exposure_var), data = d, weights = w,
      method = "weighting", un = TRUE, binary = "std",
      s.d.denom = if (binary) "pooled" else "all"
    ), error = function(e) e)
    if (inherits(bal, "error"))
      return(bind_cols(meta, tibble(covariate = NA_character_,
                                   balance_metric = NA_character_, imbalance = NA_real_,
                                   diagnostic_error = conditionMessage(bal))))
    b <- bal$Balance
    metric <- if (binary) "Diff.Adj" else "Corr.Adj"
    if (!metric %in% names(b)) stop("Missing fitted-balance metric: ", metric)
    bind_cols(meta[rep(1L, nrow(b)), ],
              tibble(covariate = rownames(b), balance_metric = metric,
                     imbalance = b[[metric]], diagnostic_error = NA_character_))
  }))
}

save_survival_fit_diagnostics <- function(fits, run_label, output_dir, nested = FALSE) {
  if (!nested) {
    warnings <- bind_rows(tibble(imputation = integer(), warning = character()),
      bind_rows(lapply(seq_along(fits), function(i) {
        w <- attr(fits[[i]], "fit_warnings")
        if (length(w)) tibble(imputation = i, warning = w) else NULL
      })))
    readr::write_csv(warnings, file.path(output_dir, paste0(run_label, "_fit_warnings.csv")))
  }
  cells <- if (nested) {
    bind_rows(lapply(seq_along(fits), function(i) {
      bind_rows(lapply(seq_along(fits[[i]]), function(k) {
        attr(fits[[i]][[k]], "weight_diagnostics") %>% mutate(imputation = i, pm_draw = k)
      }))
    }))
  } else {
    bind_rows(lapply(seq_along(fits), function(i) {
      attr(fits[[i]], "weight_diagnostics") %>% mutate(imputation = i, pm_draw = NA_integer_)
    }))
  }
  if (!nrow(cells)) stop("Fitted weight diagnostics missing for ", run_label)
  readr::write_csv(cells, file.path(output_dir, paste0(run_label, "_fitted_weight_balance.csv")))
  # Save aggregate fit cells so pooling/Monte Carlo checks can be repeated
  # without refitting or rerunning MICE. No participant records are stored.
  saveRDS(fits, file.path(output_dir, paste0(run_label, "_fit_cells.rds")))
  summary <- cells %>% group_by(interval, stage, imputation, pm_draw) %>%
    summarise(n = first(n), ess = first(ess), ess_fraction = first(ess_fraction),
              max_weight = first(max_weight), max_weight_share = first(max_weight_share),
              ess_unexposed = first(ess_unexposed), ess_exposed = first(ess_exposed),
              max_abs_imbalance = if (all(is.na(imbalance))) NA_real_ else max(abs(imbalance), na.rm = TRUE),
              n_above_0_1 = if (all(is.na(imbalance))) NA_integer_ else
                sum(abs(imbalance) > 0.1, na.rm = TRUE),
              diagnostic_error = paste(unique(na.omit(diagnostic_error)), collapse = " | "),
              .groups = "drop")
  readr::write_csv(summary, file.path(output_dir, paste0(run_label, "_fitted_weight_summary.csv")))
  invisible(summary)
}

##############################################################################
### Rubin's Rules pooling helpers                                         ###
##############################################################################
# * Degrees of freedom use the Barnard & Rubin (1999) small-sample adjustment
#   when the complete-data degrees of freedom (nu_com) are supplied;
#   nu_com = Inf gives the large-sample degrees of freedom (Rubin 1987).
# * pool_nested_mi() uses the two-stage (nested) multiple-imputation variance
#   (Shen 2000; Rubin 2003; Harel & Schafer 2003):
#       T = Ubar + (1 + 1/m) B + (1 - 1/K) W
#   where B is the between-imputation variance of the imputation means and W
#   is the mean within-imputation variance across the K PM2.5 draws. The draws
#   represent uncertainty in each person's true exposure, so W is added.
#   variance_rule = "reiter2004" gives Reiter's (2004) synthetic-data rule
#   (T = Ubar + (1 + 1/m) B - W/K), for comparison only.
#
# Complete-data degrees of freedom for svyglm() are the design degrees of
# freedom (number of clusters minus 1): the number of distinct cluster_id
# values minus 1 for svydesign(ids = ~cluster_id) (see cluster_lookup()).

barnard_rubin_df <- function(df_large_sample, lambda, nu_com) {
  if (!is.finite(nu_com)) return(df_large_sample)
  nu_obs <- ((nu_com + 1) / (nu_com + 3)) * nu_com * (1 - lambda)
  if (!is.finite(df_large_sample)) return(nu_obs)
  (df_large_sample * nu_obs) / (df_large_sample + nu_obs)
}

# Keep RD and log-RR simulation diagnostics distinct when contrast tables
# are column-bound. The original five inferential columns retain their names.
prefix_pool_diagnostics <- function(x, prefix) {
  keep <- paste0(prefix, c("estimate", "se", "df", "ci_low", "ci_high"))
  extra <- setdiff(names(x), keep)
  names(x)[match(extra, names(x))] <- paste0(prefix, extra)
  x
}

write_pool_precision <- function(rd, rr, path) {
  readr::write_csv(bind_rows(
    mutate(rd, measure = "Risk difference", simulation_scale = "risk difference", .before = 1),
    mutate(rr, measure = "Risk ratio", simulation_scale = "log risk ratio", .before = 1)
  ) %>% mutate(mcse_above_10pct_se = mcse_over_se > 0.1), path)
}

# Delete-one-outer-imputation jackknife of an analysis statistic. This measures
# simulation error of the estimated SE, separately from MCSE of the estimate.
jackknife_mcse <- function(x) {
  n <- length(x)
  if (n < 3L || any(!is.finite(x))) return(NA_real_)
  sqrt((n - 1) / n * sum((x - mean(x))^2))
}

# Standard (single-level) Rubin's Rules combination across m completed-data
# analyses, each contributing a point estimate Q and its model-based variance U.
pool_rubin <- function(Q, U, nu_com = Inf) {
  m <- length(Q)
  Qbar <- mean(Q)
  Ubar <- mean(U)
  B <- var(Q)
  Tvar <- Ubar + (1 + 1 / m) * B

  # Rubin's (1987) large-sample degrees of freedom, then the Barnard-Rubin
  # small-sample adjustment (identical to the large-sample value if
  # nu_com = Inf). lambda = proportion of total variance due to missing data.
  df_ls <- if (B > 0) (m - 1) * (1 + Ubar / ((1 + 1 / m) * B))^2 else Inf
  lambda <- ((1 + 1 / m) * B) / Tvar
  df <- barnard_rubin_df(df_ls, lambda, nu_com)

  se <- sqrt(Tvar)
  tval <- qt(0.975, df)

  # Monte Carlo error of the pooled point estimate, not its sampling SE.
  # For ratios these diagnostics remain on the log scale after back-transform.
  mcse <- sqrt(B / m)
  se_mcse <- if (m > 2L) jackknife_mcse(vapply(seq_len(m), function(i) {
    sqrt(mean(U[-i]) + (1 + 1 / (m - 1)) * var(Q[-i]))
  }, numeric(1))) else NA_real_
  tibble(estimate = Qbar, se = se, df = df,
        ci_low = Qbar - tval * se, ci_high = Qbar + tval * se,
        m = m, between_variance = B, within_variance = Ubar,
        variance_fraction_imputation = lambda,
        mcse = mcse, mcse_over_se = mcse / se,
        mcse_se = se_mcse, mcse_se_over_se = se_mcse / se)
}

# Nested (two-stage) multiple-imputation combining rule for two nested
# sources of imputation uncertainty: the m MICE imputations of missing data,
# and K PM2.5 measurement-error draws nested within each of those m
# imputations. Q_mat/U_mat are m x K matrices of point estimates and their
# model-based variances. See the header comment above for the variance rule.
pool_nested_mi <- function(Q_mat, U_mat, nu_com = Inf,
                           variance_rule = c("shen2000", "reiter2004")) {
  variance_rule <- match.arg(variance_rule)
  m <- nrow(Q_mat)
  K <- ncol(Q_mat)

  Qbar_i <- rowMeans(Q_mat)   # within-imputation mean across the K PM draws
  Ubar_i <- rowMeans(U_mat)   # within-imputation mean variance across the K PM draws
  B_i <- apply(Q_mat, 1, var) # within-imputation variance *due to* PM draws

  Qbar <- mean(Qbar_i)
  Ubar <- mean(Ubar_i)
  Bm <- var(Qbar_i)   # between-imputation variance
  Bw <- mean(B_i)     # average within-imputation PM-draw variance

  if (variance_rule == "reiter2004") {
    # Reiter (2004) rule, for comparison only.
    Tvar <- max(Ubar + (1 + 1 / m) * Bm - Bw / K, 0)
    df_ls <- if (Bm > 0) (m - 1) * (1 + Ubar / ((1 + 1 / m) * Bm))^2 else Inf
    lambda <- min(((1 + 1 / m) * Bm) / Tvar, 1)
  } else {
    Tvar <- Ubar + (1 + 1 / m) * Bm + (1 - 1 / K) * Bw
    # Shen (2000) degrees of freedom for two-stage MI
    a <- ((1 + 1 / m) * Bm) / Tvar
    b <- ((1 - 1 / K) * Bw) / Tvar
    inv_df <- (if (m > 1) a^2 / (m - 1) else 0) + (if (K > 1) b^2 / (m * (K - 1)) else 0)
    df_ls <- if (inv_df > 0) 1 / inv_df else Inf
    lambda <- a + b
  }
  df <- barnard_rubin_df(df_ls, lambda, nu_com)

  se <- sqrt(Tvar)
  tval <- qt(0.975, df)

  # Independent outer-imputation means contain both missing-data and finite-K
  # PM-draw simulation error. Do not add Bw/(m*K) to Bm/m again.
  mcse <- sqrt(Bm / m)
  se_mcse <- if (m > 2L) jackknife_mcse(vapply(seq_len(m), function(i) {
    u <- mean(Ubar_i[-i]); b_outer <- var(Qbar_i[-i]); b_pm <- mean(B_i[-i])
    total <- if (variance_rule == "reiter2004") {
      max(u + (1 + 1 / (m - 1)) * b_outer - b_pm / K, 0)
    } else u + (1 + 1 / (m - 1)) * b_outer + (1 - 1 / K) * b_pm
    sqrt(total)
  }, numeric(1))) else NA_real_
  tibble(estimate = Qbar, se = se, df = df,
        ci_low = Qbar - tval * se, ci_high = Qbar + tval * se,
        m = m, K = K, between_variance = Bm, within_variance = Ubar,
        pm_draw_variance = Bw, variance_fraction_missing_and_pm = lambda,
        mcse = mcse, mcse_over_se = mcse / se,
        mcse_se = se_mcse, mcse_se_over_se = se_mcse / se,
        mcse_pm_component = sqrt(Bw / (m * K)))
}

# Ratio measures (risk ratios) are pooled on the log scale and
# back-transformed, so confidence limits are positive and asymmetric.
# Var(log RR) = Var(RR) / RR^2 (delta method). The pooled estimate is the
# back-transformed mean log RR; the returned `se` and `df` refer to log(RR).
pool_rubin_log <- function(Q, U, nu_com = Inf) {
  if (any(!is.finite(Q)) || any(Q <= 0)) {
    stop("pool_rubin_log(): every risk-ratio estimate must be positive and finite.")
  }
  pool_rubin(log(Q), U / Q^2, nu_com = nu_com) %>%
    mutate(estimate = exp(estimate), ci_low = exp(ci_low), ci_high = exp(ci_high))
}

pool_nested_mi_log <- function(Q_mat, U_mat, nu_com = Inf,
                               variance_rule = c("shen2000", "reiter2004")) {
  if (any(!is.finite(Q_mat)) || any(Q_mat <= 0)) {
    stop("pool_nested_mi_log(): every risk-ratio estimate must be positive and finite.")
  }
  pool_nested_mi(log(Q_mat), U_mat / Q_mat^2, nu_com = nu_com,
                 variance_rule = variance_rule) %>%
    mutate(estimate = exp(estimate), ci_low = exp(ci_low), ci_high = exp(ci_high))
}

# Absolute risks (and prevalences) are pooled on the logit scale and
# back-transformed, so that confidence limits stay within (0, 1). The
# delta-method variance on the logit scale is U / (p(1 - p))^2. The returned
# se and df refer to the logit scale.
pool_rubin_logit <- function(Q, U, nu_com = Inf) {
  if (any(!is.finite(Q)) || any(Q <= 0) || any(Q >= 1)) {
    stop("pool_rubin_logit(): every absolute risk must lie strictly between 0 and 1.")
  }
  pool_rubin(qlogis(Q), U / (Q * (1 - Q))^2, nu_com = nu_com) %>%
    mutate(estimate = plogis(estimate), ci_low = plogis(ci_low), ci_high = plogis(ci_high))
}

pool_nested_mi_logit <- function(Q_mat, U_mat, nu_com = Inf,
                                 variance_rule = c("shen2000", "reiter2004")) {
  if (any(!is.finite(Q_mat)) || any(Q_mat <= 0) || any(Q_mat >= 1)) {
    stop("pool_nested_mi_logit(): every absolute risk must lie strictly between 0 and 1.")
  }
  pool_nested_mi(qlogis(Q_mat), U_mat / (Q_mat * (1 - Q_mat))^2, nu_com = nu_com,
                 variance_rule = variance_rule) %>%
    mutate(estimate = plogis(estimate), ci_low = plogis(ci_low), ci_high = plogis(ci_high))
}


##############################################################################
### Cluster identifiers for cluster-robust standard errors                ###
##############################################################################
# Outcome-model variances are cluster-robust (sandwich) estimates from
# svydesign(ids = ~cluster_id). cluster_level sets the cluster:
#   "household"  - contacts of the same index case (RELATED.INDEX.PID); a
#                  contact with no recorded index case is its own cluster.
#                  Default.
#   "person"     - each contact (allows only for a contact's repeated
#                  person-periods; no clustering in the baseline models).
#   "site_month" - study site x calendar month of the household's first Entry
#                  visit. Households are nested within these clusters, so this
#                  also allows for correlation between households enrolled in
#                  the same city and month, which share ambient PM2.5.
# Study site alone gives only two clusters, too few for a cluster-robust
# variance; site is adjusted for in the models instead.
# Cluster membership is taken from the raw Entry visits in Data/data.Rdata.
cluster_lookup <- function(level) {
  stopifnot(level %in% c("person", "household", "site_month"))
  env <- new.env()
  load(here::here("Data", "data.Rdata"), envir = env)
  entry <- env$complete %>%
    dplyr::filter(stringr::str_replace_all(VISIT, " ", "") == "Entry") %>%
    dplyr::distinct(PID, .keep_all = TRUE) %>%
    dplyr::transmute(
      PID = as.character(PID),
      household_id = dplyr::if_else(is.na(RELATED.INDEX.PID) |
                                      !nzchar(trimws(as.character(RELATED.INDEX.PID))),
                                    paste0("own_", as.character(PID)),
                                    as.character(RELATED.INDEX.PID)),
      Site = as.character(Site),
      entry_date = as.Date(Entry.Date)) %>%
    dplyr::group_by(household_id) %>%
    dplyr::mutate(first_entry = if (all(is.na(entry_date))) as.Date(NA) else min(entry_date, na.rm = TRUE),
                  site_month_id = paste(dplyr::first(Site),
                                        ifelse(is.na(first_entry), paste0("nodate_", household_id),
                                               format(first_entry, "%Y-%m")), sep = "_")) %>%
    dplyr::ungroup()
  entry$cluster_id <- switch(level, person = entry$PID, household = entry$household_id,
                             site_month = entry$site_month_id)
  message("cluster_level = '", level, "': ", dplyr::n_distinct(entry$cluster_id),
          " clusters for ", nrow(entry), " contacts with an Entry visit")
  entry %>% dplyr::select(PID, cluster_id)
}

# Attach cluster_id to a completed dataset (or person-period data) by PID.
add_cluster_id <- function(df, lookup) {
  df$cluster_id <- NULL
  out <- df %>%
    dplyr::mutate(.PID_chr = as.character(PID)) %>%
    dplyr::left_join(lookup, by = c(.PID_chr = "PID")) %>%
    dplyr::select(-.PID_chr)
  if (anyNA(out$cluster_id))
    stop(sum(is.na(out$cluster_id)), " row(s) have no cluster_id; check PIDs against Data/data.Rdata")
  out
}
