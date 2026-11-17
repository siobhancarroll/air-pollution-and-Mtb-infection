library(tidyverse)
library(here)
library(mitools)
library(tableone)
library(MatchIt)
library(cobalt)
library(mice)
library(survey)
library(WeightIt)
library(splines)
library(future.apply)
library(marginaleffects)

source(here("mi_pooling_helpers.R"))
source(here("outcome_adjustment_helpers.R"))

set.seed(123)

# Household fuel use and ambient PM2.5 are not adjusted for each other: the
# ambient PM2.5 models do not include household fuel use, and the household
# fuel (HAP) models do not include ambient PM2.5.
# adjust_pm_for_fuel <- TRUE adds baseline any polluting fuel use
# (solid_fuel_use) to the PM2.5 weighting model (output folder tagged "fuel-").
if (!exists("adjust_pm_for_fuel")) adjust_pm_for_fuel <- FALSE
pm_fuel_adjustment <- if (adjust_pm_for_fuel) "solid_fuel_use" else character(0)

# Run settings (set before source()-ing):
# weight_method: WeightIt method for the exposure weights in all models
#   (household fuel and PM2.5). "cbps" (default) = covariate balancing
#   propensity score; "glm" = stabilized logistic / normal-linear GPS;
#   "ebal" = entropy balancing.
# site_subset: NULL (both sites, default) or "NIRT" / "BJGMC" for the
#   site-stratified sensitivity analysis (data restricted to that site, Site
#   removed from the covariate lists, overall population only).
# Non-default settings are added to the output folder name.
if (!exists("weight_method")) weight_method <- "cbps"
if (!exists("site_subset")) site_subset <- NULL
if (!exists("fuel_comparator")) fuel_comparator <- "nonusers"
if (!exists("hap_only")) hap_only <- FALSE
if (!exists("outcome_adjustment_mode")) outcome_adjustment_mode <- "full"
# cluster_level: cluster for the cluster-robust standard errors of the
#   outcome models: "household" (default; contacts of the same index case),
#   "person" (no clustering), or "site_month" (site x calendar month of the
#   household's first Entry visit). See cluster_lookup() in
#   mi_pooling_helpers.R. Non-default values are added to the output folder name.
if (!exists("cluster_level")) cluster_level <- "household"
stopifnot(cluster_level %in% c("person", "household", "site_month"))
if (!exists("outcome_model_comparison")) outcome_model_comparison <- FALSE
stopifnot(outcome_adjustment_mode %in% c("current", "full"))
if (outcome_model_comparison && !exists("outcome_comparison_cache"))
  stop("Set outcome_comparison_cache to a shared Data/ directory for the comparison")
stopifnot(fuel_comparator %in% c("nonusers", "exclusive_clean"),
          is.logical(hap_only), length(hap_only) == 1L, !is.na(hap_only))
# outcome_adjust_covariates: covariates also included in every outcome model
#   (household fuel and PM2.5), with risks standardized over the population
#   (see fit_binary_exposure() / fit_continuous_exposure_draw()).
#   Default: the covariates that remained imbalanced after weighting in at
#   least one model: season and Site (strongly associated with the city-month
#   PM2.5 exposure), living area and housing type. character(0) gives
#   weighting-only estimates when outcome_adjustment_mode = "current".
#   Site is dropped in site-stratified runs, and a
#   covariate with a single level in the data is dropped from that model.
#   Full adjustment is the default: weighting covariates are added to these
#   terms. Set outcome_adjustment_mode <- "current" for reduced adjustment.
# season_detail: "collapsed" (default, primary analysis) uses Monsoon vs
#   Non-Monsoon in every weighting, censoring and outcome model; "full" uses
#   the four seasons (Winter, Pre-Monsoon, Monsoon [reference], Post-Monsoon)
#   as a sensitivity analysis (adds "season4-" to the output folder name).
#   Season was collapsed because the four-level version left sparse
#   season-by-interval cells in the smaller risk sets, and some models could
#   not be estimated.
if (!exists("season_detail")) season_detail <- "collapsed"
stopifnot(season_detail %in% c("full", "collapsed"))
SEASON <- if (season_detail == "full") "season" else "season_collapsed"
if (!exists("outcome_adjust_covariates"))
  outcome_adjust_covariates <- c(SEASON, "Site", "Living.Area_collapsed", "Type.of.House")
stopifnot(weight_method %in% c("glm", "cbps", "ebal"),
          is.null(site_subset) || site_subset %in% c("NIRT", "BJGMC"))
run_tag <- paste0(if (fuel_comparator == "exclusive_clean") "cleanref-" else "",
                  if (hap_only) "haponly-" else "",
                  if (outcome_model_comparison || outcome_adjustment_mode != "full")
                    paste0("outcome-", outcome_adjustment_mode, "-") else "",
                  if (!is.null(site_subset)) paste0("site-", site_subset, "-") else "",
                  if (weight_method != "cbps") paste0(weight_method, "-") else "",
                  if (length(outcome_adjust_covariates) == 0) "weightonly-" else "",
                  if (cluster_level != "household") paste0("cluster-", cluster_level, "-") else "",
                  if (season_detail == "full") "season4-" else "")
fit_weights <- function(formula, data) {
  args <- list(formula = formula, data = data, method = weight_method)
  if (weight_method == "glm") args$stabilize <- TRUE
  comparison_weightit(args)
}
restrict_site <- function(df) {
  if (is.null(site_subset)) df else droplevels(df[as.character(df$Site) == site_subset, ])
}

### Set up folder to store plots ##############################################
main_dir <- here("tables_and_figures")
if (!exists("output_stamp")) output_stamp <- as.character(Sys.Date())
sub_dir <- paste0("ps-", if (adjust_pm_for_fuel) "fuel-" else "", run_tag, output_stamp)
fp <- file.path(main_dir, sub_dir)  # output path for tables and figures
dir.create(fp, recursive = TRUE)  # create output folder

# Load the most recent imputed datasets
load_latest("^list_imputed_dataframes_\\d{4}-\\d{2}-\\d{2}\\.Rdata$")
load_latest("^long_imputed_dataframes_\\d{4}-\\d{2}-\\d{2}\\.Rdata$")

allImputations$imputations <- lapply(allImputations$imputations, collapse_covariates)
allImputations$imputations <- lapply(allImputations$imputations, restrict_site)  # no change unless site_subset is set
# cluster_id (see cluster_level) for the cluster-robust outcome-model variances
allImputations$imputations <- lapply(allImputations$imputations, add_cluster_id,
                                     lookup = cluster_lookup(cluster_level))

### Covariate lists per analysis population ###################################
# Lists are used as given. assess_balance() writes an observations-per-
# parameter diagnostic for each population x exposure to
# covariate_capacity_summary.csv. base_covariates_children is a shorter list
# for the smaller children subgroup.
# The baseline adjustment set is the same as the longitudinal set
# (base_covariates_survival in discrete_time_survival.R), with season at Entry
# (season_collapsed, Monsoon vs Non-Monsoon) in place of season_Entry_collapsed.
# Height, HbA1c, sleeping location, alcohol use and smoking status were
# included in earlier versions of the baseline models only; they were removed
# so that all analyses use one adjustment set.
base_covariates_overall <- c("B_Age", "B_Gender", "Site", SEASON,
                             "Income_collapsed", "Education_collapsed", "Employment_collapsed",
                             "Type.of.House", "Living.Area_collapsed",
                             "BMI", "Household.size", "B_BCG")

base_covariates_children <- c("Income_collapsed", "Living.Area_collapsed",
                              "B_Age", "B_Gender", "Site", SEASON,
                              "BMI_for_age", "Household.size", "B_BCG")

# Same as base_covariates_overall without B_Gender, for the single-sex
# (female/male) analyses
base_covariates_single_sex <- setdiff(base_covariates_overall, "B_Gender")

# Notes on the covariate lists:
# - season_collapsed is the season at the Entry visit, used for both the HAP
#   and PM2.5 models.
# - BMI (from B_Weight and B_Height) replaces B_Weight and MUAC, which were
#   highly collinear. See collapse_covariates() in mi_pooling_helpers.R.
#   Children use BMI_for_age (prepare_children()).
# - Children: Education, Employment and Type.of.House are omitted, as in
#   base_covariates_survival_children.
# - The wood_use models in the single-sex subgroups have an effective sample
#   size of about 46% of n.

# Household.size was assessed as a candidate confounder and added to the
# lists above, as in the longitudinal models. HIV (<5% prevalence) is not
# included.
#
# candidate_confounders: additional covariates to assess (capacity,
# multicollinearity, balance) without changing the fitted models; currently
# empty. See prepare_population().
candidate_confounders <- character(0)

##############################################################################
### Per-completed-dataset model fitting                                  ###
##############################################################################
# Binary HAP exposure (see hap_exposures below): one (estimate, variance)
# pair per imputation for each of the risk difference (RD) and risk ratio
# (RR). outcome_var/outcome_family allow continuous (gaussian) or binary
# (quasibinomial) outcomes. Every imputation for a population is fitted with
# the same covariate list.
fit_binary_exposure <- function(data_frame, cov_list, exposure_var, outcome_var, outcome_family) {
  hap_formula <- as.formula(paste(exposure_var, "~",
                                  paste(cov_list, collapse = " + ")))
  # Exposure weights from weight_method, plus outcome-model adjustment for
  # outcome_adjust_covariates, with risks standardized over the population.
  # Without adjustment covariates this gives the weighted marginal risks.
  dr_covs <- resolve_outcome_covariates(cov_list,
    intersect(outcome_adjust_covariates, cov_list), outcome_adjustment_mode)
  dr_covs <- dr_covs[vapply(dr_covs, function(v) length(unique(stats::na.omit(data_frame[[v]]))) > 1, logical(1))]
  outcome_formula <- as.formula(paste(outcome_var, "~", paste(c(exposure_var, dr_covs), collapse = " + ")))

  w_hap <- fit_weights(hap_formula, data_frame)
  design_hap <- svydesign(ids = ~cluster_id, weights = ~w_hap$weights, data = data_frame)
  model_hap <- fit_adjustment_outcome(outcome_formula, design = design_hap,
    family = outcome_family, prediction_data = data_frame)

  if (fuel_comparator == "exclusive_clean" &&
      (anyNA(coef(model_hap)) || any(!is.finite(coef(model_hap))) ||
       any(!is.finite(vcov(model_hap))) || !isTRUE(model_hap$converged)))
    stop("Non-converged or unidentified clean-reference outcome model")

  b_hat <- coef(model_hap)
  b_hat <- b_hat[!is.na(b_hat)]
  V_hat <- vcov(model_hap)[names(b_hat), names(b_hat), drop = FALSE]
  tt <- delete.response(terms(model_hap))
  monitor <- prediction_range_monitor()
  std_at <- function(a) {
    nd <- data_frame
    nd[[exposure_var]] <- a
    X <- model.matrix(tt, model.frame(tt, nd, xlev = model_hap$xlevels, na.action = na.pass))
    X <- X[, names(b_hat), drop = FALSE]
    X <- X[stats::complete.cases(X), , drop = FALSE]
    p <- plogis(drop(X %*% b_hat))   # quasibinomial(logit)
    if (outcome_model_comparison && (nrow(X) != nrow(nd) || anyNA(X)))
      stop("Incomplete prediction covariates; refusing to drop target members")
    monitor$add(p)
    list(m = mean(p), g = colMeans(X * (p * (1 - p))))
  }
  qform <- function(g) as.numeric(t(g) %*% V_hat %*% g)
  s0 <- std_at(0)
  s1 <- std_at(1)
  if (fuel_comparator == "exclusive_clean" &&
      (!is.finite(s0$m) || !is.finite(s1$m) || s0$m <= 0 || s1$m <= 0))
    stop("Invalid clean-reference standardized risks")
  rr <- s1$m / s0$m
  out <- list(rd_estimate = s1$m - s0$m, rd_variance = qform(s1$g - s0$g),
      rr_estimate = rr, rr_variance = rr^2 * qform(s1$g / s1$m - s0$g / s0$m),
      # Absolute standardized prevalences in the reference and exposed arms
      reference_risk = s0$m, reference_risk_variance = qform(s0$g),
      exposed_risk = s1$m, exposed_risk_variance = qform(s1$g))
  attr(out, "outcome_diagnostics") <- monitor$finish(model_hap)
  if (fuel_comparator == "exclusive_clean") {
    diagnostic_data <- data_frame
    diagnostic_data$interval <- 0L
    if (!"PID" %in% names(diagnostic_data)) diagnostic_data$PID <- seq_len(nrow(diagnostic_data))
    diagnostic_data$.hap_weight <- w_hap$weights
    attr(out, "weight_diagnostics") <- fitted_survival_balance(
      diagnostic_data, cov_list, exposure_var, ".hap_weight", "exposure")
  }
  out
}

# Continuous exposure (PM2.5, mean_pm_recentavg): one log-normal draw of
# pm_sampled, giving one (estimate, variance) pair per curve grid point and
# per contrast, for one (imputation x PM draw) cell of the nested design.
fit_continuous_exposure_draw <- function(data_frame, cov_list, pred_grid,
                                         ref_val, contrast_vals,
                                         outcome_var, outcome_family,
                                         max_attempts = 5) {
  # If the weighting fit fails for a PM draw, it is retried with a new draw,
  # so that every imputation has exactly K draws (the nested pooling needs a
  # complete m x K grid).
  if (outcome_model_comparison) max_attempts <- 1L
  last_error <- NULL
  for (attempt in seq_len(max_attempts)) {
  result <- tryCatch({
  # Draws use SE_pm_recentavg (uncertainty in the estimated site-month mean,
  # propagated across the 3 months; see data_processing.R).
  mean_col <- data_frame[["mean_pm_recentavg"]]
  se_col <- data_frame[["SE_pm_recentavg"]]

  var_log <- log(1 + (se_col^2 / mean_col^2))
  mu_log <- log(mean_col) - (0.5 * var_log)
  data_frame$pm_sampled <- comparison_lognormal_draw(mu_log, sqrt(var_log))

  pm_formula <- as.formula(paste("pm_sampled ~",
                                 paste(cov_list, collapse = " + ")))

  # Generalized propensity score weights for the continuous exposure
  w_pm <- fit_weights(pm_formula, data_frame)

  # The natural-spline basis is built once and stored, rather than with ns()
  # in the model formula, because svyglm() predictions at a few new points do
  # not reliably reuse the training knots. predict() on the stored basis
  # reuses the original knots for both the curve and the contrasts.
  ns_basis <- ns(data_frame$pm_sampled, df = 3)
  data_frame$ns1 <- ns_basis[, 1]
  data_frame$ns2 <- ns_basis[, 2]
  data_frame$ns3 <- ns_basis[, 3]

  # Additional outcome adjustment: any outcome_adjust_covariates in cov_list are
  # added to the outcome model, and risks are standardized over the
  # population (g-computation: predicted risks at a given PM2.5 level,
  # averaged over everyone), with delta-method variances from the robust
  # covariance matrix.
  dr_covs <- resolve_outcome_covariates(cov_list,
    intersect(outcome_adjust_covariates, cov_list), outcome_adjustment_mode)
  dr_covs <- dr_covs[vapply(dr_covs, function(v) length(unique(stats::na.omit(data_frame[[v]]))) > 1, logical(1))]
  outcome_formula <- as.formula(paste(outcome_var, "~", paste(c("ns1 + ns2 + ns3", dr_covs), collapse = " + ")))
  design_pm <- svydesign(ids = ~cluster_id, weights = ~w_pm$weights, data = data_frame)
  model_pm <- fit_adjustment_outcome(outcome_formula, design = design_pm,
    family = outcome_family, prediction_data = data_frame)

  b_hat <- coef(model_pm)
  keep_b <- !is.na(b_hat)
  b_hat <- b_hat[keep_b]
  V_hat <- vcov(model_pm)[names(b_hat), names(b_hat), drop = FALSE]
  tt <- delete.response(terms(model_pm))
  monitor <- prediction_range_monitor()
  std_at <- function(x, domain = "contrasts") {
    b <- predict(ns_basis, rep(x, nrow(data_frame)))
    nd <- data_frame
    nd$ns1 <- b[, 1]; nd$ns2 <- b[, 2]; nd$ns3 <- b[, 3]
    X <- model.matrix(tt, model.frame(tt, nd, xlev = model_pm$xlevels, na.action = na.pass))
    X <- X[, names(b_hat), drop = FALSE]
    X <- X[stats::complete.cases(X), , drop = FALSE]
    p <- plogis(drop(X %*% b_hat))  # quasibinomial(logit) outcome models
    if (outcome_model_comparison && (nrow(X) != nrow(nd) || anyNA(X)))
      stop("Incomplete prediction covariates; refusing to drop target members")
    monitor$add(p, domain)
    list(m = mean(p), g = colMeans(X * (p * (1 - p))))
  }
  qform <- function(g) as.numeric(t(g) %*% V_hat %*% g)

  curve_list <- lapply(pred_grid, std_at, domain = "curve")
  curve_est <- sapply(curve_list, `[[`, "m")
  curve_var <- sapply(curve_list, function(s) qform(s$g))

  ref_std <- std_at(ref_val)
  contrast_rd_est <- numeric(length(contrast_vals))
  contrast_rd_var <- numeric(length(contrast_vals))
  contrast_rr_est <- numeric(length(contrast_vals))
  contrast_rr_var <- numeric(length(contrast_vals))
  contrast_std <- vector("list", length(contrast_vals))
  for (k in seq_along(contrast_vals)) {
    s1 <- std_at(contrast_vals[k])
    contrast_std[[k]] <- s1
    contrast_rd_est[k] <- s1$m - ref_std$m
    contrast_rd_var[k] <- qform(s1$g - ref_std$g)
    rr <- s1$m / ref_std$m
    contrast_rr_est[k] <- rr
    # RR variance stored on the RR scale; pool_*_log() converts to log scale
    contrast_rr_var[k] <- rr^2 * qform(s1$g / s1$m - ref_std$g / ref_std$m)
  }

  level_std <- c(list(ref_std), contrast_std)
  out <- list(curve_est = curve_est,
      curve_var = curve_var,
      # Absolute standardized prevalences at the reference and contrast levels
      level_risk_est = vapply(level_std, `[[`, numeric(1), "m"),
      level_risk_var = vapply(level_std, function(s) qform(s$g), numeric(1)),
      contrast_rd_est = contrast_rd_est,
      contrast_rd_var = contrast_rd_var,
      contrast_rr_est = contrast_rr_est,
      contrast_rr_var = contrast_rr_var)
  attr(out, "outcome_diagnostics") <- monitor$finish(model_pm)
  out
  }, error = function(e) e)
  if (!inherits(result, "error")) return(result)
  last_error <- result
  }
  if (outcome_model_comparison) stop(last_error)
  stop("fit_continuous_exposure_draw failed after ", max_attempts,
      " attempts (last error: ", conditionMessage(last_error), ")")
}

##############################################################################
### Population definitions                                               ###
##############################################################################
prepare_overall <- function(df) df

prepare_children <- function(df) {
  child_df <- df[df$B_Age < 15, ]
  tryCatch({
    zscorer::addWGSR(data = child_df, sex = "B_Gender_numeric",
                     firstPart = "B_Weight", secondPart = "B_Height",
                     thirdPart = "B_Age_months", index = "bfa",
                     output = "BMI_for_age")
  }, error = function(e) NULL)
}

prepare_adults <- function(df) df[df$B_Age >= 15, ]

prepare_female <- function(df) df %>% filter(B_Age >= 15 & B_Gender == "Female")

prepare_male <- function(df) df %>% filter(B_Age >= 15 & B_Gender == "Male")

populations <- list(
  overall  = list(label = "overall",  covariates = base_covariates_overall,    prepare = prepare_overall),
  children = list(label = "children", covariates = base_covariates_children,  prepare = prepare_children),
  adults   = list(label = "adults",   covariates = base_covariates_overall,    prepare = prepare_adults),
  female   = list(label = "female",   covariates = base_covariates_single_sex, prepare = prepare_female),
  male     = list(label = "male",     covariates = base_covariates_single_sex, prepare = prepare_male)
)
# The children subgroup is not analysed: effective sample size and covariate
# balance were inadequate under every weighting method assessed
# (compare_pm_weighting_methods.R). prepare_children() is kept for reference.
populations$children <- NULL
if (exists("analysis_populations")) {
  stopifnot(length(analysis_populations) > 0, all(analysis_populations %in% names(populations)))
  populations <- populations[analysis_populations]
}
if (!is.null(site_subset)) {
  # Site-stratified analysis: overall population only, without Site
  populations <- populations["overall"]
  populations$overall$covariates <- setdiff(populations$overall$covariates, "Site")
}

### Outcome definitions #########################################################
# Two binary outcomes at Entry: QGIT positivity (QGIT.Result, imputed directly)
# and TST positivity (>= 5 mm). Models use a logit link, but effects are
# reported as risk differences and ratios from standardized predicted
# probabilities.
outcomes <- list(
  qgit = list(var = "qgit_positive_entry", family = quasibinomial(),
             label = "qgit",
             y_label = "Predicted probability of IGRA positivity"),
  tst  = list(var = "tst_positive_entry", family = quasibinomial(),
             label = "tst",
             y_label = "Predicted probability of TST positivity (≥5mm)")
)

### Analysis-wide constants ####################################################
# Reference PM2.5 level: 40 ug/m3 (India's annual National Ambient Air Quality
# Standard). Contrasts in 10 ug/m3 steps up to 80 ug/m3 (twice the standard,
# the highest level well within the observed exposure range). reference_in_range and
# contrast_in_range flag values outside each population's 5th-95th
# percentile exposure range.
ref_val <- 40
contrast_vals <- c(50, 60, 70, 80)
pred_grid <- seq(5, 110, length.out = 50)
K <- 10  # PM2.5 measurement-error draws nested within each imputation
writeLines(c(paste("pm_reference_ug_m3:", ref_val),
             paste("pm_contrasts_ug_m3:", paste(contrast_vals, collapse = ",")),
             paste("pm_draws_per_imputation:", K),
             paste("weight_method:", weight_method),
             paste("fuel_comparator:", fuel_comparator),
             paste("hap_only:", hap_only),
             paste("outcome_adjustment_mode:", outcome_adjustment_mode),
             paste("outcome_model_comparison:", outcome_model_comparison),
             paste("site_subset:", if (is.null(site_subset)) "both" else site_subset),
             paste("cluster_level:", cluster_level),
             paste("season_detail:", season_detail),
             "standardization: mean predicted probability over analysis population"),
           file.path(fp, "run_settings.txt"))

# Binary HAP exposures (names are labels; values are column names):
# solid_fuel_use = any polluting fuel (primary or secondary); wood_use and
# kerosene_use = wood and kerosene use (primary or secondary);
# solid_fuel_primary = polluting primary cooking fuel.
hap_exposures <- c(solid_fuel = "solid_fuel_use", wood = "wood_use",
                   kerosene = "kerosene_use", solid_primary = "solid_fuel_primary")
if (exists("hap_exposure_subset")) {
  stopifnot(length(hap_exposure_subset) > 0, all(hap_exposure_subset %in% names(hap_exposures)))
  hap_exposures <- hap_exposures[hap_exposure_subset]
}
writeLines(c(readLines(file.path(fp, "run_settings.txt")),
             paste("hap_exposures:", paste(hap_exposures, collapse = ",")),
             "clean_reference_rule: exclude contradictory clean/polluting classifications; allow fuel stacking in exposed arm",
             "clean_reference_target: selected-fuel users plus exclusive-clean users with an observed assay outcome, in each imputation"),
           file.path(fp, "run_settings.txt"))

##############################################################################
### Prepare one population's completed datasets and run its diagnostics   #
##############################################################################
# Returns the completed datasets, the covariate lists for the weighting models
# (pop$covariates, unchanged) and a capacity diagnostic per exposure. Every
# imputation is fitted with the same model.
prepare_population <- function(pop) {
  m <- length(allImputations$imputations)

  pop_datasets <- lapply(allImputations$imputations, pop$prepare)
  bad <- sapply(pop_datasets, function(d) is.null(d) || nrow(d) < 30)
  if (any(bad)) {
    stop(pop$label, ": ", sum(bad), " of ", m,
        " imputed datasets failed preparation or had fewer than 30 rows - ",
        "investigate before pooling (dropping imputations ad hoc would break Rubin's Rules).")
  }

  # Balance and capacity diagnostics use one randomly chosen imputation
  # (reproducible via set.seed(123))
  check_idx <- sample(seq_len(m), 1)
  message(pop$label, ": using imputation #", check_idx, " of ", m, " for balance diagnostics/thinning")
  check_data <- pop_datasets[[check_idx]]

  # Multicollinearity among the confounders (once per population)
  multicollinearity <- check_multicollinearity(pop$covariates, check_data, label = pop$label)

  # Balance for each binary HAP exposure model
  hap_bals <- lapply(names(hap_exposures), function(exp_name) {
    assess_balance(check_data, pop$covariates,
                   hap_exposures[[exp_name]], exposure_type = "binary", fuel_comparator = fuel_comparator, method = weight_method,
                   label = paste0(pop$label, "_", exp_name))
  })
  names(hap_bals) <- names(hap_exposures)
  hap_covariates <- lapply(hap_bals, function(b) attr(b, "cov_list_used"))

  # Balance for the continuous PM2.5 exposure model
  continuous_bal <- if (hap_only) NULL else assess_balance(check_data, c(pop$covariates, pm_fuel_adjustment),
                                   "mean_pm_recentavg", method = weight_method,
                                   exposure_type = "continuous",
                                   label = paste0(pop$label, "_pm25_recentavg"))
  continuous_covariates <- attr(continuous_bal, "cov_list_used")

  capacity <- c(lapply(hap_bals, function(b) attr(b, "capacity")),
               list(attr(continuous_bal, "capacity")))

  # Candidate-confounder diagnostics: the same checks for pop$covariates plus
  # candidate_confounders. The fitted models use pop$covariates only.
  candidate_capacity <- NULL
  candidate_multicollinearity <- NULL
  if (length(candidate_confounders) > 0) {
    extended_covariates <- union(pop$covariates, candidate_confounders)

    candidate_multicollinearity <- check_multicollinearity(
      extended_covariates, check_data, label = paste0(pop$label, "_plus_candidates"))

    candidate_hap_bals <- lapply(names(hap_exposures), function(exp_name) {
      assess_balance(check_data, extended_covariates,
                     hap_exposures[[exp_name]], exposure_type = "binary", fuel_comparator = fuel_comparator, method = weight_method,
                     label = paste0(pop$label, "_", exp_name, "_plus_candidates"))
    })

    candidate_continuous_bal <- assess_balance(check_data, extended_covariates,
                                               "mean_pm_recentavg",
                                               exposure_type = "continuous",
                                               label = paste0(pop$label, "_pm25_recentavg_plus_candidates"))

    candidate_capacity <- c(lapply(candidate_hap_bals, function(b) attr(b, "capacity")),
                            list(attr(candidate_continuous_bal, "capacity")))
  }

  list(datasets = pop_datasets, hap_covariates = hap_covariates,
      continuous_covariates = continuous_covariates, capacity = capacity,
      multicollinearity = multicollinearity,
      candidate_capacity = candidate_capacity,
      candidate_multicollinearity = candidate_multicollinearity)
}

##############################################################################
### Run one binary HAP exposure model for one population x outcome         ###
##############################################################################
run_hap_for_population <- function(pop, pop_datasets, hap_covariates, exposure_var, outcome) {
  run_label <- paste0(pop$label, "_", exposure_var, "_", outcome$label)
  message("Running HAP model for population x exposure x outcome: ", run_label)

  if (fuel_comparator == "exclusive_clean") {
    pop_datasets <- prepare_hap_comparisons(pop_datasets, exposure_var, fuel_comparator,
                                            run_label, fp, outcome_var = outcome$var)
    if (any(vapply(pop_datasets, function(d) length(unique(d[[exposure_var]])) != 2L, logical(1))))
      stop("Both exposure groups must be present in every imputation; see comparison_counts.csv")
  }

  binary_fits <- future_lapply(pop_datasets, function(d) {
    tryCatch(if (fuel_comparator == "exclusive_clean")
               capture_hap_warnings(fit_binary_exposure(d, hap_covariates, exposure_var, outcome$var, outcome$family)) else
               fit_binary_exposure(d, hap_covariates, exposure_var, outcome$var, outcome$family),
            error = function(e) structure(list(), fit_error = conditionMessage(e),
              outcome_diagnostics = attr(e, "outcome_diagnostics")))
  }, future.seed = TRUE)

  save_outcome_diagnostics(binary_fits, run_label, fp)
  fit_errs <- unlist(lapply(binary_fits, attr, "fit_error"))
  if (length(fit_errs) > 0) {
    stop(run_label, ": ", length(fit_errs), " binary-exposure model fit(s) failed - ",
        "investigate before pooling. First error: ", fit_errs[1])
  }

  # Complete-data df for Barnard-Rubin: design df of svydesign(ids = ~cluster_id)
  nu_com <- n_distinct(pop_datasets[[1]]$cluster_id) - 1
  if (fuel_comparator == "exclusive_clean") {
    nu_com <- min(vapply(pop_datasets, function(d) n_distinct(d$cluster_id), integer(1))) - 1
    save_survival_fit_diagnostics(binary_fits, run_label, fp)
  }
  rd_pooled <- pool_rubin(sapply(binary_fits, function(x) x$rd_estimate),
                          sapply(binary_fits, function(x) x$rd_variance), nu_com = nu_com)
  rr_pooled <- pool_rubin_log(sapply(binary_fits, function(x) x$rr_estimate),
                          sapply(binary_fits, function(x) x$rr_variance), nu_com = nu_com)
  saveRDS(binary_fits, file.path(fp, paste0(run_label, "_fit_cells.rds")))
  write_pool_precision(rd_pooled, rr_pooled,
    file.path(fp, paste0(run_label, "_monte_carlo_precision.csv")))
  # Absolute standardized prevalences, pooled on the logit scale
  absolute_risks <- bind_rows(lapply(c("reference", "exposed"), function(arm) {
    pooled <- pool_rubin_logit(
      vapply(binary_fits, function(x) x[[paste0(arm, "_risk")]], numeric(1)),
      vapply(binary_fits, function(x) x[[paste0(arm, "_risk_variance")]], numeric(1)),
      nu_com = nu_com)
    mutate(pooled, population = pop$label, outcome = outcome$label, exposure = exposure_var,
           arm = arm, exposure_value = as.integer(arm == "exposed"), .before = 1)
  }))
  write.csv(absolute_risks, file.path(fp, paste0(run_label, "_absolute_risks.csv")), row.names = FALSE)

  message(run_label, " HAP risk difference: ", round(rd_pooled$estimate, 3),
         " (95% CI: ", round(rd_pooled$ci_low, 3), ", ",
         round(rd_pooled$ci_high, 3), "); risk ratio: ", round(rr_pooled$estimate, 3),
         " (95% CI: ", round(rr_pooled$ci_low, 3), ", ", round(rr_pooled$ci_high, 3), ")")

  list(rd = rd_pooled, rr = rr_pooled)
}

##############################################################################
### Run the continuous PM2.5 (mean_pm_recentavg) dose-response model for   ###
### one population x outcome: pooled curve, pooled contrasts, and a plot   ###
##############################################################################
run_continuous_for_population <- function(pop, pop_datasets, continuous_covariates, outcome) {
  run_label <- paste0(pop$label, "_", outcome$label, "_recentavg")
  message("Running continuous-PM model for population x outcome: ", run_label)

  m <- length(pop_datasets)

  ## Continuous exposure: nested m x K design, nested-MI combining rule #######
  continuous_fits <- future_lapply(pop_datasets, function(d) {
    lapply(seq_len(K), function(k) {
      tryCatch(fit_continuous_exposure_draw(d, continuous_covariates, pred_grid,
                                            ref_val, contrast_vals,
                                            outcome$var, outcome$family),
              error = function(e) structure(list(), fit_error = conditionMessage(e),
              outcome_diagnostics = attr(e, "outcome_diagnostics")))
    })
  }, future.seed = TRUE)

  save_outcome_diagnostics(continuous_fits, run_label, fp, nested = TRUE)
  fit_errs <- unlist(lapply(continuous_fits, function(x) lapply(x, attr, "fit_error")))
  if (length(fit_errs) > 0) {
    stop(run_label, ": ", length(fit_errs), " continuous-exposure model fit(s) failed - ",
        "investigate before pooling. First error: ", fit_errs[1])
  }

  n_grid <- length(pred_grid)
  saveRDS(continuous_fits, file.path(fp, paste0(run_label, "_fit_cells.rds")))
  n_contrast <- length(contrast_vals)
  # Complete-data df for Barnard-Rubin: design df of svydesign(ids = ~cluster_id)
  nu_com <- n_distinct(pop_datasets[[1]]$cluster_id) - 1

  curve_Q <- array(NA_real_, dim = c(m, K, n_grid))
  curve_U <- array(NA_real_, dim = c(m, K, n_grid))
  contrast_rd_Q <- array(NA_real_, dim = c(m, K, n_contrast))
  contrast_rd_U <- array(NA_real_, dim = c(m, K, n_contrast))
  contrast_rr_Q <- array(NA_real_, dim = c(m, K, n_contrast))
  contrast_rr_U <- array(NA_real_, dim = c(m, K, n_contrast))

  for (i in seq_len(m)) {
    for (k in seq_len(K)) {
      fit <- continuous_fits[[i]][[k]]
      curve_Q[i, k, ] <- fit$curve_est
      curve_U[i, k, ] <- fit$curve_var
      contrast_rd_Q[i, k, ] <- fit$contrast_rd_est
      contrast_rd_U[i, k, ] <- fit$contrast_rd_var
      contrast_rr_Q[i, k, ] <- fit$contrast_rr_est
      contrast_rr_U[i, k, ] <- fit$contrast_rr_var
    }
  }

  curve_pooled <- map_dfr(seq_len(n_grid), function(p) {
    pool_nested_mi(curve_Q[, , p], curve_U[, , p], nu_com = nu_com)
  }) %>%
    mutate(pm_levels = pred_grid, .before = 1) %>%
    rename(lower = ci_low, upper = ci_high)

  rd_pooled <- map_dfr(seq_len(n_contrast), function(p) {
    pool_nested_mi(contrast_rd_Q[, , p], contrast_rd_U[, , p], nu_com = nu_com)
  }) %>%
    rename(rd_estimate = estimate, rd_se = se, rd_df = df, rd_ci_low = ci_low, rd_ci_high = ci_high) %>%
      prefix_pool_diagnostics("rd_")

  rr_pooled <- map_dfr(seq_len(n_contrast), function(p) {
    pool_nested_mi_log(contrast_rr_Q[, , p], contrast_rr_U[, , p], nu_com = nu_com)
  }) %>%
    rename(rr_estimate = estimate, rr_se = se, rr_df = df, rr_ci_low = ci_low, rr_ci_high = ci_high) %>%
      prefix_pool_diagnostics("rr_")

  # Observed PM2.5 range for this population, pooled across imputations. Used
  # to flag contrasts outside the observed range and to trim the plotted
  # curve.
  achieved_pm <- unlist(lapply(pop_datasets, function(d) d[["mean_pm_recentavg"]]))
  safe_bounds <- quantile(achieved_pm, probs = c(0.05, 0.95), na.rm = TRUE)

  contrast_pooled <- bind_cols(
    tibble(reference = ref_val, contrast = contrast_vals),
    rd_pooled, rr_pooled
  ) %>%
    mutate(
      pm_p05_achieved = safe_bounds[1],
      pm_p95_achieved = safe_bounds[2],
      reference_in_range = reference >= pm_p05_achieved & reference <= pm_p95_achieved,
      contrast_in_range = contrast >= pm_p05_achieved & contrast <= pm_p95_achieved
    )

  n_out_of_range <- sum(!contrast_pooled$contrast_in_range)
  if (n_out_of_range > 0) {
    message(run_label, ": ", n_out_of_range, " of ", nrow(contrast_pooled),
           " contrast_vals fall outside this population's achieved 5th-95th ",
           "percentile PM range (", round(safe_bounds[1], 1), "-", round(safe_bounds[2], 1),
           " ug/m3) - those rows are extrapolating beyond the data (see ",
           "contrast_in_range column).")
  }

  print(contrast_pooled)

  # Absolute standardized prevalences at the reference and contrast levels,
  # pooled on the logit scale
  levels_abs <- c(ref_val, contrast_vals)
  level_Q <- array(NA_real_, dim = c(m, K, length(levels_abs)))
  level_U <- array(NA_real_, dim = c(m, K, length(levels_abs)))
  for (i in seq_len(m)) for (k in seq_len(K)) {
    level_Q[i, k, ] <- continuous_fits[[i]][[k]]$level_risk_est
    level_U[i, k, ] <- continuous_fits[[i]][[k]]$level_risk_var
  }
  absolute_pooled <- map_dfr(seq_along(levels_abs), function(p) {
    pool_nested_mi_logit(level_Q[, , p], level_U[, , p], nu_com = nu_com)
  }) %>%
    mutate(population = pop$label, outcome = outcome$label, pm_levels = levels_abs,
           arm = if_else(levels_abs == ref_val, "reference", "contrast"), .before = 1)
  write.csv(absolute_pooled, file.path(fp, paste0(run_label, "_absolute_risks.csv")), row.names = FALSE)

  ## Plot ######################################################################
  plot_df_trimmed <- curve_pooled %>%
    filter(pm_levels >= safe_bounds[1] & pm_levels <= safe_bounds[2])

  plot_path <- paste0(fp, "/", run_label, "_pm25_results_plot.jpg")
  jpeg(file = plot_path, height = "7", width = "13", units = "in", res = 300)
  print(
    ggplot(plot_df_trimmed, aes(x = pm_levels, y = estimate)) +
      geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey", alpha = 0.4) +
      geom_line(color = "black", size = 1) +
      geom_rug(data = tibble(pm = achieved_pm),
              aes(x = pm), inherit.aes = FALSE, alpha = 0.2, sides = "b") +
      labs(
        x = bquote(PM[2.5] ~ "3-month pre-visit average concentration (" * mu * "g/" * m^3 * ")"),
        y = outcome$y_label) +
      theme_minimal()
  )
  dev.off()

  # Save the pooled curve for re-plotting
  write_csv(curve_pooled, paste0(fp, "/", run_label, "_pm25_curve.csv"))

  list(curve = curve_pooled, contrasts = contrast_pooled)
}

##############################################################################
### Run all populations x outcomes                                       ###
##############################################################################
plan(multisession, workers = parallel::detectCores() - 1)

pop_datasets_by_pop <- lapply(populations, prepare_population)

# Covariate-capacity report: one row per population x exposure, written
# before the main run
capacity_summary <- bind_rows(lapply(names(populations), function(pop_name) {
  bind_rows(pop_datasets_by_pop[[pop_name]]$capacity) %>%
    mutate(population = pop_name, .before = 1)
}))
cat("\n", strrep("=", 78), "\n", sep = "")
cat("COVARIATE CAPACITY (obs/param for the covariate list as currently specified)\n")
cat(strrep("=", 78), "\n", sep = "")
print(capacity_summary, n = Inf, width = Inf)
under_capacity <- capacity_summary %>% filter(!meets_threshold)
if (nrow(under_capacity) > 0) {
  cat("\n", nrow(under_capacity), " combination(s) below min_obs_per_param - see rows above.\n", sep = "")
}
write_csv(capacity_summary, file.path(fp, "covariate_capacity_summary.csv"))
message("Saved: ", file.path(fp, "covariate_capacity_summary.csv"))

# Multicollinearity report: one row per population x covariate
# (check_multicollinearity(); gvif_adj > threshold is flagged)
multicollinearity_summary <- bind_rows(lapply(names(populations), function(pop_name) {
  pop_datasets_by_pop[[pop_name]]$multicollinearity %>%
    mutate(population = pop_name, .before = 1)
}))
cat("\n", strrep("=", 78), "\n", sep = "")
cat("MULTICOLLINEARITY AMONG CONFOUNDERS (GVIF^(1/(2*Df)), threshold = 2 ~ VIF = 4)\n")
cat(strrep("=", 78), "\n", sep = "")
print(multicollinearity_summary %>% mutate(across(c(gvif, gvif_adj), ~round(., 2))), n = Inf, width = Inf)
n_flagged <- sum(multicollinearity_summary$flagged)
if (n_flagged > 0) {
  cat("\n", n_flagged, " population x covariate combination(s) flagged - see rows above.\n", sep = "")
}
write_csv(multicollinearity_summary, file.path(fp, "multicollinearity_summary.csv"))
message("Saved: ", file.path(fp, "multicollinearity_summary.csv"))

# Candidate-confounder diagnostics (no-op while candidate_confounders is
# empty). Not used by the fitted models.
if (length(candidate_confounders) > 0) {
  candidate_capacity_summary <- bind_rows(lapply(names(populations), function(pop_name) {
    bind_rows(pop_datasets_by_pop[[pop_name]]$candidate_capacity) %>%
      mutate(population = pop_name, .before = 1)
  }))
  cat("\n", strrep("=", 78), "\n", sep = "")
  cat("CANDIDATE CONFOUNDERS (", paste(candidate_confounders, collapse = ", "),
     ") - CAPACITY IF ADDED (diagnostic only, not used in the models below)\n", sep = "")
  cat(strrep("=", 78), "\n", sep = "")
  print(candidate_capacity_summary, n = Inf, width = Inf)
  write_csv(candidate_capacity_summary, file.path(fp, "candidate_covariate_capacity_summary.csv"))
  message("Saved: ", file.path(fp, "candidate_covariate_capacity_summary.csv"))

  candidate_multicollinearity_summary <- bind_rows(lapply(names(populations), function(pop_name) {
    pop_datasets_by_pop[[pop_name]]$candidate_multicollinearity %>%
      mutate(population = pop_name, .before = 1)
  }))
  cat("\n", strrep("=", 78), "\n", sep = "")
  cat("CANDIDATE CONFOUNDERS - MULTICOLLINEARITY IF ADDED (diagnostic only)\n")
  cat(strrep("=", 78), "\n", sep = "")
  print(candidate_multicollinearity_summary %>% mutate(across(c(gvif, gvif_adj), ~round(., 2))),
       n = Inf, width = Inf)
  write_csv(candidate_multicollinearity_summary, file.path(fp, "candidate_multicollinearity_summary.csv"))
  message("Saved: ", file.path(fp, "candidate_multicollinearity_summary.csv"))
  message("\nBalance/ESS for the '_plus_candidates' labels (love plots and messages above) show how ",
         "well weighting balances these characteristics if added - see the ",
         "love_plot_<label>_plus_candidates.jpg files in this run's output folder.")
}

skipped_combinations <- character(0)

results <- lapply(names(populations), function(pop_name) {
  pop <- populations[[pop_name]]
  prep <- pop_datasets_by_pop[[pop_name]]
  pd <- prep$datasets

  lapply(names(outcomes), function(outcome_name) {
    outcome <- outcomes[[outcome_name]]

    # A combination whose model fails persistently is skipped with a message,
    # so that the rest of the run continues
    hap_results <- lapply(names(hap_exposures), function(exp_name) {
      tryCatch(
        run_hap_for_population(pop, pd, prep$hap_covariates[[exp_name]],
                               hap_exposures[[exp_name]], outcome),
        error = function(e) {
          combo <- paste0(pop_name, "_", hap_exposures[[exp_name]], "_", outcome_name, "_HAP")
          message("SKIPPING ", combo, ": ", conditionMessage(e))
          if (outcome_model_comparison) record_outcome_failure(combo, e, fp) else
            if (fuel_comparator == "exclusive_clean") record_hap_failure(combo, e, fp)
          skipped_combinations[[length(skipped_combinations) + 1]] <<- combo
          NULL
        }
      )
    }) %>% setNames(names(hap_exposures))

    continuous_result <- if (hap_only) NULL else tryCatch(
      run_continuous_for_population(pop, pd, prep$continuous_covariates, outcome),
      error = function(e) {
        combo <- paste0(pop_name, "_", outcome_name, "_recentavg")
        message("SKIPPING ", combo, ": ", conditionMessage(e))
        if (outcome_model_comparison) record_outcome_failure(combo, e, fp)
        skipped_combinations[[length(skipped_combinations) + 1]] <<- combo
        NULL
      }
    )

    list(hap = hap_results, continuous = continuous_result)
  }) %>% setNames(names(outcomes))
}) %>% setNames(names(populations))

plan(sequential)

if (length(skipped_combinations) > 0) {
  message("\n", length(skipped_combinations), " combination(s) skipped due to persistent model-fitting failures:")
  message(paste(" -", skipped_combinations, collapse = "\n"))
}

##############################################################################
### Results summary tables                                                ###
##############################################################################
# HAP effects: one row per population x exposure x outcome x measure (RD, RR)
hap_summary <- map_dfr(names(populations), function(pop_name) {
  map_dfr(names(outcomes), function(outcome_name) {
    map_dfr(names(hap_exposures), function(exp_name) {
      res <- results[[pop_name]][[outcome_name]]$hap[[exp_name]]
      if (is.null(res)) return(NULL) # skipped combination (see message log above)
      bind_rows(
        mutate(res$rd, population = pop_name, exposure = hap_exposures[[exp_name]],
              outcome = outcome_name, measure = "Risk difference", .before = 1),
        mutate(res$rr, population = pop_name, exposure = hap_exposures[[exp_name]],
              outcome = outcome_name, measure = "Risk ratio", .before = 1)
      )
    })
  })
}) %>%
  bind_rows(empty_hap_summary()) %>%
  dplyr::select(population, exposure, outcome, measure, estimate, se, ci_low, ci_high, df) %>%
  mutate(across(c(estimate, se, ci_low, ci_high), ~round(., 4)),
        df = round(df, 1))

# PM2.5 (mean_pm_recentavg) dose-response contrasts: one row per population x
# outcome x contrast, with RD and RR vs the reference level (ref_val)
pm_contrast_summary <- if (hap_only) tibble() else map_dfr(names(populations), function(pop_name) {
  map_dfr(names(outcomes), function(outcome_name) {
    res <- results[[pop_name]][[outcome_name]]$continuous
    if (is.null(res)) return(NULL) # skipped combination (see message log above)
    res$contrasts %>%
      mutate(population = pop_name, outcome = outcome_name, .before = 1)
  })
}) %>%
  mutate(across(where(is.numeric) & !any_of(c("reference", "contrast")), ~round(., 4)))

cat("\n", strrep("=", 78), "\n", sep = "")
cat("HAP EFFECTS (solid_fuel_use, wood_use, kerosene_use, solid_fuel_primary): risk difference and risk ratio\n")
cat(strrep("=", 78), "\n", sep = "")
print(hap_summary, n = Inf, width = Inf)

cat("\n", strrep("=", 78), "\n", sep = "")
cat("PM2.5 DOSE-RESPONSE CONTRASTS (vs. reference = ", ref_val, " ug/m3)\n", sep = "")
cat(strrep("=", 78), "\n", sep = "")
print(pm_contrast_summary, n = Inf, width = Inf)

# Save summary tables
write_csv(hap_summary, file.path(fp, "hap_effects_summary.csv"))
if (!hap_only) write_csv(pm_contrast_summary, file.path(fp, "pm25_contrasts_summary.csv"))
message("\nSaved: ", file.path(fp, "hap_effects_summary.csv"))
if (!hap_only) message("Saved: ", file.path(fp, "pm25_contrasts_summary.csv"))
