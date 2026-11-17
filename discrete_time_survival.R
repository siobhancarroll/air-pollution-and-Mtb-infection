library(tidyverse)
library(here)
library(mitools)
library(survey)
library(WeightIt)
library(splines)
library(future.apply)
library(marginaleffects)

source(here("mi_pooling_helpers.R"))
source(here("outcome_adjustment_helpers.R"))

set.seed(123)

### Conversion definitions to analyse ########################################
# "primary"     - QGIT Negative -> Positive and TST <5 mm -> >=5 mm (cohorts
#                 "qgit" and "tst"). Default.
# "sensitivity" - stricter definitions built in multiple_imputation.R (cohorts
#                 "qgit_strict", "tst_strict10", "tst_incr10").
# All other settings are the same in both runs, so the outputs are directly
# comparable (see compare_conversion_definitions.R).
# Set conversion_definition <- "sensitivity" before source()-ing to run the
# sensitivity analysis.
if (!exists("conversion_definition")) conversion_definition <- "primary"
stopifnot(conversion_definition %in% c("primary", "sensitivity"))

# Degrees of freedom for the PM2.5 natural spline (default 3). Use 1 (linear)
# or 2 if the event counts (conversion_definition_event_counts.csv) are too
# few for a 3-df spline.
if (!exists("pm_spline_df")) pm_spline_df <- 3

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
# The clean-reference trial is opt-in; ordinary runs retain their original
# non-user comparisons. See run_fuel_clean_reference.R for a targeted run.
if (!exists("fuel_comparator")) fuel_comparator <- "nonusers"
if (!exists("hap_only")) hap_only <- FALSE
if (!exists("outcome_adjustment_mode")) outcome_adjustment_mode <- "full"
# cluster_level: cluster for the cluster-robust standard errors of the
#   outcome models: "household" (default; contacts of the same index case),
#   "person", or "site_month" (site x calendar month of the household's first
#   Entry visit). See cluster_lookup() in mi_pooling_helpers.R. Non-default
#   values are added to the output folder name.
if (!exists("cluster_level")) cluster_level <- "household"
# pm_window: PM2.5 exposure window for the conversion models. "interval"
#   (primary) is the average over each follow-up interval; the alternatives in
#   pm_windows.R (pre_entry3, early2, first_interval) are timing
#   sensitivity analyses and add "window-<name>-" to the output folder name.
if (!exists("pm_window")) pm_window <- "interval"
stopifnot(pm_window %in% c("interval", "pre_entry3", "early2", "first_interval"))
# adjust_window_season: TRUE adds pm_window_monsoon (share of the exposure
# window's months in the monsoon) to the PM2.5 weighting and outcome models,
# alongside the interval-end season. Sensitivity analysis for the timing
# windows; adds "wseason-" to the output folder name.
if (!exists("adjust_window_season")) adjust_window_season <- FALSE
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
SEASON_ENTRY <- if (season_detail == "full") "season_Entry" else "season_Entry_collapsed"
stopifnot(cluster_level %in% c("person", "household", "site_month"))
if (!exists("outcome_model_comparison")) outcome_model_comparison <- FALSE
stopifnot(outcome_adjustment_mode %in% c("current", "full"))
if (outcome_model_comparison && !exists("outcome_comparison_cache"))
  stop("Set outcome_comparison_cache to a shared Data/ directory for the comparison")
stopifnot(fuel_comparator %in% c("nonusers", "exclusive_clean"),
          is.logical(hap_only), length(hap_only) == 1L, !is.na(hap_only))
stopifnot(weight_method %in% c("glm", "cbps", "ebal"),
          is.null(site_subset) || site_subset %in% c("NIRT", "BJGMC"))
run_tag <- paste0(if (fuel_comparator == "exclusive_clean") "cleanref-" else "",
                  if (hap_only) "haponly-" else "",
                  if (outcome_model_comparison || outcome_adjustment_mode != "full")
                    paste0("outcome-", outcome_adjustment_mode, "-") else "",
                  if (!is.null(site_subset)) paste0("site-", site_subset, "-") else "",
                  if (weight_method != "cbps") paste0(weight_method, "-") else "",
                  if (cluster_level != "household") paste0("cluster-", cluster_level, "-") else "",
                  if (pm_window != "interval") paste0("window-", pm_window, "-") else "",
                  if (adjust_window_season) "wseason-" else "",
                  if (season_detail == "full") "season4-" else "")
fit_weights <- function(formula, data) {
  args <- list(formula = formula, data = data, method = weight_method)
  if (weight_method == "glm") args$stabilize <- TRUE
  comparison_weightit(args)
}
# outcome_adjust_covariates: covariates also included in every hazard model
#   (household fuel and PM2.5), with individual cumulative risks standardized over the
#   fixed baseline prediction population. Default: the covariates that
#   remained imbalanced after weighting in at least one model: season of the
#   interval, Site (dropped in site-stratified runs), living area and housing
#   type. With outcome_adjustment_mode = "current", character(0) gives
#   weighting-only estimates, except that the
#   household fuel models always include season. Full adjustment is the
#   default: the terms below are retained and weighting covariates are added.
#   Set outcome_adjustment_mode <- "current" for the reduced sensitivity set.
#   A covariate with a single
#   level in the data is dropped from that model.
if (!exists("outcome_adjust_covariates"))
  outcome_adjust_covariates <- c(SEASON, if (is.null(site_subset)) "Site",
                                 "Living.Area_collapsed", "Type.of.House")
restrict_site <- function(df) {
  if (is.null(site_subset)) df else droplevels(df[as.character(df$Site) == site_subset, ])
}

### Set up folder to store plots ##############################################
main_dir <- here("tables_and_figures")
if (!exists("output_stamp")) output_stamp <- as.character(Sys.Date())
sub_dir <- paste0("survival-",
                  if (conversion_definition == "primary") "" else paste0(conversion_definition, "-"),
                  if (pm_spline_df == 3) "" else paste0("df", pm_spline_df, "-"),
                  if (adjust_pm_for_fuel) "fuel-" else "",
                  run_tag,
                  output_stamp)
fp <- file.path(main_dir, sub_dir)
dir.create(fp, recursive = TRUE)
writeLines(c(paste("conversion_definition:", conversion_definition),
             paste("pm_spline_df:", pm_spline_df),
             paste("adjust_pm_for_fuel:", adjust_pm_for_fuel),
             paste("weight_method:", weight_method),
             paste("fuel_comparator:", fuel_comparator),
             paste("hap_only:", hap_only),
             paste("outcome_adjustment_mode:", outcome_adjustment_mode),
             paste("outcome_model_comparison:", outcome_model_comparison),
             paste("site_subset:", if (is.null(site_subset)) "both" else site_subset),
             paste("cluster_level:", cluster_level),
             paste("pm_window:", pm_window),
             paste("adjust_window_season:", adjust_window_season),
             paste("season_detail:", season_detail)),
           file.path(fp, "run_settings.txt"))

### Load the person-period imputation lists ###################################
# Most recent survival_person_period_<definition>_*.Rdata in Data/ (from
# multiple_imputation.R). collapse_covariates() derives the collapsed income,
# education and employment variables row by row.
cohort_specs <- if (conversion_definition == "primary") {
  list(
    # Plot y-axis labels; "(<horizon>)" is appended in
    # run_survival_continuous(). Cohort names keep the "qgit" prefix.
    qgit = "Cumulative risk of IGRA conversion",
    tst  = "Cumulative risk of TST conversion (≥5mm)"
  )
} else {
  list(
    qgit_strict         = "Cumulative risk of IGRA conversion (<0.20 to >0.70 IU/mL)",
    tst_strict10        = "Cumulative risk of TST conversion (<5 to ≥10mm)",
    tst_incr10          = "Cumulative risk of TST conversion (≥10mm, ≥10mm increase)"
  )
}

# cluster_id (see cluster_level) is attached to every person-period row.
cluster_ids <- cluster_lookup(cluster_level)
cohorts <- setNames(lapply(names(cohort_specs), function(def) {
  load_latest(paste0("^survival_person_period_", def, "_\\d{4}-\\d{2}-\\d{2}\\.Rdata$"))
  imp <- get(paste0("allImputations_survival_", def), envir = .GlobalEnv)
  list(label = def, imputations = lapply(lapply(lapply(imp$imputations, collapse_covariates), restrict_site),
                                         add_cluster_id, lookup = cluster_ids),
       y_label = cohort_specs[[def]])
}), names(cohort_specs))
if (pm_window != "interval" || adjust_window_season) {
  source(here("pm_windows.R"))
  pm_lookup <- pm_window_lookup(pm_window)
  for (nm in names(cohorts)) cohorts[[nm]]$imputations <- lapply(cohorts[[nm]]$imputations, apply_pm_window,
                                                                lookup = pm_lookup, replace_pm = pm_window != "interval")
  write.csv(pm_lookup %>% group_by(interval) %>%
              summarise(n = n(), n_defined = sum(!is.na(mean_pm_intervalavg)),
                        mean = mean(mean_pm_intervalavg, na.rm = TRUE), sd = sd(mean_pm_intervalavg, na.rm = TRUE), .groups = "drop"),
            file.path(fp, "pm_window_summary.csv"), row.names = FALSE)
}

### Covariates ##################################################################
# The same adjustment set is used in the baseline models (base_covariates_overall
# in propensity_scores.R, with season at Entry as season_collapsed); here
# season_Entry_collapsed (season at Entry) is used for the fixed-baseline HAP models.
# The continuous PM2.5 model uses the time-varying season of each interval's
# end visit instead (via substitute_confounder()).
#
# The list is used as given for every interval, imputation and PM draw, in
# both outcome cohorts (QGIT-negative and TST-negative at Entry). It was
# chosen to fit the smallest risk set (TST cohort, interval 2), so
# average_drinks_per_day and B_HbA1c are not included. Contact hours with the
# index case (recorded at Month 4) are not included in any model.
# BMI and B_Height replace B_Weight and MUAC, as in propensity_scores.R.
base_covariates_survival <- c("B_Age", "B_Gender", "Site", SEASON_ENTRY,
                              "Income_collapsed", "Education_collapsed", "Employment_collapsed",
                              "Type.of.House", "Living.Area_collapsed",
                              "BMI", "Household.size", "B_BCG")

# Candidate confounders assessed (capacity, multicollinearity, balance)
# without being added to the fitted models; results are written to
# candidate_covariate_capacity_summary.csv,
# candidate_multicollinearity_summary.csv and "_plus_candidates" balance
# output. HIV (<5% prevalence) is not included.
candidate_confounders <- c("smoking_status", "Household.size")
if (hap_only) candidate_confounders <- character(0)

##############################################################################
### Population subgroups (diagnostics only)                              ###
##############################################################################
# Covariate capacity, multicollinearity and balance diagnostics for each
# population subgroup x outcome cohort (prepare_population_cohort()), written
# to population_*.csv and love_plot_<population>_* files. The effect
# estimates are for the overall population only.
if (!is.null(site_subset)) base_covariates_survival <- setdiff(base_covariates_survival, "Site")
base_covariates_survival_children <- c("Income_collapsed", "Living.Area_collapsed",
                                       "B_Age", "B_Gender", "Site", SEASON_ENTRY,
                                       "BMI_for_age", "Household.size", "B_BCG")

# Same as base_covariates_survival without B_Gender, for the single-sex
# (female/male) subgroups
base_covariates_survival_single_sex <- c("B_Age", "Site", SEASON_ENTRY,
                                         "Income_collapsed", "Education_collapsed",
                                         "Employment_collapsed", "Type.of.House",
                                         "Living.Area_collapsed",
                                         "BMI", "Household.size", "B_BCG")

prepare_survival_overall <- function(df) df

prepare_survival_children <- function(df) {
  child_df <- df[df$B_Age < 15, ]
  tryCatch({
    zscorer::addWGSR(data = child_df, sex = "B_Gender_numeric",
                     firstPart = "B_Weight", secondPart = "B_Height",
                     thirdPart = "B_Age_months", index = "bfa",
                     output = "BMI_for_age")
  }, error = function(e) NULL)
}

prepare_survival_adults <- function(df) df[df$B_Age >= 15, ]

prepare_survival_female <- function(df) df %>% filter(B_Age >= 15 & B_Gender == "Female")

prepare_survival_male <- function(df) df %>% filter(B_Age >= 15 & B_Gender == "Male")

# overall and adults use base_covariates_survival
populations <- list(
  overall  = list(label = "overall",  covariates = base_covariates_survival,            prepare = prepare_survival_overall),
  children = list(label = "children", covariates = base_covariates_survival_children,   prepare = prepare_survival_children),
  adults   = list(label = "adults",   covariates = base_covariates_survival,            prepare = prepare_survival_adults),
  female   = list(label = "female",   covariates = base_covariates_survival_single_sex, prepare = prepare_survival_female),
  male     = list(label = "male",     covariates = base_covariates_survival_single_sex, prepare = prepare_survival_male)
)

### Analysis-wide constants ####################################################
# Reference PM2.5 level: 40 ug/m3 (India's annual National Ambient Air Quality
# Standard). Contrasts in 10 ug/m3 steps up to 80 ug/m3 (twice the standard,
# the highest level well within the observed exposure range). reference_in_range and
# contrast_in_range flag values outside each cohort's 5th-95th percentile
# exposure range.
ref_val <- 40
contrast_vals <- c(50, 60, 70, 80)
pred_grid <- seq(5, 110, length.out = 50)
K <- 10  # PM2.5 measurement-error draws nested within each imputation
writeLines(c(readLines(file.path(fp, "run_settings.txt")),
             paste("pm_reference_ug_m3:", ref_val),
             paste("pm_contrasts_ug_m3:", paste(contrast_vals, collapse = ",")),
             paste("pm_draws_per_imputation:", K),
             "standardization: mean of individual cumulative risks over fixed baseline cohort"),
           file.path(fp, "run_settings.txt"))

# Binary HAP exposures (as in propensity_scores.R): solid_fuel_use = any
# polluting fuel (primary or secondary); wood_use and kerosene_use = wood and
# kerosene use (primary or secondary); solid_fuel_primary = polluting primary
# cooking fuel.
# Change the exposure definition and conversion definition separately. The
# stricter outcome analyses retain the main any-polluting-fuel comparison;
# alternative HAP definitions are evaluated only with primary conversion.
hap_exposures <- resolve_hap_analysis_scope(conversion_definition,
  if (exists("hap_exposure_subset")) hap_exposure_subset else NULL)
writeLines(c(readLines(file.path(fp, "run_settings.txt")),
             paste("hap_exposures:", paste(hap_exposures, collapse = ",")),
             "hap_sensitivity_scope: alternative HAP definitions use primary conversion only; stricter conversion retains any polluting fuel only",
             "clean_reference_rule: exclude contradictory clean/polluting classifications; allow fuel stacking in exposed arm",
             "clean_reference_target: baseline-eligible selected-fuel users plus exclusive-clean users in each imputation"),
           file.path(fp, "run_settings.txt"))
# Discrete-time hazard models include factor(interval) as the baseline
# hazard (Singer & Willett). Results are cumulative risks of conversion by
# each horizon, from the discrete-time g-formula: 1 - prod(1 - hazard) over
# the horizon's intervals. month4 = interval 1 (Entry-Month 4); month12 =
# intervals 1-2 (Entry-Month 12). Month 24 is not analysed because the
# interval 3 risk sets are too small (77-133 people).
#
# Every model uses only intervals <= max(unlist(horizons)).
horizons <- list(month4 = 1:1, month12 = 1:2)

##############################################################################
### Inverse probability of censoring weights (IPCW)                      ###
##############################################################################
# build_person_period() (multiple_imputation.R) keeps each person's first
# missing/indeterminate follow-up interval as a row with event = NA,
# censored = 1, for a censoring model. IPCW corrects for censoring that is
# related to exposure or covariates.
#
# add_ipcw_weights() fits a pooled discrete-time hazard model for censoring
# (censored ~ exposure_terms + cov_list + factor(interval)) on all at-risk
# rows, including censored rows. The weight is the inverse of the cumulative
# probability of remaining uncensored, stabilized by the same quantity from
# a factor(interval)-only model.
#
# exposure_terms (the HAP exposure or the PM2.5 spline basis) lets censoring
# depend on exposure; exposure_terms = NULL fits an exposure-free model.
#
# Returns person_period_df with an `ipcw` column. Combining with the exposure
# weight, dropping event = NA rows and truncating are done in
# finalize_survival_weights() below.
add_ipcw_weights <- function(person_period_df, cov_list, exposure_terms = NULL) {
  rhs_full <- paste(c(exposure_terms, cov_list, "factor(interval)"), collapse = " + ")
  cens_formula_full <- as.formula(paste("censored ~", rhs_full))
  cens_formula_marg <- as.formula("censored ~ factor(interval)")

  cens_model <- glm(cens_formula_full, data = person_period_df, family = binomial(link = "logit"))
  marg_model <- glm(cens_formula_marg, data = person_period_df, family = binomial(link = "logit"))

  person_period_df$p_censored_full <- predict(cens_model, type = "response")
  person_period_df$p_censored_marg <- predict(marg_model, type = "response")

  person_period_df %>%
    arrange(PID, interval) %>%
    group_by(PID) %>%
    mutate(
      surv_full = cumprod(1 - p_censored_full),
      surv_marg = cumprod(1 - p_censored_marg),
      ipcw = surv_marg / surv_full
    ) %>%
    ungroup() %>%
    dplyr::select(-p_censored_full, -p_censored_marg, -surv_full, -surv_marg)
}

# Multiplies the exposure weight by the IPCW weight, drops event = NA rows
# (used only in the censoring model), and truncates the combined weight at
# the 1st/99th percentiles of the remaining rows. trunc_probs = c(0, 1)
# disables truncation.
finalize_survival_weights <- function(person_period_df, exposure_weight_col,
                                      trunc_probs = c(0.01, 0.99)) {
  person_period_df %>%
    mutate(combined_weight = .data[[exposure_weight_col]] * ipcw) %>%
    filter(!is.na(event)) %>%
    mutate(combined_weight = pmin(
      pmax(combined_weight, quantile(combined_weight, trunc_probs[1])),
      quantile(combined_weight, trunc_probs[2])
    ))
}

##############################################################################
### PM2.5 monitoring coverage gap                                        ###
##############################################################################
# PM2.5 data cover December 2014 to May 2019 at both sites, so for some early
# enrollees mean_pm_intervalavg/SE_pm_intervalavg are undefined for an
# interval (no PM2.5 data for those months). This differs from a missed visit,
# which is handled by censoring and IPCW.
#
# A person is excluded from the continuous PM2.5 model from the first
# interval with undefined PM2.5 onward, since discrete-time risk sets are
# sequential. The HAP models are unaffected. Used both in model fitting and
# in the PM2.5 diagnostics, so both use the same population.
truncate_at_first_pm_gap <- function(df) {
  first_na <- df %>%
    filter(is.na(mean_pm_intervalavg) | is.na(SE_pm_intervalavg)) %>%
    group_by(PID) %>%
    summarise(first_na = min(interval), .groups = "drop")
  df %>%
    left_join(first_na, by = "PID") %>%
    filter(is.na(first_na) | interval < first_na) %>%
    dplyr::select(-first_na)
}

##############################################################################
### Per-completed-dataset model fitting                                  ###
##############################################################################
# Binary HAP exposure (see hap_exposures): fixed at baseline, so the weighting
# model is fitted once per person (not per person-period row) and the weight
# is attached to all of that person's rows. The hazard model uses
# svydesign(ids = ~cluster_id) for cluster-robust standard errors (cluster set
# by cluster_level; household by default).
#
# The hazard model includes the exposure, factor(interval), season and
# outcome_adjust_covariates. Individual cumulative risks are averaged over
# the baseline cohort, with analytic delta-method variances.
# Build prediction rows independently of observed conversion/censoring.
# Every baseline member contributes every requested interval. The wide season
# columns must have been carried through by the revised multiple_imputation.R.
make_standardization_data <- function(df, max_interval) {
  visits <- c("Month4", "Month12", "Month24")[seq_len(max_interval)]
  required <- paste0("season_", visits)
  if (!all(required %in% names(df)))
    stop("Prediction seasons are missing. Rebuild survival person-period files ",
         "from the existing completed imputations using the revised data builder.")
  baseline <- df %>% filter(interval == 1) %>% distinct(PID, .keep_all = TRUE)
  if (!nrow(baseline)) stop("Empty baseline standardization population")
  nd <- bind_rows(lapply(seq_len(max_interval), function(iv) {
    d <- baseline
    d$interval <- iv
    d$season <- d[[paste0("season_", visits[iv])]]
    d
  }))
  nd <- collapse_covariates(nd)
  nd %>% arrange(PID, interval)
}

# Predict conditional hazards, combine survival WITHIN participant, then
# average cumulative risk. The analytic gradient retains joint covariance
# across intervals and exposure scenarios (conditional on fitted weights and
# the empirical standardization population, as in the previous variance code).
make_risk_predictor <- function(model, prediction_data, horizons, set_exposure) {
  beta <- coef(model)
  if (anyNA(beta)) stop("Aliased outcome-model coefficients: resolve before prediction")
  V <- vcov(model)[names(beta), names(beta), drop = FALSE]
  tt <- delete.response(terms(model))
  ids <- unique(prediction_data$PID)
  monitor <- prediction_range_monitor()
  evaluate <- function(value, iv, domain = "curve") {
    nd <- set_exposure(prediction_data, value)
    X <- model.matrix(tt, model.frame(tt, nd, xlev = model$xlevels,
                                     na.action = na.pass))
    X <- X[, names(beta), drop = FALSE]
    if (nrow(X) != nrow(nd) || anyNA(X) || any(!is.finite(X)))
      stop("Incomplete prediction covariates; refusing to silently drop target members")
    eta <- drop(X %*% beta)
    # For cloglog, cumulative hazard over intervals is sum(exp(eta)).
    # Individual risk = 1-exp(-sum(exp(eta))).
    selected <- nd$interval %in% iv
    group <- match(nd$PID[selected], ids)
    counts <- tabulate(group, nbins = length(ids))
    if (any(counts != length(iv))) stop("Incomplete or duplicate prediction intervals")
    lambda <- exp(eta[selected])
    if (any(!is.finite(lambda))) stop("Non-finite predicted cumulative hazard")
    totals <- rowsum(cbind(lambda, X[selected, , drop = FALSE] * lambda),
                     group, reorder = TRUE)
    survival <- exp(-totals[, 1])
    monitor$add(-expm1(-totals[, 1]), domain)
    list(r = mean(-expm1(-totals[, 1])),
         g = colMeans(totals[, -1, drop = FALSE] * survival))
  }
  qform <- function(g) as.numeric(t(g) %*% V %*% g)
  contrast <- function(value, reference, iv) {
    a <- evaluate(value, iv, "contrasts"); b <- evaluate(reference, iv, "contrasts")
    if (b$r <= 0 || a$r <= 0) stop("Non-positive risk: RR undefined")
    list(reference_risk = b$r, reference_risk_variance = qform(b$g),
         exposed_risk = a$r, exposed_risk_variance = qform(a$g),
         rd_estimate = a$r - b$r, rd_variance = qform(a$g - b$g),
         rr_estimate = a$r / b$r,
         rr_variance = (a$r / b$r)^2 * qform(a$g / a$r - b$g / b$r))
  }
  list(evaluate = evaluate, contrast = contrast, qform = qform,
       diagnostics = function() monitor$finish(model))
}

fit_survival_binary_exposure <- function(person_period_df, cov_list, exposure_var, horizons) {
  max_interval <- max(unlist(horizons))
  person_period_df <- person_period_df %>% filter(interval <= max_interval)
  prediction_data <- make_standardization_data(person_period_df, max_interval)

  person_level <- person_period_df %>% distinct(PID, .keep_all = TRUE)

  hap_formula <- as.formula(paste(exposure_var, "~", paste(cov_list, collapse = " + ")))
  w_hap <- fit_weights(hap_formula, person_level)
  person_level$w_hap <- w_hap$weights

  weighted_df <- person_period_df %>%
    left_join(person_level %>% dplyr::select(PID, w_hap), by = "PID")

  # IPCW: the censoring model includes exposure_var as well as cov_list
  weighted_df <- add_ipcw_weights(weighted_df, cov_list, exposure_terms = exposure_var)
  exposure_diagnostics <- fitted_survival_balance(
    weighted_df, cov_list, exposure_var, "w_hap", "exposure")
  weighted_df <- finalize_survival_weights(weighted_df, "w_hap")
  weight_diagnostics <- bind_rows(exposure_diagnostics, fitted_survival_balance(
    weighted_df, cov_list, exposure_var, "combined_weight", "combined_truncated"))

  design_hap <- svydesign(ids = ~cluster_id, weights = ~combined_weight, data = weighted_df)
  dr_covs <- resolve_outcome_covariates(cov_list, outcome_adjust_covariates,
    outcome_adjustment_mode, mandatory = SEASON)
  dr_covs <- dr_covs[vapply(dr_covs, function(v) length(unique(stats::na.omit(weighted_df[[v]]))) > 1, logical(1))]
  hap_outcome_formula <- as.formula(paste("event ~", paste(c(exposure_var,
    dr_covs, "factor(interval)"), collapse = " + ")))
  model_hap <- fit_adjustment_outcome(hap_outcome_formula,
    design = design_hap, family = binomial(link = "cloglog"), prediction_data = prediction_data)
  write_fit_diagnostics(model_hap, weighted_df, paste0("hap_", exposure_var))
  if (fuel_comparator == "exclusive_clean" &&
      (!isTRUE(model_hap$converged) || any(!is.finite(coef(model_hap))) ||
       any(!is.finite(vcov(model_hap))) || degf(design_hap) <= 0))
    stop("Non-converged or unidentified clean-reference hazard model")

  risk_predictor <- make_risk_predictor(model_hap, prediction_data, horizons,
    function(nd, value) { nd[[exposure_var]] <- value; nd })
  out <- setNames(lapply(names(horizons), function(hname) {
    risk_predictor$contrast(1, 0, horizons[[hname]])
  }), names(horizons))
  attr(out, "outcome_diagnostics") <- risk_predictor$diagnostics()
  attr(out, "weight_diagnostics") <- weight_diagnostics
  attr(out, "complete_data_df") <- degf(design_hap)
  out
}

# Continuous exposure (PM2.5): one log-normal draw of pm_sampled per
# person-period row (each interval has its own mean and SE), giving one
# (estimate, variance) pair per curve grid point and contrast for each
# horizon, for one (imputation x PM draw) cell of the nested design. The
# weighting model is fitted separately for each interval's risk set.
# Cumulative risk by each horizon is computed as in
# fit_survival_binary_exposure().

# diagnose_outcome_fits: TRUE writes aggregate diagnostics for every fitted
# conversion hazard model (coefficients with SEs, and events by interval x
# site x season) to <output folder>/fit_diagnostics/, to investigate models
# that cannot be estimated. No participant-level data are written.
if (!exists("diagnose_outcome_fits")) diagnose_outcome_fits <- FALSE
write_fit_diagnostics <- function(model, df, kind) {
  if (!isTRUE(diagnose_outcome_fits)) return(invisible(NULL))
  dd <- file.path(fp, "fit_diagnostics"); dir.create(dd, showWarnings = FALSE)
  tag <- paste0(kind, "_n", nrow(df), "_", Sys.getpid(), "_", sample.int(1e9, 1))
  b <- stats::coef(model); V <- tryCatch(stats::vcov(model), error = function(e) NULL)
  se <- if (is.null(V)) NA_real_ else sqrt(pmax(diag(V), 0))[names(b)]
  utils::write.csv(data.frame(kind = kind, n_rows = nrow(df), n_people = dplyr::n_distinct(df$PID),
                              converged = isTRUE(model$converged), iterations = model$iter %||% NA,
                              term = names(b), estimate = unname(b), se = unname(se)),
                   file.path(dd, paste0("coef_", tag, ".csv")), row.names = FALSE)
  cells <- df %>% dplyr::mutate(.season = as.character(.data[[SEASON]])) %>%
    dplyr::group_by(interval, Site, .season) %>%
    dplyr::summarise(n = dplyr::n(), events = sum(event == 1, na.rm = TRUE),
                     weight_sum = sum(combined_weight, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(kind = kind, n_rows = nrow(df), .before = 1)
  utils::write.csv(cells, file.path(dd, paste0("cells_", tag, ".csv")), row.names = FALSE)
  invisible(NULL)
}
`%||%` <- function(a, b) if (is.null(a)) b else a

fit_survival_continuous_draw <- function(person_period_df, cov_list,
                                         pred_grid, ref_val, contrast_vals,
                                         horizons, max_attempts = 5) {
  if (outcome_model_comparison) max_attempts <- 1L
  max_interval <- max(unlist(horizons))
  # If a weighting fit fails for a PM draw, it is retried with a new draw
  # (as in fit_continuous_exposure_draw() in propensity_scores.R).
  last_error <- NULL
  original_person_period_df <- person_period_df
  for (attempt in seq_len(max_attempts)) {
  result <- tryCatch({
  person_period_df <- original_person_period_df
  person_period_df <- person_period_df[person_period_df$interval <= max_interval, ]

  # Drop each person from the first interval with undefined PM2.5 onward
  # (see truncate_at_first_pm_gap())
  person_period_df <- truncate_at_first_pm_gap(person_period_df)
  # Target: baseline-negative contacts with first-interval PM coverage.
  # Build predictions BEFORE dropping censored rows. Later PM gaps restrict
  # fitting rows, not membership of this fixed prediction population.
  prediction_data <- make_standardization_data(person_period_df, max_interval)

  # Interval-average PM2.5: every calendar month between the interval's start
  # and end visits (see data_processing.R)
  mean_col <- person_period_df[["mean_pm_intervalavg"]]
  se_col <- person_period_df[["SE_pm_intervalavg"]]

  var_log <- log(1 + (se_col^2 / mean_col^2))
  mu_log <- log(mean_col) - (0.5 * var_log)
  person_period_df$pm_sampled <- comparison_lognormal_draw(mu_log, sqrt(var_log))

  # PM2.5 weights are fitted separately for each interval, on that interval's
  # risk set, so each at-risk person contributes once per interval and each
  # interval's confounding is adjusted separately. The same covariate list is
  # used for every interval, imputation and draw. A failed fit raises an
  # error that the retry loop above handles.
  intervals_present <- sort(unique(person_period_df$interval))
  pm_weighted_list <- lapply(intervals_present, function(iv) {
    iv_df <- person_period_df[person_period_df$interval == iv, ]
    iv_formula <- as.formula(paste("pm_sampled ~", paste(cov_list, collapse = " + ")))
    w_iv <- tryCatch(
      fit_weights(iv_formula, iv_df),
      error = function(e) {
        stop("interval ", iv, ": weightit() failed with the ",
            length(cov_list), "-covariate list as given (n = ", nrow(iv_df),
            " at risk this interval): ", conditionMessage(e), call. = FALSE)
      }
    )
    iv_df$w_pm <- w_iv$weights
    iv_df
  })
  person_period_df <- bind_rows(pm_weighted_list)

  # The natural-spline basis is built once and stored, rather than with ns()
  # in the model formula, because svyglm() predictions at a few new points do
  # not reliably reuse the training knots. predict() on the stored basis
  # reuses the original knots for both the curve and the contrasts.
  # Degrees of freedom: pm_spline_df.
  ns_basis <- ns(person_period_df$pm_sampled, df = pm_spline_df)
  ns_terms <- paste0("ns", seq_len(ncol(ns_basis)))
  for (j in seq_along(ns_terms)) person_period_df[[ns_terms[j]]] <- ns_basis[, j]

  # w_pm is already attached per row, so rows can be reordered from here on

  # IPCW: the censoring model includes the PM2.5 spline basis (ns1/ns2/...).
  # PM2.5 is defined for censored rows too, since it depends only on site and
  # calendar time.
  person_period_df <- add_ipcw_weights(person_period_df, cov_list,
                                       exposure_terms = ns_terms)
  exposure_diagnostics <- fitted_survival_balance(
    person_period_df, cov_list, "pm_sampled", "w_pm", "exposure")
  person_period_df <- finalize_survival_weights(person_period_df, "w_pm")
  weight_diagnostics <- bind_rows(exposure_diagnostics, fitted_survival_balance(
    person_period_df, cov_list, "pm_sampled", "combined_weight", "combined_truncated"))

  design_pm <- svydesign(ids = ~cluster_id, weights = ~combined_weight, data = person_period_df)
  # Additional outcome adjustment; no double-robustness claim is made.
  # Prediction below averages individual cumulative risks over baseline members.
  dr_covs <- intersect(resolve_outcome_covariates(cov_list, outcome_adjust_covariates,
    outcome_adjustment_mode), names(person_period_df))
  dr_covs <- dr_covs[vapply(dr_covs, function(v) length(unique(stats::na.omit(person_period_df[[v]]))) > 1, logical(1))]
  model_pm <- fit_adjustment_outcome(as.formula(paste("event ~", paste(c(ns_terms, dr_covs, "factor(interval)"), collapse = " + "))),
    design = design_pm, family = binomial(link = "cloglog"), prediction_data = prediction_data)
  write_fit_diagnostics(model_pm, person_period_df, "pm")

  risk_predictor <- make_risk_predictor(model_pm, prediction_data, horizons,
    function(nd, value) {
      basis <- predict(ns_basis, rep(value, nrow(nd)))
      for (j in seq_along(ns_terms)) nd[[ns_terms[j]]] <- basis[, j]
      nd
    })
  out <- setNames(lapply(names(horizons), function(hname) {
    iv <- horizons[[hname]]
    curve <- lapply(pred_grid, risk_predictor$evaluate, iv = iv)
    contrasts <- lapply(contrast_vals, risk_predictor$contrast,
                        reference = ref_val, iv = iv)
    # Include the exact reference even when it is absent from the plotting grid.
    reference <- risk_predictor$evaluate(ref_val, iv)
    list(reference_risk = reference$r,
         reference_risk_variance = risk_predictor$qform(reference$g),
         # Absolute risks at the reference and each contrast level
         level_risk_est = c(reference$r, vapply(contrasts, function(z) z$exposed_risk, numeric(1))),
         level_risk_var = c(risk_predictor$qform(reference$g),
                            vapply(contrasts, function(z) z$exposed_risk_variance, numeric(1))),
         curve_est = vapply(curve, function(z) z$r, numeric(1)),
         curve_var = vapply(curve, function(z) risk_predictor$qform(z$g), numeric(1)),
         contrast_rd_est = vapply(contrasts, function(z) z$rd_estimate, numeric(1)),
         contrast_rd_var = vapply(contrasts, function(z) z$rd_variance, numeric(1)),
         contrast_rr_est = vapply(contrasts, function(z) z$rr_estimate, numeric(1)),
         contrast_rr_var = vapply(contrasts, function(z) z$rr_variance, numeric(1)))
  }), names(horizons))
  attr(out, "outcome_diagnostics") <- risk_predictor$diagnostics()
  attr(out, "weight_diagnostics") <- weight_diagnostics
  out
  }, error = function(e) e)
  if (!inherits(result, "error")) return(result)
  last_error <- result
  }
  if (outcome_model_comparison) stop(last_error)
  stop("fit_survival_continuous_draw failed after ", max_attempts,
      " attempts (last error: ", conditionMessage(last_error), ")")
}

##############################################################################
### Prepare one cohort and run its diagnostics                           ###
##############################################################################
# Returns the imputations, the covariate lists for the weighting models
# (unchanged) and capacity diagnostics: per exposure, and per interval for
# the continuous PM2.5 model, since the risk set shrinks over follow-up.
prepare_cohort <- function(cohort) {
  m <- length(cohort$imputations)
  bad <- sapply(cohort$imputations, function(d) is.null(d) || length(unique(d$PID)) < 30)
  if (any(bad)) {
    stop(cohort$label, ": ", sum(bad), " of ", m,
        " imputed person-period datasets had fewer than 30 people at risk - ",
        "investigate before pooling.")
  }

  # Diagnostics use one randomly chosen imputation (reproducible via
  # set.seed(123))
  check_idx <- sample(seq_along(cohort$imputations), 1)
  message(cohort$label, ": using imputation #", check_idx, " of ",
         length(cohort$imputations), " for balance diagnostics/thinning")
  check_data <- cohort$imputations[[check_idx]]

  # Restrict to intervals <= max(unlist(horizons)), as in model fitting
  max_interval <- max(unlist(horizons))
  check_data <- check_data %>% filter(interval <= max_interval)
  person_level_check <- check_data %>% distinct(PID, .keep_all = TRUE)

  # Censoring rates by interval (descriptive; the IPCW weights are refitted in
  # each imputation and draw)
  message("[", cohort$label, "] censoring (event = NA) rows by interval, imputation #", check_idx, ":")
  print(check_data %>%
         group_by(interval) %>%
         summarise(n_at_risk = n(), n_censored = sum(is.na(event)),
                   pct_censored = round(100 * mean(is.na(event)), 1), .groups = "drop"))

  hap_bals <- lapply(names(hap_exposures), function(exp_name) {
    assess_balance(person_level_check, base_covariates_survival, hap_exposures[[exp_name]],
                   exposure_type = "binary", fuel_comparator = fuel_comparator, method = weight_method, label = paste0(cohort$label, "_", exp_name))
  })
  names(hap_bals) <- names(hap_exposures)
  hap_covariates <- lapply(hap_bals, function(b) attr(b, "cov_list_used"))

  # Continuous PM2.5 model: time-varying season of each interval's end visit
  # in place of season_Entry_collapsed
  continuous_covs <- substitute_confounder(base_covariates_survival, SEASON, SEASON_ENTRY)
  continuous_covs <- c(continuous_covs, pm_fuel_adjustment)  # see adjust_pm_for_fuel
  if (adjust_window_season) continuous_covs <- c(continuous_covs, window_season_covariates(season_detail))

  # Same PM2.5 coverage-gap truncation as in model fitting
  check_data_continuous <- truncate_at_first_pm_gap(check_data)
  continuous_bal <- if (hap_only) NULL else assess_balance(check_data_continuous, continuous_covs, "mean_pm_intervalavg", method = weight_method,
                                   exposure_type = "continuous", label = paste0(cohort$label, "_pm25_intervalavg"))
  continuous_covariates <- attr(continuous_bal, "cov_list_used")

  # Multicollinearity among the confounders, for both versions of the list
  # (Entry season for HAP models; time-varying season for the PM2.5 model)
  multicollinearity <- bind_rows(
    check_multicollinearity(base_covariates_survival, person_level_check,
                            label = paste0(cohort$label, "_hap")),
    check_multicollinearity(continuous_covs, check_data,
                            label = paste0(cohort$label, "_pm25_intervalavg"))
  )

  # Per-interval capacity check for the continuous PM2.5 model (diagnostic
  # only)
  intervals_present <- sort(unique(check_data_continuous$interval))
  interval_capacity <- bind_rows(lapply(intervals_present, function(iv) {
    iv_df <- check_data_continuous[check_data_continuous$interval == iv, ]
    covariate_capacity(continuous_covs, iv_df,
                       label = paste0(cohort$label, "_pm25_intervalavg_interval", iv)) %>%
      mutate(interval = iv, .after = "label")
  }))

  capacity <- c(lapply(hap_bals, function(b) attr(b, "capacity")),
               list(attr(continuous_bal, "capacity")),
               list(interval_capacity))

  # Candidate-confounder diagnostics: the same checks with
  # candidate_confounders added. Not used by the fitted models.
  candidate_capacity <- NULL
  candidate_multicollinearity <- NULL
  if (length(candidate_confounders) > 0) {
    extended_hap_covs <- union(base_covariates_survival, candidate_confounders)
    extended_continuous_covs <- union(continuous_covs, candidate_confounders)

    candidate_multicollinearity <- bind_rows(
      check_multicollinearity(extended_hap_covs, person_level_check,
                              label = paste0(cohort$label, "_hap_plus_candidates")),
      check_multicollinearity(extended_continuous_covs, check_data,
                              label = paste0(cohort$label, "_pm25_intervalavg_plus_candidates"))
    )

    candidate_hap_bals <- lapply(names(hap_exposures), function(exp_name) {
      assess_balance(person_level_check, extended_hap_covs, hap_exposures[[exp_name]],
                     exposure_type = "binary", fuel_comparator = fuel_comparator, method = weight_method,
                     label = paste0(cohort$label, "_", exp_name, "_plus_candidates"))
    })

    candidate_continuous_bal <- assess_balance(check_data_continuous, extended_continuous_covs,
                                               "mean_pm_intervalavg",
                                               exposure_type = "continuous",
                                               label = paste0(cohort$label, "_pm25_intervalavg_plus_candidates"))

    candidate_capacity <- c(lapply(candidate_hap_bals, function(b) attr(b, "capacity")),
                            list(attr(candidate_continuous_bal, "capacity")))
  }

  list(imputations = cohort$imputations, hap_covariates = hap_covariates,
      continuous_covariates = continuous_covariates, capacity = capacity,
      multicollinearity = multicollinearity,
      candidate_capacity = candidate_capacity,
      candidate_multicollinearity = candidate_multicollinearity)
}

##############################################################################
### Prepare one population x cohort combination (diagnostics only)        ###
##############################################################################
# As prepare_cohort(), but restricted to one population subgroup
# (pop$prepare()) and using that population's covariate list. Combinations
# with fewer than 30 people at risk are skipped with a warning.
prepare_population_cohort <- function(pop, cohort) {
  pop_imputations <- lapply(cohort$imputations, function(d) {
    if (is.null(d)) return(NULL)
    tryCatch(pop$prepare(d), error = function(e) NULL)
  })

  m <- length(pop_imputations)
  bad <- sapply(pop_imputations, function(d) is.null(d) || length(unique(d$PID)) < 30)
  if (any(bad)) {
    warning(pop$label, " x ", cohort$label, ": ", sum(bad), " of ", m,
           " imputed person-period datasets had fewer than 30 people at risk after",
           " restricting to this population - skipping this combination's diagnostics.")
    return(NULL)
  }

  check_idx <- sample(seq_along(pop_imputations), 1)
  check_data <- pop_imputations[[check_idx]] %>% filter(interval <= max(unlist(horizons)))
  person_level_check <- check_data %>% distinct(PID, .keep_all = TRUE)

  label_prefix <- paste0(pop$label, "_", cohort$label)

  hap_bals <- lapply(names(hap_exposures), function(exp_name) {
    assess_balance(person_level_check, pop$covariates, hap_exposures[[exp_name]],
                   exposure_type = "binary", fuel_comparator = fuel_comparator, method = weight_method, label = paste0(label_prefix, "_", exp_name))
  })
  names(hap_bals) <- names(hap_exposures)

  continuous_covs <- substitute_confounder(pop$covariates, SEASON, SEASON_ENTRY)
  continuous_covs <- c(continuous_covs, pm_fuel_adjustment)  # see adjust_pm_for_fuel
  if (adjust_window_season) continuous_covs <- c(continuous_covs, window_season_covariates(season_detail))

  # Same PM2.5 coverage-gap truncation as in model fitting
  check_data_continuous <- truncate_at_first_pm_gap(check_data)
  continuous_bal <- if (hap_only) NULL else assess_balance(check_data_continuous, continuous_covs, "mean_pm_intervalavg", method = weight_method,
                                   exposure_type = "continuous",
                                   label = paste0(label_prefix, "_pm25_intervalavg"))

  multicollinearity <- bind_rows(
    check_multicollinearity(pop$covariates, person_level_check,
                            label = paste0(label_prefix, "_hap")),
    check_multicollinearity(continuous_covs, check_data,
                            label = paste0(label_prefix, "_pm25_intervalavg"))
  )

  intervals_present <- sort(unique(check_data_continuous$interval))
  interval_capacity <- bind_rows(lapply(intervals_present, function(iv) {
    iv_df <- check_data_continuous[check_data_continuous$interval == iv, ]
    covariate_capacity(continuous_covs, iv_df,
                       label = paste0(label_prefix, "_pm25_intervalavg_interval", iv)) %>%
      mutate(interval = iv, .after = "label")
  }))

  capacity <- c(lapply(hap_bals, function(b) attr(b, "capacity")),
               list(attr(continuous_bal, "capacity")),
               list(interval_capacity))

  # Covariate lists and population-restricted imputations
  hap_covariates <- lapply(hap_bals, function(b) attr(b, "cov_list_used"))
  continuous_covariates <- attr(continuous_bal, "cov_list_used")

  list(population = pop$label, cohort = cohort$label,
      capacity = capacity, multicollinearity = multicollinearity,
      imputations = pop_imputations, hap_covariates = hap_covariates,
      continuous_covariates = continuous_covariates)
}

##############################################################################
### Run one binary HAP hazard model for one cohort. Returns a list keyed   ###
### by horizon (see `horizons`), each with $rd/$rr.                        ###
##############################################################################
run_survival_hap <- function(cohort_label, imputations, hap_covariates, exposure_var, horizons) {
  run_label <- paste0(cohort_label, "_", exposure_var)
  message("Running survival HAP model for cohort x exposure: ", run_label)

  if (fuel_comparator == "exclusive_clean") {
    imputations <- lapply(imputations, function(d) d[d$interval <= max(unlist(horizons)), , drop = FALSE])
    imputations <- prepare_hap_comparisons(imputations, exposure_var, fuel_comparator, run_label, fp)
    if (any(vapply(imputations, function(d) length(unique(d[[exposure_var]])) != 2L, logical(1))))
      stop("Both exposure groups must be present in every imputation; see comparison_counts.csv")
  }

  binary_fits <- future_lapply(imputations, function(d) {
    tryCatch(if (fuel_comparator == "exclusive_clean")
               capture_hap_warnings(fit_survival_binary_exposure(d, hap_covariates, exposure_var, horizons)) else
               fit_survival_binary_exposure(d, hap_covariates, exposure_var, horizons),
            error = function(e) e)
  }, future.seed = TRUE)

  save_outcome_diagnostics(binary_fits, run_label, fp)
  failed <- sapply(binary_fits, inherits, "error")
  if (any(failed)) {
    failure_msgs <- unique(sapply(binary_fits[failed], conditionMessage))
    stop(run_label, ": ", sum(failed), " of ", length(binary_fits),
        " binary-exposure hazard model fits failed - investigate before pooling. ",
        "Distinct error(s): ", paste(failure_msgs, collapse = " | "))
  }

  # Complete-data df for Barnard-Rubin: design df of svydesign(ids = ~cluster_id)
  save_survival_fit_diagnostics(binary_fits, run_label, fp)
  nu_com <- n_distinct(imputations[[1]]$cluster_id) - 1
  if (fuel_comparator == "exclusive_clean")
    nu_com <- min(vapply(binary_fits, attr, numeric(1), "complete_data_df"))

  setNames(lapply(names(horizons), function(hname) {
    rd_pooled <- pool_rubin(sapply(binary_fits, function(x) x[[hname]]$rd_estimate),
                            sapply(binary_fits, function(x) x[[hname]]$rd_variance), nu_com = nu_com)
    rr_pooled <- pool_rubin_log(sapply(binary_fits, function(x) x[[hname]]$rr_estimate),
                            sapply(binary_fits, function(x) x[[hname]]$rr_variance), nu_com = nu_com)
    # Exact standardized absolute risks from the same fits as the contrasts,
    # pooled on the logit scale (pool_rubin_logit()). Never reconstruct these
    # from the pooled RD and RR.
    absolute_risks <- bind_rows(lapply(c("reference", "exposed"), function(arm) {
      pooled <- pool_rubin_logit(
        vapply(binary_fits, function(x) x[[hname]][[paste0(arm, "_risk")]], numeric(1)),
        vapply(binary_fits, function(x) x[[hname]][[paste0(arm, "_risk_variance")]], numeric(1)),
        nu_com = nu_com)
      mutate(pooled, cohort = cohort_label, exposure = exposure_var,
             horizon = hname, arm = arm, exposure_value = as.integer(arm == "exposed"),
             .before = 1)
    }))
    write.csv(absolute_risks,
      file.path(fp, paste0(run_label, "_", hname, "_absolute_risks.csv")), row.names = FALSE)
    write_pool_precision(rd_pooled, rr_pooled,
      file.path(fp, paste0(run_label, "_", hname, "_monte_carlo_precision.csv")))

    message(run_label, " HAP cumulative-risk difference (", hname, "): ", round(rd_pooled$estimate, 3),
           " (95% CI: ", round(rd_pooled$ci_low, 3), ", ", round(rd_pooled$ci_high, 3),
           "); risk ratio: ", round(rr_pooled$estimate, 3),
           " (95% CI: ", round(rr_pooled$ci_low, 3), ", ", round(rr_pooled$ci_high, 3), ")")

    list(rd = rd_pooled, rr = rr_pooled)
  }), names(horizons))
}

##############################################################################
### Run the continuous PM2.5 (mean_pm_intervalavg) dose-response hazard    ###
### model for one cohort: pooled curve, contrasts and plot for each        ###
### horizon. Returns a list keyed by horizon, each with $curve/$contrasts. ###
##############################################################################
run_survival_continuous <- function(cohort_label, imputations, continuous_covariates, y_label, horizons) {
  run_label <- paste0(cohort_label, "_intervalavg")
  message("Running survival continuous-PM model for cohort: ", run_label)

  m <- length(imputations)

  continuous_fits <- future_lapply(imputations, function(d) {
    lapply(seq_len(K), function(k) {
      tryCatch(fit_survival_continuous_draw(d, continuous_covariates, pred_grid,
                                            ref_val, contrast_vals, horizons),
              error = function(e) e)
    })
  }, future.seed = TRUE)

  save_outcome_diagnostics(continuous_fits, run_label, fp, nested = TRUE)
  failed_flat <- unlist(lapply(continuous_fits, function(x) sapply(x, inherits, "error")))
  if (any(failed_flat)) {
    all_errors <- unlist(lapply(continuous_fits, function(x) x[sapply(x, inherits, "error")]),
                         recursive = FALSE)
    failure_msgs <- unique(sapply(all_errors, conditionMessage))
    stop(run_label, ": ", sum(failed_flat), " of ", length(failed_flat),
        " (imputation x PM-draw) continuous-exposure hazard model fits failed ",
        "(each already retried up to max_attempts times inside ",
        "fit_survival_continuous_draw()) - investigate before pooling. ",
        "Distinct error(s): ", paste(failure_msgs, collapse = " | "))
  }

  n_grid <- length(pred_grid)
  save_survival_fit_diagnostics(continuous_fits, run_label, fp, nested = TRUE)
  n_contrast <- length(contrast_vals)
  # Complete-data df for Barnard-Rubin: design df of svydesign(ids = ~cluster_id)
  nu_com <- n_distinct(imputations[[1]]$cluster_id) - 1

  # Observed PM2.5 range for this cohort, pooled across imputations. Used to
  # flag contrasts outside the observed range and to trim the plotted curve.
  achieved_pm <- unlist(lapply(imputations, function(d) d[["mean_pm_intervalavg"]]))
  safe_bounds <- quantile(achieved_pm, probs = c(0.05, 0.95), na.rm = TRUE)

  setNames(lapply(names(horizons), function(hname) {
    reference_Q <- matrix(NA_real_, nrow = m, ncol = K)
    reference_U <- matrix(NA_real_, nrow = m, ncol = K)
    curve_Q <- array(NA_real_, dim = c(m, K, n_grid))
    curve_U <- array(NA_real_, dim = c(m, K, n_grid))
    contrast_rd_Q <- array(NA_real_, dim = c(m, K, n_contrast))
    contrast_rd_U <- array(NA_real_, dim = c(m, K, n_contrast))
    contrast_rr_Q <- array(NA_real_, dim = c(m, K, n_contrast))
    contrast_rr_U <- array(NA_real_, dim = c(m, K, n_contrast))

    for (i in seq_len(m)) {
      for (k in seq_len(K)) {
        fit <- continuous_fits[[i]][[k]][[hname]]
        reference_Q[i, k] <- fit$reference_risk
        reference_U[i, k] <- fit$reference_risk_variance
        curve_Q[i, k, ] <- fit$curve_est
        curve_U[i, k, ] <- fit$curve_var
        contrast_rd_Q[i, k, ] <- fit$contrast_rd_est
        contrast_rd_U[i, k, ] <- fit$contrast_rd_var
        contrast_rr_Q[i, k, ] <- fit$contrast_rr_est
        contrast_rr_U[i, k, ] <- fit$contrast_rr_var
      }
    }

    # Absolute cumulative risks at the reference and contrast PM2.5 levels,
    # from the same fits as the contrasts, pooled on the logit scale.
    levels_abs <- c(ref_val, contrast_vals)
    level_Q <- array(NA_real_, dim = c(m, K, length(levels_abs)))
    level_U <- array(NA_real_, dim = c(m, K, length(levels_abs)))
    for (i in seq_len(m)) for (k in seq_len(K)) {
      fit <- continuous_fits[[i]][[k]][[hname]]
      level_Q[i, k, ] <- fit$level_risk_est
      level_U[i, k, ] <- fit$level_risk_var
    }
    absolute_pooled <- map_dfr(seq_along(levels_abs), function(p) {
      pool_nested_mi_logit(level_Q[, , p], level_U[, , p], nu_com = nu_com)
    }) %>%
      mutate(cohort = cohort_label, horizon = hname, pm_levels = levels_abs,
             arm = if_else(levels_abs == ref_val, "reference", "contrast"), .before = 1)
    write.csv(absolute_pooled,
      file.path(fp, paste0(run_label, "_", hname, "_absolute_risks.csv")), row.names = FALSE)
    reference_pooled <- filter(absolute_pooled, arm == "reference")
    write.csv(reference_pooled,
      file.path(fp, paste0(run_label, "_", hname, "_absolute_reference_risk.csv")),
      row.names = FALSE)

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
      message(run_label, " (", hname, "): ", n_out_of_range, " of ", nrow(contrast_pooled),
             " contrast_vals fall outside this cohort's achieved 5th-95th ",
             "percentile PM range (", round(safe_bounds[1], 1), "-", round(safe_bounds[2], 1),
             " ug/m3) - those rows are extrapolating beyond the data (see ",
             "contrast_in_range column).")
    }

    print(contrast_pooled)

    ## Plot ####################################################################
    plot_df_trimmed <- curve_pooled %>%
      filter(pm_levels >= safe_bounds[1] & pm_levels <= safe_bounds[2])

    plot_path <- paste0(fp, "/", run_label, "_", hname, "_survival_results_plot.jpg")
    jpeg(file = plot_path, height = "7", width = "13", units = "in", res = 300)
    print(
      ggplot(plot_df_trimmed, aes(x = pm_levels, y = estimate)) +
        geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey", alpha = 0.4) +
        geom_line(color = "black", size = 1) +
        geom_rug(data = tibble(pm = achieved_pm),
                aes(x = pm), inherit.aes = FALSE, alpha = 0.2, sides = "b") +
        labs(
          x = bquote(PM[2.5] ~ "interval-average concentration (" * mu * "g/" * m^3 * ")"),
          y = paste0(y_label, " by ", c(month4 = "Month 4", month12 = "Month 12")[[hname]])) +
        theme_minimal()
    )
    dev.off()

    # Save the pooled curve for re-plotting
    write_csv(curve_pooled %>% mutate(cohort = cohort_label, horizon = hname, .before = 1),
              paste0(fp, "/", run_label, "_", hname, "_survival_curve.csv"))

    list(curve = curve_pooled, contrasts = contrast_pooled)
  }), names(horizons))
}

##############################################################################
### Run both cohorts                                                      ###
##############################################################################
plan(multisession, workers = parallel::detectCores() - 1)

cohort_imputations <- lapply(cohorts, prepare_cohort)

# Covariate-capacity report: one row per cohort x exposure, plus one row per
# cohort x interval for the continuous PM2.5 model. Written before the main
# run.
capacity_summary <- bind_rows(lapply(names(cohorts), function(cohort_name) {
  bind_rows(cohort_imputations[[cohort_name]]$capacity) %>%
    mutate(cohort = cohort_name, .before = 1)
}))
cat("\n", strrep("=", 78), "\n", sep = "")
cat("COVARIATE CAPACITY (obs/param for the covariate list as currently specified)\n")
cat(strrep("=", 78), "\n", sep = "")
print(capacity_summary, n = Inf, width = Inf)
under_capacity <- capacity_summary %>% filter(!meets_threshold)
if (nrow(under_capacity) > 0) {
  cat("\n", nrow(under_capacity), " combination(s) below min_obs_per_param - see rows above",
     " (interval-level rows flag which lag/interval within a cohort is underpowered).\n", sep = "")
}
write_csv(capacity_summary, file.path(fp, "covariate_capacity_summary.csv"))
message("Saved: ", file.path(fp, "covariate_capacity_summary.csv"))

# Multicollinearity report: one row per cohort x list version (hap: Entry
# season; pm25_intervalavg: time-varying season) x covariate
# (check_multicollinearity(); gvif_adj > threshold is flagged)
multicollinearity_summary <- bind_rows(lapply(names(cohorts), function(cohort_name) {
  cohort_imputations[[cohort_name]]$multicollinearity %>%
    mutate(cohort = cohort_name, .before = 1)
}))
cat("\n", strrep("=", 78), "\n", sep = "")
cat("MULTICOLLINEARITY AMONG CONFOUNDERS (GVIF^(1/(2*Df)), threshold = 2 ~ VIF = 4)\n")
cat(strrep("=", 78), "\n", sep = "")
print(multicollinearity_summary %>% mutate(across(c(gvif, gvif_adj), ~round(., 2))), n = Inf, width = Inf)
n_flagged <- sum(multicollinearity_summary$flagged)
if (n_flagged > 0) {
  cat("\n", n_flagged, " cohort x covariate combination(s) flagged - see rows above.\n", sep = "")
}
write_csv(multicollinearity_summary, file.path(fp, "multicollinearity_summary.csv"))
message("Saved: ", file.path(fp, "multicollinearity_summary.csv"))

# Candidate-confounder diagnostics (no-op while candidate_confounders is
# empty). Not used by the fitted models.
if (length(candidate_confounders) > 0) {
  candidate_capacity_summary <- bind_rows(lapply(names(cohorts), function(cohort_name) {
    bind_rows(cohort_imputations[[cohort_name]]$candidate_capacity) %>%
      mutate(cohort = cohort_name, .before = 1)
  }))
  cat("\n", strrep("=", 78), "\n", sep = "")
  cat("CANDIDATE CONFOUNDERS (", paste(candidate_confounders, collapse = ", "),
     ") - CAPACITY IF ADDED (diagnostic only, not used in the models below)\n", sep = "")
  cat(strrep("=", 78), "\n", sep = "")
  print(candidate_capacity_summary, n = Inf, width = Inf)
  write_csv(candidate_capacity_summary, file.path(fp, "candidate_covariate_capacity_summary.csv"))
  message("Saved: ", file.path(fp, "candidate_covariate_capacity_summary.csv"))

  candidate_multicollinearity_summary <- bind_rows(lapply(names(cohorts), function(cohort_name) {
    cohort_imputations[[cohort_name]]$candidate_multicollinearity %>%
      mutate(cohort = cohort_name, .before = 1)
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

##############################################################################
### Population x cohort diagnostics                                       ###
##############################################################################
# Runs prepare_population_cohort() for each population x cohort and writes
# population_covariate_capacity_summary.csv,
# population_multicollinearity_summary.csv and
# love_plot_<population>_<cohort>_*.jpg. Combinations with fewer than 30
# people at risk are skipped with a warning.
population_cohort_diagnostics <- list()
# Run for the primary definitions, both sites, only. The children subgroup
# is not analysed (see propensity_scores.R).
populations$children <- NULL
for (pop_name in if (!hap_only && conversion_definition == "primary" && is.null(site_subset)) names(populations) else character(0)) {
  for (cohort_name in names(cohorts)) {
    key <- paste0(pop_name, "_", cohort_name)
    population_cohort_diagnostics[[key]] <- prepare_population_cohort(
      populations[[pop_name]], cohorts[[cohort_name]])
  }
}
population_cohort_diagnostics <- Filter(Negate(is.null), population_cohort_diagnostics)

if (length(population_cohort_diagnostics) > 0) {
  population_capacity_summary <- bind_rows(lapply(population_cohort_diagnostics, function(x) {
    bind_rows(x$capacity) %>%
      mutate(population = x$population, cohort = x$cohort, .before = 1)
  }))
  cat("\n", strrep("=", 78), "\n", sep = "")
  cat("POPULATION x COHORT COVARIATE CAPACITY (obs/param, diagnostics-only)\n")
  cat(strrep("=", 78), "\n", sep = "")
  print(population_capacity_summary, n = Inf, width = Inf)
  pop_under_capacity <- population_capacity_summary %>% filter(!meets_threshold)
  if (nrow(pop_under_capacity) > 0) {
    cat("\n", nrow(pop_under_capacity), " population x cohort combination(s) below ",
       "min_obs_per_param - see rows above.\n", sep = "")
  }
  write_csv(population_capacity_summary, file.path(fp, "population_covariate_capacity_summary.csv"))
  message("Saved: ", file.path(fp, "population_covariate_capacity_summary.csv"))

  population_multicollinearity_summary <- bind_rows(lapply(population_cohort_diagnostics, function(x) {
    x$multicollinearity %>%
      mutate(population = x$population, cohort = x$cohort, .before = 1)
  }))
  cat("\n", strrep("=", 78), "\n", sep = "")
  cat("POPULATION x COHORT MULTICOLLINEARITY (GVIF^(1/(2*Df)), threshold = 2 ~ VIF = 4)\n")
  cat(strrep("=", 78), "\n", sep = "")
  print(population_multicollinearity_summary %>% mutate(across(c(gvif, gvif_adj), ~round(., 2))),
       n = Inf, width = Inf)
  pop_n_flagged <- sum(population_multicollinearity_summary$flagged)
  if (pop_n_flagged > 0) {
    cat("\n", pop_n_flagged, " population x cohort x covariate combination(s) flagged - see rows above.\n", sep = "")
  }
  write_csv(population_multicollinearity_summary, file.path(fp, "population_multicollinearity_summary.csv"))
  message("Saved: ", file.path(fp, "population_multicollinearity_summary.csv"))
} else if (conversion_definition == "primary" && is.null(site_subset)) {
  warning("Every population x cohort combination was skipped (see warnings above) - ",
         "no population_covariate_capacity_summary.csv/population_multicollinearity_summary.csv written.")
}

results <- lapply(names(cohorts), function(cohort_name) {
  cohort <- cohorts[[cohort_name]]
  prep <- cohort_imputations[[cohort_name]]
  imputations <- prep$imputations

  hap_results <- lapply(names(hap_exposures), function(exp_name) {
    tryCatch(run_survival_hap(cohort_name, imputations, prep$hap_covariates[[exp_name]],
                              hap_exposures[[exp_name]], horizons), error = function(e) {
      if (fuel_comparator != "exclusive_clean" && !outcome_model_comparison) stop(e)
      label <- paste0(cohort_name, "_", hap_exposures[[exp_name]])
      if (outcome_model_comparison) record_outcome_failure(label, e, fp) else record_hap_failure(label, e, fp)
      message("SKIPPING ", label, ": ", conditionMessage(e))
      NULL
    })
  }) %>% setNames(names(hap_exposures))

  continuous_result <- if (hap_only) NULL else tryCatch(
    run_survival_continuous(cohort_name, imputations, prep$continuous_covariates, cohort$y_label, horizons),
    error = function(e) {
      if (!outcome_model_comparison) stop(e)
      record_outcome_failure(paste0(cohort_name, "_intervalavg"), e, fp)
      message("SKIPPING ", cohort_name, " PM2.5: ", conditionMessage(e))
      NULL
    })

  list(hap = hap_results, continuous = continuous_result)
}) %>% setNames(names(cohorts))

plan(sequential)

##############################################################################
### Results summary tables                                                ###
##############################################################################
# HAP cumulative risk by horizon: one row per cohort x exposure x horizon x
# measure (RD, RR)
hap_summary <- map_dfr(names(cohorts), function(cohort_name) {
  map_dfr(names(hap_exposures), function(exp_name) {
    map_dfr(names(horizons), function(hname) {
      res <- results[[cohort_name]]$hap[[exp_name]][[hname]]
      if (is.null(res)) return(NULL)
      bind_rows(
        mutate(res$rd, cohort = cohort_name, exposure = hap_exposures[[exp_name]],
              horizon = hname, measure = "Risk difference", .before = 1),
        mutate(res$rr, cohort = cohort_name, exposure = hap_exposures[[exp_name]],
              horizon = hname, measure = "Risk ratio", .before = 1)
      )
    })
  })
}) %>%
  bind_rows(empty_hap_summary(survival = TRUE)) %>%
  dplyr::select(cohort, exposure, horizon, measure, estimate, se, ci_low, ci_high, df) %>%
  mutate(across(c(estimate, se, ci_low, ci_high), ~round(., 4)),
        df = round(df, 1))

# PM2.5 (mean_pm_intervalavg) cumulative-risk contrasts by horizon: one row
# per cohort x horizon x contrast, with RD and RR vs the reference level
# (ref_val)
pm_contrast_summary <- if (hap_only) tibble() else map_dfr(names(cohorts), function(cohort_name) {
  if (is.null(results[[cohort_name]]$continuous)) return(NULL)
  map_dfr(names(horizons), function(hname) {
    results[[cohort_name]]$continuous[[hname]]$contrasts %>%
      mutate(cohort = cohort_name, horizon = hname, .before = 1)
  })
}) %>%
  mutate(across(where(is.numeric) & !any_of(c("reference", "contrast")), ~round(., 4)))

cat("\n", strrep("=", 78), "\n", sep = "")
cat("HAP CUMULATIVE RISK BY HORIZON (solid_fuel_use, wood_use, kerosene_use, solid_fuel_primary): risk difference and risk ratio\n")
cat(strrep("=", 78), "\n", sep = "")
print(hap_summary, n = Inf, width = Inf)

cat("\n", strrep("=", 78), "\n", sep = "")
cat("PM2.5 CUMULATIVE-RISK CONTRASTS BY HORIZON (vs. reference = ", ref_val, " ug/m3)\n", sep = "")
cat(strrep("=", 78), "\n", sep = "")
print(pm_contrast_summary, n = Inf, width = Inf)

# Save summary tables
write_csv(hap_summary, file.path(fp, "survival_hap_effects_summary.csv"))
if (!hap_only) write_csv(pm_contrast_summary, file.path(fp, "survival_pm25_contrasts_summary.csv"))
message("\nSaved: ", file.path(fp, "survival_hap_effects_summary.csv"))
if (!hap_only) message("Saved: ", file.path(fp, "survival_pm25_contrasts_summary.csv"))

##############################################################################
### Population x cohort effect estimation (not run)                      ###
##############################################################################
# Commented out. Would fit the binary HAP hazard models for each population
# subgroup except overall (whose results are those above). Continuous PM2.5
# models are not fitted by subgroup because the interval 2 risk sets are too
# small (about 1.7-3.7 observations per parameter).
plan(multisession, workers = parallel::detectCores() - 1)

# population_results <- list()
# for (pop_name in setdiff(names(populations), "overall")) {
#   for (cohort_name in names(cohorts)) {
#     key <- paste0(pop_name, "_", cohort_name)
#     prep <- population_cohort_diagnostics[[key]]
#     if (is.null(prep)) {
#       warning(key, ": no diagnostics available (skipped above) - skipping model fit too.")
#       next
#     }
#     cohort <- cohorts[[cohort_name]]
#     imputations <- prep$imputations
#     run_label <- paste0(pop_name, "_", cohort_name)

#     hap_results <- lapply(names(hap_exposures), function(exp_name) {
#       run_survival_hap(run_label, imputations, prep$hap_covariates[[exp_name]],
#                        hap_exposures[[exp_name]], horizons)
#     }) %>% setNames(names(hap_exposures))

#     # No run_survival_continuous() call here - see comment above the loop.
#     population_results[[key]] <- list(population = pop_name, cohort = cohort_name,
#                                       hap = hap_results, continuous = NULL)
#   }
# }

plan(sequential)

# ##############################################################################
# ### Population x cohort results summary                                   ###
# ##############################################################################
# population_hap_summary <- map_dfr(population_results, function(res) {
#   map_dfr(names(hap_exposures), function(exp_name) {
#     map_dfr(names(horizons), function(hname) {
#       r <- res$hap[[exp_name]][[hname]]
#       bind_rows(
#         mutate(r$rd, population = res$population, cohort = res$cohort,
#               exposure = hap_exposures[[exp_name]], horizon = hname,
#               measure = "Risk difference", .before = 1),
#         mutate(r$rr, population = res$population, cohort = res$cohort,
#               exposure = hap_exposures[[exp_name]], horizon = hname,
#               measure = "Risk ratio", .before = 1)
#       )
#     })
#   })
# }) %>%
#   dplyr::select(population, cohort, exposure, horizon, measure, estimate, se, ci_low, ci_high, df) %>%
#   mutate(across(c(estimate, se, ci_low, ci_high), ~round(., 4)), df = round(df, 1))
# 
# # "overall" rows: same estimates as hap_summary above, just relabeled with
# # population = "overall" rather than refit (see comment above this section).
# population_hap_summary <- bind_rows(
#   hap_summary %>% mutate(population = "overall", .before = 1),
#   population_hap_summary
# )
# 
# # Continuous-PM (longitudinal dose-response) results are reported for
# # "overall" only - see comment above the population loop. population_results
# # entries all have continuous = NULL now, so there's nothing to map over;
# # this table is just pm_contrast_summary (the pooled/overall run above)
# # relabeled, kept as its own object/CSV for a consistent population x
# # cohort x horizon x contrast shape alongside population_hap_summary.
# population_pm_contrast_summary <- pm_contrast_summary %>%
#   mutate(population = "overall", .before = 1)
# 
# cat("\n", strrep("=", 78), "\n", sep = "")
# cat("POPULATION x COHORT HAP CUMULATIVE RISK BY HORIZON: risk difference and risk ratio\n")
# cat(strrep("=", 78), "\n", sep = "")
# print(population_hap_summary, n = Inf, width = Inf)
# 
# cat("\n", strrep("=", 78), "\n", sep = "")
# cat("POPULATION x COHORT PM2.5 CUMULATIVE-RISK CONTRASTS BY HORIZON (vs. reference = ", ref_val, " ug/m3)\n", sep = "")
# cat(strrep("=", 78), "\n", sep = "")
# print(population_pm_contrast_summary, n = Inf, width = Inf)
# 
# write_csv(population_hap_summary, file.path(fp, "population_survival_hap_effects_summary.csv"))
# write_csv(population_pm_contrast_summary, file.path(fp, "population_survival_pm25_contrasts_summary.csv"))
# message("\nSaved: ", file.path(fp, "population_survival_hap_effects_summary.csv"))
# message("Saved: ", file.path(fp, "population_survival_pm25_contrasts_summary.csv"))
