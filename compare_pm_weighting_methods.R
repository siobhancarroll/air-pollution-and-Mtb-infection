### Compare weighting strategies for the continuous PM2.5 exposure #############
# PM2.5 is assigned by city and calendar month, so season and site predict it
# strongly and a normal-linear GPS model may not balance them. This script
# compares, on the first `n_imp` imputations:
#   glm              - stabilized normal-linear GPS
#   glm_site_season  - glm plus a site x season interaction in the exposure model
#   cbps             - covariate balancing propensity score
#   ebal             - entropy balancing
#   glm_temp         - glm with temperature in place of season
#   cbps_temp        - cbps with temperature in place of season
#   cbps_season_temp - cbps with both season and temperature
# and, for the baseline analysis, each method restricted to the PM2.5 range
# observed at both sites (common support).
#
# Balance is reported on season, temperature and site whether or not they are
# in the exposure model. Baseline estimates (RD at 60 and 80 vs 40 ug/m3,
# observed city-month PM2.5, point estimates averaged over imputations) are
# shown for a 3-df spline and a linear term, to show how much the choice of
# method matters.
#
# Covariate lists are copied from propensity_scores.R and
# discrete_time_survival.R (plus any polluting fuel use); keep them in sync.
#
# Outputs, tables_and_figures/pm-weighting-<date>/:
#   baseline_balance_by_method.csv, baseline_estimates_by_method.csv,
#   survival_balance_by_method.csv, common_support.csv

library(tidyverse)
library(here)
library(mitools)
library(WeightIt)
library(cobalt)
library(splines)

source(here("mi_pooling_helpers.R"))
set.seed(123)
if (!exists("n_imp")) n_imp <- 5
# run_survival <- FALSE runs the baseline comparison only
if (!exists("run_survival")) run_survival <- TRUE
err_rows <- list()

out_dir <- file.path(here("tables_and_figures"), paste0("pm-weighting-", Sys.Date()))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

### Baseline data ###############################################################
load_latest("^list_imputed_dataframes_\\d{4}-\\d{2}-\\d{2}\\.Rdata$")
base_imps <- lapply(allImputations$imputations[seq_len(n_imp)], collapse_covariates)

cov_overall <- c("Income_collapsed", "Education_collapsed", "Employment_collapsed",
                 "Living.Area", "B_Age", "B_Gender", "Site", "season",
                 "alcohol_disorder", "BMI", "B_Height",
                 "B_HbA1c", "B_BCG", "B_Sleepindex", "Type.of.House",
                 "smoking_status", "Household.size", "solid_fuel_use")
cov_children <- c("Income_collapsed", "Living.Area_collapsed", "B_Age", "B_Gender",
                  "Site", "season", "BMI_for_age", "Household.size", "solid_fuel_use")
prepare_children <- function(df) {
  zscorer::addWGSR(data = df[df$B_Age < 15, ], sex = "B_Gender_numeric",
                   firstPart = "B_Weight", secondPart = "B_Height",
                   thirdPart = "B_Age_months", index = "bfa", output = "BMI_for_age")
}
base_pops <- list(
  overall  = list(prep = function(d) d,                  covs = cov_overall),
  children = list(prep = prepare_children,               covs = cov_children),
  adults   = list(prep = function(d) d[d$B_Age >= 15, ], covs = cov_overall)
)

### Helpers #####################################################################
make_specs <- function(covs, season_var, temp_var, has_temp) {
  sp <- list(
    glm             = list(covs = covs, method = "glm"),
    glm_site_season = list(covs = c(covs, paste0("Site:", season_var)), method = "glm"),
    cbps            = list(covs = covs, method = "cbps"),
    ebal            = list(covs = covs, method = "ebal")
  )
  if (has_temp) {
    no_season <- setdiff(covs, season_var)
    sp$glm_temp         <- list(covs = c(no_season, temp_var), method = "glm")
    sp$cbps_temp        <- list(covs = c(no_season, temp_var), method = "cbps")
    sp$cbps_season_temp <- list(covs = c(covs, temp_var), method = "cbps")
  }
  sp
}

fit_w <- function(d, exposure, spec) {
  f <- as.formula(paste(exposure, "~", paste(spec$covs, collapse = " + ")))
  args <- list(formula = f, data = d, method = spec$method)
  if (spec$method == "glm") args$stabilize <- TRUE
  suppressWarnings(do.call(weightit, args))
}

balance_row <- function(w, d, season_var, temp_var, has_temp) {
  addl_vars <- c(season_var, if (has_temp) temp_var, "Site")
  addl_vars <- addl_vars[sapply(addl_vars, function(v) length(unique(na.omit(d[[v]]))) > 1)]
  bt <- bal.tab(w, un = TRUE, stats = "correlations",
                addl = as.formula(paste("~", paste(addl_vars, collapse = " + "))), data = d)
  b <- bt$Balance
  adj <- abs(b$Corr.Adj); un <- abs(b$Corr.Un); nm <- rownames(b)
  mx <- function(v, pat) { k <- grepl(pat, nm); if (any(k)) max(v[k], na.rm = TRUE) else NA_real_ }
  ww <- w$weights
  tibble(n = nrow(d), ess_pct = 100 * sum(ww)^2 / sum(ww^2) / length(ww), max_weight = max(ww),
         max_abs_corr_all = max(adj, na.rm = TRUE),
         n_covariates_over_0.1 = sum(adj > 0.1, na.rm = TRUE),
         season_before = mx(un, "^season"), season_after = mx(adj, "^season"),
         temp_before = mx(un, "^temp"), temp_after = mx(adj, "^temp"),
         site_after = mx(adj, "^Site"))
}

contrast_rd <- function(d, y, w, df_spline) {
  pm <- d$mean_pm_recentavg
  basis <- if (df_spline == 1) NULL else ns(pm, df = df_spline)
  X <- function(x) if (is.null(basis)) data.frame(b1 = x) else setNames(as.data.frame(predict(basis, x)), paste0("b", 1:df_spline))
  dd <- bind_cols(d, X(pm)); keep <- !is.na(dd[[y]])
  f <- as.formula(paste(y, "~", paste(names(X(40)), collapse = " + ")))
  fit <- suppressWarnings(glm(f, data = dd[keep, ], weights = w[keep], family = quasibinomial()))
  r <- function(x) as.numeric(predict(fit, newdata = X(x), type = "response"))
  tibble(contrast = c(60, 80), rd = c(r(60), r(80)) - r(40))
}

### Common support for the baseline exposure ####################################
support <- function(d) {
  s <- d %>% group_by(Site) %>% summarise(lo = quantile(mean_pm_recentavg, .05, na.rm = TRUE),
                                          hi = quantile(mean_pm_recentavg, .95, na.rm = TRUE), .groups = "drop")
  c(lo = max(s$lo), hi = min(s$hi))
}

### Baseline run ################################################################
bal_rows <- list(); est_rows <- list(); sup_rows <- list()
for (pn in names(base_pops)) {
  for (i in seq_len(n_imp)) {
    d0 <- tryCatch(base_pops[[pn]]$prep(base_imps[[i]]), error = function(e) e)
    if (inherits(d0, "error")) {
      err_rows[[length(err_rows) + 1]] <- tibble(part = "baseline prepare", population = pn, method = NA, imputation = i, error = conditionMessage(d0))
      next
    }
    if (i == 1) message(pn, ": temperature available = ", "temp_recentavg" %in% names(d0))
    has_temp <- "temp_recentavg" %in% names(d0)
    cs <- support(d0)
    if (i == 1) sup_rows[[pn]] <- tibble(population = pn, common_lo = cs[["lo"]], common_hi = cs[["hi"]],
                                         n = nrow(d0), n_in_common = sum(d0$mean_pm_recentavg >= cs[["lo"]] & d0$mean_pm_recentavg <= cs[["hi"]]))
    specs <- make_specs(base_pops[[pn]]$covs, "season", "temp_recentavg", has_temp)
    for (restrict in c(FALSE, TRUE)) {
      d <- if (restrict) droplevels(d0[d0$mean_pm_recentavg >= cs[["lo"]] & d0$mean_pm_recentavg <= cs[["hi"]], ]) else d0
      for (sn in names(specs)) {
        lab <- paste0(sn, if (restrict) " (common support)" else "")
        w <- tryCatch(fit_w(d, "mean_pm_recentavg", specs[[sn]]), error = function(e) e)
        if (inherits(w, "error")) {
          err_rows[[length(err_rows) + 1]] <- tibble(part = "baseline weights", population = pn, method = lab, imputation = i, error = conditionMessage(w))
          next
        }
        br <- tryCatch(balance_row(w, d, "season", "temp_recentavg", has_temp), error = function(e) e)
        if (inherits(br, "error")) {
          err_rows[[length(err_rows) + 1]] <- tibble(part = "baseline balance", population = pn, method = lab, imputation = i, error = conditionMessage(br))
          next
        }
        bal_rows[[length(bal_rows) + 1]] <- br %>% mutate(population = pn, method = lab, imputation = i, .before = 1)
        for (y in c(IGRA = "qgit_positive_entry", TST = "tst_positive_entry")) for (dfs in c(3, 1)) {
          r <- tryCatch(contrast_rd(d, y, w$weights, dfs), error = function(e) {
            err_rows[[length(err_rows) + 1]] <<- tibble(part = "baseline estimate", population = pn, method = lab, imputation = i, error = conditionMessage(e)); NULL })
          if (!is.null(r)) est_rows[[length(est_rows) + 1]] <- r %>%
            mutate(population = pn, method = lab, outcome = if (y == "qgit_positive_entry") "IGRA" else "TST",
                   pm_term = if (dfs == 3) "spline (3 df)" else "linear", imputation = i)
        }
      }
    }
    message("baseline: ", pn, " imputation ", i, " done")
  }
}
summ <- function(x, by) x %>% group_by(across(all_of(by))) %>%
  summarise(across(where(is.numeric) & !any_of("imputation"), ~ round(mean(.x, na.rm = TRUE), 3)),
            n_imputations = n(), .groups = "drop")
baseline_balance <- summ(bind_rows(bal_rows), c("population", "method"))
write_csv(baseline_balance, file.path(out_dir, "baseline_balance_by_method.csv"))
print(baseline_balance, n = Inf, width = Inf)
baseline_est <- bind_rows(est_rows) %>% group_by(population, outcome, pm_term, method, contrast) %>%
  summarise(rd = round(mean(rd), 3), .groups = "drop") %>%
  pivot_wider(names_from = contrast, values_from = rd, names_prefix = "rd_at_")
write_csv(baseline_est, file.path(out_dir, "baseline_estimates_by_method.csv"))
print(baseline_est, n = Inf, width = Inf)
write_csv(bind_rows(sup_rows), file.path(out_dir, "common_support.csv"))
print(bind_rows(sup_rows))

### Survival (conversion) exposure models: balance per interval ################
if (run_survival) {
cov_surv <- c("B_Age", "B_Gender", "Site", "season_collapsed",
              "Income_collapsed", "Education_collapsed", "Employment_collapsed",
              "Type.of.House", "Living.Area_collapsed",
              "BMI", "Household.size", "B_BCG", "solid_fuel_use")
truncate_pm_gap <- function(df) {
  first_na <- df %>% filter(is.na(mean_pm_intervalavg)) %>% group_by(PID) %>%
    summarise(first_na = min(interval), .groups = "drop")
  df %>% left_join(first_na, by = "PID") %>% filter(is.na(first_na) | interval < first_na) %>% select(-first_na)
}
surv_rows <- list()
for (coh in c("qgit", "tst")) {
  load_latest(paste0("^survival_person_period_", coh, "_\\d{4}-\\d{2}-\\d{2}\\.Rdata$"))
  imps <- get(paste0("allImputations_survival_", coh))$imputations
  for (i in seq_len(n_imp)) {
    pp <- truncate_pm_gap(collapse_covariates(imps[[i]]) %>% filter(interval <= 2))
    has_temp <- "temp_intervalavg" %in% names(pp)
    specs <- make_specs(cov_surv, "season_collapsed", "temp_intervalavg", has_temp)
    for (iv in 1:2) {
      d <- droplevels(pp[pp$interval == iv, ])
      for (sn in names(specs)) {
        w <- tryCatch(fit_w(d, "mean_pm_intervalavg", specs[[sn]]), error = function(e) e)
        if (inherits(w, "error")) {
          err_rows[[length(err_rows) + 1]] <- tibble(part = "survival weights", population = paste(coh, iv), method = sn, imputation = i, error = conditionMessage(w))
          next
        }
        surv_rows[[length(surv_rows) + 1]] <- balance_row(w, d, "season_collapsed", "temp_intervalavg", has_temp) %>%
          mutate(cohort = coh, interval = iv, method = sn, imputation = i, .before = 1)
      }
    }
    message("survival: ", coh, " imputation ", i, " done")
  }
}
surv_balance <- summ(bind_rows(surv_rows), c("cohort", "interval", "method"))
write_csv(surv_balance, file.path(out_dir, "survival_balance_by_method.csv"))
print(surv_balance, n = Inf, width = Inf)
}
errs <- bind_rows(err_rows)
write_csv(errs, file.path(out_dir, "errors.csv"))
if (nrow(errs)) { message(nrow(errs), " fits failed - see errors.csv"); print(distinct(errs, part, population, method, error), n = 50) }
message("Saved to ", out_dir)
