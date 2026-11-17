### Household fuel use vs ambient PM2.5 #######################################
# Describes how household fuel use relates to assigned ambient (city-month)
# PM2.5, to inform whether the two exposures overlap or confound each other.
# No outcome models are fitted.
#
# Exposures: solid_fuel_use (any polluting fuel), wood_use, kerosene_use,
#   solid_fuel_primary, each vs its manuscript comparator (select_hap_comparison()
#   with fuel_comparator = "exclusive_clean": wood and kerosene users versus
#   exclusively clean fuel users; any polluting fuel is unaffected; polluting
#   primary fuel versus all other primary fuels). Set
#   fuel_comparator <- "nonusers" before source()-ing for users versus non-users.
# PM2.5: (1) baseline 3-month average (mean_pm_recentavg; baseline analysis
#   dataset) and (2) first-interval average, Entry to Month 4
#   (mean_pm_intervalavg; IGRA and TST person-period data, interval 1).
# Strata: overall, and within each site.
#
# Per imputation, a linear model of PM2.5 on the fuel indicator (svyglm with
# household-clustered SEs, matching the main analyses) gives the mean PM2.5
# difference (users - comparator); within-site differences come from the
# site-stratified fits, and "Overall, site-adjusted" adds Site to the model.
# Differences are pooled with Rubin's rules (Barnard-Rubin df, clusters - 1).
# Also reported, averaged over imputations: arm sizes, mean (SD) PM2.5 by arm,
# the standardized mean difference (pooled-SD denominator) and the
# point-biserial correlation (averaged on the Fisher-z scale).
# PM2.5 here is the assigned point value; the exposure draws used in the
# models are not propagated (descriptive purpose).
#
# From project root: source("scripts/audit_fuel_ambient_relationship.R")
# Output (aggregate only): tables_and_figures/fuel-ambient-<date>/

library(tidyverse)
library(here)
library(survey)
source(here("mi_pooling_helpers.R"))   # load_latest, select_hap_comparison, pool_rubin, cluster_lookup

if (!exists("fuel_comparator")) fuel_comparator <- "exclusive_clean"
stopifnot(fuel_comparator %in% c("nonusers", "exclusive_clean"))
out_dir <- file.path(here("tables_and_figures"),
                     paste0("fuel-ambient-", if (fuel_comparator == "nonusers") "nonusers-" else "", get0("output_stamp", ifnotfound = as.character(Sys.Date()))))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fuel_vars <- c(solid_fuel_use = "Any polluting fuel use", wood_use = "Wood use",
               kerosene_use = "Kerosene use", solid_fuel_primary = "Polluting primary cooking fuel")
clusters <- cluster_lookup("household")

## Data sources: list of imputations, each one row per person ################
load_latest("^list_imputed_dataframes_\\d{4}-\\d{2}-\\d{2}\\.Rdata$")
sources <- list(baseline_recentavg = list(
  label = "Baseline 3-month average (baseline analysis population)", pm = "mean_pm_recentavg",
  imps = allImputations$imputations))
for (cohort in c("qgit", "tst")) {
  obj <- paste0("allImputations_survival_", cohort)
  ok <- tryCatch({ load_latest(paste0("^survival_person_period_", cohort, "_\\d{4}-\\d{2}-\\d{2}\\.Rdata$")); TRUE },
                 error = function(e) { message("Skipping ", cohort, " person-period data: ", conditionMessage(e)); FALSE })
  if (!ok || !exists(obj)) next
  sources[[paste0(cohort, "_interval1")]] <- list(
    label = paste0("Entry to Month 4 average (", toupper(sub("qgit", "IGRA", cohort)), " conversion cohort)"),
    pm = "mean_pm_intervalavg",
    imps = lapply(get(obj)$imputations, function(d) d %>% filter(interval == 1) %>% distinct(PID, .keep_all = TRUE)))
}

## One imputation x stratum x exposure #######################################
fit_one <- function(d, fuel, pm, adjust_site = FALSE) {
  d <- d[is.finite(d[[pm]]) & !is.na(d[[fuel]]), , drop = FALSE]
  x <- as.integer(d[[fuel]] == 1); y <- d[[pm]]
  if (length(unique(x)) < 2 || nrow(d) < 10) return(NULL)
  d$.x <- x; d$.y <- y
  des <- svydesign(ids = ~cluster_id, data = d, weights = ~1)
  f <- if (adjust_site && n_distinct(d$Site) > 1) .y ~ .x + Site else .y ~ .x
  fit <- svyglm(f, design = des)
  sd_pool <- sqrt(((sum(x == 1) - 1) * var(y[x == 1]) + (sum(x == 0) - 1) * var(y[x == 0])) / (length(y) - 2))
  tibble(n_exposed = sum(x == 1), n_comparator = sum(x == 0),
         mean_exposed = mean(y[x == 1]), sd_exposed = sd(y[x == 1]),
         mean_comparator = mean(y[x == 0]), sd_comparator = sd(y[x == 0]),
         Q = unname(coef(fit)[".x"]), U = unname(vcov(fit)[".x", ".x"]),
         smd = (mean(y[x == 1]) - mean(y[x == 0])) / sd_pool,
         z = atanh(cor(x, y)), n_clusters = n_distinct(d$cluster_id))
}

results <- list(); arms <- list()
for (src_name in names(sources)) {
  src <- sources[[src_name]]
  imps <- lapply(src$imps, function(d) add_cluster_id(d %>% mutate(PID = as.character(PID)), clusters))
  sites <- sort(unique(as.character(imps[[1]]$Site)))
  strata <- c(list(Overall = list(site = NULL, adj = FALSE), `Overall, site-adjusted` = list(site = NULL, adj = TRUE)),
              setNames(lapply(sites, function(s) list(site = s, adj = FALSE)), sites))
  for (fuel in names(fuel_vars)) for (st in names(strata)) {
    per_imp <- map_dfr(seq_along(imps), function(i) {
      d <- select_hap_comparison(imps[[i]], fuel, fuel_comparator)
      if (!is.null(strata[[st]]$site)) d <- d[as.character(d$Site) == strata[[st]]$site, , drop = FALSE]
      r <- fit_one(d, fuel, src$pm, strata[[st]]$adj)
      if (is.null(r)) NULL else mutate(r, imputation = i)
    })
    if (nrow(per_imp) < length(imps)) {
      message(src_name, " / ", fuel, " / ", st, ": not estimable in ", length(imps) - nrow(per_imp), " imputation(s); skipped")
      next
    }
    pooled <- pool_rubin(per_imp$Q, per_imp$U, nu_com = min(per_imp$n_clusters) - 1)
    results[[length(results) + 1]] <- tibble(
      pm_measure = src$label, exposure = fuel_vars[[fuel]], stratum = st,
      n_exposed = mean(per_imp$n_exposed), n_comparator = mean(per_imp$n_comparator),
      mean_pm_exposed = mean(per_imp$mean_exposed), sd_pm_exposed = mean(per_imp$sd_exposed),
      mean_pm_comparator = mean(per_imp$mean_comparator), sd_pm_comparator = mean(per_imp$sd_comparator),
      mean_difference = pooled$estimate, md_ci_low = pooled$ci_low, md_ci_high = pooled$ci_high,
      md_df = pooled$df, smd = mean(per_imp$smd),
      smd_min = min(per_imp$smd), smd_max = max(per_imp$smd),
      point_biserial_r = tanh(mean(per_imp$z)), m = length(imps))
  }
  # Fuel-use prevalence by site (any polluting fuel etc.), averaged over imputations
  arms[[src_name]] <- map_dfr(seq_along(imps), function(i) imps[[i]] %>%
      filter(is.finite(.data[[src$pm]])) %>%
      group_by(Site) %>%
      summarise(n = n(), across(all_of(names(fuel_vars)), ~ mean(.x == 1)),
                pm_mean = mean(.data[[src$pm]]), pm_sd = sd(.data[[src$pm]]),
                pm_p05 = quantile(.data[[src$pm]], 0.05), pm_p95 = quantile(.data[[src$pm]], 0.95),
                .groups = "drop")) %>%
    group_by(Site) %>% summarise(across(everything(), mean), .groups = "drop") %>%
    mutate(pm_measure = src$label, .before = 1)
}

res <- bind_rows(results) %>%
  mutate(across(where(is.double), ~ round(., 3)),
         mean_difference_fmt = sprintf("%.1f (%.1f to %.1f)", mean_difference, md_ci_low, md_ci_high),
         exposed_fmt = sprintf("%.1f (%.1f)", mean_pm_exposed, sd_pm_exposed),
         comparator_fmt = sprintf("%.1f (%.1f)", mean_pm_comparator, sd_pm_comparator))
write_csv(res, file.path(out_dir, "fuel_ambient_pm_differences.csv"))
write_csv(bind_rows(arms) %>% mutate(across(where(is.double), ~ round(., 3))),
          file.path(out_dir, "fuel_prevalence_and_pm_by_site.csv"))
writeLines(c(paste("fuel_comparator:", fuel_comparator),
             paste("m:", length(allImputations$imputations)),
             "mean_difference: PM2.5 (ug/m3) among fuel users minus comparator; Rubin-pooled, household-clustered SEs.",
             "'Overall, site-adjusted' = within-site difference averaged across sites (Site in the model).",
             "smd uses the pooled-SD denominator; smd_min/max = range across imputations (not a CI).",
             "point_biserial_r averaged on the Fisher-z scale across imputations.",
             "PM2.5 is the assigned city-month value; exposure draws are not propagated."),
           file.path(out_dir, "README.txt"))
print(res %>% select(pm_measure, exposure, stratum, n_exposed, n_comparator, exposed_fmt, comparator_fmt,
                     mean_difference_fmt, smd, point_biserial_r), n = Inf, width = Inf)
message("Saved to ", out_dir)
