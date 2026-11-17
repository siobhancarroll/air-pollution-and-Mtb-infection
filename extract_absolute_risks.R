### Absolute risks alongside risk contrasts: table-ready extract ##############
# Gathers the pooled standardized absolute risks written by
# discrete_time_survival.R (cumulative conversion risk) and propensity_scores.R
# (baseline prevalence of test positivity) and places them next to the pooled
# RD and RR from the same run.
#
# Absolute risks are pooled with Rubin's rules (nested rules for PM2.5) on the
# logit scale and back-transformed, so CIs stay within 0-1. They come from the
# same fits as the contrasts; RD is NOT exactly exposed minus reference after
# pooling, because each quantity is pooled separately.
#
# Reads the most recent primary folders (survival-<stamp>, ps-<stamp>) unless
# `survival_dir` / `ps_dir` are set before source()-ing. Only aggregate output.
#
# Output: tables_and_figures/absolute-risks-<date>/
#   longitudinal_hap_absolute_risks.csv   cohort x exposure x horizon
#   longitudinal_pm25_absolute_risks.csv  cohort x horizon x PM2.5 level
#   baseline_hap_absolute_risks.csv       population x outcome x exposure
#   baseline_pm25_absolute_risks.csv      population x outcome x PM2.5 level
# Each has numeric columns plus formatted "x.x% (a.a-b.b)" columns.

library(tidyverse)
library(here)

fig_root <- here("tables_and_figures")
latest_dir <- function(prefix) {
  dirs <- list.files(fig_root, pattern = paste0("^", prefix, "-\\d{4}-\\d{2}-\\d{2}(-\\d{6})?$"),
                     full.names = TRUE)
  dirs <- dirs[vapply(dirs, function(d) length(list.files(d, "_absolute_risks\\.csv$")) > 0, logical(1))]
  if (!length(dirs)) stop("No ", prefix, "- folder with *_absolute_risks.csv files; rerun the analysis first.")
  dirs[order(basename(dirs), decreasing = TRUE)][1]
}
if (!exists("survival_dir")) survival_dir <- latest_dir("survival")
if (!exists("ps_dir")) ps_dir <- latest_dir("ps")
# Wood and kerosene are reported versus exclusively clean fuel use, from the
# clean-reference runs (run_all.R steps 02b and 04b). Their rows replace the
# primary-run (versus non-users) rows when a clean-reference folder exists.
latest_clean <- function(prefix) {
  dirs <- list.files(fig_root, pattern = paste0("^", prefix, "-cleanref-haponly-\\d{4}-\\d{2}-\\d{2}(-\\d{6})?$"),
                     full.names = TRUE)
  if (!length(dirs)) return(NA_character_)
  dirs[order(basename(dirs), decreasing = TRUE)][1]
}
if (!exists("survival_clean_dir")) survival_clean_dir <- latest_clean("survival")
if (!exists("ps_clean_dir")) ps_clean_dir <- latest_clean("ps")
clean_fuels <- c("wood_use", "kerosene_use")
swap_clean <- function(primary, clean_dir, file) {
  if (is.na(clean_dir) || !file.exists(file.path(clean_dir, file))) {
    warning("No clean-reference ", file, "; wood and kerosene rows are versus non-users")
    return(primary)
  }
  clean <- read_csv(file.path(clean_dir, file), show_col_types = FALSE) %>% filter(exposure %in% clean_fuels)
  bind_rows(filter(primary, !exposure %in% clean_fuels), clean)
}
message("Clean-reference estimates: ", survival_clean_dir, " | ", ps_clean_dir)
message("Longitudinal estimates: ", survival_dir, "\nBaseline estimates:     ", ps_dir)
out_dir <- file.path(fig_root, paste0("absolute-risks-", get0("output_stamp", ifnotfound = as.character(Sys.Date()))))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_abs <- function(dir) {
  files <- list.files(dir, "_absolute_risks\\.csv$", full.names = TRUE)
  out <- map_dfr(files, ~ read_csv(.x, show_col_types = FALSE) %>% mutate(source_file = basename(.x)))
  for (col in c("exposure", "pm_levels")) if (!col %in% names(out)) out[[col]] <- NA
  out
}
pct <- function(est, lo, hi, digits = 1) {
  sprintf(paste0("%.", digits, "f%% (%.", digits, "f–%.", digits, "f)"),
          100 * est, 100 * lo, 100 * hi)
}
num <- function(est, lo, hi, digits = 2) {
  sprintf(paste0("%.", digits, "f (%.", digits, "f–%.", digits, "f)"), est, lo, hi)
}
# RD/RR summaries -> one row per comparison, RD and RR side by side
wide_effects <- function(eff, keys) {
  eff %>%
    mutate(measure = recode(measure, "Risk difference" = "rd", "Risk ratio" = "rr")) %>%
    select(all_of(keys), measure, estimate, ci_low, ci_high) %>%
    pivot_wider(names_from = measure, values_from = c(estimate, ci_low, ci_high),
                names_glue = "{measure}_{.value}")
}

## Longitudinal (cumulative conversion risk) ##################################
surv_abs <- read_abs(survival_dir)
surv_hap_abs <- surv_abs %>% filter(!is.na(exposure))
if (!is.na(survival_clean_dir)) surv_hap_abs <- bind_rows(filter(surv_hap_abs, !exposure %in% clean_fuels),
  read_abs(survival_clean_dir) %>% filter(exposure %in% clean_fuels))
surv_pm_abs <- surv_abs %>% filter(!is.na(pm_levels))

hap_eff <- read_csv(file.path(survival_dir, "survival_hap_effects_summary.csv"), show_col_types = FALSE) %>%
  swap_clean(survival_clean_dir, "survival_hap_effects_summary.csv") %>%
  wide_effects(c("cohort", "exposure", "horizon"))
long_hap <- surv_hap_abs %>%
  select(cohort, exposure, horizon, arm, estimate, ci_low, ci_high) %>%
  pivot_wider(names_from = arm, values_from = c(estimate, ci_low, ci_high),
              names_glue = "{arm}_risk_{.value}") %>%
  left_join(hap_eff, by = c("cohort", "exposure", "horizon")) %>%
  mutate(reference_risk = pct(reference_risk_estimate, reference_risk_ci_low, reference_risk_ci_high),
         exposed_risk = pct(exposed_risk_estimate, exposed_risk_ci_low, exposed_risk_ci_high),
         risk_difference = num(rd_estimate, rd_ci_low, rd_ci_high),
         risk_ratio = num(rr_estimate, rr_ci_low, rr_ci_high)) %>%
  arrange(cohort, horizon, exposure)
write_csv(long_hap, file.path(out_dir, "longitudinal_hap_absolute_risks.csv"))

pm_file <- file.path(survival_dir, "survival_pm25_contrasts_summary.csv")
long_pm <- surv_pm_abs %>%
  select(cohort, horizon, pm_levels, arm, estimate, ci_low, ci_high) %>%
  mutate(absolute_risk = pct(estimate, ci_low, ci_high))
if (file.exists(pm_file)) {
  pm_eff <- read_csv(pm_file, show_col_types = FALSE)
  long_pm <- long_pm %>% left_join(pm_eff, by = c("cohort", "horizon", "pm_levels" = "contrast"),
                                   suffix = c("", "_contrast"))
}
long_pm <- arrange(long_pm, cohort, horizon, pm_levels)
write_csv(long_pm, file.path(out_dir, "longitudinal_pm25_absolute_risks.csv"))

## Baseline (standardized prevalence of test positivity) ######################
ps_abs <- read_abs(ps_dir)
ps_hap_abs <- ps_abs %>% filter(!is.na(exposure))
if (!is.na(ps_clean_dir)) ps_hap_abs <- bind_rows(filter(ps_hap_abs, !exposure %in% clean_fuels),
  read_abs(ps_clean_dir) %>% filter(exposure %in% clean_fuels))
ps_pm_abs <- ps_abs %>% filter(!is.na(pm_levels))

bhap_eff <- read_csv(file.path(ps_dir, "hap_effects_summary.csv"), show_col_types = FALSE) %>%
  swap_clean(ps_clean_dir, "hap_effects_summary.csv") %>%
  wide_effects(c("population", "exposure", "outcome"))
base_hap <- ps_hap_abs %>%
  select(population, outcome, exposure, arm, estimate, ci_low, ci_high) %>%
  pivot_wider(names_from = arm, values_from = c(estimate, ci_low, ci_high),
              names_glue = "{arm}_risk_{.value}") %>%
  left_join(bhap_eff, by = c("population", "exposure", "outcome")) %>%
  mutate(reference_prevalence = pct(reference_risk_estimate, reference_risk_ci_low, reference_risk_ci_high),
         exposed_prevalence = pct(exposed_risk_estimate, exposed_risk_ci_low, exposed_risk_ci_high),
         risk_difference = num(rd_estimate, rd_ci_low, rd_ci_high),
         risk_ratio = num(rr_estimate, rr_ci_low, rr_ci_high)) %>%
  arrange(population, outcome, exposure)
write_csv(base_hap, file.path(out_dir, "baseline_hap_absolute_risks.csv"))

bpm_file <- file.path(ps_dir, "pm25_contrasts_summary.csv")
base_pm <- ps_pm_abs %>%
  select(population, outcome, pm_levels, arm, estimate, ci_low, ci_high) %>%
  mutate(absolute_prevalence = pct(estimate, ci_low, ci_high))
if (file.exists(bpm_file)) {
  bpm_eff <- read_csv(bpm_file, show_col_types = FALSE)
  base_pm <- base_pm %>% left_join(bpm_eff, by = c("population", "outcome", "pm_levels" = "contrast"),
                                   suffix = c("", "_contrast"))
}
base_pm <- arrange(base_pm, population, outcome, pm_levels)
write_csv(base_pm, file.path(out_dir, "baseline_pm25_absolute_risks.csv"))

writeLines(c(paste("survival_dir:", survival_dir), paste("ps_dir:", ps_dir),
             paste("survival_clean_dir (wood, kerosene):", survival_clean_dir),
             paste("ps_clean_dir (wood, kerosene):", ps_clean_dir),
             "Absolute risks pooled on the logit scale (Rubin; nested rules for PM2.5) and back-transformed.",
             "Longitudinal: standardized cumulative conversion risk by Month 4 / Month 12.",
             "Baseline: standardized prevalence of test positivity at enrolment.",
             "Pooled RD need not equal pooled exposed minus pooled reference risk exactly."),
           file.path(out_dir, "README.txt"))
message("Saved to ", out_dir)
