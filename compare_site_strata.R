### Site-stratified sensitivity analysis: combine results ###################
# Tabulates PM2.5 results from the site-specific runs alongside the pooled
# (both-site) results.
#
# Requires (each in a fresh R session):
#   site_subset <- "NIRT";  source("propensity_scores.R")
#   site_subset <- "BJGMC"; source("propensity_scores.R")
# Site-specific conversion models are not fitted because the per-site risk
# sets are too small; survival-site-* folders are included only if they exist.
# method_tag must match the run tag used for the site runs, if any.
#
# Site-specific estimates use each site's own 5th-95th percentile exposure
# range, so some contrasts may fall outside the data (see in_range).
# RD_difference is the Chennai (NIRT) minus Pune (BJGMC) risk difference, with a
# 95% CI from the two site-specific standard errors.
#
# Output: tables_and_figures/site-strata-<date>/pm25_by_site.csv

library(tidyverse)
library(here)

if (!exists("method_tag")) method_tag <- ""
tf <- here("tables_and_figures")
latest <- function(prefix) {
  dirs <- list.files(tf, pattern = paste0("^", prefix, "\\d{4}-\\d{2}-\\d{2}$"), full.names = TRUE)
  if (!length(dirs)) stop("No folder matching ", prefix, "<date>")
  dirs[which.max(as.Date(str_extract(dirs, "\\d{4}-\\d{2}-\\d{2}$")))]
}
sites <- c(both = "", NIRT = "site-NIRT-", BJGMC = "site-BJGMC-")
site_names <- c(both = "Both sites (pooled)", NIRT = "Chennai (NIRT)", BJGMC = "Pune (BJGMC)")

has_surv_sites <- all(sapply(c("NIRT", "BJGMC"), function(s)
  length(list.files(tf, pattern = paste0("^survival-site-", s, "-", method_tag, "\\d{4}"))) > 0))
surv <- if (!has_surv_sites) NULL else map_dfr(names(sites), function(s) {
  d <- latest(paste0("survival-", sites[[s]], method_tag))
  message(s, ": ", d)
  read_csv(file.path(d, "survival_pm25_contrasts_summary.csv"), show_col_types = FALSE) %>%
    mutate(site = s, analysis = paste0("Conversion, ", if_else(horizon == "month4", "Month 4", "Month 12")),
           assay = if_else(cohort == "qgit", "IGRA", "TST"))
})
base <- map_dfr(names(sites), function(s) {
  d <- latest(paste0("ps-", sites[[s]], method_tag))
  message(s, ": ", d)
  read_csv(file.path(d, "pm25_contrasts_summary.csv"), show_col_types = FALSE) %>%
    filter(population == "overall") %>%
    mutate(site = s, analysis = "Baseline prevalence", assay = if_else(outcome == "qgit", "IGRA", "TST"))
})
all <- bind_rows(surv, base) %>% filter(contrast %in% c(50, 60, 70, 80))

fmt <- function(e, l, h, n) sprintf(paste0("%.", n, "f (%.", n, "f, %.", n, "f)"), e, l, h)
het <- all %>% filter(site != "both") %>%
  select(analysis, assay, contrast, site, rd_estimate, rd_se, rr_estimate, rr_se) %>%
  pivot_wider(names_from = site, values_from = c(rd_estimate, rd_se, rr_estimate, rr_se)) %>%
  transmute(analysis, assay, contrast,
            rd_diff = rd_estimate_NIRT - rd_estimate_BJGMC,
            rd_diff_se = sqrt(rd_se_NIRT^2 + rd_se_BJGMC^2)) %>%
  transmute(analysis, assay, contrast,
            RD_difference = fmt(rd_diff, rd_diff - 1.96 * rd_diff_se, rd_diff + 1.96 * rd_diff_se, 3))

out <- all %>%
  transmute(analysis, assay, site = site_names[site], contrast,
            RD = fmt(rd_estimate, rd_ci_low, rd_ci_high, 3),
            RR = fmt(rr_estimate, rr_ci_low, rr_ci_high, 2),
            in_range = contrast_in_range,
            pm_p05 = round(pm_p05_achieved, 1), pm_p95 = round(pm_p95_achieved, 1)) %>%
  left_join(het, by = c("analysis", "assay", "contrast")) %>%
  mutate(analysis = factor(analysis, c("Conversion, Month 4", "Conversion, Month 12", "Baseline prevalence")),
         site = factor(site, site_names)) %>%
  arrange(analysis, assay, contrast, site)

out_dir <- file.path(tf, paste0("site-strata-", method_tag, get0("output_stamp", ifnotfound = as.character(Sys.Date()))))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write_csv(out, file.path(out_dir, "pm25_by_site.csv"))
print(out, n = Inf, width = Inf)
message("Saved to ", out_dir)
