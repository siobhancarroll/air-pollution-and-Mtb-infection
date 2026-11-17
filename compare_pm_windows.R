### PM2.5 exposure-timing sensitivity analyses: compiled table ################
# Collects the PM2.5 contrasts (60 and 80 vs 40 ug/m3) and absolute
# cumulative risks from the primary run (interval average) and each
# alternative-window run (survival-window-<name>-<date>; see pm_windows.R).
# Output: tables_and_figures/pm-windows-<date>/pm_windows_summary.csv
library(tidyverse); library(here)
fig_root <- here("tables_and_figures")
newest <- function(pattern) {
  d <- list.files(fig_root, pattern = pattern, full.names = TRUE)
  d <- d[file.exists(file.path(d, "survival_pm25_contrasts_summary.csv"))]
  if (!length(d)) return(NA_character_)
  d[order(basename(d), decreasing = TRUE)][1]
}
windows <- c(interval = "^survival-\\d{4}-\\d{2}-\\d{2}$",
             pre_entry3 = "^survival-window-pre_entry3-\\d{4}-\\d{2}-\\d{2}$",
             early2 = "^survival-window-early2-\\d{4}-\\d{2}-\\d{2}$",
             first_interval = "^survival-window-first_interval-\\d{4}-\\d{2}-\\d{2}$")
dirs <- vapply(windows, newest, character(1))
print(dirs)
out <- map_dfr(names(dirs)[!is.na(dirs)], function(w) {
  d <- dirs[[w]]
  con <- read_csv(file.path(d, "survival_pm25_contrasts_summary.csv"), show_col_types = FALSE) %>%
    filter(contrast %in% c(60, 80)) %>%
    select(cohort, horizon, contrast, starts_with("rd_estimate"), rd_ci_low, rd_ci_high, rr_estimate, rr_ci_low, rr_ci_high,
           pm_p05_achieved, pm_p95_achieved, contrast_in_range, reference_in_range)
  abs <- map_dfr(list.files(d, "_intervalavg_month(4|12)_absolute_risks\\.csv$", full.names = TRUE), read_csv, show_col_types = FALSE) %>%
    select(cohort, horizon, pm_levels, risk = estimate, risk_low = ci_low, risk_high = ci_high)
  ref <- abs %>% filter(pm_levels == 40) %>% select(cohort, horizon, ref_risk = risk, ref_low = risk_low, ref_high = risk_high)
  con %>% left_join(abs, by = c("cohort", "horizon", "contrast" = "pm_levels")) %>%
    left_join(ref, by = c("cohort", "horizon")) %>% mutate(window = w, source = basename(d), .before = 1)
})
od <- file.path(fig_root, paste0("pm-windows-", get0("output_stamp", ifnotfound = as.character(Sys.Date())))); dir.create(od, showWarnings = FALSE, recursive = TRUE)
write_csv(out, file.path(od, "pm_windows_summary.csv"))
message("Saved to ", od)
