### Compare survival results across conversion definitions ##################
# Tabulates the primary and stricter-definition runs of
# discrete_time_survival.R side by side, for the PM2.5 contrasts
# (50/60/70/80 vs 40 ug/m3) and the binary HAP exposures.
#
# Requires:
#   1. multiple_imputation.R (primary and sensitivity person-period files)
#   2. discrete_time_survival.R with conversion_definition <- "primary"
#   3. discrete_time_survival.R with conversion_definition <- "sensitivity"
#
# Uses the most recent output folder of each kind unless primary_dir /
# sensitivity_dir are set.

library(tidyverse)
library(here)

latest_dir <- function(pattern) {
  dirs <- list.files(here("tables_and_figures"), pattern = pattern, full.names = TRUE)
  if (length(dirs) == 0) stop("No output folder matching ", pattern)
  dates <- as.Date(sub(".*(\\d{4}-\\d{2}-\\d{2})$", "\\1", dirs))
  dirs[which.max(dates)]
}

if (!exists("primary_dir")) primary_dir <- latest_dir("^survival-\\d{4}-\\d{2}-\\d{2}$")
if (!exists("sensitivity_dir")) sensitivity_dir <- latest_dir("^survival-sensitivity-(df\\d-)?\\d{4}-\\d{2}-\\d{2}$")
message("Primary results:     ", primary_dir)
message("Sensitivity results: ", sensitivity_dir)

definition_labels <- c(
  qgit                = "QGIT: Negative to Positive (primary)",
  qgit_strict         = "QGIT: <0.2 to >0.7 IU/mL",
  tst                 = "TST: <5 to ≥5 mm (primary)",
  tst_strict10        = "TST: <5 to ≥10 mm",
  tst_incr10          = "TST: ≥10 mm and ≥10 mm increase"
)

read_both <- function(file) {
  bind_rows(read_csv(file.path(primary_dir, file), show_col_types = FALSE),
            read_csv(file.path(sensitivity_dir, file), show_col_types = FALSE)) %>%
    mutate(assay = if_else(str_detect(cohort, "^qgit"), "QGIT", "TST"),
           definition = factor(definition_labels[cohort], levels = definition_labels)) %>%
    arrange(assay, definition)
}

pm_comparison <- read_both("survival_pm25_contrasts_summary.csv") %>%
  filter(contrast %in% c(50, 60, 70, 80)) %>%
  transmute(assay, definition, horizon, reference, contrast,
            RD = sprintf("%.3f (%.3f, %.3f)", rd_estimate, rd_ci_low, rd_ci_high),
            RR = sprintf("%.2f (%.2f, %.2f)", rr_estimate, rr_ci_low, rr_ci_high)) %>%
  arrange(horizon, assay, definition, contrast)

hap_comparison <- read_both("survival_hap_effects_summary.csv") %>%
  # Compare outcome definitions only for the main HAP exposure. This also
  # prevents historical combined sensitivity runs entering the review table.
  filter(exposure == "solid_fuel_use") %>%
  transmute(assay, definition, exposure, horizon, measure,
            estimate_ci = sprintf("%.3f (%.3f, %.3f)", estimate, ci_low, ci_high)) %>%
  pivot_wider(names_from = measure, values_from = estimate_ci) %>%
  arrange(horizon, exposure, assay, definition)

print(pm_comparison, n = Inf, width = Inf)
print(hap_comparison, n = Inf, width = Inf)

out_dir <- file.path(here("tables_and_figures"), paste0("conversion-definitions-", get0("output_stamp", ifnotfound = as.character(Sys.Date()))))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write_csv(pm_comparison, file.path(out_dir, "pm25_contrasts_by_conversion_definition.csv"))
write_csv(hap_comparison, file.path(out_dir, "hap_effects_by_conversion_definition.csv"))
writeLines(c(paste("primary_dir:", primary_dir), paste("sensitivity_dir:", sensitivity_dir)),
           file.path(out_dir, "source_runs.txt"))
message("Saved comparison tables to ", out_dir)
