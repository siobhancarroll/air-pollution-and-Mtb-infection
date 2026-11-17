# Trial: wood users / kerosene users versus exclusively clean-fuel users.
# Usage from a fresh R session in this project:
#   source("run_fuel_clean_reference.R")
# Reuses the completed imputations. Does not run imputation or PM2.5 models.
# Runs primary conversion definitions, then baseline positivity
# in the overall population. Normal run_all.R outputs are kept separate.

library(here)
library(dplyr)
library(readr)

trial_stamp <- paste0(Sys.Date(), "-", format(Sys.time(), "%H%M%S"))
trial_log_dir <- here("logs", paste0("fuel-clean-reference-", trial_stamp))
dir.create(trial_log_dir, recursive = TRUE, showWarnings = FALSE)
trial_rscript <- file.path(R.home("bin"), "Rscript")

trial_settings <- paste(
  'fuel_comparator <- "exclusive_clean";',
  'hap_exposure_subset <- c("wood", "kerosene");',
  'hap_only <- TRUE;',
  'analysis_populations <- "overall";',
  paste0('output_stamp <- "', trial_stamp, '";')
)
trial_steps <- list(
  list(name = "01_survival_primary", script = "discrete_time_survival.R",
       settings = 'conversion_definition <- "primary";'),
  list(name = "03_baseline", script = "propensity_scores.R", settings = "")
)

trial_status <- lapply(trial_steps, function(step) {
  message(format(Sys.time(), "%H:%M"), "  running ", step$name)
  code <- paste(trial_settings, step$settings, paste0('source("', step$script, '")'))
  logfile <- file.path(trial_log_dir, paste0(step$name, ".log"))
  status <- system2(trial_rscript, c("-e", shQuote(code)), stdout = logfile, stderr = logfile)
  message("Finished ", step$name, ": exit status ", status)
  data.frame(step = step$name, exit_status = status, log_file = logfile)
})
trial_status <- do.call(rbind, trial_status)
write.csv(trial_status, file.path(trial_log_dir, "step_status.csv"), row.names = FALSE)

trial_dirs <- here("tables_and_figures", paste0(c("survival-", "ps-"),
                                               "cleanref-haponly-", trial_stamp))
writeLines(trial_dirs, file.path(trial_log_dir, "output_directories.txt"))

# Compact review tables, retaining each fit's full imputation-level exports.
trial_counts <- list()
trial_weights <- list()
trial_precision <- list()
trial_failures <- list()
for (folder in trial_dirs) {
  if (!dir.exists(folder)) next
  read_cells <- function(pattern) {
    files <- list.files(folder, pattern = pattern, full.names = TRUE)
    dplyr::bind_rows(lapply(files, function(f) {
      readr::read_csv(f, show_col_types = FALSE) %>%
        dplyr::mutate(analysis = basename(folder), source_file = basename(f), .before = 1)
    }))
  }
  trial_counts[[folder]] <- read_cells("_comparison_counts\\.csv$")
  trial_weights[[folder]] <- read_cells("_fitted_weight_summary\\.csv$")
  trial_precision[[folder]] <- read_cells("_monte_carlo_precision\\.csv$")
  trial_failures[[folder]] <- read_cells("^hap_fit_failures\\.csv$")
}

trial_counts <- dplyr::bind_rows(trial_counts)
if (nrow(trial_counts)) {
  counts_review <- trial_counts %>%
    dplyr::group_by(analysis, source_file, exposure, interval, comparison_arm) %>%
    dplyr::summarise(
      n_people_min = min(n_people), n_people_max = max(n_people),
      events_min = min(n_events), events_max = max(n_events),
      observed_min = min(n_observed), observed_max = max(n_observed),
      prediction_n_min = min(n_prediction_people), prediction_n_max = max(n_prediction_people),
      inconsistent_excluded_min = min(n_inconsistent_excluded),
      inconsistent_excluded_max = max(n_inconsistent_excluded), .groups = "drop")
  readr::write_csv(counts_review, file.path(trial_log_dir, "comparison_counts_review.csv"))
}
trial_weights <- dplyr::bind_rows(trial_weights)
if (nrow(trial_weights)) {
  weight_review <- trial_weights %>%
    dplyr::group_by(analysis, source_file, interval, stage) %>%
    dplyr::summarise(
      ess_min = min(ess), ess_max = max(ess),
      ess_fraction_min = min(ess_fraction),
      clean_arm_ess_min = min(ess_unexposed), exposed_arm_ess_min = min(ess_exposed),
      largest_weight_share = max(max_weight_share),
      largest_abs_imbalance = if (all(is.na(max_abs_imbalance))) NA_real_ else
        max(max_abs_imbalance, na.rm = TRUE),
      diagnostic_errors = sum(!is.na(diagnostic_error) & nzchar(diagnostic_error)),
      .groups = "drop")
  readr::write_csv(weight_review, file.path(trial_log_dir, "fitted_weights_review.csv"))
}
trial_precision <- dplyr::bind_rows(trial_precision)
if (nrow(trial_precision))
  readr::write_csv(trial_precision, file.path(trial_log_dir, "monte_carlo_precision_review.csv"))
trial_failures <- dplyr::bind_rows(trial_failures)
if (nrow(trial_failures))
  readr::write_csv(trial_failures, file.path(trial_log_dir, "fit_failures_review.csv"))

message("Trial logs and review tables: ", trial_log_dir)
message("Trial result folders:\n", paste(trial_dirs, collapse = "\n"))
if (any(trial_status$exit_status != 0) || nrow(trial_failures)) {
  warning("Some trial analyses failed. Inspect step_status.csv, fit_failures_review.csv (if present), and logs. Failed imputations were not dropped for pooling.")
} else {
  message("All trial models completed; assess counts, balance, ESS and precision before choosing these comparisons.")
}
