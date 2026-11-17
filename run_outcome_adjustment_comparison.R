# Planned current-versus-full outcome-adjustment comparison. No imputation.
# source("run_outcome_adjustment_comparison.R") from a fresh project R session.
# Current/full specifications share cached weights and PM draws; hash exports
# verify the weighted fitting inputs and prediction populations actually match.
library(here)
library(dplyr)
library(readr)
if (!requireNamespace("digest", quietly = TRUE)) stop("Install digest before running the comparison")
source(here("tests", "test_outcome_adjustment_comparison.R"), local = new.env(parent = globalenv()))
adjustment_resuming <- exists("adjustment_resume_stamp")
adjustment_stamp <- if (adjustment_resuming) as.character(adjustment_resume_stamp) else
  paste0(Sys.Date(), "-", format(Sys.time(), "%H%M%S"))
stopifnot(length(adjustment_stamp) == 1L,
  grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}-[0-9]{6}$", adjustment_stamp))
adjustment_log_dir <- here("logs", paste0("outcome-adjustment-comparison-", adjustment_stamp))
adjustment_cache <- here("Data", "outcome_adjustment_comparison", adjustment_stamp, "inputs")
dir.create(adjustment_log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(adjustment_cache, recursive = TRUE, showWarnings = FALSE)
adjustment_rscript <- file.path(R.home("bin"), "Rscript")

# Record exactly which dated input objects were available at the start.
adjustment_input_files <- list.files(here("Data"), pattern = "\\.Rdata$", full.names = TRUE)
adjustment_input_files <- adjustment_input_files[grepl(
  "^(list_imputed_dataframes_|long_imputed_dataframes_|survival_person_period_)", basename(adjustment_input_files))]
adjustment_input_groups <- split(adjustment_input_files,
  sub("_[0-9]{4}-[0-9]{2}-[0-9]{2}\\.Rdata$", "", basename(adjustment_input_files)))
adjustment_inputs <- vapply(adjustment_input_groups, function(x) tail(sort(x), 1), character(1))
adjustment_manifest_path <- file.path(adjustment_log_dir, "input_manifest.csv")
if (adjustment_resuming) {
  if (!file.exists(adjustment_manifest_path)) stop("Cannot resume without the original input manifest")
  adjustment_manifest <- read_csv(adjustment_manifest_path,
    col_types = cols(.default = col_character()), show_col_types = FALSE)
  # Completed steps are reusable only with the original analysis code.
  analysis_sources <- here(c("discrete_time_survival.R", "propensity_scores.R",
    "outcome_adjustment_helpers.R", "mi_pooling_helpers.R", "fuel_processing_helpers.R"))
  if (any(file.info(analysis_sources)$mtime > file.info(adjustment_manifest_path)$mtime))
    stop("Analysis source files changed since this run started; start a fresh comparison instead")
} else {
  adjustment_input_hash <- tools::md5sum(adjustment_inputs)
  if (anyNA(adjustment_input_hash)) stop("Could not checksum all comparison inputs")
  adjustment_manifest <- tibble(file = names(adjustment_input_hash), md5 = unname(adjustment_input_hash))
  write_csv(adjustment_manifest, adjustment_manifest_path)
}

verify_adjustment_inputs <- function() {
  # Compare checksums by explicit path, rather than identical() on a named
  # vector. Retain evidence for every unreadable, missing or changed input.
  current <- unname(tools::md5sum(adjustment_manifest$file))
  unreadable <- which(is.na(current))
  if (length(unreadable)) current[unreadable] <-
    unname(tools::md5sum(adjustment_manifest$file[unreadable]))
  changed <- is.na(current) | current != adjustment_manifest$md5
  if (any(changed)) {
    report <- file.path(adjustment_log_dir, "input_checksum_mismatches.csv")
    write_csv(tibble(file = adjustment_manifest$file[changed],
      expected_md5 = adjustment_manifest$md5[changed], current_md5 = current[changed]), report)
    stop("Comparison input verification failed; see ", report,
      ". No affected step was started and no datasets were mixed for pooling.")
  }
  invisible(TRUE)
}
verify_adjustment_inputs()

adjustment_scenarios <- list(
  list(label = "survival_primary", script = "discrete_time_survival.R",
       settings = 'conversion_definition <- "primary";', prefix = "survival-"),
  list(label = "survival_sensitivity", script = "discrete_time_survival.R",
       settings = 'conversion_definition <- "sensitivity";', prefix = "survival-sensitivity-"),
  list(label = "baseline", script = "propensity_scores.R",
       settings = 'analysis_populations <- "overall";', prefix = "ps-"),
  list(label = "survival_primary_clean", script = "discrete_time_survival.R",
       settings = 'conversion_definition <- "primary"; fuel_comparator <- "exclusive_clean"; hap_only <- TRUE; hap_exposure_subset <- c("wood", "kerosene");',
       prefix = "survival-cleanref-haponly-"),
  list(label = "baseline_clean", script = "propensity_scores.R",
       settings = 'analysis_populations <- "overall"; fuel_comparator <- "exclusive_clean"; hap_only <- TRUE; hap_exposure_subset <- c("wood", "kerosene");',
       prefix = "ps-cleanref-haponly-")
)

# Quote paths as R strings; system2() receives one shell-quoted expression.
adjustment_records <- list()
adjustment_previous_records <- if (adjustment_resuming &&
    file.exists(file.path(adjustment_log_dir, "step_status.csv")))
  read_csv(file.path(adjustment_log_dir, "step_status.csv"), show_col_types = FALSE) else tibble()
for (scenario in adjustment_scenarios) for (mode in c("current", "full")) {
  verify_adjustment_inputs()
  previous_step <- if (nrow(adjustment_previous_records)) adjustment_previous_records %>%
    filter(analysis == scenario$label, outcome_adjustment == mode, exit_status == 0) else tibble()
  if (nrow(previous_step) == 1L && dir.exists(previous_step$output_dir) &&
      length(list.files(previous_step$output_dir, "_outcome_model_diagnostics\\.csv$")) > 0L) {
    message("Reusing completed step: ", scenario$label, " / ", mode)
    adjustment_records[[length(adjustment_records) + 1L]] <- previous_step
    next
  }
  code <- paste('outcome_model_comparison <- TRUE;',
    paste0('outcome_adjustment_mode <- ', encodeString(mode, quote = '"'), ';'),
    paste0('outcome_comparison_cache <- ', encodeString(adjustment_cache, quote = '"'), ';'),
    paste0('output_stamp <- ', encodeString(adjustment_stamp, quote = '"'), ';'),
    scenario$settings, paste0('source(', encodeString(scenario$script, quote = '"'), ')'))
  logfile <- file.path(adjustment_log_dir, paste0(scenario$label, "_", mode, ".log"))
  message(format(Sys.time(), "%H:%M"), "  ", scenario$label, " / ", mode)
  status <- system2(adjustment_rscript, c("-e", shQuote(code)), stdout = logfile, stderr = logfile)
  folder <- here("tables_and_figures", paste0(scenario$prefix, "outcome-", mode, "-", adjustment_stamp))
  adjustment_records[[length(adjustment_records) + 1L]] <- tibble(
    analysis = scenario$label, outcome_adjustment = mode, exit_status = status,
    output_dir = folder, log_file = logfile)
  write_csv(bind_rows(adjustment_records), file.path(adjustment_log_dir, "step_status.csv"))
}
adjustment_records <- bind_rows(adjustment_records)
verify_adjustment_inputs()

# Collect model-level diagnostics, including failed cells before pooling.
adjustment_diagnostics <- list(); adjustment_failures <- list(); adjustment_effects <- list()
for (i in seq_len(nrow(adjustment_records))) {
  run <- adjustment_records[i, ];folder <- run$output_dir
  if (!dir.exists(folder)) next
  for (file in list.files(folder, "_outcome_model_diagnostics\\.csv$", full.names = TRUE)) {
    d <- read_csv(file, col_types = cols(.default = col_character()), show_col_types = FALSE)
    d$outcome_adjustment <- run$outcome_adjustment
    adjustment_diagnostics[[length(adjustment_diagnostics) + 1L]] <- d %>% mutate(
      analysis = run$analysis, run_label = sub("_outcome_model_diagnostics\\.csv$", "", basename(file)), .before = 1)
  }
  file <- file.path(folder, "model_fit_failures.csv")
  if (file.exists(file)) adjustment_failures[[length(adjustment_failures) + 1L]] <-
    read_csv(file, col_types = cols(.default = col_character()), show_col_types = FALSE) %>%
      mutate(analysis = run$analysis, outcome_adjustment = run$outcome_adjustment, .before = 1)
  is_survival <- startsWith(run$analysis, "survival")
  file <- file.path(folder, if (is_survival) "survival_hap_effects_summary.csv" else "hap_effects_summary.csv")
  if (file.exists(file) && file.info(file)$size > 1) {
    d <- read_csv(file, show_col_types = FALSE)
    if (nrow(d)) {
      d <- d %>% mutate(cohort = if (is_survival) cohort else outcome,
        horizon = if (is_survival) horizon else "baseline",
        run_label = if (is_survival) paste0(cohort, "_", exposure) else paste0(population, "_", exposure, "_", outcome),
        reference = NA_real_, contrast = NA_real_) %>%
        select(cohort, exposure, horizon, run_label, reference, contrast, measure, estimate, se, ci_low, ci_high)
      adjustment_effects[[length(adjustment_effects) + 1L]] <- d %>% mutate(
        analysis = run$analysis, outcome_adjustment = run$outcome_adjustment, .before = 1)
    }
  }
  file <- file.path(folder, if (is_survival) "survival_pm25_contrasts_summary.csv" else "pm25_contrasts_summary.csv")
  if (file.exists(file) && file.info(file)$size > 1) {
    d <- read_csv(file, show_col_types = FALSE)
    if (nrow(d)) for (effect_prefix in c("rd", "rr")) {
      e <- d %>% mutate(cohort = if (is_survival) cohort else outcome,
        horizon = if (is_survival) horizon else "baseline", exposure = "pm25",
        run_label = if (is_survival) paste0(cohort, "_intervalavg") else paste0(population, "_", outcome, "_recentavg"),
        measure = if (effect_prefix == "rd") "Risk difference" else "Risk ratio",
        estimate = .data[[paste0(effect_prefix, "_estimate")]], se = .data[[paste0(effect_prefix, "_se")]],
        ci_low = .data[[paste0(effect_prefix, "_ci_low")]], ci_high = .data[[paste0(effect_prefix, "_ci_high")]]) %>%
        select(cohort, exposure, horizon, run_label, reference, contrast, measure, estimate, se, ci_low, ci_high)
      adjustment_effects[[length(adjustment_effects) + 1L]] <- e %>% mutate(
        analysis = run$analysis, outcome_adjustment = run$outcome_adjustment, .before = 1)
    }
  }
}
adjustment_diagnostics <- bind_rows(tibble(analysis = character(), run_label = character(),
  outcome_adjustment = character(), imputation = character(), pm_draw = character(),
  status = character(), error = character(), fit_input_hash = character(), prediction_population_hash = character()),
  bind_rows(adjustment_diagnostics))
adjustment_failures <- bind_rows(adjustment_failures)
adjustment_effects <- bind_rows(adjustment_effects)
if (nrow(adjustment_diagnostics)) {
  write_csv(adjustment_diagnostics, file.path(adjustment_log_dir, "outcome_model_diagnostics.csv"))
  # One compact table for the survival algorithm checks; all failed checks
  # also trigger the strict fit guard, so affected cells cannot be pooled.
  verification_columns <- grep("^(original_|alternate_start_|tight_tolerance_|initialization$|max_iterations$|iterations$)",
    names(adjustment_diagnostics), value = TRUE)
  if (length(verification_columns)) write_csv(
    adjustment_diagnostics %>% filter(grepl("^survival", analysis)) %>%
      select(analysis, run_label, outcome_adjustment, imputation, pm_draw, status, error,
        all_of(verification_columns)),
    file.path(adjustment_log_dir, "survival_numerical_verification.csv"))
  paired_fields <- intersect(c("analysis", "run_label", "imputation", "pm_draw", "status", "error",
    "fit_input_hash", "prediction_population_hash"), names(adjustment_diagnostics))
  cur <- adjustment_diagnostics %>% filter(outcome_adjustment == "current") %>% select(all_of(paired_fields))
  full <- adjustment_diagnostics %>% filter(outcome_adjustment == "full") %>% select(all_of(paired_fields))
  pairing <- full_join(cur, full, by = c("analysis", "run_label", "imputation", "pm_draw"), suffix = c("_current", "_full")) %>%
    mutate(fitting_inputs_match = !is.na(fit_input_hash_current) & !is.na(fit_input_hash_full) &
             fit_input_hash_current == fit_input_hash_full,
      prediction_populations_match = !is.na(prediction_population_hash_current) & !is.na(prediction_population_hash_full) &
             prediction_population_hash_current == prediction_population_hash_full,
      both_models_ok = !is.na(status_current) & !is.na(status_full) & status_current == "ok" & status_full == "ok")
  write_csv(pairing, file.path(adjustment_log_dir, "paired_input_checks.csv"))
  verification <- pairing %>% group_by(analysis, run_label) %>% summarise(
    n_cells = n(), all_inputs_match = all(fitting_inputs_match & prediction_populations_match),
    all_models_ok = all(both_models_ok), .groups = "drop")
  write_csv(verification, file.path(adjustment_log_dir, "paired_model_verification.csv"))
}
if (nrow(adjustment_effects)) {
  write_csv(adjustment_effects, file.path(adjustment_log_dir, "pooled_effects_long.csv"))
  keys <- c("analysis", "run_label", "cohort", "exposure", "horizon", "reference", "contrast", "measure")
  cur <- adjustment_effects %>% filter(outcome_adjustment == "current") %>% select(-outcome_adjustment)
  full <- adjustment_effects %>% filter(outcome_adjustment == "full") %>% select(-outcome_adjustment)
  comparison <- full_join(cur, full, by = keys, suffix = c("_current", "_full")) %>%
    mutate(estimate_change = estimate_full - estimate_current,
      se_ratio_full_over_current = se_full / se_current,
      ci_width_ratio_full_over_current = (ci_high_full - ci_low_full) / (ci_high_current - ci_low_current),
      se_scale = if_else(measure == "Risk ratio", "log risk ratio", "risk difference"))
  if (exists("verification")) comparison <- left_join(comparison, verification, by = c("analysis", "run_label"))
  write_csv(comparison, file.path(adjustment_log_dir, "outcome_adjustment_comparison.csv"))
}
if (nrow(adjustment_failures)) write_csv(adjustment_failures, file.path(adjustment_log_dir, "model_fit_failures.csv"))
message("Comparison logs and review tables: ", adjustment_log_dir)
if (any(adjustment_records$exit_status != 0) || nrow(adjustment_failures))
  warning("Some models failed. Inspect diagnostics and step_status.csv; failed cells were not dropped for pooling.")
if (exists("verification") && any(!verification$all_inputs_match))
  warning("Some paired inputs differ or could not be verified. Do not interpret those rows as an isolated outcome-model comparison.")
