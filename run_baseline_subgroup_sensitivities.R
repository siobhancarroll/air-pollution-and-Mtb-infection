# Refresh the baseline site and adult analyses with full outcome adjustment.
# Run in a fresh R session. No imputation, conversion models,
# or combinations of adult restriction with site stratification are run.
library(here)
library(dplyr)
library(readr)
subgroup_resuming <- exists("subgroup_resume_stamp")
subgroup_stamp <- if (subgroup_resuming) as.character(subgroup_resume_stamp) else
  paste0(Sys.Date(), "-", format(Sys.time(), "%H%M%S"))
stopifnot(length(subgroup_stamp) == 1L,
  grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}-[0-9]{6}$", subgroup_stamp))
subgroup_output_stamp <- paste0("subgroups-", subgroup_stamp)
subgroup_log <- here("logs", paste0("baseline-subgroups-", subgroup_stamp))
subgroup_cache <- here("Data", "baseline_subgroups", subgroup_stamp, "inputs")
dir.create(subgroup_log, recursive = TRUE, showWarnings = FALSE)
dir.create(subgroup_cache, recursive = TRUE, showWarnings = FALSE)
subgroup_inputs <- vapply(c("list_imputed_dataframes", "long_imputed_dataframes"), function(prefix) {
  files <- list.files(here("Data"), paste0("^", prefix, "_[0-9]{4}-[0-9]{2}-[0-9]{2}\\.Rdata$"), full.names = TRUE)
  if (!length(files)) stop("Missing completed input: ", prefix)
  tail(sort(files), 1)
}, character(1))
subgroup_manifest_path <- file.path(subgroup_log, "input_manifest.csv")
if (subgroup_resuming) {
  if (!file.exists(subgroup_manifest_path)) stop("Cannot resume without the original subgroup manifest")
  subgroup_manifest <- read_csv(subgroup_manifest_path,
    col_types = cols(.default = col_character()), show_col_types = FALSE)
  if (!setequal(unname(subgroup_inputs), subgroup_manifest$file))
    stop("Selected completed-data files differ from the original subgroup run; start a fresh run")
  subgroup_sources <- here(c("propensity_scores.R", "outcome_adjustment_helpers.R",
    "mi_pooling_helpers.R", "fuel_processing_helpers.R"))
  if (any(file.info(subgroup_sources)$mtime > file.info(subgroup_manifest_path)$mtime))
    stop("Analysis source files changed since this run started; start a fresh subgroup run")
} else {
  subgroup_hash <- unname(tools::md5sum(subgroup_inputs))
  if (anyNA(subgroup_hash)) stop("Could not checksum subgroup inputs")
  subgroup_manifest <- tibble(file = unname(subgroup_inputs), md5 = subgroup_hash)
  write_csv(subgroup_manifest, subgroup_manifest_path)
}
if (anyNA(subgroup_manifest$md5)) stop("Original subgroup manifest contains missing checksums")
verify_subgroup_inputs <- function() {
  audit_path <- file.path(subgroup_log, "input_checksum_checks.csv")
  for (attempt in 1:2) {
    read_warnings <- character()
    now <- withCallingHandlers(unname(tools::md5sum(subgroup_manifest$file)),
      warning = function(w) {
        read_warnings <<- c(read_warnings, conditionMessage(w)); invokeRestart("muffleWarning")
      })
    bad <- is.na(now) | now != subgroup_manifest$md5
    audit <- tibble(file = subgroup_manifest$file, expected_md5 = subgroup_manifest$md5,
      current_md5 = now, matches = !bad, attempt = attempt,
      checked_at = as.character(Sys.time()),
      read_warnings = paste(unique(read_warnings), collapse = " | "))
    write_csv(audit, audit_path, append = file.exists(audit_path))
    if (!any(bad)) return(invisible(TRUE))
    if (attempt == 1L) message("Input checksum check failed; retrying once. Details: ", audit_path)
  }
  report <- file.path(subgroup_log, "input_checksum_mismatches.csv")
  write_csv(audit[bad, ], report)
  stop("Subgroup input verification failed twice; see ", report,
    ". No next analysis step was started.")
}
subgroup_steps <- list(
  list(label = "adults", settings = 'analysis_populations <- c("adults", "female", "male");',
       prefix = "ps-outcome-full-"),
  list(label = "adults_clean", settings = paste(
       'analysis_populations <- c("adults", "female", "male"); hap_only <- TRUE;',
       'fuel_comparator <- "exclusive_clean"; hap_exposure_subset <- c("wood", "kerosene");'),
       prefix = "ps-cleanref-haponly-outcome-full-"),
  list(label = "site_NIRT", settings = paste(
       'analysis_populations <- "overall"; site_subset <- "NIRT";',
       'hap_exposure_subset <- "solid_fuel";'), prefix = "ps-outcome-full-site-NIRT-"),
  list(label = "site_BJGMC", settings = paste(
       'analysis_populations <- "overall"; site_subset <- "BJGMC";',
       'hap_exposure_subset <- "solid_fuel";'), prefix = "ps-outcome-full-site-BJGMC-")
)
# Check explicit filenames: directory listings and exit status alone do not
# establish that every model was run or every export was collected.
expected_subgroup_diagnostics <- function(label) {
  pops <- if (startsWith(label, "adults")) c("adults", "female", "male") else "overall"
  fuels <- if (label == "adults_clean") c("wood_use", "kerosene_use") else
    if (label == "adults") c("solid_fuel_use", "wood_use", "kerosene_use", "solid_fuel_primary") else "solid_fuel_use"
  rows <- list()
  for (pop in pops) for (assay in c("qgit", "tst")) {
    for (fuel in fuels) rows[[length(rows) + 1L]] <- tibble(
      run_label = paste(pop, fuel, assay, sep = "_"), expected_cells = 30L)
    if (label != "adults_clean") rows[[length(rows) + 1L]] <- tibble(
      run_label = paste(pop, assay, "recentavg", sep = "_"), expected_cells = 300L)
  }
  bind_rows(rows)
}
read_subgroup_diagnostic <- function(file, expected_cells) {
  for (attempt in 1:2) {
    value <- tryCatch(read_csv(file, col_types = cols(.default = col_character()),
      show_col_types = FALSE), error = function(e) NULL)
    if (!is.null(value) && nrow(value) == expected_cells &&
        all(c("imputation", "pm_draw", "status") %in% names(value))) {
      counts <- table(value$imputation)
      if (setequal(names(counts), as.character(1:30)) &&
          all(counts == expected_cells / 30L) &&
          !anyDuplicated(paste(value$imputation, value$pm_draw))) return(value)
    }
  }
  NULL
}
subgroup_outputs_complete <- function(folder, label) {
  expected <- expected_subgroup_diagnostics(label)
  all(vapply(seq_len(nrow(expected)), function(i) {
    !is.null(read_subgroup_diagnostic(file.path(folder,
      paste0(expected$run_label[i], "_outcome_model_diagnostics.csv")), expected$expected_cells[i]))
  }, logical(1)))
}
subgroup_records <- list()
subgroup_previous_records <- if (subgroup_resuming &&
    file.exists(file.path(subgroup_log, "step_status.csv")))
  read_csv(file.path(subgroup_log, "step_status.csv"), show_col_types = FALSE) else tibble()
for (step in subgroup_steps) {
  verify_subgroup_inputs()
  previous <- if (nrow(subgroup_previous_records)) subgroup_previous_records %>%
    filter(analysis == step$label, exit_status == 0) else tibble()
  if (nrow(previous) == 1L && subgroup_outputs_complete(previous$output_dir, step$label)) {
    message("Reusing completed subgroup step: ", step$label)
    subgroup_records[[length(subgroup_records) + 1L]] <- mutate(previous, outputs_complete = TRUE)
    next
  }
  code <- paste('outcome_adjustment_mode <- "full"; outcome_model_comparison <- TRUE;',
    paste0('outcome_comparison_cache <- ', encodeString(subgroup_cache, quote = '"'), ';'),
    paste0('output_stamp <- ', encodeString(subgroup_output_stamp, quote = '"'), ';'),
    step$settings)
  # Freeze the analysis body into a standalone script. This avoids relying
  # on a relative source() call returning the complete file in the child.
  analysis_path <- here("propensity_scores.R")
  analysis_lines <- character()
  for (attempt in 1:2) {
    analysis_lines <- readLines(analysis_path, warn = FALSE)
    if (length(analysis_lines) > 500L &&
        any(grepl("base_covariates_overall <-", analysis_lines, fixed = TRUE)) &&
        any(grepl("write_csv(hap_summary", analysis_lines, fixed = TRUE))) break
  }
  if (length(analysis_lines) <= 500L ||
      !any(grepl("write_csv(hap_summary", analysis_lines, fixed = TRUE)))
    stop("Could not read the complete baseline analysis source; no child analysis was launched: ", analysis_path)
  # Include dependencies as well: a prior child reached the analysis body
  # but its source() of mi_pooling_helpers.R did not define load_latest().
  helper_specs <- list(
    list(file = "fuel_processing_helpers.R", marker = "complete_fuel_imputation <- function"),
    list(file = "mi_pooling_helpers.R", marker = "load_latest <- function"),
    list(file = "outcome_adjustment_helpers.R", marker = "fit_adjustment_outcome <- function")
  )
  helper_lines <- character()
  for (helper in helper_specs) {
    lines <- character()
    for (attempt in 1:3) {
      lines <- readLines(here(helper$file), warn = FALSE)
      if (length(lines) > 100L && any(grepl(helper$marker, lines, fixed = TRUE)) &&
          !inherits(tryCatch(parse(text = lines), error = function(e) e), "error")) break
    }
    if (length(lines) <= 100L || !any(grepl(helper$marker, lines, fixed = TRUE)) ||
        inherits(tryCatch(parse(text = lines), error = function(e) e), "error"))
      stop("Could not capture complete helper source; no child analysis was launched: ", helper$file)
    helper_lines <- c(helper_lines, paste0("# Embedded dependency: ", helper$file), lines, "")
  }
  helper_source_calls <- c('source(here("mi_pooling_helpers.R"))',
    'source(here("outcome_adjustment_helpers.R"))')
  if (!all(vapply(helper_source_calls, function(call)
      sum(trimws(analysis_lines) == call) == 1L, logical(1))))
    stop("Baseline helper imports changed; review the standalone driver before continuing")
  analysis_lines <- analysis_lines[!trimws(analysis_lines) %in% helper_source_calls]
  driver <- file.path(subgroup_log, paste0(step$label, "_driver.R"))
  writeLines(c(
    paste0('setwd(', encodeString(here(), quote = '"'), ')'),
    paste0('.libPaths(c(', paste(encodeString(.libPaths(), quote = '"'), collapse = ", "), '))'),
    paste0('message("Starting standalone subgroup analysis: ', step$label, '")'),
    'message("R: ", R.version.string, "; working directory: ", getwd())',
    helper_lines,
    'stopifnot(exists("load_latest", mode = "function"), exists("fit_adjustment_outcome", mode = "function"), exists("complete_fuel_imputation", mode = "function"))',
    'message("Validated baseline helper functions are loaded")',
    code, analysis_lines,
    'if (!exists("fp") || !file.exists(file.path(fp, "hap_effects_summary.csv"))) stop("Analysis did not produce its expected output")'
  ), driver, useBytes = TRUE)
  logfile <- file.path(subgroup_log, paste0(step$label, ".log"))
  message("Running ", step$label, "; log: ", logfile)
  status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(driver)),
    stdout = logfile, stderr = logfile)
  output_dir <- here("tables_and_figures", paste0(step$prefix, subgroup_output_stamp))
  complete <- subgroup_outputs_complete(output_dir, step$label)
  if (status == 0 && !complete) {
    status <- 1L
    warning(step$label, ": process returned zero but expected diagnostics are missing or incomplete")
  }
  subgroup_records[[length(subgroup_records) + 1L]] <- tibble(analysis = step$label,
    exit_status = status, outputs_complete = complete, output_dir = output_dir,
    log_file = logfile)
  write_csv(bind_rows(subgroup_records), file.path(subgroup_log, "step_status.csv"))
}
verify_subgroup_inputs()
subgroup_records <- bind_rows(subgroup_records)
subgroup_diagnostics <- list(); subgroup_failures <- list(); subgroup_pm <- list(); subgroup_hap <- list(); subgroup_coverage <- list()
for (i in seq_len(nrow(subgroup_records))) {
  step <- subgroup_records[i, ]; folder <- step$output_dir
  annotate <- function(file) read_csv(file, col_types = cols(.default = col_character()),
    show_col_types = FALSE) %>% mutate(analysis = step$analysis, .before = 1)
  expected <- expected_subgroup_diagnostics(step$analysis)
  for (j in seq_len(nrow(expected))) {
    file <- file.path(folder, paste0(expected$run_label[j], "_outcome_model_diagnostics.csv"))
    value <- read_subgroup_diagnostic(file, expected$expected_cells[j])
    subgroup_coverage[[length(subgroup_coverage) + 1L]] <- tibble(analysis = step$analysis,
      run_label = expected$run_label[j], expected_cells = expected$expected_cells[j],
      complete = !is.null(value))
    if (!is.null(value)) subgroup_diagnostics[[length(subgroup_diagnostics) + 1L]] <- value %>%
      mutate(analysis = step$analysis, run_label = expected$run_label[j], .before = 1) else
      subgroup_failures[[length(subgroup_failures) + 1L]] <- tibble(analysis = step$analysis,
        run_label = expected$run_label[j], error = "Expected model diagnostics missing or incomplete")
  }
  file <- file.path(folder, "model_fit_failures.csv")
  if (file.exists(file)) subgroup_failures[[length(subgroup_failures) + 1L]] <- annotate(file)
  file <- file.path(folder, "pm25_contrasts_summary.csv")
  if (file.exists(file) && file.info(file)$size > 1L)
    subgroup_pm[[length(subgroup_pm) + 1L]] <- annotate(file)
  file <- file.path(folder, "hap_effects_summary.csv")
  if (file.exists(file) && file.info(file)$size > 1L) {
    values <- annotate(file)
    # Adult wood/kerosene rows use their clean-reference exports, matching
    # the manuscript's overall-population comparisons.
    if (step$analysis == "adults") values <- values %>% filter(!exposure %in% c("wood_use", "kerosene_use"))
    subgroup_hap[[length(subgroup_hap) + 1L]] <- values
  }
}
subgroup_diagnostics <- bind_rows(subgroup_diagnostics)
subgroup_failures <- bind_rows(subgroup_failures)
write_csv(bind_rows(subgroup_coverage), file.path(subgroup_log, "output_coverage.csv"))
write_csv(subgroup_diagnostics, file.path(subgroup_log, "outcome_model_diagnostics.csv"))
write_csv(bind_rows(subgroup_pm), file.path(subgroup_log, "pm25_subgroup_effects.csv"))
write_csv(bind_rows(subgroup_hap), file.path(subgroup_log, "hap_subgroup_effects.csv"))
if (nrow(subgroup_failures)) {
  write_csv(subgroup_failures, file.path(subgroup_log, "model_fit_failures.csv"))
} else {
  # A successful resumed run must not retain an earlier failure report.
  write_csv(tibble(analysis = character(), run_label = character(), error = character()),
    file.path(subgroup_log, "model_fit_failures.csv"))
}
if (nrow(subgroup_diagnostics)) write_csv(subgroup_diagnostics %>%
    group_by(analysis, run_label) %>% summarise(n_cells = n(),
      all_models_ok = all(status == "ok"), .groups = "drop"),
    file.path(subgroup_log, "model_verification.csv"))
writeLines(c(
  "Full outcome adjustment; adult restriction is age >=15 years, matching the previous subgroup definition.",
  "The adult models use the overall baseline weighting covariate set; site models omit Site because it is constant.",
  "Adult female and male models use the same full baseline set except B_Gender, which is constant within each stratum.",
  "Weights and standardization populations are re-estimated within each subgroup.",
  "Site analyses include PM2.5 and the main any-polluting-fuel definition only; adult-by-site and stricter-conversion combinations are not run.",
  "Interpret site PM2.5 contrasts only when both reference_in_range and contrast_in_range are TRUE; comparisons between sites require common support.",
  "Wide confidence intervals and low ESS limit interpretation; failed or unidentified fits must not be reported as valid estimates with a caveat.",
  "Do not infer subgroup differences by comparing statistical significance; these are exploratory subgroup estimates."
), file.path(subgroup_log, "interpretation_notes.txt"))
if (any(subgroup_records$exit_status != 0) || nrow(subgroup_failures) ||
    !nrow(subgroup_diagnostics) || any(subgroup_diagnostics$status != "ok"))
  warning("Some subgroup steps have incomplete outputs or failed models; inspect step_status.csv, output_coverage.csv and the analysis logs before reporting")
message("Subgroup outputs and review tables: ", subgroup_log)
