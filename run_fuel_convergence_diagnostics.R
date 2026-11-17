# Run in a fresh R session. Reuses completed imputations and cached weights.
# Diagnostic refits are exported separately and never substituted for pooling.
local({
  library(here)
  previous <- "2026-10-26-220919"
  manifest <- readr::read_csv(here("logs", paste0("outcome-adjustment-comparison-", previous),
    "input_manifest.csv"), col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE)
  if (any(!file.exists(manifest$file)) ||
      any(unname(tools::md5sum(manifest$file)) != manifest$md5))
    stop("Original comparison inputs are unavailable or changed")
  stamp <- paste0("convergence-", Sys.Date(), "-", format(Sys.time(), "%H%M%S"))
  directory <- here("logs", stamp)
  dir.create(directory, recursive = TRUE)
  old <- Sys.getenv("FUEL_CONVERGENCE_DIAGNOSTICS", unset = NA_character_)
  on.exit(if (is.na(old)) Sys.unsetenv("FUEL_CONVERGENCE_DIAGNOSTICS") else
    Sys.setenv(FUEL_CONVERGENCE_DIAGNOSTICS = old), add = TRUE)
  Sys.setenv(FUEL_CONVERGENCE_DIAGNOSTICS = directory)
  # The analysis and its sourced helpers expect settings (including fp) in
  # the global environment. Use a fresh child process, as the comparison
  # runner does, rather than sourcing into this local wrapper environment.
  quote_r <- function(x) encodeString(x, quote = '"')
  code <- paste(
    'outcome_model_comparison <- TRUE; outcome_adjustment_mode <- "full";',
    paste0('outcome_comparison_cache <- ', quote_r(here("Data", "outcome_adjustment_comparison", previous, "inputs")), ';'),
    paste0('output_stamp <- ', quote_r(stamp), ';'),
    'conversion_definition <- "primary"; fuel_comparator <- "exclusive_clean"; hap_only <- TRUE;',
    # Preserve the original exposure order and RNG stream for matched inputs.
    'hap_exposure_subset <- c("wood", "kerosene");',
    paste0('source(', quote_r(here("discrete_time_survival.R")), ')'))
  logfile <- file.path(directory, "analysis.log")
  message("Running convergence exploration; log: ", logfile)
  status <- system2(file.path(R.home("bin"), "Rscript"),
    c("-e", shQuote(code)), stdout = logfile, stderr = logfile)
  if (status != 0L) stop("Convergence exploration failed; inspect ", logfile)
  original <- readr::read_csv(here("logs", paste0("outcome-adjustment-comparison-", previous),
    "outcome_model_diagnostics.csv"),
    col_types = readr::cols(.default = readr::col_character()), show_col_types = FALSE)
  readr::write_csv(dplyr::filter(original, analysis == "survival_primary_clean",
    outcome_adjustment == "full", run_label == "tst_wood_use"),
    file.path(directory, "original_tst_wood_imputation_lookup.csv"))
  message("Convergence diagnostics saved to: ", directory)
})
