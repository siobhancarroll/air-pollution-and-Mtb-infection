### Run the full Manuscript 2 pipeline ######################################
# Runs each step in a separate Rscript process so that settings do not carry
# over between steps. Each step's console output is written to
# logs/<date>/<step>.log, and the run stops at the first failed step.
#
# Usage (from the project folder): source("run_all.R")
# run_imputation <- FALSE skips multiple_imputation.R and uses the most recent
# (start_from <- "<step name>" resumes from that step.)
# imputed datasets in Data/.
#
# Also saves sessionInfo() and package versions, and compares the new summary
# tables with those from a reference run (compare_to).

library(here)
if (!exists("run_imputation")) run_imputation <- TRUE
# diagnose_outcome_fits: TRUE writes aggregate coefficient and cell-count
# diagnostics for each survival outcome model to <output>/fit_diagnostics/.
if (!exists("diagnose_outcome_fits")) diagnose_outcome_fits <- FALSE
# compare_to: date stamp of the reference run. Default: the most recent
# earlier primary run (newest survival-YYYY-MM-DD folder before today).
if (!exists("compare_to")) {
  prev <- sub("^survival-", "", list.files(here("tables_and_figures"), pattern = "^survival-\\d{4}-\\d{2}-\\d{2}$"))
  prev <- sort(prev[prev < get0("output_stamp", ifnotfound = as.character(Sys.Date()))], decreasing = TRUE)
  compare_to <- if (length(prev)) prev[1] else NA_character_
}

log_dir <- here("logs", as.character(Sys.Date()))
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
rscript <- file.path(R.home("bin"), "Rscript")

steps <- list(
  list(name = "01_multiple_imputation",      code = 'source("multiple_imputation.R")'),
  list(name = "02_survival_primary",         code = 'source("discrete_time_survival.R")'),
  # Wood and kerosene versus exclusively clean fuel use (the comparisons
  # reported in the manuscript); output folders *-cleanref-haponly-<date>
  list(name = "02b_survival_clean_reference", code = 'fuel_comparator <- "exclusive_clean"; hap_only <- TRUE; hap_exposure_subset <- c("wood", "kerosene"); source("discrete_time_survival.R")'),
  list(name = "03_survival_sensitivity",     code = 'conversion_definition <- "sensitivity"; source("discrete_time_survival.R")'),
  # PM2.5 exposure-timing sensitivity analyses (pm_windows.R); the
  # any-polluting-fuel models are repeated only because the script requires a
  # fuel exposure
  list(name = "03c_window_pre_entry3",     code = 'pm_window <- "pre_entry3"; hap_exposure_subset <- "solid_fuel"; source("discrete_time_survival.R")', optional = TRUE),
  list(name = "03e_window_early2",         code = 'pm_window <- "early2"; hap_exposure_subset <- "solid_fuel"; source("discrete_time_survival.R")', optional = TRUE),
  list(name = "03f_window_first_interval", code = 'pm_window <- "first_interval"; hap_exposure_subset <- "solid_fuel"; source("discrete_time_survival.R")', optional = TRUE),
  list(name = "04_baseline",                 code = 'source("propensity_scores.R")'),
  list(name = "04b_baseline_clean_reference", code = 'fuel_comparator <- "exclusive_clean"; hap_only <- TRUE; hap_exposure_subset <- c("wood", "kerosene"); analysis_populations <- "overall"; source("propensity_scores.R")'),
  list(name = "05_baseline_site_NIRT",       code = 'site_subset <- "NIRT"; source("propensity_scores.R")'),
  list(name = "06_baseline_site_BJGMC",      code = 'site_subset <- "BJGMC"; source("propensity_scores.R")'),
  list(name = "07_compare_definitions",      code = 'source("compare_conversion_definitions.R")'),
  list(name = "08_compare_site_strata",      code = 'source("compare_site_strata.R")'),
  list(name = "09_compare_assay_cohorts",    code = 'source("compare_assay_cohorts.R")'),
  list(name = "10_check_counts",             code = 'source("check_counts.R")'),
  # Reads the newest primary survival-<date> and ps-<date> folders, i.e. those
  # written by steps 02 and 04 of this run
  list(name = "11_absolute_risks",           code = 'source("extract_absolute_risks.R")'),
  list(name = "12_fuel_ambient_pm",          code = 'source("scripts/audit_fuel_ambient_relationship.R")'),
  # Entry versus Month 4 fuel classification (descriptive; Methods)
  list(name = "12b_fuel_switching",          code = 'source("scripts/fuel_switching_summary.R")'),
  list(name = "13_compare_pm_windows",       code = 'source("compare_pm_windows.R")')
)
if (!run_imputation) steps <- steps[-1]
# start_from: name of a step to resume from (e.g. after a failure), skipping
# the steps before it
if (exists("start_from")) {
  k <- match(start_from, vapply(steps, `[[`, character(1), "name"))
  if (is.na(k)) stop("start_from is not a step name: ", start_from)
  steps <- steps[k:length(steps)]
}

# One date stamp for every output folder of this run, even if it runs past
# midnight
if (!exists("output_stamp")) output_stamp <- as.character(Sys.Date())
step_settings <- sprintf('diagnose_outcome_fits <- %s; output_stamp <- "%s"; ',
                         if (isTRUE(diagnose_outcome_fits)) "TRUE" else "FALSE", output_stamp)
for (s in steps) {
  s$code <- paste0(step_settings, s$code)
  log_file <- file.path(log_dir, paste0(s$name, ".log"))
  message(format(Sys.time(), "%H:%M"), "  running ", s$name, " ...")
  t0 <- Sys.time()
  status <- system2(rscript, c("-e", shQuote(s$code)), stdout = log_file, stderr = log_file)
  message("      finished in ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min (exit status ", status, ")")
  # Exploratory sensitivity steps (optional = TRUE) are reported and skipped
  # on failure; every other step stops the run.
  if (status != 0 && isTRUE(s$optional)) {
    warning(s$name, " failed (optional sensitivity step; continuing) - see ", log_file, call. = FALSE)
    failed_optional <- c(get0("failed_optional", ifnotfound = character(0)), s$name)
  } else if (status != 0) stop(s$name, " failed - see ", log_file)
}

### Package versions ##########################################################
pkgs <- c("mice", "WeightIt", "cobalt", "survey", "marginaleffects", "mitools",
          "splines", "future.apply", "zscorer", "tidyverse")
versions <- data.frame(package = pkgs,
                       version = sapply(pkgs, function(p) tryCatch(as.character(packageVersion(p)), error = function(e) NA)))
write.csv(versions, file.path(log_dir, "package_versions.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo.txt"))
print(versions, row.names = FALSE)

### Compare with the reference run #############################################
tf <- here("tables_and_figures")
today <- output_stamp
files <- c("survival-%s/survival_pm25_contrasts_summary.csv", "survival-%s/survival_hap_effects_summary.csv",
           "survival-sensitivity-%s/survival_pm25_contrasts_summary.csv",
           "ps-%s/pm25_contrasts_summary.csv", "ps-%s/hap_effects_summary.csv",
           "site-strata-%s/pm25_by_site.csv",
           "absolute-risks-%s/longitudinal_hap_absolute_risks.csv",
           "absolute-risks-%s/longitudinal_pm25_absolute_risks.csv",
           "absolute-risks-%s/baseline_hap_absolute_risks.csv",
           "absolute-risks-%s/baseline_pm25_absolute_risks.csv",
           "fuel-ambient-%s/fuel_ambient_pm_differences.csv")
# Reference files that are stored only in iCloud (not downloaded) can read as
# empty; these are reported and skipped rather than stopping the run.
read_ref <- function(path) tryCatch(read.csv(path), error = function(e) NULL)
if (is.na(compare_to)) message("No earlier run to compare with.")
if (!is.na(compare_to) && today != compare_to) {
  message("Comparing with the run of ", compare_to)
  for (f in files) {
    new <- file.path(tf, sprintf(f, today)); old <- file.path(tf, sprintf(f, compare_to))
    if (!file.exists(new) || !file.exists(old)) { message("Not compared (missing): ", sprintf(f, "<date>")); next }
    a <- read_ref(new); b <- read_ref(old)
    if (is.null(a) || is.null(b)) {
      message("Not compared (could not be read; if stored only in iCloud, download it in Finder): ",
              if (is.null(a)) new else old)
      next
    }
    num <- intersect(names(a)[sapply(a, is.numeric)], names(b)[sapply(b, is.numeric)])
    same_shape <- nrow(a) == nrow(b)
    diff <- if (same_shape) max(abs(as.matrix(a[num]) - as.matrix(b[num])), na.rm = TRUE) else NA
    message(sprintf("%-55s rows %s/%s, largest numeric difference: %s",
                    sprintf(f, "<date>"), nrow(a), nrow(b), if (is.na(diff)) "n/a" else signif(diff, 3)))
  }
}
if (length(get0("failed_optional", ifnotfound = character(0))))
  message("Optional sensitivity steps that failed: ", paste(failed_optional, collapse = ", "))
message("Done. Logs and versions in ", log_dir)
