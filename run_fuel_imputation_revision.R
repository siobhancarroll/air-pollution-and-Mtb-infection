# Rebuild corrected imputations and the full pipeline, then the clean-reference
# trial. Start with source("audit_fuel_questionnaire.R") and review its outputs.
# No full rerun of data_processing.R is needed: cached raw questionnaire fields
# in Data/data.Rdata are reprocessed inside multiple_imputation.R.
library(here)
source(here("audit_fuel_questionnaire.R"))
if (any(entry_flags, na.rm = TRUE))
  stop("Resolve within-person Entry none/positive contradictions before rerunning imputation; differing household-member reports are permitted")
revision_stamp <- format(Sys.time(), "%Y-%m-%d-%H%M%S")
revision_backup <- here("Data", paste0("before-fuel-revision-", revision_stamp))
dir.create(revision_backup, recursive = TRUE, showWarnings = FALSE)
# Preserve same-day objects and output folders which ordinary scripts replace.
revision_objects <- list.files(here("Data"), full.names = TRUE)
revision_objects <- revision_objects[grepl(paste0("_", Sys.Date(), "\\.(Rdata|rds)$"), basename(revision_objects))]
if (length(revision_objects) && !all(file.copy(revision_objects, revision_backup)))
  stop("Could not preserve existing same-day data objects")
revision_dirs <- c(paste0("imputation_", Sys.Date()),
  paste0(c("survival-", "survival-sensitivity-", "ps-", "ps-site-NIRT-", "ps-site-BJGMC-"), Sys.Date()))
for (name in revision_dirs) {
  path <- here("tables_and_figures", name)
  if (dir.exists(path) && !file.copy(path, revision_backup, recursive = TRUE))
    stop("Could not preserve ", path)
}
if (dir.exists(here("logs", as.character(Sys.Date()))) &&
    !file.copy(here("logs", as.character(Sys.Date())), revision_backup, recursive = TRUE))
  stop("Could not preserve existing pipeline logs")
writeLines(c(paste("Backup:", revision_backup),
  "Fuel rule: all named secondary fuels No and secondary_fuel_none No => secondary answers missing; primary answers retained",
  "Fuel indicators are derived from completed source responses, not imputed separately"),
  file.path(revision_backup, "revision_notes.txt"))
run_imputation <- TRUE
source(here("run_all.R"))
source(here("run_fuel_clean_reference.R"))
