# Audit cached source responses only. No R models or imputation are fitted.
library(here)
library(dplyr)
library(readr)
source(here("tests", "test_fuel_processing.R"))
fuel_audit_input <- new.env(parent = emptyenv())
load(here("Data", "data.Rdata"), envir = fuel_audit_input)
fuel_audit_data <- process_fuel_questionnaire(fuel_audit_input$complete)
fuel_audit_dir <- here("tables_and_figures", paste0("fuel-questionnaire-audit-", Sys.Date()))
write_fuel_questionnaire_audit(fuel_audit_data, fuel_audit_dir)
message("Aggregate audit: ", fuel_audit_dir)
message("Identifiable flagged records: Data/fuel_questionnaire_audit/")
entry_flags <- gsub(" ", "", as.character(fuel_audit_data$VISIT)) == "Entry" &
  fuel_audit_data$fuel_pre_imputation_conflict
if (any(entry_flags, na.rm = TRUE))
  warning("Within-person Entry none/positive contradictions need review before imputation. Differences between household members are retained and do not block the run.")
