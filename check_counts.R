### Check analysis counts ###################################################
# 1. Contacts enrolled vs contacts in the analysis dataset, and which enrolled
#    contacts are not in the analysis dataset.
# 2. Contacts and conversions in each interval of the PM2.5 conversion models,
#    after dropping person-intervals without PM2.5 data (before December 2014).
#
# Participant IDs are printed to the console only; saved files hold counts.

library(tidyverse)
library(here)
library(mitools)
source(here("mi_pooling_helpers.R"))   # load_latest()

out_dir <- file.path(here("tables_and_figures"), paste0("count-checks-", get0("output_stamp", ifnotfound = as.character(Sys.Date()))))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

load(here("Data", "data.Rdata"))   # `complete`
load_latest("^list_imputed_dataframes_\\d{4}-\\d{2}-\\d{2}\\.Rdata$")
load_latest("^survival_person_period_qgit_\\d{4}-\\d{2}-\\d{2}\\.Rdata$")
load_latest("^survival_person_period_tst_\\d{4}-\\d{2}-\\d{2}\\.Rdata$")

complete$VISIT <- str_replace_all(complete$VISIT, " ", "")
all_ids   <- unique(as.character(complete$PID))
entry_ids <- unique(as.character(complete$PID[complete$VISIT == "Entry"]))
analysis_ids <- unique(as.character(allImputations$imputations[[1]]$PID))

entry <- complete %>% filter(VISIT == "Entry") %>% distinct(PID, .keep_all = TRUE) %>%
  mutate(PID = as.character(PID))
not_in_analysis <- entry %>% filter(!PID %in% analysis_ids)
no_entry <- setdiff(all_ids, entry_ids)

cat("\n== Enrolment and analysis dataset ==\n")
cat("Contacts in data.Rdata (any visit):", length(all_ids), "\n")
cat("Contacts with an Entry visit:      ", length(entry_ids), "\n")
cat("Contacts in the analysis dataset:  ", length(analysis_ids), "\n")
cat("Entry contacts NOT in the analysis dataset:", nrow(not_in_analysis), "\n")
print(not_in_analysis %>% count(baseline_IGRA = QGIT.Result, TST_recorded = !is.na(TST)))
cat("PIDs (for you to look up):", paste(not_in_analysis$PID, collapse = ", "), "\n")
if (length(no_entry)) cat("Contacts with no Entry visit row:", length(no_entry), "-", paste(no_entry, collapse = ", "), "\n")
cat("Visits recorded for contacts with no Entry row:\n")
print(complete %>% filter(PID %in% no_entry) %>% count(VISIT))

### PM2.5 conversion-model counts ##############################################
pm_counts <- function(pp, cohort) {
  pp <- pp %>% filter(interval <= 2)
  first_na <- pp %>% filter(is.na(mean_pm_intervalavg) | is.na(SE_pm_intervalavg)) %>%
    group_by(PID) %>% summarise(first_na = min(interval), .groups = "drop")
  kept <- pp %>% left_join(first_na, by = "PID") %>% filter(is.na(first_na) | interval < first_na)
  pp %>% group_by(interval) %>%
    summarise(at_risk = n_distinct(PID), conversions = sum(event == 1, na.rm = TRUE), .groups = "drop") %>%
    left_join(kept %>% group_by(interval) %>%
                summarise(in_pm_model = n_distinct(PID), conversions_in_pm_model = sum(event == 1, na.rm = TRUE),
                          .groups = "drop"), by = "interval") %>%
    mutate(cohort = cohort, excluded_no_pm = at_risk - in_pm_model, .before = 1)
}
counts <- bind_rows(pm_counts(allImputations_survival_qgit$imputations[[1]], "IGRA"),
                    pm_counts(allImputations_survival_tst$imputations[[1]], "TST"))
cat("\n== PM2.5 conversion models: contacts per interval ==\n")
print(counts)

write_csv(tibble(item = c("contacts_any_visit", "contacts_with_entry", "analysis_dataset",
                          "entry_not_in_analysis", "no_entry_row"),
                 n = c(length(all_ids), length(entry_ids), length(analysis_ids),
                       nrow(not_in_analysis), length(no_entry))),
          file.path(out_dir, "enrolment_counts.csv"))
write_csv(counts, file.path(out_dir, "pm_model_counts.csv"))
message("Saved to ", out_dir)
