### Household fuel classification at Entry versus Month 4 #####################
# How often a contact's fuel classification differed between the Entry and
# Month 4 questionnaires, among contacts with both classifications observed.
# Uses the same corrected fuel processing as multiple_imputation.R
# (process_fuel_questionnaire()), and the first record of each visit, as in
# that script. Writes aggregate counts only:
#   tables_and_figures/fuel-switching-<date>/fuel_switching_summary.csv
library(tidyverse)
library(here)

load(here("Data", "data.Rdata"))
source(here("fuel_processing_helpers.R"))
complete <- process_fuel_questionnaire(complete)
complete$VISIT <- gsub(" ", "", as.character(complete$VISIT))

indicators <- c("Any polluting fuel use" = "exclusive_clean", "Wood use" = "wood_use",
                "Kerosene use" = "kerosene_use", "Polluting primary cooking fuel" = "solid_fuel_primary")

wide <- complete %>%
  filter(VISIT %in% c("Entry", "Month4")) %>%
  mutate(VISIT = factor(VISIT, levels = c("Entry", "Month4"))) %>%
  arrange(PID, VISIT) %>%
  distinct(PID, VISIT, .keep_all = TRUE) %>%
  dplyr::select(PID, RELATED.INDEX.PID, VISIT, all_of(unname(indicators))) %>%
  mutate(across(all_of(unname(indicators)), as.character)) %>%
  pivot_wider(id_cols = c(PID, RELATED.INDEX.PID), names_from = VISIT,
              values_from = all_of(unname(indicators)))

summary_tbl <- imap_dfr(indicators, function(var, label) {
  e <- wide[[paste0(var, "_Entry")]]; m <- wide[[paste0(var, "_Month4")]]
  # exclusive_clean is coded YES for clean-only households; flip it so that
  # "exposed" means polluting fuel use for every indicator
  exposed <- function(x) if (var == "exclusive_clean") x == "No" else x == "YES"
  both <- !is.na(e) & !is.na(m)
  changed <- both & e != m
  tibble(
    indicator = label,
    contacts_entry_observed = sum(!is.na(e)),
    contacts_both_observed = sum(both),
    contacts_changed = sum(changed),
    pct_changed = round(100 * mean(changed[both]), 1),
    unexposed_to_exposed = sum(both & !exposed(e) & exposed(m)),
    exposed_to_unexposed = sum(both & exposed(e) & !exposed(m)),
    households_both_observed = n_distinct(wide$RELATED.INDEX.PID[both]),
    households_any_change = n_distinct(wide$RELATED.INDEX.PID[changed])
  )
})
print(summary_tbl, width = Inf)

out_dir <- here("tables_and_figures",
                paste0("fuel-switching-", get0("output_stamp", ifnotfound = as.character(Sys.Date()))))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write_csv(summary_tbl, file.path(out_dir, "fuel_switching_summary.csv"))
message("Saved ", file.path(out_dir, "fuel_switching_summary.csv"))
