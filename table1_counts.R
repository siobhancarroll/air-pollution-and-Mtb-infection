### Table 1 summaries: baseline characteristics by fuel use and PM2.5 group ###
# Summarises observed baseline values for Table 1a (by household fuel use) and
# Table 1b (by ambient PM2.5 group: fixed cut-points <32, 32-52, 52-72,
# >=72 ug/m3, and quartiles of the observed values in each source), from two
# sources:
#   raw      - Entry-visit rows of data.Rdata (Month 4 row for contact hours),
#              for all enrolled contacts and for the analysis dataset only
#   analysis - the pre-imputation (observed) values in the long imputed
#              dataset (.imp == 0), i.e. after the cleaning steps and
#              exclusions in multiple_imputation.R
# Categorical variables are summarised at their raw levels (n per level,
# missing as NA), continuous variables as n, mean and SD.
#
# Output (aggregate only): tables_and_figures/table1-<date>/
#   table1_categorical.csv, table1_continuous.csv

library(tidyverse)
library(here)
source(here("mi_pooling_helpers.R"))   # load_latest()
source(here("fuel_processing_helpers.R"))

out_dir <- file.path(here("tables_and_figures"), paste0("table1-", Sys.Date()))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

load(here("Data", "data.Rdata"))   # `complete`
# Apply the same questionnaire correction used before imputation. Cached
# data.Rdata may contain the older dropdown-default fuel classifications.
complete <- process_fuel_questionnaire(complete)
load_latest("^list_imputed_dataframes_\\d{4}-\\d{2}-\\d{2}\\.Rdata$")
load_latest("^long_imputed_dataframes_\\d{4}-\\d{2}-\\d{2}\\.Rdata$")
analysis_ids <- unique(as.character(allImputations$imputations[[1]]$PID))

cat_vars <- c("B_Gender", "Site", "B_HH_Income.Categories", "B_Education", "B_Employment",
              "Living.Area", "Type.of.House", "alcohol_disorder", "B_BCG", "B_Sleepindex")
cont_vars <- c("B_Age", "B_Height", "B_Weight", "MUAC", "Household.size",
               "contact_hours_per_week_Month4", "B_HbA1c", "mean_pm_recentavg")

### Raw Entry-visit data ########################################################
complete$VISIT <- str_replace_all(complete$VISIT, " ", "")
raw <- complete %>% filter(VISIT == "Entry") %>% distinct(PID, .keep_all = TRUE) %>%
  mutate(PID = as.character(PID), MUAC = Mid.upper.arm.circumference..cm.,
         fuel = exclusive_clean) %>%
  left_join(complete %>% filter(VISIT == "Month4") %>% distinct(PID, .keep_all = TRUE) %>%
              transmute(PID = as.character(PID), contact_hours_per_week_Month4 = contact_hours_per_week),
            by = "PID")

### Observed values in the analysis dataset (.imp == 0) ##########################
obs <- impdata %>% filter(.imp == 0) %>%
  mutate(PID = as.character(PID),
         fuel = if ("exclusive_clean_Entry" %in% names(.)) exclusive_clean_Entry else exclusive_clean,
         mean_pm_recentavg = if ("mean_pm_recentavg" %in% names(.)) mean_pm_recentavg else mean_pm_recentavg_Entry)

missing_vars <- list(raw = setdiff(c(cat_vars, cont_vars), names(raw)),
                     analysis = setdiff(c(cat_vars, cont_vars), names(obs)))
print(missing_vars)

### Group definitions ##########################################################
add_groups <- function(df) {
  df %>% mutate(
    fuel_group = case_when(toupper(as.character(fuel)) == "YES" ~ "Exclusive clean fuel",
                           toupper(as.character(fuel)) == "NO" ~ "Any polluting fuel use",
                           TRUE ~ "Missing fuel use information"),
    pm_group = case_when(is.na(mean_pm_recentavg) ~ NA_character_,
                         mean_pm_recentavg < 32 ~ "Q1 <32",
                         mean_pm_recentavg < 52 ~ "Q2 32-52",
                         mean_pm_recentavg < 72 ~ "Q3 52-72",
                         TRUE ~ "Q4 >=72"),
    pm_quartile = {
      qs <- quantile(mean_pm_recentavg, c(0.25, 0.5, 0.75), na.rm = TRUE)
      lab <- sprintf(c("Q1 <%.1f", "Q2 %.1f-%.1f", "Q3 %.1f-%.1f", "Q4 >=%.1f"),
                     c(qs[1], qs[1], qs[2], qs[3]), c(NA, qs[2], qs[3], NA))
      lab <- sub("-NA", "", lab)
      as.character(cut(mean_pm_recentavg, c(-Inf, qs, Inf), labels = lab, right = FALSE))
    })
}
sources <- list(
  raw_all      = add_groups(raw),
  raw_analysis = add_groups(raw %>% filter(PID %in% analysis_ids)),
  analysis     = add_groups(obs)
)

summarise_source <- function(df, src) {
  long <- bind_rows(mutate(df, table = "1a", group = "Overall"),
                    mutate(df, table = "1a", group = fuel_group),
                    mutate(df, table = "1b", group = "Overall"),
                    df %>% filter(!is.na(pm_group)) %>% mutate(table = "1b", group = pm_group),
                    df %>% filter(!is.na(pm_quartile)) %>% mutate(table = "1b_quartile", group = pm_quartile))
  cats <- map_dfr(intersect(c(cat_vars, "fuel_group"), names(df)), function(v)
    long %>% mutate(level = as.character(.data[[v]])) %>% count(table, group, level) %>%
      mutate(variable = v, .before = 1))
  conts <- map_dfr(intersect(cont_vars, names(df)), function(v)
    long %>% group_by(table, group) %>%
      summarise(n = sum(!is.na(.data[[v]])), mean = mean(.data[[v]], na.rm = TRUE),
                sd = sd(.data[[v]], na.rm = TRUE), .groups = "drop") %>%
      mutate(variable = v, .before = 1))
  ns <- long %>% count(table, group, name = "N")
  list(cat = mutate(cats, source = src, .before = 1),
       cont = mutate(conts, source = src, .before = 1),
       n = mutate(ns, source = src, .before = 1))
}
res <- imap(sources, summarise_source)
write_csv(bind_rows(map(res, "cat")), file.path(out_dir, "table1_categorical.csv"))
write_csv(bind_rows(map(res, "cont")), file.path(out_dir, "table1_continuous.csv"))
write_csv(bind_rows(map(res, "n")), file.path(out_dir, "table1_group_sizes.csv"))
print(bind_rows(map(res, "n")), n = Inf)
message("Saved to ", out_dir)
