
### Load packages #############################################################
library(tidyverse)
library(mice)
library(here)
library(knitr)
library(flextable)
library(lubridate)
library(corrr)
library(mitools)
library(naniar)
library(gtsummary)
library(reformulas)

### Load data (created by data_processing.R) ###############################
load(here("Data", "data.Rdata"))
source(here("fuel_processing_helpers.R"))

### Outcome variables as imputation predictors ##############################
# TRUE  - multiple imputation then deletion (MID; von Hippel 2007): the
#         continuous outcome measurements (TST induration in mm and QGIT
#         antigen minus nil at Entry and Months 4, 12 and 24) are used as
#         predictors when imputing incomplete baseline covariates and household
#         fuel exposures (Moons et al. 2006; White, Royston & Wood 2011).
#         Imputed outcome values are deleted afterwards; events are defined
#         from observed values only (see build_person_period()).
# FALSE - baseline covariates and exposures are predicted by other baseline
#         covariates only, and follow-up antigen minus nil is not included.
# Output files are named by date only, so two runs on the same day overwrite
# each other.
mid_outcome_predictors <- TRUE

### Set up folder to store plots ##############################################
main_dir <- here("tables_and_figures")  # main output directory
sub_dir <- paste0("imputation_", Sys.Date())  # sub-directory named with the current date
fp <- file.path(main_dir, sub_dir)  # output path for tables and figures
dir.create(fp, recursive = TRUE)  # create output folder

# Reprocess preserved source fields, including cached data produced by the old
# processing script. No need to repeat PM2.5/calendar construction for this fix.
complete <- process_fuel_questionnaire(complete)
write_fuel_questionnaire_audit(complete, file.path(fp, "fuel_questionnaire_audit"))
entry_fuel_conflict <- gsub(" ", "", as.character(complete$VISIT)) == "Entry" &
  complete$fuel_pre_imputation_conflict
if (any(entry_fuel_conflict, na.rm = TRUE))
  stop("Within-person Entry responses report both no secondary fuel and a positive secondary fuel; review Data/fuel_questionnaire_audit/flagged_records_*.csv. Differences between household members do not block imputation.")


### Set up the data for imputation ############################################
# Remove spaces in the variable (otherwise, causes problems with MICE)
complete$VISIT <- str_replace_all(complete$VISIT, " ", "")

# Remove all spaces and characters from the Living Area variable
complete$Living.Area <- str_replace_all(complete$Living.Area, " ", "")
complete$Living.Area <- str_replace_all(complete$Living.Area, "/", "")
complete$Living.Area <- str_replace_all(complete$Living.Area, "-", "")

# Remove spaces from sleep index 
complete$B_Sleepindex <- str_replace_all(complete$B_Sleepindex, " ", "")
complete$B_Sleepindex <- str_replace_all(complete$B_Sleepindex, "/", "")

# Remove spaces and characters from education
complete$B_Education <- str_replace_all(complete$B_Education, " ", "")
complete$B_Education <- str_replace_all(complete$B_Education, "\\(", "")
complete$B_Education <- str_replace_all(complete$B_Education, "\\)", "")
complete$B_Education <- str_replace_all(complete$B_Education, "\\+", "plus")

# Remove spaces and characters from employment
complete$B_Employment <- str_replace_all(complete$B_Employment, " ", "")
complete$B_Employment <- str_replace_all(complete$B_Employment, "-", "")
complete$B_Employment <- str_replace_all(complete$B_Employment, "/", "")

# Remove spaces and characters from HH income
complete$B_HH_Income.Categories <- 
  str_replace_all(complete$B_HH_Income.Categories, ",", "")
complete$B_HH_Income.Categories <- 
  str_replace_all(complete$B_HH_Income.Categories, " – ", "to")
complete$B_HH_Income.Categories <- 
  str_replace_all(complete$B_HH_Income.Categories, "≥ ", "morethan")

# Remove spaces and characters from chronic medical conditions
complete$Chronic.medical.conditions <- 
  str_replace_all(complete$Chronic.medical.conditions, "0. ", "")
complete$Chronic.medical.conditions <- 
  str_replace_all(complete$Chronic.medical.conditions, "1. ", "")

data_imp <- complete %>%
  dplyr::select(PID, RELATED.INDEX.PID, VISIT, Site, all_of(fuel_source_columns), exclusive_clean,
         wood_use, kerosene_use, solid_fuel_primary, QGIT.Result, TST,
         B_Age, B_Gender, B_Education, B_Employment, B_HH_Income.Categories,
         Living.Area, B_Sleepindex, B_HbA1c, B_BCG,
         B_Weight, B_Height, Mid.upper.arm.circumference..cm., Household.size,
         smoking_status, pack_years, alcohol_use_status, average_drinks_per_day,
         alcohol_disorder, Chronic.medical.conditions, Type.of.House,
         delay, contact_hours_per_week, b08fu0028, TB.Antigen...Nil,
         mean_pm_recentavg, SE_pm_recentavg, temp_recentavg,
         mean_pm_intervalavg, SE_pm_intervalavg, temp_intervalavg,
         temp_lag1, season_lag1, season) %>%
  mutate(
    across(c(all_of(fuel_source_columns), exclusive_clean, wood_use, kerosene_use, solid_fuel_primary,
             B_Gender, B_Employment, Site,
             Living.Area, B_Sleepindex, B_BCG,
             QGIT.Result, B_Education, PID, RELATED.INDEX.PID, 
             smoking_status, alcohol_use_status, Chronic.medical.conditions,
             Type.of.House, b08fu0028, season_lag1,
             alcohol_disorder),
           ~as.factor(.)),
    
    MUAC = Mid.upper.arm.circumference..cm.
  ) %>%
  dplyr::select(-Mid.upper.arm.circumference..cm.)

# Chronological factor levels for VISIT, so that arrange(PID, VISIT) puts
# Entry first (alphabetical order would put Month12 before Month4).
# Unscheduled events follow the scheduled visits, in alphabetical order.
scheduled_visit_levels <- c("ActiveTb", "Entry", "Month4", "Month6", 
                            "Month12", "Month18", "Month24")
visit_levels <- c(scheduled_visit_levels,
                  sort(setdiff(unique(data_imp$VISIT), scheduled_visit_levels)))
data_imp$VISIT <- factor(data_imp$VISIT, levels = visit_levels)

# Variables imputed with polr must be ordered factors
hh_income_levels <- c("1000to4999", "5000to9999", "10000to14999", 
                      "15000to19999", "20000to24999", "25000to29999", 
                      "morethan30000")
data_imp$B_HH_Income.Categories <- factor(data_imp$B_HH_Income.Categories,
                                          levels = hh_income_levels,
                                          ordered = TRUE)

meal_levels <- c("None", "One", "Two", "Three or More")

data_imp$b08fu0028 <- factor(data_imp$b08fu0028,
                             levels = meal_levels,
                             ordered = TRUE)

alcohol_disorder_levels <- c("Low Risk", 
                                 "Harmful Consumption", 
                                 "Alcohol Dependence")
data_imp$alcohol_disorder <- factor(data_imp$alcohol_disorder,
                                    levels = alcohol_disorder_levels,
                                    ordered = TRUE)

education_levels <- c("NotApplicable", "Illiterate", "PrimarySchool0to5years",
                      "HighSchool6to10years", "JuniorCollege11to12years",
                      "College12plusyears", "Professional", "Graduate")
data_imp$B_Education <- factor(data_imp$B_Education,
                               levels = education_levels,
                               ordered = TRUE)

data_imp <- data_imp %>%
  filter(VISIT != "ActiveTb")

# Create a unique VISIT_ID per person for pivot_wider 
data_imp <- data_imp %>%
  arrange(PID, VISIT) %>%  # sort by VISIT's chronological factor levels
  group_by(PID, VISIT) %>%
  mutate(visit_n = row_number()) %>%  # sequence number within each visit type for each person
  ungroup() %>%
  # Unique visit identifier: the VISIT name for the first occurrence of a
  # visit type, with _2, _3, etc. appended for later occurrences
  mutate(VISIT_ID = ifelse(visit_n > 1, paste0(VISIT, "_", visit_n), as.character(VISIT))) %>%
  dplyr::select(-visit_n) # Remove the temporary visit_n column

# Review the new VISIT_ID values 
table(data_imp$VISIT_ID)

# Reshape to wide format (one row per person) for MICE
data_imp_wide <- data_imp %>%
  # Non-time-varying variables
  dplyr::select(PID, RELATED.INDEX.PID, Site, B_Age, B_Gender, B_Education,
         B_Employment, B_HH_Income.Categories, Living.Area, B_Sleepindex,
         B_BCG, Household.size,
         B_Weight, B_Height, MUAC, smoking_status,
         pack_years, alcohol_use_status, average_drinks_per_day,
         alcohol_disorder, Chronic.medical.conditions, Type.of.House, delay,
         B_HbA1c, mean_pm_recentavg, SE_pm_recentavg, temp_recentavg,
         temp_lag1, season_lag1, season) %>%
  # Keep the first row for each PID (visit level to person level)
  distinct(PID, .keep_all = TRUE)

time_varying <- data_imp %>%
  dplyr::select(PID, VISIT_ID, all_of(fuel_source_columns), exclusive_clean, wood_use, kerosene_use,
         solid_fuel_primary, QGIT.Result, TST,
         contact_hours_per_week, b08fu0028, TB.Antigen...Nil) %>%
  pivot_wider(id_cols = PID,
              names_from = VISIT_ID,
              names_sep = "_",
              values_from = c(all_of(fuel_source_columns), "exclusive_clean", "wood_use", "kerosene_use",
                              "solid_fuel_primary", "QGIT.Result", "TST",
                              "contact_hours_per_week", "b08fu0028",
                              "TB.Antigen...Nil"))

data_imp_wide <- left_join(data_imp_wide, 
                           time_varying,
                           by = "PID")

data_imp_wide <- data_imp_wide %>% dplyr::select(PID, RELATED.INDEX.PID, Site,
                                          B_Age, B_Gender, B_Education, B_Employment,
                                          B_HH_Income.Categories,
                                          Living.Area, B_Sleepindex,
                                          B_HbA1c, B_BCG,
                                          B_Weight, B_Height, MUAC,
                                          Household.size, smoking_status,
                                          average_drinks_per_day, alcohol_disorder,
                                          contact_hours_per_week_Month4,
                                          Chronic.medical.conditions,
                                          all_of(paste0(fuel_source_columns, "_Entry")),
                                          exclusive_clean_Entry,
                                          exclusive_clean_Month4,
                                          wood_use_Entry, wood_use_Month4,
                                          kerosene_use_Entry, kerosene_use_Month4,
                                          solid_fuel_primary_Entry,
                                          solid_fuel_primary_Month4,
                                          TB.Antigen...Nil_Entry,
                                          QGIT.Result_Entry, TST_Entry,
                                          mean_pm_recentavg, SE_pm_recentavg,
                                          temp_recentavg, temp_lag1, season_lag1,
                                          season,
                                          Type.of.House)

# An indeterminate Entry QGIT result (3 participants) is kept as an observed
# result. These participants are not eligible for the IGRA conversion cohort
# (which requires a negative Entry result) and are excluded from the baseline
# IGRA model (qgit_positive_entry is set to NA below), but remain in the TST
# analyses if their baseline TST was observed.
message("Indeterminate Entry QGIT results: ",
        sum(as.character(data_imp_wide$QGIT.Result_Entry) %in% "Indeterminate"))

# CLEANING THE DATA PRE-IMPUTATION ############################################

# Remove values for individuals whose weight/height/MUAC is listed as 999 - 
# set it to NA so it is imputed
data_imp_wide[data_imp_wide$PID == "PUNB991013101", "B_Weight"] <- NA

data_imp_wide[data_imp_wide$PID == "PUNB991013101", "MUAC"] <- NA

data_imp_wide[data_imp_wide$PID == "PUNB991013101", "B_Height"] <- NA
data_imp_wide[data_imp_wide$PID == "PUNB991004501", "B_Height"] <- NA
data_imp_wide[data_imp_wide$PID == "PUNB991006603", "B_Height"] <- NA
data_imp_wide[data_imp_wide$PID == "PUNB991006604", "B_Height"] <- NA
data_imp_wide[data_imp_wide$PID == "PUNB991006605", "B_Height"] <- NA
data_imp_wide[data_imp_wide$PID == "PUNB991006606", "B_Height"] <- NA
data_imp_wide[data_imp_wide$PID == "PUNB991006607", "B_Height"] <- NA
data_imp_wide[data_imp_wide$PID == "PUNB991006608", "B_Height"] <- NA
data_imp_wide[data_imp_wide$PID == "PUNB991008601", "B_Height"] <- NA
data_imp_wide[data_imp_wide$PID == "PUNB991008602", "B_Height"] <- NA
data_imp_wide[data_imp_wide$PID == "PUNB991009001", "B_Height"] <- NA

# Remove values for some participants whose weight and height is basically 
# impossible for their age
data_imp_wide[data_imp_wide$PID == "CHEB173011205", "B_Weight"] <- NA
data_imp_wide[data_imp_wide$PID == "CHEB153030403", "B_Weight"] <- NA
data_imp_wide[data_imp_wide$PID == "CHEB173016002", "B_Weight"] <- NA

data_imp_wide[data_imp_wide$PID == "CHEB153030403", "B_Height"] <- NA
data_imp_wide[data_imp_wide$PID == "CHEB173011303", "B_Height"] <- NA
data_imp_wide[data_imp_wide$PID == "CHEB173016002", "B_Height"] <- NA

# Remove one individual with every variable missing and another with
# implausible values throughout
data_imp_wide <- data_imp_wide %>%
  filter(PID != "CHEB173003802" & PID != "PUNB991014602")

# Chronic medical conditions: set the 3 missing values to "No" (the most
# common value)
data_imp_wide$Chronic.medical.conditions <- 
  as.character(data_imp_wide$Chronic.medical.conditions)
data_imp_wide[is.na(data_imp_wide$Chronic.medical.conditions), 
              "Chronic.medical.conditions"] <- "No"
data_imp_wide$Chronic.medical.conditions <-
  factor(data_imp_wide$Chronic.medical.conditions)

# Assign WHO 50th-percentile height, weight and MUAC by age and sex for very
# young children, which MICE imputes poorly
data_imp_wide[data_imp_wide$PID == "PUNB991044306", "B_Weight"] <- 9.5
data_imp_wide[data_imp_wide$PID == "PUNB991044306", "B_Height"] <- 75
data_imp_wide[data_imp_wide$PID == "PUNB991044306", "MUAC"] <- 14.5

# Set one infant's biologically implausible weight and height to the 50th
# percentile
data_imp_wide[data_imp_wide$PID == "PUNB991033302", "B_Weight"] <- 3.5
data_imp_wide[data_imp_wide$PID == "PUNB991033302", "B_Height"] <- 50

# How often fuel use changes between Entry and Month 4/6, among people with
# both, to decide whether Entry values can be filled from Month 4/6
data_imp_wide %>%
  filter(!is.na(exclusive_clean_Entry) & !is.na(exclusive_clean_Month4)) %>%
  count(exclusive_clean_Entry, exclusive_clean_Month4) %>%
  mutate(pct = n / sum(n) * 100)

# About 20% of people change fuel use between Entry and Month 4/6, so Entry
# fuel use is imputed with MICE rather than copied from Month 4/6. The same
# applies to all household fuel variables.
data_imp_wide %>%
  filter(!is.na(wood_use_Entry) & !is.na(wood_use_Month4)) %>%
  count(wood_use_Entry, wood_use_Month4) %>%
  mutate(pct = n / sum(n) * 100)

data_imp_wide %>%
  filter(!is.na(kerosene_use_Entry) & !is.na(kerosene_use_Month4)) %>%
  count(kerosene_use_Entry, kerosene_use_Month4) %>%
  mutate(pct = n / sum(n) * 100)

data_imp_wide %>%
  filter(!is.na(solid_fuel_primary_Entry) & !is.na(solid_fuel_primary_Month4)) %>%
  count(solid_fuel_primary_Entry, solid_fuel_primary_Month4) %>%
  mutate(pct = n / sum(n) * 100)

# Rename columns
data_imp_wide <- rename(data_imp_wide,
                        antigen_minus_nil_Entry = TB.Antigen...Nil_Entry)

# Amount of missingness in each variable
vars_to_check <- c("B_HH_Income.Categories", "B_Education", "B_Employment", 
                   "Living.Area", "B_Age", "B_Gender", 
                   "Site", "B_Weight", "B_BCG", "B_Height", "MUAC", 
                   "B_HbA1c", "Type.of.House",
                   "Chronic.medical.conditions", "B_Sleepindex", 
                   "contact_hours_per_week_Month4", 
                   "antigen_minus_nil_Entry", "QGIT.Result_Entry", "TST_Entry",
                   "exclusive_clean_Entry", "wood_use_Entry", "kerosene_use_Entry",
                   "solid_fuel_primary_Entry", "Household.size",
                   "average_drinks_per_day", "mean_pm_recentavg",
                   "temp_lag1", "season_lag1")

data_imp_wide %>%
  dplyr::select(all_of(vars_to_check)) %>%
  gg_miss_var(show_pct = TRUE) + 
  labs(title = "Percentage of Missing Data by Variable")

data_imp_wide %>%
  dplyr::select(all_of(vars_to_check)) %>%
  tbl_summary(
    missing = "always",  # show the count of NAs
    missing_text = "Missing"
  ) %>%
  bold_labels()

# Flag missing Entry outcomes before imputation, for multiple imputation then
# deletion (MID): continuous QGIT antigen minus nil, QGIT positivity and TST
# positivity
data_imp_wide <- data_imp_wide %>%
  mutate(missing_outcome = ifelse(
    is.na(antigen_minus_nil_Entry),
    1,
    0
  ))

# Examine correlations
data_imp_wide_corr <- as.data.frame(lapply(data_imp_wide,as.numeric))
corr <- data_imp_wide_corr %>% correlate()

B_Age_preds <- arrange(corr, desc(abs(B_Age)))
B_HH_Income.Categories_preds <- arrange(corr, desc(abs(B_HH_Income.Categories)))
B_Sleepindex_preds <- arrange(corr, desc(abs(B_Sleepindex)))
B_BCG_preds <- arrange(corr, desc(abs(B_BCG)))
B_Weight_preds <- arrange(corr, desc(abs(B_Weight)))
B_Height_preds <- arrange(corr, desc(abs(B_Height)))
MUAC_preds <- arrange(corr, desc(abs(MUAC)))
avg_drinks_preds <- arrange(corr, desc(abs(average_drinks_per_day)))
chronic_preds <- arrange(corr, desc(abs(Chronic.medical.conditions)))
fuel_entry_preds <- arrange(corr, desc(abs(exclusive_clean_Entry)))
wood_entry_preds <- arrange(corr, desc(abs(wood_use_Entry)))
kerosene_entry_preds <- arrange(corr, desc(abs(kerosene_use_Entry)))
solid_fuel_primary_entry_preds <- arrange(corr, desc(abs(solid_fuel_primary_Entry)))
contact_hours_preds <- arrange(corr, desc(abs(contact_hours_per_week_Month4)))
HbA1c_preds <- arrange(corr, desc(abs(B_HbA1c)))
antigen_minus_nil_Entry_preds <- arrange(corr, desc(abs(antigen_minus_nil_Entry)))
mean_pm_preds <- arrange(corr, desc(abs(mean_pm_recentavg)))

# Check variable structure and types
str(data_imp_wide)

################################################################################
### Single MICE model: Entry and Month 4/12/24                              ###
################################################################################
# One MICE run on the wide dataset (Entry and all follow-up variables) is
# used for both the baseline analysis (Entry columns) and the survival
# analysis (all visits, reshaped to person-period format). Baseline
# covariates are therefore imputed once for both analyses, and follow-up
# QGIT/TST/PM2.5 data are used as auxiliary information.

# mean_pm_recentavg (baseline exposure) is defined at Entry only and
# mean_pm_intervalavg (survival exposure) at Month 4/12/24 only (see
# data_processing.R), so each is pivoted only for those visits.
pm_recentavg_vars <- c("mean_pm_recentavg", "SE_pm_recentavg", "temp_recentavg")
pm_intervalavg_vars <- c("mean_pm_intervalavg", "SE_pm_intervalavg", "temp_intervalavg")

# season (at each visit's own date) is needed at Entry (baseline and HAP
# models) and at Month 4/12/24 (survival model, at each interval's end visit).
followup_season_vars <- c("season")

# Follow-up (Month 4/12/24) QGIT/TST results and time-varying PM2.5 and
# climate variables, pivoted to one row per PID with a "_<visit>" suffix.
# Interval-average variables are defined only at these visits.
followup_source <- complete %>%
  filter(VISIT %in% c("Month4", "Month12", "Month24")) %>%
  distinct(PID, VISIT, .keep_all = TRUE) %>%
  dplyr::select(PID, VISIT, QGIT.Result, TST,
                any_of(if (mid_outcome_predictors) "TB.Antigen...Nil" else character(0)),
                all_of(pm_intervalavg_vars), all_of(followup_season_vars)) %>%
  mutate(QGIT.Result = as.factor(QGIT.Result))
if (mid_outcome_predictors) {
  followup_source <- followup_source %>%
    rename(antigen_minus_nil = TB.Antigen...Nil) %>%
    mutate(antigen_minus_nil = suppressWarnings(as.numeric(as.character(antigen_minus_nil))))
}

followup_wide <- followup_source %>%
  pivot_wider(id_cols = PID, names_from = VISIT, names_sep = "_",
              values_from = c("QGIT.Result", "TST",
                              any_of(if (mid_outcome_predictors) "antigen_minus_nil" else character(0)),
                              all_of(pm_intervalavg_vars), all_of(followup_season_vars)))

# Replace mean_pm_intervalavg_MonthX/SE_pm_intervalavg_MonthX, which are NA
# when a visit row is missing, with interval_pm_by_pid from data_processing.R:
# the same values where a visit occurred, otherwise the window based on the
# protocol schedule (Entry.Date + 4/12/24 months).
followup_wide <- followup_wide %>%
  dplyr::select(-any_of(c("mean_pm_intervalavg_Month4", "SE_pm_intervalavg_Month4",
                          "mean_pm_intervalavg_Month12", "SE_pm_intervalavg_Month12",
                          "mean_pm_intervalavg_Month24", "SE_pm_intervalavg_Month24"))) %>%
  left_join(interval_pm_by_pid, by = "PID")

# Entry-visit variables from data_imp_wide, with the PM2.5 and climate
# columns given an "_Entry" suffix to match the Month 4/12/24 columns.
# antigen_minus_nil_Entry and missing_outcome are carried through.
# season_Entry is the Entry-visit season used in the HAP and baseline
# PM2.5 models; temp_lag1 and season_lag1 are auxiliary variables only.
entry_fixed_vars <- c("temp_lag1", "season_lag1", "season")

# Every variable to be imputed or analysed (e.g. alcohol_disorder) must be
# included in this select(), in cross_sectional_cols below, and in the
# person-period allow-list, or it is silently dropped. Check names() of the
# saved list_imputed_dataframes_*.Rdata after adding a variable.
survival_wide <- data_imp_wide %>%
  dplyr::select(PID, RELATED.INDEX.PID, Site, B_Age, B_Gender, B_Education,
                B_Employment, B_HH_Income.Categories, Living.Area,
                B_Sleepindex, B_HbA1c, B_BCG, B_Weight, B_Height, MUAC,
                Household.size, smoking_status, average_drinks_per_day, alcohol_disorder,
                contact_hours_per_week_Month4,
                Chronic.medical.conditions, Type.of.House,
                all_of(paste0(fuel_source_columns, "_Entry")),
                exclusive_clean_Entry, wood_use_Entry, kerosene_use_Entry,
                solid_fuel_primary_Entry, QGIT.Result_Entry, TST_Entry,
                antigen_minus_nil_Entry, missing_outcome,
                all_of(pm_recentavg_vars), all_of(entry_fixed_vars)) %>%
  rename_with(~ paste0(.x, "_Entry"), .cols = all_of(c(pm_recentavg_vars, entry_fixed_vars))) %>%
  left_join(followup_wide, by = "PID")

visit_suffixes <- c("Entry", "Month4", "Month12", "Month24")

# Flag missing follow-up outcomes before imputation; these visits are
# censored rather than filled with imputed values
for (vis in visit_suffixes) {
  qgit_col <- paste0("QGIT.Result_", vis)
  tst_col <- paste0("TST_", vis)
  survival_wide[[paste0("missing_qgit_", vis)]] <- as.integer(is.na(survival_wide[[qgit_col]]))
  survival_wide[[paste0("missing_tst_", vis)]] <- as.integer(is.na(survival_wide[[tst_col]]))
  amn_col <- paste0("antigen_minus_nil_", vis)
  if (mid_outcome_predictors && amn_col %in% names(survival_wide)) {
    survival_wide[[paste0("missing_amn_", vis)]] <- as.integer(is.na(survival_wide[[amn_col]]))
  }
}

### Set up the unified imputation model #######################################
init_surv <- mice(survival_wide, maxit = 0)
init_surv$nmis

meth_surv <- init_surv$method
pred_surv <- init_surv$predictorMatrix
pred_surv[, ] <- 0

baseline_covs_for_pred <- c("B_HH_Income.Categories", "B_Education",
                            "B_Employment", "Living.Area",
                            "average_drinks_per_day", "alcohol_disorder", "B_Age",
                            "Type.of.House", "B_Gender", "Site", "B_Weight",
                            "B_BCG", "B_Height", "MUAC", "B_HbA1c",
                            "Chronic.medical.conditions", "B_Sleepindex",
                            "smoking_status",
                            paste0(fuel_source_columns, "_Entry"))

# Predictor matrix: each visit's variables are predicted by the baseline
# covariates and that visit's other measured variables, rather than every
# visit's PM2.5 and climate variables (which would add many collinear
# predictors).
for (vis in visit_suffixes) {
  qgit_col <- paste0("QGIT.Result_", vis)
  tst_col <- paste0("TST_", vis)

  # Recent-average PM2.5/temperature, temp_lag1 and season_lag1 exist only at
  # Entry; interval-average PM2.5/temperature only at Month 4/12/24.
  if (vis == "Entry") {
    pm_cols <- paste0(c(pm_recentavg_vars, entry_fixed_vars), "_Entry")
  } else {
    pm_cols <- paste0(c(pm_intervalavg_vars, followup_season_vars), "_", vis)
  }

  # antigen_minus_nil_Entry is a mutual predictor for the Entry visit only
  entry_only_extra <- if (vis == "Entry") "antigen_minus_nil_Entry" else character(0)
  # With mid_outcome_predictors, follow-up antigen minus nil is imputed (pmm)
  # with that visit's QGIT/TST, and deleted before analysis.
  if (mid_outcome_predictors && vis != "Entry") {
    entry_only_extra <- paste0("antigen_minus_nil_", vis)
  }

  vis_vars <- c(qgit_col, tst_col, pm_cols, entry_only_extra)
  vis_vars <- vis_vars[vis_vars %in% names(meth_surv)]

  # Imputation methods: polyreg for QGIT results at every visit (Negative,
  # Positive or Indeterminate); pmm for continuous variables.
  if (qgit_col %in% names(meth_surv)) meth_surv[qgit_col] <- "polyreg"
  if (tst_col %in% names(meth_surv)) meth_surv[tst_col] <- "pmm"
  if (vis == "Entry" && "antigen_minus_nil_Entry" %in% names(meth_surv)) {
    meth_surv["antigen_minus_nil_Entry"] <- "pmm"
  }
  amn_vis_col <- paste0("antigen_minus_nil_", vis)
  if (mid_outcome_predictors && vis != "Entry" && amn_vis_col %in% names(meth_surv)) {
    meth_surv[amn_vis_col] <- "pmm"
  }

  # Entry mean_pm_recentavg/SE_pm_recentavg are imputed with pmm; follow-up
  # interval averages are not imputed (see below).
  mean_se_cols <- pm_cols[grepl("^(mean|SE)_pm_recentavg", pm_cols)]
  mean_se_cols <- mean_se_cols[mean_se_cols %in% names(meth_surv)]
  meth_surv[mean_se_cols] <- "pmm"

  # Temperature and season variables are determined by site and visit date
  # (ERA5 reanalysis temperature; season from the visit month) and are
  # essentially complete, so they are not imputed. Follow-up
  # mean_pm_intervalavg/SE_pm_intervalavg (site x calendar-month lookup, from
  # interval_pm_by_pid) are also not imputed.
  temp_season_cols <- pm_cols[grepl(
    "^(temp_recentavg|temp_intervalavg|temp_lag|season|mean_pm_intervalavg|SE_pm_intervalavg)",
    pm_cols)]
  temp_season_cols <- temp_season_cols[temp_season_cols %in% names(meth_surv)]
  meth_surv[temp_season_cols] <- ""

  imputed_targets <- setdiff(vis_vars, temp_season_cols)

  for (v in imputed_targets) {
    others_same_visit <- setdiff(vis_vars, v)
    pred_surv[v, c(baseline_covs_for_pred, others_same_visit)] <- 1
  }
}

# Baseline and ID variables that are not imputed
vars_none_surv <- c("PID", "RELATED.INDEX.PID", "Site", "B_Gender", 
                    "B_Employment", "B_Education", "B_Age", "Type.of.House",
                    "Household.size", "Living.Area", 
                    "Chronic.medical.conditions", "missing_outcome",
                    grep("^missing_(qgit|tst|amn)_", names(meth_surv), value = TRUE))
meth_surv[vars_none_surv[vars_none_surv %in% names(meth_surv)]] <- ""

# Baseline covariates imputed once, for both the baseline and survival
# analyses
meth_surv["B_HH_Income.Categories"] <- "polr"
meth_surv["alcohol_disorder"] <- "polr"
meth_surv["B_Sleepindex"] <- "polyreg"
# smoking_status (never/former/current): unordered factor, imputed with
# polyreg
meth_surv["smoking_status"] <- "polyreg"
meth_surv[c("MUAC", "B_Height", "B_Weight", "B_HbA1c",
           "average_drinks_per_day", "contact_hours_per_week_Month4")] <- "pmm"
meth_surv["B_BCG"] <- "logreg"
fuel_entry_sources <- paste0(fuel_source_columns, "_Entry")
fuel_derived_entry <- paste0(c("exclusive_clean", "wood_use", "kerosene_use", "solid_fuel_primary"), "_Entry")
# Derived indicators are reconstructed after imputation, not fitted separately.
meth_surv[fuel_derived_entry] <- ""
pred_surv[, fuel_derived_entry] <- 0
meth_surv["primary_fuel_Entry"] <- "polyreg"
# Logistic coefficient draws failed with a non-positive-definite covariance
# matrix in the questionnaire-aware run. CART samples observed donor responses
# within leaves, retaining the same covariate and MID outcome predictors.
meth_surv[paste0(fuel_secondary_columns, "_Entry")] <- "cart"

pred_surv["B_HH_Income.Categories", baseline_covs_for_pred] <- 1
pred_surv["B_Sleepindex", baseline_covs_for_pred] <- 1
pred_surv["smoking_status", baseline_covs_for_pred] <- 1
pred_surv["alcohol_disorder", baseline_covs_for_pred] <- 1
pred_surv[c("MUAC", "B_Height", "B_Weight", "B_HbA1c",
           "average_drinks_per_day", "contact_hours_per_week_Month4",
           "B_BCG", fuel_entry_sources),
         baseline_covs_for_pred] <- 1

### Outcomes as predictors of incomplete baseline variables (MID) ############
# With mid_outcome_predictors, each incomplete baseline covariate/exposure is
# also predicted by the continuous outcome measurements at every visit: TST
# induration (mm) and QGIT antigen minus nil (IU/mL). The categorical QGIT
# result is not added because it is nearly determined by antigen minus nil.
# Missing outcome values are imputed in the chained equations and deleted
# afterwards.
if (mid_outcome_predictors) {
  outcome_predictors <- intersect(
    c(paste0("TST_", visit_suffixes), paste0("antigen_minus_nil_", visit_suffixes)),
    colnames(pred_surv))
  incomplete_baseline <- intersect(
    c(baseline_covs_for_pred, "contact_hours_per_week_Month4"),
    names(meth_surv)[meth_surv != ""])
  pred_surv[incomplete_baseline, outcome_predictors] <- 1
  message("MID: outcome predictors (", paste(outcome_predictors, collapse = ", "),
          ") added to the imputation models for: ",
          paste(incomplete_baseline, collapse = ", "))

  # mice() returns NA imputations for rows where a predictor is missing and
  # not imputed. Follow-up interval-average PM2.5 is not imputed and is NA
  # before PM2.5 monitoring began, so incomplete, non-imputed variables are
  # removed as predictors when mid_outcome_predictors is TRUE. They are still
  # used in the analyses.
  incomplete_unimputed <- names(meth_surv)[
    meth_surv == "" & colSums(is.na(survival_wide[names(meth_surv)])) > 0]
  incomplete_unimputed <- intersect(incomplete_unimputed, colnames(pred_surv))
  if (length(incomplete_unimputed) > 0) {
    message("MID: removing incomplete, non-imputed predictor(s): ",
            paste(incomplete_unimputed, collapse = ", "))
    pred_surv[, incomplete_unimputed] <- 0
  }
}

# Fuel sources with no observed variation cannot support a conditional
# categorical model. Use observed-value sampling explicitly for those
# sources, report that assumption, and never silently leave missing cells.
fuel_source_methods <- bind_rows(lapply(fuel_entry_sources, function(v) {
  observed <- survival_wide[[v]][!is.na(survival_wide[[v]])]
  n_values <- length(unique(observed))
  if (!length(observed)) stop("No observed answers for ", v, "; review the questionnaire audit before imputing")
  if (n_values == 1L) {
    meth_surv[v] <<- "sample"
    pred_surv[v, ] <<- 0
  }
  tibble(variable = v, n_observed = length(observed), n_missing = sum(is.na(survival_wide[[v]])),
         observed_categories = paste(unique(as.character(observed)), collapse = " | "),
         method = meth_surv[v],
         assumption = if (n_values == 1L) "Observed-value sampling has no category variation; unobserved categories cannot be learned" else
           if (meth_surv[v] == "cart") "CART donor imputation (minbucket 5, cp 0.0001)" else "Conditional categorical imputation")
}))
write_csv(fuel_source_methods, file.path(fp, "fuel_source_imputation_methods.csv"))
# Keep constant source targets in the sampler (sample method), while excluding
# constant predictors from every regression. Other derived fields stay unused.
constant_predictors <- names(survival_wide)[vapply(survival_wide,
  function(x) length(unique(stats::na.omit(x))) <= 1L, logical(1))]
pred_surv[, constant_predictors] <- 0
pred_surv[, fuel_derived_entry] <- 0
diag(pred_surv) <- 0

### Run the MICE model ########################################################
# m = 30 imputations, maxit = 50 iterations. Check mice_all$loggedEvents and
# the trace plots below, and increase maxit if the chains have not stabilized.
# Save every five completed iterations. Resume only identical model inputs.
fuel_checkpoint_path <- here("Data", paste0("mice_fuel_sources_checkpoint_", Sys.Date(), ".rds"))
fuel_checkpoint_signature <- list(data = survival_wide, method = meth_surv,
  predictorMatrix = pred_surv, m = 30L, seed = 123L,
  mice_version = as.character(packageVersion("mice")),
  processing_hash = unname(tools::md5sum(here("fuel_processing_helpers.R"))))
fuel_checkpoint <- if (file.exists(fuel_checkpoint_path)) readRDS(fuel_checkpoint_path) else NULL
if (!is.null(fuel_checkpoint) && identical(fuel_checkpoint$signature, fuel_checkpoint_signature) &&
    inherits(fuel_checkpoint$model, "mids")) {
  mice_all <- fuel_checkpoint$model
  message("Resuming fuel imputation checkpoint at iteration ", mice_all$iteration)
} else {
  if (file.exists(fuel_checkpoint_path)) {
    preserved_path <- paste0(fuel_checkpoint_path, ".incompatible-", format(Sys.time(), "%H%M%S"))
    if (!file.copy(fuel_checkpoint_path, preserved_path)) stop("Could not preserve incompatible checkpoint")
    message("Model/data differ; preserved previous checkpoint at ", preserved_path)
  }
  mice_all <- mice(data = survival_wide, maxit = 0, m = 30, seed = 123,
                   method = meth_surv, predictorMatrix = pred_surv,
                   remove.constant = FALSE)
}
save_fuel_checkpoint <- function(model) {
  temp_path <- paste0(fuel_checkpoint_path, ".tmp")
  saveRDS(list(signature = fuel_checkpoint_signature, model = model), temp_path)
  if (!file.rename(temp_path, fuel_checkpoint_path)) stop("Could not save fuel imputation checkpoint")
}
save_fuel_checkpoint(mice_all)
while (mice_all$iteration < 50L) {
  mice_all <- mice::mice.mids(mice_all, maxit = min(5L, 50L - mice_all$iteration),
                             printFlag = TRUE)
  save_fuel_checkpoint(mice_all)
  message("Saved fuel imputation checkpoint: iteration ", mice_all$iteration, " / 50")
}

# Logged events (predictors removed for collinearity or constancy), one row
# per variable and removed predictor
if (!is.null(mice_all$loggedEvents) && nrow(mice_all$loggedEvents)) {
  print(dplyr::count(mice_all$loggedEvents, dep, meth, out))
} else {
  message("No logged imputation events")
}

# Check that every imputed variable is complete in every completed dataset.
# Remaining NAs indicate an incomplete predictor (see above).
imputed_targets_all <- names(mice_all$method)[mice_all$method != ""]
na_left <- Reduce(`+`, lapply(seq_len(mice_all$m), function(i) {
  colSums(is.na(mice::complete(mice_all, i)[imputed_targets_all]))
}))
if (any(na_left > 0)) {
  print(na_left[na_left > 0])
  stop("mice() left missing values in ", sum(na_left > 0), " imputed variable(s) (counts summed ",
       "across imputations, printed above). Check for incomplete predictors with method = \"\".")
}
message("All imputed variables are complete in all ", mice_all$m, " imputations.")
# Preserve chains and the original missingness masks for diagnostics without
# repeating the expensive imputation. This participant-level object stays in Data/.
saveRDS(mice_all, here("Data", paste0("mice_all_", Sys.Date(), ".rds")))

# Validate source preservation and reconstruct all indicators before any
# analysis files are written. Also describe completed within-household differences;
# the existing MICE model remains participant-level and does not itself enforce
# a common fuel draw for households with no observed questionnaire response.
# Completed differences can also reflect differing observed reports, which
# are retained by design and are not logical contradictions.
fuel_validation <- bind_rows(lapply(seq_len(mice_all$m), function(i) {
  d <- complete_fuel_imputation(mice_all, i)
  for (v in fuel_entry_sources) {
    observed <- !is.na(survival_wide[[v]])
    if (any(as.character(d[[v]][observed]) != as.character(survival_wide[[v]][observed])))
      stop("Observed fuel source response changed in imputation ", i, ": ", v)
  }
  conflicts <- d %>% filter(!is.na(RELATED.INDEX.PID), nzchar(as.character(RELATED.INDEX.PID))) %>%
    group_by(RELATED.INDEX.PID) %>%
    summarise(n_contacts = n(),
      different_analysis_classifications = any(vapply(pick(all_of(fuel_derived_entry)),
        function(x) n_distinct(x) > 1L, logical(1))), .groups = "drop")
  tibble(imputation = i, n_people = nrow(d),
    n_clean = sum(d$exclusive_clean_Entry == "YES"),
    n_wood = sum(d$wood_use_Entry == "YES"),
    n_kerosene = sum(d$kerosene_use_Entry == "YES"),
    n_contradictions = sum(d$exclusive_clean_Entry == "YES" &
      (d$wood_use_Entry == "YES" | d$kerosene_use_Entry == "YES" | d$solid_fuel_primary_Entry == "YES")),
    households_with_different_completed_classifications = sum(conflicts$different_analysis_classifications))
}))
stopifnot(all(fuel_validation$n_contradictions == 0))
write_csv(fuel_validation, file.path(fp, "fuel_completed_data_validation.csv"))

tp_path <- paste0(fp, "/traceplots.pdf")
# Trace plots built from mice_all$chainMean/chainVar with base graphics (one
# page per imputed variable, mean and SD), since plot.mids() renders blank
# panels in some mice/lattice versions
imputed_vars <- names(mice_all$method)[mice_all$method != ""]

pdf(file = tp_path, onefile = TRUE, width = 11, height = 4.5)
for (v in imputed_vars) {
  cm <- mice_all$chainMean[v, , ]
  cv <- mice_all$chainVar[v, , ]
  # Variables with very few missing values have undefined chain variances;
  # a placeholder page is printed for these.
  if (all(is.na(cm)) || all(is.na(cv))) {
    plot.new()
    text(0.5, 0.5, paste0(v, "\n(too few missing values for trace diagnostics: n = ",
                          sum(is.na(mice_all$data[[v]])), ")"))
    next
  }
  par(mfrow = c(1, 2))
  matplot(cm, type = "l", lty = 1, xlab = "Iteration", ylab = "mean",
         main = paste0(v, " - chain means"))
  matplot(sqrt(cv), type = "l", lty = 1, xlab = "Iteration", ylab = "SD",
         main = paste0(v, " - chain SDs"))
}
dev.off()
par(mfrow = c(1, 1))

# Density plots (print() is needed for lattice plots to render inside
# source()). Follow-up interval-average PM2.5 is not imputed, so it has no
# density plot.
dp_path <- paste0(fp, "/density-plot-weight.jpg")
jpeg(file = dp_path, 
     height = "7", width = "13", units = "in", res = 300)
print(densityplot(mice_all, ~ B_Weight))
dev.off()

dp_path <- paste0(fp, "/density-plot-height.jpg")
jpeg(file = dp_path, 
     height = "7", width = "13", units = "in", res = 300)
print(densityplot(mice_all, ~ B_Height))
dev.off()

dp_path <- paste0(fp, "/density-plot-muac.jpg")
jpeg(file = dp_path, 
     height = "7", width = "13", units = "in", res = 300)
print(densityplot(mice_all, ~ MUAC))
dev.off()

dp_path <- paste0(fp, "/density-plot-drinks-per-day.jpg")
jpeg(file = dp_path, 
     height = "7", width = "13", units = "in", res = 300)
print(densityplot(mice_all, ~ average_drinks_per_day))
dev.off()

dp_path <- paste0(fp, "/density-plot-contact-hours.jpg")
jpeg(file = dp_path, 
     height = "7", width = "13", units = "in", res = 300)
print(densityplot(mice_all, ~ contact_hours_per_week_Month4))
dev.off()

dp_path <- paste0(fp, "/density-plot-HbA1c.jpg")
jpeg(file = dp_path, 
     height = "7", width = "13", units = "in", res = 300)
print(densityplot(mice_all, ~ B_HbA1c))
dev.off()

dp_path <- paste0(fp, "/density-plot-mean_pm_recentavg_Entry.jpg")
jpeg(file = dp_path,
     height = "7", width = "13", units = "in", res = 300)
print(densityplot(mice_all, ~ mean_pm_recentavg_Entry))
dev.off()

dp_path <- paste0(fp, "/density-plot-SE_pm_recentavg_Entry.jpg")
jpeg(file = dp_path,
     height = "7", width = "13", units = "in", res = 300)
print(densityplot(mice_all, ~ SE_pm_recentavg_Entry))
dev.off()

dp_path <- paste0(fp, "/density-plot-TST_Entry.jpg")
jpeg(file = dp_path, 
     height = "7", width = "13", units = "in", res = 300)
print(densityplot(mice_all, ~ TST_Entry))
dev.off()

# Follow-up outcome density plots (MID run only). Imputed outcome values are
# deleted before analysis but used as predictors during imputation.
if (mid_outcome_predictors) {
  for (dp_var in c(paste0("antigen_minus_nil_", c("Month4", "Month12", "Month24")),
                   paste0("TST_", c("Month4", "Month12", "Month24")))) {
    if (!dp_var %in% names(mice_all$data) || sum(is.na(mice_all$data[[dp_var]])) < 2) next
    jpeg(file = paste0(fp, "/density-plot-", dp_var, ".jpg"),
         height = "7", width = "13", units = "in", res = 300)
    print(densityplot(mice_all, as.formula(paste("~", dp_var))))
    dev.off()
  }
}

dp_path <- paste0(fp, "/density-plot-antigen-minus-nil-entry.jpg")
jpeg(file = dp_path, 
     height = "7", width = "13", units = "in", res = 300)
print(densityplot(mice_all, ~ antigen_minus_nil_Entry))
dev.off()

# Strip plots
sp_path <- paste0(fp, "/strip-plot-income.jpg")
jpeg(file = sp_path,
     height = "7", width = "13", units = "in", res = 300)
print(stripplot(mice_all, B_HH_Income.Categories))
dev.off()

sp_path <- paste0(fp, "/strip-plot-sleepindex.jpg")
jpeg(file = sp_path,
     height = "7", width = "13", units = "in", res = 300)
print(stripplot(mice_all, B_Sleepindex))
dev.off()

sp_path <- paste0(fp, "/strip-plot-BCG.jpg")
jpeg(file = sp_path,
     height = "7", width = "13", units = "in", res = 300)
print(stripplot(mice_all, B_BCG))
dev.off()

sp_path <- paste0(fp, "/strip-plot-primary_fuel_entry.jpg")
jpeg(file = sp_path,
     height = "7", width = "13", units = "in", res = 300)
print(stripplot(mice_all, primary_fuel_Entry))
dev.off()

sp_path <- paste0(fp, "/strip-plot-qgit-result-entry.jpg")
jpeg(file = sp_path,
     height = "7", width = "13", units = "in", res = 300)
print(stripplot(mice_all, QGIT.Result_Entry))
dev.off()

# Check the correct number of data sets were imputed
imputed_data <- complete_fuel_imputation(mice_all, "long")
length(unique(imputed_data$.imp)) # should = m

# Get the first completed dataset
imputed_data_check <- complete_fuel_imputation(mice_all, 1)

# Check for any remaining missing values in this imputed dataset
missing_counts <- colSums(is.na(imputed_data_check))
print(missing_counts)

# Extract the original data with missing values to get observed statistics
original_data <- mice_all$data

# Variables to compare (list all continuous/scaled variables)
vars_to_compare <- c("MUAC", "B_Height", "B_Weight", "contact_hours_per_week_Month4",
                     "B_HbA1c", "antigen_minus_nil_Entry", "TST_Entry",
                     "average_drinks_per_day", "mean_pm_recentavg_Entry",
                     "SE_pm_recentavg_Entry", "mean_pm_intervalavg_Month4",
                     "SE_pm_intervalavg_Month4", "mean_pm_intervalavg_Month12",
                     "SE_pm_intervalavg_Month12", "mean_pm_intervalavg_Month24",
                     "SE_pm_intervalavg_Month24")

# Calculate means and SDs from the observed data
observed_stats <- data.frame(
  Variable = vars_to_compare,
  Observed_Mean = sapply(original_data[vars_to_compare], function(x) mean(x, na.rm = TRUE)),
  Observed_SD = sapply(original_data[vars_to_compare], function(x) sd(x, na.rm = TRUE))
)

# Extract means and SDs from each imputed dataset and then average them
imputed_means <- sapply(1:mice_all$m, function(i) {
  imputed_data <- complete_fuel_imputation(mice_all, i)
  sapply(imputed_data[vars_to_compare], mean, na.rm = TRUE)
})
imputed_mean_avg <- rowMeans(imputed_means)

imputed_sds <- sapply(1:mice_all$m, function(i) {
  imputed_data <- complete_fuel_imputation(mice_all, i)
  sapply(imputed_data[vars_to_compare], sd, na.rm = TRUE)
})
imputed_sd_avg <- rowMeans(imputed_sds)

# Create the final comparison table
comparison_table <- data.frame(
  Variable = vars_to_compare,
  Observed_Mean = observed_stats$Observed_Mean,
  Imputed_Mean = imputed_mean_avg,
  Observed_SD = observed_stats$Observed_SD,
  Imputed_SD = imputed_sd_avg
)

# Format table
ft_comparison <- flextable(comparison_table) %>%
  set_header_labels(Imputed_Mean = "Completed Mean", Imputed_SD = "Completed SD") %>%
  set_caption("Comparison of Observed vs. Completed Dataset Means and Standard Deviations") %>%
  colformat_double(digits = 3) %>%
  autofit()

# Save table
save_as_image(ft_comparison, path = file.path(fp, "stats_comparison.png"))
# These legacy tables compare observed data with the full completed datasets.
# Label them accurately; imputed-only comparisons are exported below.
readr::write_csv(comparison_table %>% rename(
  Completed_Mean = Imputed_Mean, Completed_SD = Imputed_SD),
  file.path(fp, "observed_vs_completed_continuous.csv"))

# 1. Identify the categorical variables in the mice_all object
cat_vars_to_compare <- c("B_BCG", "exclusive_clean_Entry","B_Sleepindex",
                         "B_HH_Income.Categories", "QGIT.Result_Entry")
fuel_diagnostic_vars <- paste0(fuel_source_columns, "_Entry")

# 2. Prepare empty data frame to store results
comparison_list <- list()

# 3. Loop through each categorical variable to calculate observed and imputed proportions
for (v in cat_vars_to_compare) {
  # Get observed proportions
  observed_proportions <- mice_all$data %>%
    filter(!is.na(!!sym(v))) %>%
    count(!!sym(v)) %>%
    mutate(Proportion = n / sum(n)) %>%
    dplyr::select(Category = !!sym(v), Observed_Proportion = Proportion) %>%
    mutate(Category = as.character(Category)) 
  
  # Get imputed proportions from all datasets
  imputed_proportions_list <- list()
  for (i in 1:mice_all$m) {
    imputed_data <- complete_fuel_imputation(mice_all, i)
    imputed_proportions_list[[i]] <- imputed_data %>%
      count(!!sym(v)) %>%
      mutate(Proportion = n / sum(n)) %>%
      dplyr::select(Category = !!sym(v), Imputed_Proportion = Proportion) %>%
      mutate(Category = as.character(Category)) 
  }
  
  # Average the imputed proportions
  avg_imputed_proportions <- imputed_proportions_list %>%
    bind_rows() %>%
    group_by(Category) %>%
    summarise(Imputed_Proportion = mean(Imputed_Proportion, na.rm = TRUE))
  
  # Join observed and imputed proportions
  comparison_table <- full_join(observed_proportions, avg_imputed_proportions, by = "Category") %>%
    mutate(Variable = v) %>%
    dplyr::select(Variable, Category, everything())
  
  # Add to the list
  comparison_list[[v]] <- comparison_table
}

# 4. Combine all tables into a single flextable
final_comparison_table <- bind_rows(comparison_list)

# 5. Format table with flextable
ft_proportions <- flextable(final_comparison_table) %>%
  # Set the column headers
  set_header_labels(Variable = "Variable",
                    Category = "Category",
                    Observed_Proportion = "Observed Proportion",
                    Imputed_Proportion = "Completed Dataset Proportion") %>%
  # Merge rows in the "Variable" column to group categories
  merge_v(j = "Variable") %>%
  # Set a simple table theme
  theme_booktabs() %>%
  # Format numeric columns to 3 decimal places
  colformat_double(digits = 3) %>%
  # Autofit column widths
  autofit() %>%
  # Add a bottom border to the last row of each variable group
  border(
    i = which(final_comparison_table$Variable != lead(final_comparison_table$Variable, default = "")),
    border.bottom = officer::fp_border(color = "black", style = "solid", width = 1.5),
    part = "body"
  ) %>%
  # Set the table caption
  set_caption("Comparison of Observed vs. Completed Dataset Proportions for Categorical Variables")

# Save table
save_as_image(ft_proportions, path = file.path(fp, "cat_proportions.png"))

# Compare observed values with only the cells that were originally missing.
# Keep one row per imputation so between-imputation variation remains visible.
imputed_only_continuous <- map_dfr(intersect(vars_to_compare, imputed_targets_all), function(v) {
  missing <- is.na(mice_all$data[[v]])
  observed <- mice_all$data[[v]][!missing]
  map_dfr(seq_len(mice_all$m), function(i) {
    values <- complete_fuel_imputation(mice_all, i)[[v]][missing]
    tibble(variable = v, imputation = i, n_observed = length(observed),
           n_imputed = length(values), observed_mean = mean(observed),
           observed_sd = sd(observed), imputed_only_mean = mean(values),
           imputed_only_sd = sd(values))
  })
})
readr::write_csv(imputed_only_continuous, file.path(fp, "observed_vs_imputed_only_continuous.csv"))
imputed_only_categorical <- map_dfr(intersect(union(cat_vars_to_compare, fuel_diagnostic_vars), imputed_targets_all), function(v) {
  missing <- is.na(mice_all$data[[v]])
  observed <- as.character(mice_all$data[[v]][!missing])
  categories <- levels(mice_all$data[[v]])
  map_dfr(seq_len(mice_all$m), function(i) {
    values <- as.character(complete_fuel_imputation(mice_all, i)[[v]][missing])
    tibble(variable = v, imputation = i, category = categories,
           n_observed = length(observed), n_imputed = length(values),
           observed_proportion = vapply(categories, function(z) mean(observed == z), numeric(1)),
           imputed_only_proportion = vapply(categories, function(z) mean(values == z), numeric(1)))
  })
})
readr::write_csv(imputed_only_categorical, file.path(fp, "observed_vs_imputed_only_categorical.csv"))

################################################################################
### Baseline (Entry-only) imputed datasets                                  ###
################################################################################
# Selects the Entry columns from each completed dataset and removes the
# "_Entry" suffix from the PM2.5 and climate columns, to match the names used
# in propensity_scores.R.
cross_sectional_cols <- c("PID", "RELATED.INDEX.PID", "Site", "B_Age", "B_Gender",
                          "B_Education", "B_Employment", "B_HH_Income.Categories",
                          "Living.Area", "B_Sleepindex", "B_HbA1c", "B_BCG",
                          "B_Weight", "B_Height", "MUAC", "Household.size",
                          "average_drinks_per_day", "alcohol_disorder", "Chronic.medical.conditions",
                          "smoking_status",
                          paste0(fuel_source_columns, "_Entry"),
                          "contact_hours_per_week_Month4", "exclusive_clean_Entry",
                          "wood_use_Entry", "kerosene_use_Entry",
                          "solid_fuel_primary_Entry",
                          "antigen_minus_nil_Entry",
                          "QGIT.Result_Entry", "TST_Entry", "missing_outcome",
                          "missing_qgit_Entry", "missing_tst_Entry",
                          paste0(c(pm_recentavg_vars, entry_fixed_vars), "_Entry"),
                          "Type.of.House")

impdata <- complete_fuel_imputation(mice_all, action = "long", include = TRUE) %>%
  dplyr::select(.imp, .id, all_of(cross_sectional_cols)) %>%
  rename_with(~ sub("_Entry$", "", .x), .cols = matches(paste0(
    "^(", paste(c(pm_recentavg_vars, entry_fixed_vars), collapse = "|"), ")_Entry$"
  ))) %>%
  rename(missing_qgit_outcome = missing_qgit_Entry,
        missing_tst_outcome = missing_tst_Entry)

# exclusive_clean_Entry must be numeric 0/1 for downstream analyses
impdata$solid_fuel_use <- ifelse(impdata$exclusive_clean_Entry == "YES",
                                        0,
                                        1)

# wood_use, kerosene_use and solid_fuel_primary_Entry are coded so that
# YES = uses this fuel; convert YES/No to 1/0
impdata$wood_use <- ifelse(impdata$wood_use_Entry == "YES", 1, 0)
impdata$kerosene_use <- ifelse(impdata$kerosene_use_Entry == "YES", 1, 0)
impdata$solid_fuel_primary <- ifelse(impdata$solid_fuel_primary_Entry == "YES", 1, 0)

# Binary outcome 1: QGIT positivity at Entry. An observed Indeterminate result
# is neither positive nor negative, so it is set to NA (excluded from the
# baseline IGRA model).
impdata$qgit_positive_entry <- case_when(
  as.character(impdata$QGIT.Result_Entry) == "Positive" ~ 1,
  as.character(impdata$QGIT.Result_Entry) == "Negative" ~ 0,
  TRUE ~ NA_real_)

# Binary outcome 2: TST positivity at Entry, from the imputed induration
# using a >= 5 mm cutoff
impdata$tst_positive_entry <- ifelse(impdata$TST_Entry >= 5, 1, 0)

# Delete the imputed outcome for people originally missing it, separately for
# each outcome (multiple imputation then deletion)
impdata <- impdata %>%
  mutate(
    antigen_minus_nil_Entry = ifelse(
      missing_outcome == 1, NA_real_, antigen_minus_nil_Entry
    ),
    qgit_positive_entry = ifelse(
      missing_qgit_outcome == 1, NA_real_, qgit_positive_entry
    ),
    tst_positive_entry = ifelse(
      missing_tst_outcome == 1, NA_real_, tst_positive_entry
    )
  )

impdata$.id <- NULL
m <- mice_all$m
allImputations <- imputationList(lapply(1:m,
                                        function(n)
                                          subset(impdata,
                                                 subset = .imp == n)))

# Save imputed data sets

name_longdata <- paste0("long_imputed_dataframes_", Sys.Date(), ".Rdata")
save(impdata, file = here("Data", name_longdata))

name_listdata <- paste0("list_imputed_dataframes_", Sys.Date(), ".Rdata")
save(allImputations, file = here("Data", name_listdata))

################################################################################
### Discrete-time survival data: person-period datasets for QGIT/TST        ###
### conversion, with time-varying PM2.5 at each scheduled visit             ###
################################################################################
### Reshape each completed dataset into person-period format ##################
# One row per PID x interval (Entry-M4, M4-M12, M12-M24), separately for the
# QGIT and TST conversion cohorts. Events are defined from observed values
# only, using the missingness flags captured before imputation; a missing or
# indeterminate result is treated as unobserved.
#
# Follow-up is right-censored at the last confirmed negative assessment. The
# first interval with a missing/unknown result is kept (event = NA,
# censored = 1) for the censoring model used to build inverse probability of
# censoring weights (IPCW) in discrete_time_survival.R; event = NA rows are
# dropped before the outcome model is fitted.
build_person_period <- function(completed_df, test) {
  if (test == "qgit") {
    # Cohort membership requires an observed baseline result
    # (missing_qgit_Entry == 0), not an imputed one, so the cohort is the same
    # in every imputation.
    negative_at_entry <- completed_df$QGIT.Result_Entry == "Negative" &
      completed_df$missing_qgit_Entry == 0
    result_cols <- c(Month4 = "QGIT.Result_Month4", Month12 = "QGIT.Result_Month12", 
                     Month24 = "QGIT.Result_Month24")
    missing_cols <- c(Month4 = "missing_qgit_Month4", Month12 = "missing_qgit_Month12", 
                      Month24 = "missing_qgit_Month24")
    is_positive <- function(x) x == "Positive"
  } else {
    # Same for TST: require an observed baseline TST (missing_tst_Entry == 0).
    negative_at_entry <- completed_df$TST_Entry < 5 &
      completed_df$missing_tst_Entry == 0
    result_cols <- c(Month4 = "TST_Month4", Month12 = "TST_Month12", 
                     Month24 = "TST_Month24")
    missing_cols <- c(Month4 = "missing_tst_Month4", Month12 = "missing_tst_Month12", 
                      Month24 = "missing_tst_Month24")
    is_positive <- function(x) x >= 5
  }
  
  cohort <- completed_df[negative_at_entry & !is.na(negative_at_entry), ]

  # An indeterminate follow-up QGIT is treated as missing (censored).
  unusable_result <- function(df, vis) {
    miss <- df[[missing_cols[vis]]] == 1
    if (test == "qgit") {
      miss <- miss | as.character(df[[result_cols[vis]]]) %in% "Indeterminate"
    }
    miss
  }

  # Participants with no outcome information at any follow-up visit are kept:
  # they contribute a censored (event = NA) interval-1 row to the censoring
  # model, so the censoring weights account for loss to follow-up after the
  # baseline visit.
  entry_only_no_outcome <- rowSums(sapply(names(missing_cols), function(v) unusable_result(cohort, v))) ==
    length(missing_cols)
  message(test, ": ", sum(entry_only_no_outcome),
          " participant(s) with no follow-up outcome, censored at interval 1")

  intervals <- list(
    list(number = 1, visit = "Month4"),
    list(number = 2, visit = "Month12"),
    list(number = 3, visit = "Month24")
  )
  
  person_period <- map_dfr(intervals, function(iv) {
    vis <- iv$visit
    was_missing <- unusable_result(cohort, vis)
    tibble(
      PID = cohort$PID,
      interval = iv$number,
      event = ifelse(was_missing, NA_real_, as.numeric(is_positive(cohort[[result_cols[vis]]]))),
      mean_pm_intervalavg = cohort[[paste0("mean_pm_intervalavg_", vis)]],
      SE_pm_intervalavg = cohort[[paste0("SE_pm_intervalavg_", vis)]],
      # temp_intervalavg is carried through but not used as a confounder
      temp_intervalavg = cohort[[paste0("temp_intervalavg_", vis)]],
      # season at each interval's end visit: time-varying confounder for the
      # continuous PM2.5 model
      season = cohort[[paste0("season_", vis)]]
    ) %>%
      # Raw education/employment/income columns, season_Entry,
      # average_drinks_per_day and contact_hours_per_week_Month4, so that
      # discrete_time_survival.R can apply collapse_covariates(). season_Entry
      # is the fixed-baseline confounder for the HAP exposure models.
      # temp_lag1_Entry is carried through but not used.
      bind_cols(cohort %>% dplyr::select(-PID) %>%
                 dplyr::select(all_of(intersect(c(
                   "B_HH_Income.Categories", "B_Education", "B_Employment",
                   "Living.Area", "B_Age", "B_Gender", "Site", "B_Height", "B_Weight",
                   "MUAC", "B_HbA1c", "B_BCG", "B_Sleepindex", "Type.of.House",
                   "average_drinks_per_day", "contact_hours_per_week_Month4",
                   "temp_lag1_Entry", "season_Entry", "solid_fuel_use", "wood_use",
                   "kerosene_use", "solid_fuel_primary",
                   # Household.size and smoking_status are needed for the
                   # candidate_confounders diagnostic in
                   # discrete_time_survival.R
                   "Household.size", "smoking_status",
                    "season_Month4", "season_Month12", "season_Month24"
                 ), names(cohort)))))
  }) %>%
    # Right-censoring: rows are kept through the first positive result, or
    # through the first missing/unknown assessment (kept as a censored = 1,
    # event = NA row), whichever comes first. Later observed results are not
    # used after a gap, since conversion during the gap is unknown.
    arrange(PID, interval) %>%
    group_by(PID) %>%
    mutate(
      first_event_interval = suppressWarnings(min(interval[event == 1], na.rm = TRUE)),
      first_missing_interval = suppressWarnings(min(interval[is.na(event)], na.rm = TRUE)),
      keep_upto = case_when(
        is.finite(first_event_interval) &
          (!is.finite(first_missing_interval) | first_event_interval <= first_missing_interval) ~
          first_event_interval,
        # Keep the censoring interval itself (event = NA) for the censoring
        # model
        is.finite(first_missing_interval) ~ first_missing_interval,
        TRUE ~ max(interval)
      ),
      # A retained row with an unknown outcome is the interval at which this
      # person was censored
      censored = as.numeric(is.na(event))
    ) %>%
    filter(interval <= keep_upto) %>%
    dplyr::select(-first_event_interval, -first_missing_interval, -keep_upto) %>%
    ungroup()

  person_period
}

# Person-period row counts per interval, including censored (event = NA)
# rows, as a check on the censoring logic

survival_completed_list <- lapply(1:mice_all$m, function(i) complete_fuel_imputation(mice_all, i))

# Follow-up PM2.5 and season for every contact, including those with no row
# for a follow-up visit (e.g. contacts seen at baseline only, who are censored
# at interval 1): interval-average PM2.5 from interval_pm_by_pid, and season at
# the actual visit date or, without a visit, at the protocol date
# (Entry.Date + 4/12/24 months).
season_of <- function(d) {
  m <- month(d)
  case_when(m %in% c(12, 1, 2) ~ "Winter", m %in% c(3, 4, 5) ~ "Pre-Monsoon",
            m %in% c(6, 7, 8, 9) ~ "Monsoon", m %in% c(10, 11) ~ "Post-Monsoon")
}
entry_dates <- complete %>%
  filter(str_replace_all(VISIT, " ", "") == "Entry") %>%
  distinct(PID, .keep_all = TRUE) %>%
  transmute(PID_chr = as.character(PID), entry_date = as.Date(Entry.Date))
pm_det <- interval_pm_by_pid %>% mutate(PID_chr = as.character(PID)) %>% dplyr::select(-PID)
survival_completed_list <- lapply(survival_completed_list, function(df) {
  season_levels <- levels(factor(df$season_Entry))
  df <- df %>%
    mutate(PID_chr = as.character(PID)) %>%
    dplyr::select(-any_of(setdiff(names(pm_det), "PID_chr"))) %>%
    left_join(pm_det, by = "PID_chr") %>%
    left_join(entry_dates, by = "PID_chr")
  for (vis in c("Month4", "Month12", "Month24")) {
    col <- paste0("season_", vis)
    nominal <- season_of(df$entry_date %m+% months(c(Month4 = 4, Month12 = 12, Month24 = 24)[[vis]]))
    df[[col]] <- factor(coalesce(as.character(df[[col]]), nominal), levels = season_levels)
  }
  df %>% dplyr::select(-PID_chr, -entry_date)
})

# Add the household fuel exposures (from their imputed _Entry columns) and
# the collapsed covariates
survival_completed_list <- lapply(survival_completed_list, function(df) {
  df$solid_fuel_use <- ifelse(df$exclusive_clean_Entry == "YES", 0, 1)
  df$wood_use <- ifelse(df$wood_use_Entry == "YES", 1, 0)
  df$kerosene_use <- ifelse(df$kerosene_use_Entry == "YES", 1, 0)
  df$solid_fuel_primary <- ifelse(df$solid_fuel_primary_Entry == "YES", 1, 0)
  df
})

person_period_qgit <- lapply(survival_completed_list, build_person_period, test = "qgit")
person_period_tst <- lapply(survival_completed_list, build_person_period, test = "tst")

message("QGIT person-period rows per interval (imputation 1): ")
print(table(person_period_qgit[[1]]$interval, useNA = "always"))
message("TST person-period rows per interval (imputation 1): ")
print(table(person_period_tst[[1]]$interval, useNA = "always"))

allImputations_survival_qgit <- imputationList(person_period_qgit)
allImputations_survival_tst <- imputationList(person_period_tst)

name_survival_qgit <- paste0("survival_person_period_qgit_", Sys.Date(), ".Rdata")
save(allImputations_survival_qgit, file = here("Data", name_survival_qgit))

name_survival_tst <- paste0("survival_person_period_tst_", Sys.Date(), ".Rdata")
save(allImputations_survival_tst, file = here("Data", name_survival_tst))

################################################################################
### Sensitivity analysis: alternative conversion definitions                ###
################################################################################
# The primary definitions are threshold crossings (QGIT: Negative to
# Positive; TST: <5 mm to >=5 mm). Changes near a threshold can reflect assay
# or reader variability, or TST boosting, so stricter definitions are built:
#
#   qgit_strict   Entry: QGIT Negative and antigen minus nil < 0.20 IU/mL
#                 Event: QGIT Positive and antigen minus nil > 0.70 IU/mL
#                 (Nemes et al. 2017). A follow-up Positive between 0.35 and
#                 0.70 IU/mL is not an event; the person stays at risk.
#   tst_strict10  Entry: TST < 5 mm. Event: TST >= 10 mm at follow-up
#                 (implies an increase of >= 6 mm; Pai et al. 2006).
#   tst_incr10    Entry: TST < 5 mm. Event: TST >= 10 mm and an increase of
#                 >= 10 mm from Entry (Pai et al. 2006).
#
# Events are defined from observed values only, as in build_person_period():
# the missingness flags captured before mice() define what was observed, and
# follow-up antigen minus nil is taken from the raw (pre-imputation) data.
#
# build_person_period_def() uses the same censoring rules as
# build_person_period(). A check below confirms that, given the primary
# definitions, it reproduces person_period_qgit/person_period_tst exactly.

# Observed (pre-imputation) antigen minus nil at each scheduled visit, from
# the same source rows as the corresponding QGIT results
amn_entry <- time_varying %>%
  transmute(PID = as.character(PID),
            amn_obs_Entry = suppressWarnings(as.numeric(as.character(TB.Antigen...Nil_Entry))))
amn_followup <- complete %>%
  filter(VISIT %in% c("Month4", "Month12", "Month24")) %>%
  distinct(PID, VISIT, .keep_all = TRUE) %>%
  transmute(PID = as.character(PID), VISIT,
            amn = suppressWarnings(as.numeric(as.character(TB.Antigen...Nil)))) %>%
  pivot_wider(id_cols = PID, names_from = VISIT, values_from = amn, names_prefix = "amn_obs_")
amn_observed <- full_join(amn_entry, amn_followup, by = "PID")
for (vis in c("Month4", "Month12", "Month24")) {
  col <- paste0("amn_obs_", vis)
  if (!col %in% names(amn_observed)) amn_observed[[col]] <- NA_real_
}

# MID deletion check (mid_outcome_predictors = TRUE): observed antigen minus
# nil values in the completed datasets must equal the raw values, and every
# raw missing value must be flagged as missing.
if (mid_outcome_predictors) {
  for (i in seq_along(survival_completed_list)) {
    chk <- survival_completed_list[[i]] %>%
      mutate(PID = as.character(PID)) %>%
      inner_join(amn_observed, by = "PID")
    for (vis in c("Entry", "Month4", "Month12", "Month24")) {
      flag <- chk[[paste0("missing_amn_", vis)]]
      completed_val <- chk[[paste0("antigen_minus_nil_", vis)]]
      raw_val <- chk[[paste0("amn_obs_", vis)]]
      bad <- (flag == 0 & !(abs(completed_val - raw_val) < 1e-8)) | (is.na(raw_val) & flag == 0)
      if (any(bad, na.rm = TRUE)) {
        stop("MID check failed: antigen_minus_nil_", vis, " in completed dataset ", i,
             " does not match the observed raw values for ", sum(bad, na.rm = TRUE),
             " participant(s).")
      }
    }
  }
  message("MID check passed: observed antigen minus nil values are unchanged by imputation ",
          "in all ", length(survival_completed_list), " completed datasets.")
}

# definition: a list with
#   eligible(df)        -> logical vector, TRUE = in the risk set at Entry
#   observed(df, vis)   -> logical vector, TRUE = outcome usable at this visit
#   event(df, vis)      -> logical vector, TRUE = conversion at this visit
#                          (only consulted where observed() is TRUE)
build_person_period_def <- function(completed_df, definition) {
  completed_df <- completed_df %>%
    mutate(PID_chr = as.character(PID)) %>%
    left_join(amn_observed, by = c("PID_chr" = "PID")) %>%
    dplyr::select(-PID_chr)

  visits <- c("Month4", "Month12", "Month24")
  eligible <- definition$eligible(completed_df)
  cohort <- completed_df[eligible & !is.na(eligible), ]

  # Participants with no observed follow-up outcome are kept (censored at
  # interval 1), as in build_person_period()

  person_period <- map_dfr(seq_along(visits), function(i) {
    vis <- visits[i]
    obs <- definition$observed(cohort, vis)
    tibble(
      PID = cohort$PID,
      interval = i,
      event = ifelse(!obs, NA_real_, as.numeric(definition$event(cohort, vis))),
      mean_pm_intervalavg = cohort[[paste0("mean_pm_intervalavg_", vis)]],
      SE_pm_intervalavg = cohort[[paste0("SE_pm_intervalavg_", vis)]],
      temp_intervalavg = cohort[[paste0("temp_intervalavg_", vis)]],
      season = cohort[[paste0("season_", vis)]]
    ) %>%
      bind_cols(cohort %>% dplyr::select(-PID) %>%
                  dplyr::select(all_of(intersect(c(
                    "B_HH_Income.Categories", "B_Education", "B_Employment",
                    "Living.Area", "B_Age", "B_Gender", "Site", "B_Height", "B_Weight",
                    "MUAC", "B_HbA1c", "B_BCG", "B_Sleepindex", "Type.of.House",
                    "average_drinks_per_day", "contact_hours_per_week_Month4",
                    "temp_lag1_Entry", "season_Entry", "solid_fuel_use", "wood_use",
                    "kerosene_use", "solid_fuel_primary",
                    "Household.size", "smoking_status",
                    "season_Month4", "season_Month12", "season_Month24"
                  ), names(cohort)))))
  }) %>%
    # Censoring logic as in build_person_period()
    arrange(PID, interval) %>%
    group_by(PID) %>%
    mutate(
      first_event_interval = suppressWarnings(min(interval[event == 1], na.rm = TRUE)),
      first_missing_interval = suppressWarnings(min(interval[is.na(event)], na.rm = TRUE)),
      keep_upto = case_when(
        is.finite(first_event_interval) &
          (!is.finite(first_missing_interval) | first_event_interval <= first_missing_interval) ~
          first_event_interval,
        is.finite(first_missing_interval) ~ first_missing_interval,
        TRUE ~ max(interval)
      ),
      censored = as.numeric(is.na(event))
    ) %>%
    filter(interval <= keep_upto) %>%
    dplyr::select(-first_event_interval, -first_missing_interval, -keep_upto) %>%
    ungroup()

  person_period
}

qgit_result <- function(df, vis) as.character(df[[paste0("QGIT.Result_", vis)]])
tst_mm <- function(df, vis) df[[paste0("TST_", vis)]]
amn_mm <- function(df, vis) df[[paste0("amn_obs_", vis)]]

qgit_eligible_primary <- function(df) {
  as.character(df$QGIT.Result_Entry) == "Negative" & df$missing_qgit_Entry == 0
}
tst_eligible_primary <- function(df) df$TST_Entry < 5 & df$missing_tst_Entry == 0

conversion_definitions <- list(
  # Primary definitions in this framework, used only for the check below
  qgit_primary_check = list(
    eligible = qgit_eligible_primary,
    observed = function(df, vis) {
      df[[paste0("missing_qgit_", vis)]] == 0 & qgit_result(df, vis) != "Indeterminate"
    },
    event = function(df, vis) qgit_result(df, vis) == "Positive"
  ),
  tst_primary_check = list(
    eligible = tst_eligible_primary,
    observed = function(df, vis) df[[paste0("missing_tst_", vis)]] == 0,
    event = function(df, vis) tst_mm(df, vis) >= 5
  ),
  # Sensitivity definitions (see header above)
  qgit_strict = list(
    eligible = function(df) {
      qgit_eligible_primary(df) & !is.na(df$amn_obs_Entry) & df$amn_obs_Entry < 0.20
    },
    observed = function(df, vis) {
      df[[paste0("missing_qgit_", vis)]] == 0 &
        qgit_result(df, vis) != "Indeterminate" & !is.na(amn_mm(df, vis))
    },
    event = function(df, vis) qgit_result(df, vis) == "Positive" & amn_mm(df, vis) > 0.70
  ),
  tst_strict10 = list(
    eligible = tst_eligible_primary,
    observed = function(df, vis) df[[paste0("missing_tst_", vis)]] == 0,
    event = function(df, vis) tst_mm(df, vis) >= 10
  ),
  tst_incr10 = list(
    eligible = tst_eligible_primary,
    observed = function(df, vis) df[[paste0("missing_tst_", vis)]] == 0,
    event = function(df, vis) tst_mm(df, vis) >= 10 & (tst_mm(df, vis) - df$TST_Entry) >= 10
  )
)

# Check that the generalized builder reproduces the primary person-period
# data exactly (rows, events, censoring) in every imputation
same_person_period <- function(a, b) {
  key <- c("PID", "interval", "event", "censored")
  a <- a %>% dplyr::select(all_of(key)) %>% mutate(PID = as.character(PID)) %>% arrange(PID, interval)
  b <- b %>% dplyr::select(all_of(key)) %>% mutate(PID = as.character(PID)) %>% arrange(PID, interval)
  isTRUE(all.equal(as.data.frame(a), as.data.frame(b), check.attributes = FALSE))
}
for (i in seq_along(survival_completed_list)) {
  chk_q <- build_person_period_def(survival_completed_list[[i]], conversion_definitions$qgit_primary_check)
  chk_t <- build_person_period_def(survival_completed_list[[i]], conversion_definitions$tst_primary_check)
  if (!same_person_period(chk_q, person_period_qgit[[i]]) ||
      !same_person_period(chk_t, person_period_tst[[i]])) {
    stop("build_person_period_def() did not reproduce the primary person-period ",
         "data for imputation ", i, " - do not use the sensitivity datasets until ",
         "this is resolved.")
  }
}
message("Regression check passed: build_person_period_def() reproduces the primary ",
        "QGIT and TST person-period data in all ", length(survival_completed_list), " imputations.")

sensitivity_definitions <- setdiff(names(conversion_definitions),
                                   c("qgit_primary_check", "tst_primary_check"))

# Event counts by definition and interval (imputation 1; identical across
# imputations because events use observed data only). Includes the primary
# definitions. Few events may not support the 3-df PM2.5 spline.
event_count_table <- map_dfr(c(list(qgit_primary = person_period_qgit[[1]],
                                    tst_primary = person_period_tst[[1]]),
                               setNames(lapply(sensitivity_definitions, function(def) {
                                 build_person_period_def(survival_completed_list[[1]],
                                                         conversion_definitions[[def]])
                               }), sensitivity_definitions)),
                             function(pp) {
  pp %>%
    filter(interval <= 2) %>%
    group_by(interval) %>%
    summarise(n_at_risk = n(),
              n_events = sum(event == 1, na.rm = TRUE),
              n_censored = sum(is.na(event)), .groups = "drop") %>%
    mutate(n_people_in_cohort = n_distinct(pp$PID), .before = 1)
}, .id = "definition")
message("Conversion events by definition (interval 1 = Entry->Month4, 2 = Month4->Month12):")
print(event_count_table, n = Inf)
write_csv(event_count_table, file.path(fp, "conversion_definition_event_counts.csv"))

for (def in sensitivity_definitions) {
  pp_list <- lapply(survival_completed_list, build_person_period_def,
                    definition = conversion_definitions[[def]])
  obj_name <- paste0("allImputations_survival_", def)
  assign(obj_name, imputationList(pp_list))
  save(list = obj_name,
       file = here("Data", paste0("survival_person_period_", def, "_", Sys.Date(), ".Rdata")))
  message("Saved: survival_person_period_", def, "_", Sys.Date(), ".Rdata")
}
