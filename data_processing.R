library(tidyverse)
library(readxl)
library(stringr)
library(kableExtra)
library(zscorer)
library(here)
library(lubridate)
library(ecmwfr)
library(terra)

setwd(here())

# Read in the data
baseline <- read_xlsx(here("Data", "baseline.xlsx"))
IAP <- read_xlsx(here("Data", "IAP.xlsx"))
IAP <- IAP[-1, ] # Remove the row with the descriptions of variables
infection <- read_xlsx(here("Data", "TST_IGRA.xlsx"))
nutrition <- read_xlsx(here("Data", "nutritional.xlsx"))
nutrition <- nutrition[-1,]
psychosocial <- read_xlsx(here("Data", "psychosocial.xlsx"))
psychosocial <- psychosocial[-1,]
ATT_initiation <- read_xlsx(here("Data", "ATT_initiation.xlsx"))
clinical <- read_xlsx(here("Data", "followup_clinical.xlsx"))
clinical <- clinical[-1, ]
PM2.5_Chennai <- read_xlsx(here("Data", "Statistics_Chennai.xlsx"))
PM2.5_Pune <- read_xlsx(here("Data", "Statistics_pune.xlsx"))

# Reduce issues with column names by removing spaces
names(baseline) <- make.names(colnames(baseline), unique = TRUE)
names(IAP) <- make.names(colnames(IAP), unique = TRUE)
names(infection) <- make.names(colnames(infection), unique = TRUE)
names(nutrition) <- make.names(colnames(nutrition), unique = TRUE)
names(psychosocial) <- make.names(colnames(psychosocial), unique = TRUE)
names(ATT_initiation) <- make.names(colnames(ATT_initiation), unique = TRUE)
names(clinical) <- make.names(colnames(clinical), unique = TRUE)

# IAP, nutrition, psychosocial and clinical data each have a per-visit
# "visdat" column (Excel serial date). Each is renamed per source and converted
# to Date before joining, to avoid ambiguous visdat.x/visdat.y columns. The
# infection data (TST_IGRA.xlsx) has no date column of its own.
excel_serial_to_date <- function(x) as.Date(as.numeric(x), origin = "1899-12-30")

IAP <- IAP %>% mutate(visdat_iap = excel_serial_to_date(visdat)) %>% select(-visdat)
nutrition <- nutrition %>% mutate(visdat_nutrition = excel_serial_to_date(visdat)) %>% select(-visdat)
psychosocial <- psychosocial %>% mutate(visdat_psychosocial = excel_serial_to_date(visdat)) %>% select(-visdat)
clinical <- clinical %>% mutate(visdat_clinical = excel_serial_to_date(visdat)) %>% select(-visdat)

# Pivot the infection data from wide (one block of columns per visit:
# Entry.___, M4.___, M12.___, M24.___, Uns.___) to long format, giving one row
# per PID per VISIT.
infection <- infection %>%
  pivot_longer(
    cols = starts_with("Entry") | starts_with("M4") | starts_with("M12") | starts_with("M24") | starts_with("Uns"),
    names_to = c("VISIT", ".value"),
    names_pattern = "^(Entry|M\\d+|Uns)\\.(.*)$"
  )

# Relabel visit codes
infection$VISIT <- case_match(infection$VISIT,
                              "Entry" ~ "Entry",
                              "M4" ~ "Month 4",
                              "M12" ~ "Month 12",
                              "M24" ~ "Month 24",
                              "Uns" ~ "Unscheduled Event")

# Step 1: Identify participants who were QGIT negative at baseline
baseline_negative <- infection %>%
  filter(VISIT == "Entry", QGIT.Result == "Negative") %>%
  pull(PID)

# Step 2: For those participants, check follow-up visits
# Indeterminate results are treated as non-conversion
followup <- infection %>%
  filter(PID %in% baseline_negative,
         VISIT != "Entry")

# Step 3: Define conversion at each follow-up visit
conversion <- followup %>%
  mutate(converted = QGIT.Result == "Positive") %>%
  group_by(PID) %>%
  summarise(
    converted_month4  = any(VISIT == "Month 4"  & QGIT.Result == "Positive", na.rm = TRUE),
    converted_month12 = any(VISIT == "Month 12" & QGIT.Result == "Positive", na.rm = TRUE),
    converted_month24 = any(VISIT == "Month 24" & QGIT.Result == "Positive", na.rm = TRUE),
    converted_any     = any(QGIT.Result == "Positive", na.rm = TRUE),
    .groups = "drop"
  )

# Step 4: Check counts
conversion %>%
  summarise(
    n_baseline_negative = n(),
    n_converted_month4  = sum(converted_month4),
    n_converted_month12 = sum(converted_month12),
    n_converted_month24 = sum(converted_month24),
    n_converted_any     = sum(converted_any)
  )

# Step 1: Check follow-up visit attendance among baseline-negative participants
infection %>%
  filter(PID %in% baseline_negative) %>%
  group_by(VISIT) %>%
  summarise(
    n = n(),
    n_with_qgit = sum(!is.na(QGIT.Result)),
    .groups = "drop"
  )

# Step 2: Check how many participants have each combination of visits
infection %>%
  filter(PID %in% baseline_negative) %>%
  group_by(PID) %>%
  summarise(
    has_month4  = any(VISIT == "Month 4",  na.rm = TRUE),
    has_month12 = any(VISIT == "Month 12", na.rm = TRUE),
    has_month24 = any(VISIT == "Month 24", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  count(has_month4, has_month12, has_month24) %>%
  arrange(desc(n))

# Unscheduled events: if Unscheduled.Event.Name gives a specific month (not
# "Unscheduled" or empty), use it as the VISIT

process_visit_column <- function(df) {
  df %>%
    mutate(
      VISIT = if_else(
        VISIT == "Unscheduled Event" &
          !is.na(Unscheduled.Event.Name) &
          Unscheduled.Event.Name != "" &
          Unscheduled.Event.Name != "Unscheduled",
        Unscheduled.Event.Name,
        VISIT
      )
    )
}

# Rename HH.PID in ATT_initiation to PID for joining
ATT_initiation <- rename(ATT_initiation,
                         PID = HH.PID)

# Ensure consistent naming for visits
standardize_visits <- function(df) {
  df %>%
    mutate(VISIT = VISIT %>%
             str_trim() %>% # Remove trailing/leading spaces
             str_replace_all("\\s+", " ") %>% # Turn multiple spaces into one space
             str_to_title() %>% # Change, e.g., month 24 into Month 24
             str_replace_all(" ", "")) %>% # Remove all spaces
    mutate(VISIT = case_when(
      str_detect(VISIT, "Activetb") ~ "Active TB",
      str_detect(VISIT, "Entry") ~ "Entry",
      str_detect(VISIT, "Month4") ~ "Month 4",
      str_detect(VISIT, "Month6") ~ "Month 6",
      str_detect(VISIT, "Month12") ~ "Month 12",
      str_detect(VISIT, "Month18|M18|M-18") ~ "Month 18",
      str_detect(VISIT, "Month24|M24|Monthe24") ~ "Month 24",
      str_detect(VISIT, "Unschedule") ~ "Unscheduled Event",
      TRUE ~ VISIT 
    ))
}

# Unscheduled events are specific to each data collection instrument and a
# participant can have several. Because they share the VISIT label
# "Unscheduled Event", joining by (PID, VISIT) would cross-match unrelated
# events, so each is given a unique VISIT label per PID and data set.
disambiguate_unscheduled <- function(df, source_name) {
  df %>%
    group_by(PID) %>%
    mutate(VISIT = if_else(
      VISIT == "Unscheduled Event",
      paste0("Unscheduled Event (", source_name, " #", row_number(), ")"),
      VISIT
    )) %>%
    ungroup()
}

# Apply to all relevant data frames

# IAP
IAP <- IAP %>% process_visit_column() %>%
  standardize_visits() 

scheduled <- IAP %>%
  filter(VISIT != "Unscheduled Event") %>%
  distinct(PID, VISIT, .keep_all = TRUE) # This removes duplicates for M4, M12, etc.

unscheduled <- IAP %>%
  filter(VISIT == "Unscheduled Event") %>% # This keeps all unscheduled rows
  disambiguate_unscheduled("IAP")

IAP <- bind_rows(scheduled, unscheduled)

# Infection
infection <- infection %>% 
  standardize_visits() 

scheduled <- infection %>%
  filter(VISIT != "Unscheduled Event") %>%
  distinct(PID, VISIT, .keep_all = TRUE) # This removes duplicates for M4, M12, etc.

unscheduled <- infection %>%
  filter(VISIT == "Unscheduled Event") %>% # This keeps all unscheduled rows
  disambiguate_unscheduled("infection")

infection <- bind_rows(scheduled, unscheduled)

# Nutrition
nutrition <- nutrition %>% process_visit_column() %>%
  standardize_visits() 

scheduled <- nutrition %>%
  filter(VISIT != "Unscheduled Event") %>%
  distinct(PID, VISIT, .keep_all = TRUE) # This removes duplicates for M4, M12, etc.

unscheduled <- nutrition %>%
  filter(VISIT == "Unscheduled Event") %>% # This keeps all unscheduled rows
  disambiguate_unscheduled("nutrition")

nutrition <- bind_rows(scheduled, unscheduled)

# Psychosocial
psychosocial <- psychosocial %>% process_visit_column() %>%
  standardize_visits() 

scheduled <- psychosocial %>%
  filter(VISIT != "Unscheduled Event") %>%
  distinct(PID, VISIT, .keep_all = TRUE) # This removes duplicates for M4, M12, etc.

unscheduled <- psychosocial %>%
  filter(VISIT == "Unscheduled Event") %>% # This keeps all unscheduled rows
  disambiguate_unscheduled("psychosocial")

psychosocial <- bind_rows(scheduled, unscheduled)

# Clinical
clinical <- clinical %>% process_visit_column() %>%
  standardize_visits() 

scheduled <- clinical %>%
  filter(VISIT != "Unscheduled Event") %>%
  distinct(PID, VISIT, .keep_all = TRUE) # This removes duplicates for M4, M12, etc.

unscheduled <- clinical %>%
  filter(VISIT == "Unscheduled Event") %>% # This keeps all unscheduled rows
  disambiguate_unscheduled("clinical")

clinical <- bind_rows(scheduled, unscheduled)

# Join the data sets together
# First, build the list of all PID-VISIT combinations across the long-format
# data sets (IAP, infection, psychosocial, nutrition)
all_pid_visits <- bind_rows(
  ATT_initiation %>% mutate(VISIT = "Entry") %>% select(PID, VISIT),
  IAP %>% select(PID, VISIT),
  infection %>% select(PID, VISIT),
  nutrition %>% select(PID, VISIT),
  psychosocial %>% select(PID, VISIT),
  clinical %>% select(PID, VISIT)
)

# Distinct PID-VISIT pairs
full_time_points <- all_pid_visits %>%
  distinct(PID, VISIT)

# Left join each data set onto the PID-VISIT list; visits absent from a data
# set get NA for its columns
data <- full_time_points %>%
  left_join(baseline, by = "PID") %>%
  left_join(IAP, by = c("PID", "VISIT")) %>%
  left_join(infection, by = c("PID", "VISIT")) %>%
  left_join(nutrition, by = c("PID", "VISIT")) %>% 
  left_join(psychosocial, by = c("PID", "VISIT")) %>%
  left_join(ATT_initiation, by = "PID") %>%
  left_join(clinical, by = c("PID", "VISIT"))

# baseline and ATT_initiation both have Site and RELATED.INDEX.PID, giving
# .x/.y pairs after joining. Check they agree before keeping the baseline
# (.x) version.
site_mismatches <- with(data, sum(Site.x != Site.y | is.na(Site.x) != is.na(Site.y), na.rm = TRUE))
related_pid_mismatches <- with(data, sum(RELATED.INDEX.PID.x != RELATED.INDEX.PID.y |
                                            is.na(RELATED.INDEX.PID.x) != is.na(RELATED.INDEX.PID.y),
                                          na.rm = TRUE))
if (site_mismatches > 0 || related_pid_mismatches > 0) {
  warning(sprintf(
    "baseline and ATT_initiation disagree on Site (%d rows) or RELATED.INDEX.PID (%d rows); keeping baseline's (.x) version - investigate before trusting downstream results.",
    site_mismatches, related_pid_mismatches
  ))
}

# Visit date for each row, so that PM2.5 and temperature exposure windows can
# be defined relative to each visit (needed for the time-varying survival
# analysis). Entry rows use Entry.Date, falling back to the questionnaire
# visit dates if it is missing; other visits use the questionnaire dates.
data <- data %>%
  mutate(
    visit_date = case_when(
      VISIT == "Entry" & !is.na(Entry.Date) ~ Entry.Date,
      TRUE ~ coalesce(visdat_iap, visdat_nutrition, visdat_psychosocial, visdat_clinical)
    )
  )

# Check coverage and consistency of visit_date
visit_date_missing <- data %>%
  filter(VISIT %in% c("Entry", "Month 4", "Month 12", "Month 24")) %>%
  group_by(VISIT) %>%
  summarise(n = n(), n_missing_visit_date = sum(is.na(visit_date)), .groups = "drop")
print(visit_date_missing)

entry_date_mismatches <- data %>%
  filter(VISIT == "Entry", !is.na(Entry.Date)) %>%
  mutate(
    questionnaire_date = coalesce(visdat_iap, visdat_nutrition, visdat_psychosocial, visdat_clinical),
    # Convert both to Date before subtracting (Entry.Date is POSIXct)
    mismatch_days = abs(as.numeric(as.Date(questionnaire_date) - as.Date(Entry.Date)))
  ) %>%
  filter(!is.na(mismatch_days) & mismatch_days > 30)
if (nrow(entry_date_mismatches) > 0) {
  warning(sprintf(
    "%d Entry rows have a questionnaire-recorded visit date more than 30 days from Entry.Date - Entry.Date is being used for visit_date regardless, but this is worth spot-checking.",
    nrow(entry_date_mismatches)
  ))
}

# Scheduled Month 4/12/24 rows with no recorded visit date use the protocol
# schedule (Entry.Date + 4/12/24 months), so that observed outcomes at those
# visits still have an exposure window.
n_filled_before <- sum(is.na(data$visit_date))
data <- data %>%
  mutate(
    nominal_offset_months = case_when(
      VISIT == "Month 4"  ~ 4,
      VISIT == "Month 12" ~ 12,
      VISIT == "Month 24" ~ 24,
      TRUE ~ NA_real_
    ),
    visit_date = if_else(
      is.na(visit_date) & !is.na(nominal_offset_months) & !is.na(Entry.Date),
      Entry.Date %m+% months(nominal_offset_months),
      visit_date
    )
  ) %>%
  select(-nominal_offset_months)
message(sprintf(
  "Filled %d of %d missing visit_date values using the nominal Entry.Date + 4/12/24-month schedule.",
  n_filled_before - sum(is.na(data$visit_date)), n_filled_before
))

# Combine with the PM2.5 data
# Prepare the PM2.5 master dataset
PM2.5_Chennai <- PM2.5_Chennai %>%
  mutate(Site = "NIRT") %>%
  rename(mean_pm = Chennai,
         SD_pm = `Standard Deviation`,
         SE_pm = `Standard error`,
         lower_ci_pm = `CI lower`,
         upper_ci_pm = `CI upper`)
PM2.5_Pune <- PM2.5_Pune %>%
  mutate(Site = "BJGMC") %>%
  rename(mean_pm = Pune,
         SD_pm = SD,
         SE_pm = SE,
         lower_ci_pm = `lower CI`,
         upper_ci_pm = `Upper CI`)

PM2.5_data <- rbind(PM2.5_Chennai, PM2.5_Pune)

data <- data %>%
  mutate(
    # Exposure dates are relative to each row's own visit_date
    date = floor_date(visit_date, "month"),
    # date_lag1 (one month before the visit) is used only to derive temp_lag1
    # and season_lag1
    date_lag1 = floor_date(visit_date %m-% months(1), "month")
  )

################################################################################
### Multi-month average PM2.5 exposure windows                              ###
### Both windows include their endpoint months. The SE of a K-month average ###
### is sqrt(mean(SE_i^2) / K), propagating within-month sampling            ###
### uncertainty only.                                                       ###
###                                                                         ###
### - mean_pm_recentavg: the 3 calendar months before each row's            ###
###   visit_date (offsets 1-3); exposure for the baseline analysis.         ###
### - mean_pm_intervalavg: every calendar month from the previous scheduled ###
###   visit to the current visit; exposure for the survival analysis.       ###
###   Defined only for Month 4/12/24 rows.                                  ###
################################################################################
data <- data %>% mutate(.row_id = row_number())

## mean_pm_recentavg: fixed 3-month pre-visit window, every row ##############
recentavg_months <- data %>%
  dplyr::select(.row_id, Site.x, visit_date) %>%
  filter(!is.na(visit_date)) %>%
  slice(rep(seq_len(n()), each = 3)) %>%
  group_by(.row_id) %>%
  mutate(month_offset = 1:3) %>%
  ungroup() %>%
  mutate(target_month = floor_date(visit_date %m-% months(month_offset), "month")) %>%
  left_join(PM2.5_data, by = c("Site.x" = "Site", "target_month" = "Months"))

recentavg_summary <- recentavg_months %>%
  group_by(.row_id) %>%
  summarise(
    # No na.rm: the average is NA if any of the 3 months has no PM2.5 data
    mean_pm_recentavg = mean(mean_pm),
    SE_pm_recentavg = sqrt(mean(SE_pm^2) / n()),
    n_months_recentavg = n(),
    .groups = "drop"
  )

data <- data %>% left_join(recentavg_summary, by = ".row_id")

## mean_pm_intervalavg: variable-width window using ACTUAL visit dates ######
visit_dates_by_pid <- data %>%
  filter(VISIT %in% c("Entry", "Month 4", "Month 12", "Month 24")) %>%
  distinct(PID, VISIT, .keep_all = TRUE) %>%
  mutate(VISIT_clean = str_replace_all(VISIT, " ", "")) %>%
  dplyr::select(PID, VISIT_clean, visit_date) %>%
  pivot_wider(id_cols = PID, names_from = VISIT_clean, values_from = visit_date,
             names_prefix = "visit_date_")

data <- data %>%
  left_join(visit_dates_by_pid, by = "PID") %>%
  mutate(
    # Date of the preceding scheduled visit, i.e. the start of the interval
    # ending at this row's visit_date. NA for Entry and unscheduled rows.
    interval_start_date = case_when(
      VISIT == "Month 4"  ~ visit_date_Entry,
      VISIT == "Month 12" ~ visit_date_Month4,
      VISIT == "Month 24" ~ visit_date_Month12,
      TRUE ~ as.Date(NA)
    )
  ) %>%
  dplyr::select(-visit_date_Entry, -visit_date_Month4, -visit_date_Month12, -visit_date_Month24)

intervalavg_base <- data %>%
  dplyr::select(.row_id, Site.x, interval_start_date, visit_date) %>%
  filter(!is.na(interval_start_date) & !is.na(visit_date))

# Rows where the previous scheduled visit is dated after the current one
# (data irregularities) are excluded from mean_pm_intervalavg and reported.
n_interval_reversed <- sum(intervalavg_base$interval_start_date > intervalavg_base$visit_date)
if (n_interval_reversed > 0) {
  message(sprintf(
    "%d of %d rows eligible for mean_pm_intervalavg have interval_start_date > visit_date (previous visit recorded after current visit) - excluding these from mean_pm_intervalavg; this is worth spot-checking in the raw visit dates.",
    n_interval_reversed, nrow(intervalavg_base)
  ))
}
intervalavg_base <- intervalavg_base %>%
  filter(interval_start_date <= visit_date)

intervalavg_months <- intervalavg_base %>%
  mutate(
    month_seq = purrr::map2(
      floor_date(interval_start_date, "month"),
      floor_date(visit_date, "month"),
      ~ seq(.x, .y, by = "month")
    )
  ) %>%
  tidyr::unnest(month_seq) %>%
  rename(target_month = month_seq) %>%
  left_join(PM2.5_data, by = c("Site.x" = "Site", "target_month" = "Months"))

intervalavg_summary <- intervalavg_months %>%
  group_by(.row_id) %>%
  summarise(
    mean_pm_intervalavg = mean(mean_pm),
    SE_pm_intervalavg = sqrt(mean(SE_pm^2) / n()),
    n_months_intervalavg = n(),
    .groups = "drop"
  )

# Check the number of months per interval and how often the averages are NA
# because of a missing month
message("mean_pm_intervalavg: n_months_intervalavg distribution among defined rows:")
print(table(intervalavg_summary$n_months_intervalavg, useNA = "always"))

n_intervalavg_eligible <- nrow(intervalavg_base)
n_intervalavg_missing <- n_intervalavg_eligible - nrow(intervalavg_summary %>% filter(!is.na(mean_pm_intervalavg)))
message(sprintf(
  "mean_pm_recentavg missing (NA) for %d of %d rows with a visit_date; mean_pm_intervalavg missing for %d of %d eligible (Month 4/12/24) rows.",
  sum(is.na(data$mean_pm_recentavg) & !is.na(data$visit_date)), sum(!is.na(data$visit_date)),
  n_intervalavg_missing, n_intervalavg_eligible
))

data <- data %>%
  left_join(intervalavg_summary, by = ".row_id") %>%
  dplyr::select(-interval_start_date)
# .row_id is kept to join the matching temperature averages later in the
# script; it is dropped in the ERA5 temperature section.

## mean_pm_intervalavg_det: deterministic PID x interval exposure, defined ##
## even when a participant has no row for a Month 4/12/24 visit. PM2.5 is a ##
## site x calendar-month lookup, so the window only needs start and end     ##
## dates: the actual visit dates where available, otherwise the protocol    ##
## schedule (Entry.Date + 4/12/24 months). Nominal dates are always taken   ##
## from Entry.Date, so consecutive missing visits do not compound. Used in  ##
## place of the visit-row-based columns in multiple_imputation.R, where     ##
## PM2.5 is not imputed. NA remains only where the PM2.5 series itself does ##
## not cover a needed month.                                                ##
pid_entry_dates <- data %>%
  filter(VISIT == "Entry") %>%
  distinct(PID, .keep_all = TRUE) %>%
  dplyr::select(PID, Site.x, Entry.Date) %>%
  filter(!is.na(Entry.Date))

# One row per PID and scheduled visit, including visits with no row in
# `data` (filled by the nominal schedule below)
visit_dates_det <- pid_entry_dates %>%
  dplyr::select(PID, Entry.Date) %>%
  left_join(
    data %>%
      filter(VISIT %in% c("Month 4", "Month 12", "Month 24")) %>%
      distinct(PID, VISIT, .keep_all = TRUE) %>%
      mutate(VISIT_clean = str_replace_all(VISIT, " ", "")) %>%
      dplyr::select(PID, VISIT_clean, visit_date) %>%
      pivot_wider(id_cols = PID, names_from = VISIT_clean, values_from = visit_date,
                 names_prefix = "visit_date_"),
    by = "PID"
  ) %>%
  mutate(
    visit_date_Month4  = coalesce(visit_date_Month4,  Entry.Date %m+% months(4)),
    visit_date_Month12 = coalesce(visit_date_Month12, Entry.Date %m+% months(12)),
    visit_date_Month24 = coalesce(visit_date_Month24, Entry.Date %m+% months(24))
  )

visit_dates_det <- visit_dates_det %>%
  left_join(pid_entry_dates %>% dplyr::select(PID, Site.x), by = "PID")

# One data frame per interval, stacked into a long PID x interval frame
interval_bounds_det <- bind_rows(
  visit_dates_det %>% transmute(PID, Site.x, vis = "Month4",
                                start_date = as.Date(Entry.Date),
                                end_date = as.Date(visit_date_Month4)),
  visit_dates_det %>% transmute(PID, Site.x, vis = "Month12",
                                start_date = as.Date(visit_date_Month4),
                                end_date = as.Date(visit_date_Month12)),
  visit_dates_det %>% transmute(PID, Site.x, vis = "Month24",
                                start_date = as.Date(visit_date_Month12),
                                end_date = as.Date(visit_date_Month24))
) %>%
  filter(!is.na(start_date) & !is.na(end_date) & start_date <= end_date)

interval_pm_months_det <- interval_bounds_det %>%
  mutate(
    month_seq = purrr::map2(
      floor_date(start_date, "month"), floor_date(end_date, "month"),
      ~ seq(.x, .y, by = "month")
    )
  ) %>%
  tidyr::unnest(month_seq) %>%
  rename(target_month = month_seq) %>%
  left_join(PM2.5_data, by = c("Site.x" = "Site", "target_month" = "Months"))

# Column names match those used in followup_wide in multiple_imputation.R,
# which these replace
interval_pm_by_pid <- interval_pm_months_det %>%
  group_by(PID, vis) %>%
  summarise(mean_pm_intervalavg = mean(mean_pm),
           SE_pm_intervalavg = sqrt(mean(SE_pm^2) / n()),
           .groups = "drop") %>%
  pivot_wider(id_cols = PID, names_from = vis,
             values_from = c(mean_pm_intervalavg, SE_pm_intervalavg),
             names_glue = "{.value}_{vis}")

# Diagnostic: PID x interval x month combinations with a valid window but no
# matching site-month in PM2.5_data
pm_coverage_by_site <- PM2.5_data %>%
  group_by(Site) %>%
  summarise(pm_min_month = min(Months, na.rm = TRUE),
           pm_max_month = max(Months, na.rm = TRUE),
           n_months_covered = n_distinct(Months), .groups = "drop")
cat("\nPM2.5_data coverage by site:\n")
print(pm_coverage_by_site)

months_needed_vs_covered <- interval_pm_months_det %>%
  mutate(pm_missing_this_month = is.na(mean_pm)) %>%
  left_join(pm_coverage_by_site, by = c("Site.x" = "Site"))

n_month_gaps <- sum(months_needed_vs_covered$pm_missing_this_month)
cat(sprintf(
  "\n%d of %d (PID x interval x calendar-month) rows across valid windows have no matching PM2.5_data row (target_month outside/absent from that site's coverage).\n",
  n_month_gaps, nrow(months_needed_vs_covered)
))
if (n_month_gaps > 0) {
  cat("Distribution of the missing target_months (out-of-coverage months driving mean_pm_intervalavg's residual NAs):\n")
  print(months_needed_vs_covered %>%
         filter(pm_missing_this_month) %>%
         count(Site.x, target_month, pm_min_month, pm_max_month) %>%
         arrange(Site.x, target_month), n = Inf)
  cat(sprintf(
    "\n%d distinct PIDs affected by at least one out-of-coverage month.\n",
    n_distinct(months_needed_vs_covered$PID[months_needed_vs_covered$pm_missing_this_month])
  ))
}

message(sprintf(
  "mean_pm_intervalavg (deterministic, real-or-nominal window): %d of %d PIDs have a defined value at all three of Month4/12/24.",
  sum(complete.cases(interval_pm_by_pid %>% dplyr::select(starts_with("mean_pm_intervalavg_")))),
  nrow(interval_pm_by_pid)
))

# Remove one participant who appears to have no data on any of the variables
data <- data %>% filter(PID != "CHEB173003802")

# Processing fuel use data ####################################################

# Rename fuel columns
data <- rename(data,
               primary_fuel = b05ia001a,
               secondary_fuel_kerosene = b05ia01b1,
               secondary_fuel_gas = b05ia01b2,
               secondary_fuel_coal = b05ia01b3,
               secondary_fuel_electricity = b05ia01b4,
               secondary_fuel_wood = b05ia01b5,
               secondary_fuel_dung = b05ia01b6,
               secondary_fuel_crop_waste = b05ia01b7,
               secondary_fuel_none = b05ia01b8,
               secondary_fuel_other = b05ia01b9,
               RELATED.INDEX.PID = RELATED.INDEX.PID.x,
               Site = Site.x,
               # baseline.xlsx has its own "TST" field (constant per person), distinct
               # from the per-visit TST induration in TST_IGRA.xlsx. After joining, the
               # per-visit induration is kept as TST and the baseline field as
               # B_TST_baseline_raw.
               TST = TST.y,
               B_TST_baseline_raw = TST.x)

# Normalize questionnaire responses, identify likely default No values, share
# compatible source responses within household/visit, then derive indicators.
# Retains *_raw fields and flags conflicting reports for review.
source(here("fuel_processing_helpers.R"))
complete <- process_fuel_questionnaire(data)
write_fuel_questionnaire_audit(complete,
  here("tables_and_figures", paste0("fuel-questionnaire-audit-", Sys.Date())))

# Some participants report exact income and others an income category;
# convert everyone to categories
complete$B_HH_Income.Categories <- factor(complete$B_HH_Income.Categories)
levels(complete$B_HH_Income.Categories)

complete <- complete %>%
  mutate(B_HH_Income.Categories = 
           case_when(B_HH_Income >= 30000 ~ "≥ 30,000",
                     B_HH_Income >= 25000 & B_HH_Income <= 29999 ~ "25,000 – 29,999",
                     B_HH_Income >= 20000 & B_HH_Income <= 24999 ~ "20,000 – 24,999",
                     B_HH_Income >= 15000 & B_HH_Income <= 19999 ~ "15,000 – 19,999",
                     B_HH_Income >= 10000 & B_HH_Income <= 14999 ~ "10,000 – 14,999",
                     B_HH_Income >= 5000 & B_HH_Income <= 9999 ~ "5,000 – 9,999",
                     B_HH_Income >= 1000 & B_HH_Income <= 4999 ~ "1,000 – 4,999",
                     TRUE ~ B_HH_Income.Categories))

# Indicator for HbA1c > 7 (suggests uncontrolled diabetes)
complete <- complete %>%
  mutate(uncontrolled_diabetes = case_when(
    is.na(B_HbA1c) ~ NA_character_,
    B_HbA1c >= 7 ~ "YES",
    B_HbA1c < 7 ~ "No"))

names(complete)

# Smoking information #########################################################
complete$b06ps001b <- as.numeric(complete$b06ps001b)
complete$b06ps002c <- as.numeric(complete$b06ps002c)
complete$b06ps03a1 <- as.numeric(complete$b06ps03a1)

complete <- complete %>%
  mutate(smoking_status =
           case_when(B_curr_smk == "Yes" ~ "current_smoker",
                     (B_curr_smk == "No" & B_ever_smk == "Yes") ~ "former_smoker",
                     B_Smoking_yes.no == "No" ~ "never_smoker",
                     TRUE ~ NA_character_)) %>%
  mutate(years_since_quit = 
           case_when(b06ps01b1 == "4. Years" ~ b06ps001b,
                     b06ps01b1 == "3. Months" ~ (b06ps001b / 12),
                     b06ps01b1 == "1. Days" ~ (b06ps001b / 365),
                     TRUE ~ NA_real_)) %>%
  mutate(years_smoked =
           case_when(
             smoking_status == "former_smoker" ~ B_Age - years_since_quit - b06ps002c,
             smoking_status == "current_smoker" ~ B_Age - b06ps002c,
             smoking_status == "never_smoker" ~ 0,
             TRUE ~ NA_real_)) %>%
  mutate(average_cigarettes_per_day = case_when(
    smoking_status != "never_smoker" ~ b06ps03a1,
    smoking_status == "never_smoker" ~ 0)) %>%
  mutate(pack_years = 
           case_when(smoking_status == "never_smoker" ~ 0,
                     TRUE ~ (average_cigarettes_per_day / 20) * years_smoked)) 

# Set smoking and drinking variables to 0 for children aged 12 and under
complete <- complete %>%
  mutate(smoking_status = ifelse(
    (B_Age <= 12) & is.na(smoking_status),
    "never_smoker",
    smoking_status)) %>%
  mutate(years_smoked = ifelse(
    (B_Age <= 12) & smoking_status == "never_smoker",
    0,
    years_smoked)) %>%
  mutate(average_cigarettes_per_day = ifelse(
    (B_Age <= 12) & smoking_status == "never_smoker",
    0,
    average_cigarettes_per_day)) %>%
  mutate(pack_years = ifelse(
    (B_Age <= 12) & smoking_status == "never_smoker",
    0,
    pack_years))

# Alcohol information #########################################################
complete$b06ps0018 <- factor(complete$b06ps0018)
complete$b06ps0018 <- fct_recode(complete$b06ps0018,
                                 "1" = "0. 1 or 2",
                                 "3" = "1. 3 or 4",
                                 "5" = "2. 5 or 6",
                                 "8" = "3. 7  8  or 9")
complete$b06ps0018 <- as.numeric(as.character(complete$b06ps0018))


complete <- complete %>%
  mutate(alcohol_use_status = 
           case_when(
             b06ps0016 == "Yes" & b06ps0017 != "0.Never" ~ "current_user",
             b06ps0016 == "Yes" & b06ps0017 == "0.Never" ~ "former_user",
             b06ps0016 == "No" ~ "never_user",
             TRUE ~ NA_character_),
         average_drinks_per_day = 
           case_when(
             alcohol_use_status == "never_user" ~ 0,
             TRUE ~ b06ps0018),
         alcohol_disorder = 
           case_when(
            B_AUDIT < 8 ~ "Low Risk",
            B_AUDIT >= 8 & B_AUDIT <= 14 ~ "Harmful Consumption",
            B_AUDIT > 14 ~ "Alcohol Dependence"))

# Set smoking and drinking variables to 0 for children aged 12 and under
complete <- complete %>%
  mutate(alcohol_use_status = ifelse(
    (B_Age <= 12) & is.na(alcohol_use_status),
    "never_user",
    alcohol_use_status)) %>%
  mutate(average_drinks_per_day = ifelse(
    (B_Age <= 12) & alcohol_use_status == "never_user",
    0,
    average_drinks_per_day)) %>%
  mutate(B_AUDIT = ifelse(
    (B_Age <= 12) & is.na(B_AUDIT),
    0,
    B_AUDIT)) %>%
  mutate(alcohol_disorder = ifelse(
    (B_Age <= 12) & alcohol_use_status == "never_user",
    "Low Risk",
    alcohol_disorder))

# Fill covariates reported at another visit, so that values are missing only
# when they are missing for the individual
covariates <- c("B_Age", "B_Gender", "B_Education", "B_Employment", 
                "B_HH_Income.Categories", "Living.Area", "B_BCG")

complete <- complete %>%
  group_by(PID) %>%
  arrange(VISIT) %>%
  fill(all_of(covariates), .direction = "downup") %>%
  ungroup()

# Code HIV consistently (recorded as yes/no or positive/negative)
complete$HIV <- case_match(complete$HIV,
                           "No" ~ "Negative",
                           "Yes" ~ "Positive",
                           .default = complete$HIV)

# Code education and employment consistently for children
complete <- complete %>%
  mutate(B_Education = case_when(
    B_Age < 4 ~ "NotApplicable",
    (B_Age < 15) & (B_Education == "Illiterate") ~ "NotApplicable",
    TRUE ~ B_Education
  )) %>%
  mutate(B_Employment = ifelse(
    (B_Age < 15) & (B_Employment == "NotApplicable" | B_Employment == "Nonworkingunemployed" | B_Employment == "Other"),
    "NotApplicable",
    B_Employment))

# Days between participant entry and the index case's ATT start
complete$delay <- interval(ymd(complete$Index.PID.ATT.Start.Date), 
                           ymd(complete$Entry.Date))
complete$delay <- complete$delay %/% days(1)

# Hours per week the participant spends with the index case
complete$b08fu0026 <- as.numeric(complete$b08fu0026)

complete <- complete %>%
  mutate(b08fu0027 = case_when(
    b08fu0027 == "<1 hr" ~ 0.5,
    b08fu0027 == "1 to 6 hrs" ~ 3.5,
    b08fu0027 == ">6 but <12 hrs" ~ 9,
    b08fu0027 == "> 12 but <18" ~ 15,
    b08fu0027 == "> 18 hrs" ~ 18,
    b08fu0027 == "Dont know" ~ NA_real_
  ))

complete <- complete %>%
  mutate(contact_hours_per_week = b08fu0026 * b08fu0027)

# Season at each row's visit date (`season`; the confounder used in the
# models) and one month before (`season_lag1`; MICE auxiliary variable only)
complete <- complete %>%
  mutate(
    month_of_visit = month(date),
    season = case_when(
      month_of_visit %in% c(12, 1, 2) ~ "Winter",
      month_of_visit %in% c(3, 4, 5) ~ "Pre-Monsoon",
      month_of_visit %in% c(6, 7, 8, 9) ~ "Monsoon",
      month_of_visit %in% c(10, 11) ~ "Post-Monsoon"),
    season = factor(season),
    month_lag1 = month(date_lag1),
    season_lag1 = case_when(
      month_lag1 %in% c(12, 1, 2) ~ "Winter",
      month_lag1 %in% c(3, 4, 5) ~ "Pre-Monsoon",
      month_lag1 %in% c(6, 7, 8, 9) ~ "Monsoon",
      month_lag1 %in% c(10, 11) ~ "Post-Monsoon"),
    season_lag1 =
      factor(season_lag1))
# Monsoon is the reference level (lowest air pollution)

# Add the ERA5 temperature data
# Centre points for Pune and Chennai
sites <- data.frame(
  Site = c("BJGMC", "NIRT"),
  lat = c(18.52, 13.08),
  lon = c(73.85, 80.20)
)

unzip_files <- unzip("india_temp_era5.zip")
nc_path <- unzip_files[grep("\\.nc$", unzip_files)]
temp_raster <- rast(nc_path)
temp_raster <- temp_raster - 273.15 # Convert to Celsius

# Assign dates
start_date <- as.Date("2013-01-01")
date_seq <- seq(start_date, by = "month", length.out = nlyr(temp_raster))
tail(date_seq, 1)
time(temp_raster) <- date_seq

site_points <- vect(sites,
                    geom = c("lon", "lat"),
                    crs = "EPSG:4326")

city_extracted <- extract(temp_raster, site_points,
                          buffer = 5000,
                          fun = mean,
                          na.rm = TRUE)

temp_lookup <- city_extracted %>%
  mutate(Site = sites$Site) %>%
  select(-ID) %>%
  pivot_longer(-Site, names_to = "layer", values_to = "avg_temp") %>%
  group_by(Site) %>%
  mutate(match_date = date_seq) %>%
  ungroup() %>%
  select(Site, match_date, avg_temp)

complete <- complete %>%
  left_join(temp_lookup, by = c("Site", "date" = "match_date")) %>%
  rename(temp = avg_temp) %>%
  # temp_lag1: temperature one month before the visit
  left_join(temp_lookup, by = c("Site", "date_lag1" = "match_date")) %>%
  rename(temp_lag1 = avg_temp)

## temp_recentavg / temp_intervalavg: mean temperature over the same months ##
## as mean_pm_recentavg / mean_pm_intervalavg (uses recentavg_months and   ##
## intervalavg_months from the PM2.5 section above).                        ##
temp_recentavg_summary <- recentavg_months %>%
  left_join(temp_lookup, by = c("Site.x" = "Site", "target_month" = "match_date")) %>%
  group_by(.row_id) %>%
  summarise(temp_recentavg = mean(avg_temp), .groups = "drop")

temp_intervalavg_summary <- intervalavg_months %>%
  left_join(temp_lookup, by = c("Site.x" = "Site", "target_month" = "match_date")) %>%
  group_by(.row_id) %>%
  summarise(temp_intervalavg = mean(avg_temp), .groups = "drop")

complete <- complete %>%
  left_join(temp_recentavg_summary, by = ".row_id") %>%
  left_join(temp_intervalavg_summary, by = ".row_id") %>%
  dplyr::select(-.row_id)

# interval_pm_by_pid is saved with `complete` for use in
# multiple_imputation.R
save(complete, interval_pm_by_pid, file = here("Data", "data.Rdata"))
