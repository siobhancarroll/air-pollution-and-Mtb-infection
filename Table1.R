library(tidyverse)
library(here)
library(mice)
library(gtsummary)
library(gt)

load(here("Data", "long_imputed_dataframes_2026-03-12.Rdata"))
load(here("Data", "list_imputed_dataframes_2026-03-12.Rdata"))

fp <- file.path(here("tables_and_figures")) # specify file path to save tables and figures
dir.create(fp, recursive = TRUE) # create new folder for plots

mids_obj <- as.mids(impdata)

impdata <- impdata %>%
  mutate(average_drinks_per_day = case_when(
    average_drinks_per_day == "0" ~ "0",
    average_drinks_per_day == "1" ~ "1",
    average_drinks_per_day %in% c("3", "5", "8") ~ "2 or more"
  )) %>%
  mutate(Site = case_when(
    Site == "BJGMC" ~ "BJGMC (Pune)",
    Site == "NIRT" ~ "NIRT (Chennai)"
  ),
  B_HH_Income.Categories = case_when(
    B_HH_Income.Categories == "1000to4999" ~ "1,000 - 4,999",
    B_HH_Income.Categories == "5000to9999" ~ "5,000 - 9,999",
    B_HH_Income.Categories == "10000to14999" ~ "10,000 - 14,999",
    B_HH_Income.Categories == "15000to19999" ~ "15,000 - 19,999",
    B_HH_Income.Categories == "20000to24999" ~ "20,000 - 24,999",
    B_HH_Income.Categories == "25000to29999" ~ "25,000 - 29,999",
    B_HH_Income.Categories == "morethan30000" ~ ">30,000"
  ),
  B_Education = case_when(
    B_Education == "NotApplicable" ~ "Not Applicable (Pre-School Age)",
    B_Education == "Illiterate" ~ "Illiterate",
    B_Education == "PrimarySchool0to5years" ~ "Primary School",
    B_Education == "HighSchool6to10years" ~ "High School",
    B_Education == "JuniorCollege11to12years" ~ "Junior College",
    B_Education == "College12plusyears" ~ "College",
    B_Education == "Professional" ~ "Professional School"
  ),
  B_Employment = case_when(
    B_Employment == "Housewife" ~ "Housewife",
    B_Employment == "Nonworkingunemployed" ~ "Non-Working/Unemployed",
    B_Employment == "NotApplicable" ~ "Not Applicable (Child)",
    B_Employment == "Other" ~ "Other",
    B_Employment == "Professional" ~ "Professional",
    B_Employment == "Semiprofessionclerk" ~ "Semi-Professional/Clerk",
    B_Employment == "Skilledmanualworkerfarmer" ~ "Skilled Manual Worker/Farmer",
    B_Employment == "Tradesbusinessshopowner" ~ "Trades/Businessowner",
    B_Employment == "Unskilledmanualworker" ~ "Unskilled Manual Worker"
  ),
  Living.Area = case_when(
    Living.Area == "PeriUrbanTown" ~ "Peri-Urban Town",
    Living.Area == "Rural" ~ "Rural",
    Living.Area == "UrbanLargeCity" ~ "Urban (Large City)",
    Living.Area == "UrbanSmallCity" ~ "Urban (Small City)"
  ),
  solid_fuel_use = case_when(
    solid_fuel_use == 0 ~ "Exclusive Clean Fuel Use",
    solid_fuel_use == 1 ~ "Any Solid Fuel Use"
  )) 

impdata$B_HH_Income.Categories <- factor(impdata$B_HH_Income.Categories,
                                         levels = c("1,000 - 4,999",
                                                    "5,000 - 9,999",
                                                    "10,000 - 14,999",
                                                    "15,000 - 19,999",
                                                    "20,000 - 24,999",
                                                    "25,000 - 29,999",
                                                    ">30,000"),
                                         ordered = TRUE)

impdata$B_Education <- factor(impdata$B_Education,
                                         levels = c("Not Applicable (Pre-School Age)",
                                                    "Illiterate",
                                                    "Primary School",
                                                    "High School",
                                                    "Junior College",
                                                    "College",
                                                    "Professional School"),
                                         ordered = TRUE)

impdata$B_Employment <- factor(impdata$B_Employment,
                                         levels = c("Not Applicable (Child)",
                                                    "Non-Working/Unemployed",
                                                    "Unskilled Manual Worker",
                                                    "Housewife",
                                                    "Skilled Manual Worker/Farmer",
                                                    "Trades/Businessowner",
                                                    "Semi-Professional/Clerk",
                                                    "Professional",
                                                    "Other"),
                                         ordered = TRUE)

impdata$solid_fuel_use <- factor(impdata$solid_fuel_use,
                                         levels = c("Exclusive Clean Fuel Use",
                                                    "Any Solid Fuel Use"),
                                         ordered = TRUE)

impdata1 <- impdata %>% filter(.imp == 1)

quantile_25 <- quantile(impdata1$mean_pm_lag3, probs = 0.25)
quantile_50 <- quantile(impdata1$mean_pm_lag3, probs = 0.5)
quantile_75 <- quantile(impdata1$mean_pm_lag3, probs = 0.75)

impdata <- impdata %>%
  mutate(pm_quartile = case_when(
    mean_pm_lag3 < quantile_25 ~ "Q1",
    mean_pm_lag3 >= quantile_25 & mean_pm_lag3 < quantile_50 ~ "Q2",
    mean_pm_lag3 >= quantile_50 & mean_pm_lag3 < quantile_75 ~ "Q3",
    mean_pm_lag3 >= quantile_75 ~ "Q4"
  ))

continuous_vars <- c("B_Age", "B_Height", "B_Weight", "MUAC", "mean_pm_lag3",
                     "antigen_minus_nil_Entry", "contact_hours_per_week_Month4")

categorical_vars <- c("Site", "B_Gender", "B_HH_Income.Categories",
                      "B_Education", "B_Employment", "Living.Area",
                      "average_drinks_per_day")

table1_hap <- impdata %>%
  filter(.imp == 1) %>%
  mutate(across(all_of(categorical_vars), as.factor)) %>%
  tbl_summary(
    by = solid_fuel_use,
    include = c(Site, B_Gender, B_Age, B_HH_Income.Categories, B_Education,
                B_Employment, B_Height, B_Weight, MUAC, Living.Area, 
                average_drinks_per_day, contact_hours_per_week_Month4,
                mean_pm_lag3, antigen_minus_nil_Entry),
    label = list(
      Site ~ "Study Site",
      B_Age ~ "Age (years)",
      B_Gender ~ "Sex",
      B_HH_Income.Categories ~ "Household Income (annual, in USD)",
      B_Education ~ "Highest education",
      B_Employment ~ "Occupation",
      Living.Area ~ "Residential Location",
      average_drinks_per_day ~ "Average alcoholic drinks per day",
      B_Height ~ "Height (cm)",
      B_Weight ~ "Weight (kg)",
      MUAC ~ "Mid-upper arm circumference (cm)",
      mean_pm_lag3 ~ "Citywide monthly average PM2.5 3 months before infection testing (ug/m3)",
      antigen_minus_nil_Entry ~ "QGIT antigen-minus-nil at baseline",
      contact_hours_per_week_Month4 ~ "Contact hours per week with TB index case"
    ),
    statistic = list(
      all_continuous()  ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    missing = "no"
  ) %>%
  add_overall() %>%
  bold_labels() 

table1_hap

table1_hap %>%
  as_gt() %>%
  tab_options(
    table.width = pct(100), 
    data_row.padding = px(2) 
  ) %>%
  gtsave(file.path(fp, "Table1_HAP.docx"))

table1_oap <- impdata %>%
  filter(.imp == 1) %>%
  mutate(across(all_of(categorical_vars), as.factor)) %>%
  tbl_summary(
    by = pm_quartile,
    include = c(Site, B_Gender, B_Age, B_HH_Income.Categories, B_Education,
                B_Employment, B_Height, B_Weight, MUAC, Living.Area, 
                average_drinks_per_day, contact_hours_per_week_Month4,
                solid_fuel_use, antigen_minus_nil_Entry),
    label = list(
      Site ~ "Study Site",
      B_Age ~ "Age (years)",
      B_Gender ~ "Sex",
      B_HH_Income.Categories ~ "Household Income (annual, in USD)",
      B_Education ~ "Highest education",
      B_Employment ~ "Occupation",
      Living.Area ~ "Residential Location",
      average_drinks_per_day ~ "Average alcoholic drinks per day",
      B_Height ~ "Height (cm)",
      B_Weight ~ "Weight (kg)",
      MUAC ~ "Mid-upper arm circumference (cm)",
      solid_fuel_use ~ "Household Fuel Use",
      antigen_minus_nil_Entry ~ "QGIT antigen-minus-nil at baseline",
      contact_hours_per_week_Month4 ~ "Contact hours per week with TB index case"
    ),
    statistic = list(
      all_continuous()  ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    missing = "no"
  ) %>%
  add_overall() %>%
  bold_labels() 

table1_oap

table1_oap %>%
  as_gt() %>%
  tab_options(
    table.width = pct(100), 
    data_row.padding = px(2) 
  ) %>%
  gtsave(file.path(fp, "Table1_OAP.docx"))

raw_df <- complete %>%
  filter(VISIT == "Entry")

raw_data_summary <- raw_df %>%
  select(B_HH_Income.Categories, B_Education, B_Employment, 
         Living.Area, B_Age, B_Gender, Site, temp_lag3, 
         average_drinks_per_day, B_Height, B_Weight, MUAC, 
         HbA1c, B_BCG, B_Sleepindex, Type.of.House, 
         contact_hours_per_week_Month4, mean_pm_lag3) %>%
  tbl_summary(
    missing_text = "(Missing)"
  ) %>%
  modify_caption("Original missing data frequency prior to Multiple Imputation.")
