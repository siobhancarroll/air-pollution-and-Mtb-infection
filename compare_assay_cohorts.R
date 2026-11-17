### Compare contacts who received each infection test at baseline ##########
# Describes differences between contacts who received an IGRA (QFT-GIT) and
# those who received a TST at the Entry visit, regardless of result.
#
# Uses observed Entry-visit values only (no imputed values), with missing
# counts shown.
#
# Population: contacts in the analysis dataset, plus those with an
# indeterminate baseline QFT-GIT (tested, but excluded from the analyses).
#
# Received = an observed Entry result: for the IGRA, a QFT-GIT result
# (positive, negative or indeterminate) or antigen-minus-nil value; for the
# TST, a recorded induration.
#
# Groups: Received IGRA | Received TST | Both | IGRA only | TST only | Neither
# SMDs compare everyone who received each test (overlapping groups; descriptive
# only) and IGRA-only vs TST-only (non-overlapping). |SMD| > 0.1 is flagged.
#
# Outputs (aggregate only), tables_and_figures/assay-tested-<date>/:
#   tested_baseline_comparison.csv - main table
#   baseline_assay_crosstab.csv    - baseline IGRA x TST, contacts tested with both
#   baseline_assay_agreement.csv   - percent agreement and kappa (IGRA vs TST >=5 mm)
#   assay_results_by_site.csv      - baseline positivity and agreement by site
#   conversions_by_site.csv        - conversions by assay, site and interval
#   assay_tested_tables.docx       - the same tables, formatted

library(tidyverse)
library(here)
library(mitools)
library(flextable)

source(here("mi_pooling_helpers.R"))   # load_latest(), collapse_covariates()

out_dir <- file.path(here("tables_and_figures"), paste0("assay-tested-", get0("output_stamp", ifnotfound = as.character(Sys.Date()))))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

### Data #######################################################################
load(here("Data", "data.Rdata"))  # `complete`: long format, one row per visit
load_latest("^list_imputed_dataframes_\\d{4}-\\d{2}-\\d{2}\\.Rdata$")

ids_analysis <- unique(as.character(allImputations$imputations[[1]]$PID))

# Same string clean-up as multiple_imputation.R, so that collapse_covariates()
# recognises the levels
clean <- function(df) {
  df %>% mutate(
    VISIT = str_replace_all(VISIT, " ", ""),
    Living.Area = str_replace_all(Living.Area, "[ /-]", ""),
    B_Sleepindex = str_replace_all(B_Sleepindex, "[ /]", ""),
    B_Education = B_Education %>% str_replace_all(" ", "") %>%
      str_replace_all("[()]", "") %>% str_replace_all("\\+", "plus"),
    B_Employment = str_replace_all(B_Employment, "[ /-]", ""),
    B_HH_Income.Categories = B_HH_Income.Categories %>% str_replace_all(",", "") %>%
      str_replace_all(" – ", "to") %>% str_replace_all("≥ ", "morethan")
  )
}
complete_c <- clean(complete)

month4 <- complete_c %>%
  filter(VISIT == "Month4") %>%
  distinct(PID, .keep_all = TRUE) %>%
  transmute(PID = as.character(PID), contact_hours_per_week_Month4 = contact_hours_per_week)

base <- complete_c %>%
  filter(VISIT == "Entry") %>%
  distinct(PID, .keep_all = TRUE) %>%
  mutate(PID = as.character(PID)) %>%
  filter(PID %in% ids_analysis | QGIT.Result %in% "Indeterminate") %>%
  left_join(month4, by = "PID") %>%
  mutate(across(c(B_Education, B_Employment, B_HH_Income.Categories), as.factor)) %>%
  collapse_covariates() %>%
  mutate(
    Income_collapsed = fct_recode(Income_collapsed, `<10,000` = "Low", `10,000–19,999` = "Medium", `≥20,000` = "High"),
    Education_collapsed = fct_recode(Education_collapsed, `Primary school or less` = "Low",
                                     `High school or junior college` = "Medium", `College or higher` = "High"),
    Employment_collapsed = fct_recode(Employment_collapsed, `Not in labour force` = "Not_in_Labor",
                                      `Manual labour/trade` = "Manual_Trade", Professional = "Professional"),
    Income_collapsed = fct_relevel(Income_collapsed, "<10,000", "10,000–19,999", "≥20,000"),
    Education_collapsed = fct_relevel(Education_collapsed, "Primary school or less",
                                      "High school or junior college", "College or higher"),
    Employment_collapsed = fct_relevel(Employment_collapsed, "Not in labour force",
                                       "Manual labour/trade", "Professional"),
    child = if_else(B_Age < 15, "Child (<15 years)", "Adult (≥15 years)"),
    BMI_adult = if_else(B_Age >= 15, BMI, NA_real_),
    # fuel flags are coded "YES"/"No"; compared case-insensitively
    solid_fuel_use = case_when(toupper(exclusive_clean) == "YES" ~ "No", toupper(exclusive_clean) == "NO" ~ "Yes"),
    across(c(wood_use, kerosene_use, solid_fuel_primary),
           ~ case_when(toupper(.x) == "YES" ~ "Yes", toupper(.x) == "NO" ~ "No")),
    igra_entry = as.character(QGIT.Result),
    amn_entry = TB.Antigen...Nil,
    tst_entry = TST,
    tst_entry_cat = cut(TST, c(-Inf, 4.999, 9.999, Inf), labels = c("<5 mm", "5–9 mm", "≥10 mm")),
    got_igra = !is.na(QGIT.Result) | !is.na(TB.Antigen...Nil),
    got_tst = !is.na(TST)
  )
message("N = ", nrow(base), " contacts; received IGRA = ", sum(base$got_igra),
        "; received TST = ", sum(base$got_tst), "; both = ", sum(base$got_igra & base$got_tst),
        "; neither = ", sum(!base$got_igra & !base$got_tst))

groups <- list(
  `Received IGRA` = base %>% filter(got_igra),
  `Received TST`  = base %>% filter(got_tst),
  `Both`          = base %>% filter(got_igra & got_tst),
  `IGRA only`     = base %>% filter(got_igra & !got_tst),
  `TST only`      = base %>% filter(got_tst & !got_igra),
  `Neither`       = base %>% filter(!got_igra & !got_tst)
)

### Fuel coding check ###########################################################
# Tabulates the raw secondary-fuel codes at Entry, to confirm that the
# "YES"/"No" coding assumed by wood_use, kerosene_use and exclusive_clean in
# data_processing.R matches the data.
fuel_codes <- complete_c %>% filter(VISIT == "Entry") %>%
  select(starts_with("secondary_fuel_"), primary_fuel) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
  count(variable, value)
write_csv(fuel_codes, file.path(out_dir, "fuel_coding_check.csv"))
print(fuel_codes, n = Inf)

### Variables ##################################################################
vars <- tribble(
  ~var,                            ~label,                                        ~type,
  "B_Age",                         "Age (years)",                                 "mean",
  "child",                         "Age group",                                   "cat",
  "B_Gender",                      "Sex",                                         "cat",
  "Site",                          "Study site",                                  "cat",
  "Income_collapsed",              "Household income (INR/month)",                "cat",
  "Education_collapsed",           "Education",                                   "cat",
  "Employment_collapsed",          "Employment",                                  "cat",
  "Living.Area",                   "Residential location",                        "cat",
  "Type.of.House",                 "House type",                                  "cat",
  "Household.size",                "Household size",                              "mean",
  "contact_hours_per_week_Month4", "Contact hours/week with index case",          "mean",
  "B_Sleepindex",                  "Sleeping arrangement (relative to index case)", "cat",
  "BMI_adult",                     "BMI, adults (kg/m²)",                         "mean",
  "B_BCG",                         "BCG-vaccinated",                              "cat",
  "B_HbA1c",                       "HbA1c (%)",                                   "mean",
  "smoking_status",                "Smoking status",                              "cat",
  "alcohol_disorder",              "Alcohol use (AUDIT category)",                "cat",
  "solid_fuel_use",                "Any polluting fuel use",                          "cat",
  "wood_use",                      "Any wood use",                                "cat",
  "kerosene_use",                  "Any kerosene use",                            "cat",
  "solid_fuel_primary",            "Polluting fuel as primary cooking fuel",          "cat",
  "mean_pm_recentavg",             "Ambient PM2.5, 3 months before entry (µg/m³)", "mean",
  "igra_entry",                    "Baseline IGRA (QFT-GIT) result",              "cat",
  "amn_entry",                     "Baseline TB antigen minus nil (IU/mL)",       "median",
  "tst_entry",                     "Baseline TST induration (mm)",                "median",
  "tst_entry_cat",                 "Baseline TST induration",                     "cat"
)

### Summary helpers ############################################################
fmt_cont <- function(x, type) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("–")
  if (type == "mean") sprintf("%.1f (%.1f)", mean(x), sd(x))
  else sprintf("%.2f [%.2f, %.2f]", median(x), quantile(x, .25), quantile(x, .75))
}

smd_cont <- function(a, b) {
  a <- a[!is.na(a)]; b <- b[!is.na(b)]
  if (length(a) < 2 || length(b) < 2) return(NA_real_)
  (mean(a) - mean(b)) / sqrt((var(a) + var(b)) / 2)
}

# Multi-category SMD (Yang & Dalton 2012); equals the usual binary SMD for 2 levels
smd_cat <- function(a, b) {
  a <- as.character(a[!is.na(a)]); b <- as.character(b[!is.na(b)])
  lev <- sort(union(a, b))
  if (length(lev) < 2 || !length(a) || !length(b)) return(NA_real_)
  p1 <- as.numeric(table(factor(a, lev))) / length(a)
  p2 <- as.numeric(table(factor(b, lev))) / length(b)
  d <- (p1 - p2)[-1]
  q1 <- p1[-1]; q2 <- p2[-1]
  S <- (diag(q1 * (1 - q1), length(q1)) - outer(q1, q1) * (1 - diag(length(q1))) +
        diag(q2 * (1 - q2), length(q2)) - outer(q2, q2) * (1 - diag(length(q2)))) / 2
  out <- tryCatch(sqrt(as.numeric(t(d) %*% solve(S) %*% d)), error = function(e) NA_real_)
  out
}

rows <- list()
for (i in seq_len(nrow(vars))) {
  v <- vars$var[i]; lab <- vars$label[i]; type <- vars$type[i]
  if (!v %in% names(base)) { message("Skipping missing variable: ", v); next }
  smd_full <- if (type == "cat") smd_cat(groups[[1]][[v]], groups[[2]][[v]]) else smd_cont(groups[[1]][[v]], groups[[2]][[v]])
  smd_excl <- if (type == "cat") smd_cat(groups[[4]][[v]], groups[[5]][[v]]) else smd_cont(groups[[4]][[v]], groups[[5]][[v]])
  n_missing <- sapply(groups, function(g) sum(is.na(g[[v]])))
  miss_txt <- if (any(n_missing > 0)) sprintf(" (missing: %s)", paste(n_missing, collapse = "/")) else ""
  if (type == "cat") {
    rows[[length(rows) + 1]] <- tibble(Characteristic = paste0(lab, miss_txt), Level = "",
      !!!setNames(as.list(rep("", length(groups))), names(groups)),
      `SMD, received IGRA vs received TST` = sprintf("%.2f", abs(smd_full)),
      `SMD, IGRA only vs TST only` = sprintf("%.2f", abs(smd_excl)))
    lev <- if (is.factor(base[[v]])) levels(droplevels(base[[v]])) else sort(unique(na.omit(as.character(base[[v]]))))
    for (l in lev) {
      cells <- sapply(groups, function(g) {
        x <- as.character(g[[v]]); nn <- sum(!is.na(x)); k <- sum(x == l, na.rm = TRUE)
        sprintf("%d (%.0f%%)", k, if (nn) 100 * k / nn else 0)
      })
      rows[[length(rows) + 1]] <- tibble(Characteristic = "", Level = l,
        !!!as.list(cells), `SMD, received IGRA vs received TST` = "", `SMD, IGRA only vs TST only` = "")
    }
  } else {
    cells <- sapply(groups, function(g) fmt_cont(g[[v]], type))
    stat <- if (type == "mean") "mean (SD)" else "median [IQR]"
    rows[[length(rows) + 1]] <- tibble(Characteristic = paste0(lab, ", ", stat, miss_txt), Level = "",
      !!!as.list(cells),
      `SMD, received IGRA vs received TST` = sprintf("%.2f", abs(smd_full)),
      `SMD, IGRA only vs TST only` = sprintf("%.2f", abs(smd_excl)))
  }
}
n_row <- tibble(Characteristic = "N", Level = "",
                !!!as.list(sapply(groups, function(g) as.character(nrow(g)))),
                `SMD, received IGRA vs received TST` = "", `SMD, IGRA only vs TST only` = "")
comparison <- bind_rows(n_row, rows) %>% mutate(across(everything(), ~ replace_na(.x, "")),
                                                across(starts_with("SMD"), ~ if_else(.x == "NA", "", .x)))
write_csv(comparison, file.path(out_dir, "tested_baseline_comparison.csv"))
print(comparison, n = Inf, width = Inf)

### Baseline agreement between assays (everyone with both results) ###########
both_tested <- base %>% filter(!is.na(igra_entry), !is.na(tst_entry_cat))
crosstab <- both_tested %>% count(igra_entry, tst_entry_cat) %>%
  pivot_wider(names_from = tst_entry_cat, values_from = n, values_fill = 0)
write_csv(crosstab, file.path(out_dir, "baseline_assay_crosstab.csv"))

kap_df <- both_tested %>% filter(igra_entry %in% c("Positive", "Negative")) %>%
  transmute(igra = igra_entry == "Positive", tst = tst_entry >= 5)
po <- mean(kap_df$igra == kap_df$tst)
pe <- mean(kap_df$igra) * mean(kap_df$tst) + mean(!kap_df$igra) * mean(!kap_df$tst)
kappa <- (po - pe) / (1 - pe)
se_k <- sqrt(po * (1 - po) / (nrow(kap_df) * (1 - pe)^2))
agree <- tibble(n = nrow(kap_df), percent_agreement = round(100 * po, 1),
                kappa = round(kappa, 2), kappa_ci_low = round(kappa - 1.96 * se_k, 2),
                kappa_ci_high = round(kappa + 1.96 * se_k, 2),
                igra_pos_tst_neg = sum(kap_df$igra & !kap_df$tst),
                igra_neg_tst_pos = sum(!kap_df$igra & kap_df$tst))
write_csv(agree, file.path(out_dir, "baseline_assay_agreement.csv"))
print(crosstab); print(agree)

### Assay results by site ######################################################
# Baseline positivity and IGRA-TST agreement by site, and conversion by site
# and interval in the baseline-negative cohorts.
kappa_fun <- function(igra_pos, tst_pos) {
  po <- mean(igra_pos == tst_pos)
  pe <- mean(igra_pos) * mean(tst_pos) + mean(!igra_pos) * mean(!tst_pos)
  (po - pe) / (1 - pe)
}
site_label <- function(x) dplyr::recode(as.character(x), BJGMC = "Pune (BJGMC)", NIRT = "Chennai (NIRT)")
by_site_baseline <- bind_rows(base %>% mutate(age_group = "All ages"),
                              base %>% mutate(age_group = if_else(B_Age < 15, "Children", "Adults"))) %>%
  group_by(Site = site_label(Site), age_group) %>%
  summarise(
    n_igra = sum(igra_entry %in% c("Positive", "Negative")),
    igra_positive_pct = round(100 * mean(igra_entry[igra_entry %in% c("Positive", "Negative")] == "Positive"), 1),
    n_tst = sum(!is.na(tst_entry)),
    tst_5mm_pct = round(100 * mean(tst_entry[!is.na(tst_entry)] >= 5), 1),
    tst_10mm_pct = round(100 * mean(tst_entry[!is.na(tst_entry)] >= 10), 1),
    n_both = sum(igra_entry %in% c("Positive", "Negative") & !is.na(tst_entry)),
    agreement_pct = {k <- igra_entry %in% c("Positive", "Negative") & !is.na(tst_entry)
                     round(100 * mean((igra_entry[k] == "Positive") == (tst_entry[k] >= 5)), 1)},
    kappa = {k <- igra_entry %in% c("Positive", "Negative") & !is.na(tst_entry)
             round(kappa_fun(igra_entry[k] == "Positive", tst_entry[k] >= 5), 2)},
    igra_pos_tst_neg = sum(igra_entry %in% "Positive" & tst_entry < 5, na.rm = TRUE),
    igra_neg_tst_pos = sum(igra_entry %in% "Negative" & tst_entry >= 5, na.rm = TRUE),
    .groups = "drop") %>%
  arrange(age_group, Site)
write_csv(by_site_baseline, file.path(out_dir, "assay_results_by_site.csv"))
print(by_site_baseline, width = Inf)

load_latest("^survival_person_period_qgit_\\d{4}-\\d{2}-\\d{2}\\.Rdata$")
load_latest("^survival_person_period_tst_\\d{4}-\\d{2}-\\d{2}\\.Rdata$")
conv_by_site <- bind_rows(
  allImputations_survival_qgit$imputations[[1]] %>% mutate(assay = "IGRA"),
  allImputations_survival_tst$imputations[[1]] %>% mutate(assay = "TST")) %>%
  mutate(Site = site_label(Site)) %>%
  group_by(assay, Site, interval) %>%
  summarise(at_risk = n_distinct(PID), events = sum(event == 1, na.rm = TRUE),
            censored = sum(is.na(event)), .groups = "drop") %>%
  mutate(events_per_100_observed = round(100 * events / (at_risk - censored), 1))
write_csv(conv_by_site, file.path(out_dir, "conversions_by_site.csv"))
print(conv_by_site)

### Word tables ################################################################
ft <- function(df, caption) flextable(df) %>% set_caption(caption) %>% autofit() %>% fontsize(size = 9, part = "all")
save_as_docx(
  `Table A` = ft(comparison, "Baseline characteristics of contacts who received each infection test at entry (observed values)"),
  `Table B` = ft(crosstab, "Baseline IGRA result by baseline TST induration, contacts tested with both assays"),
  `Table C` = ft(agree, "Agreement between baseline IGRA and TST (≥5 mm)"),
  `Table D` = ft(by_site_baseline, "Baseline assay positivity and agreement by study site"),
  `Table E` = ft(conv_by_site, "Conversions by assay, site and follow-up interval (baseline-negative cohorts)"),
  path = file.path(out_dir, "assay_tested_tables.docx"))
message("Saved to ", out_dir)
