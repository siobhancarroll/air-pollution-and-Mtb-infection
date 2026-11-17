### Dual-assay conversion cohort: counts ########################################
# Cohort: contacts in the analysis dataset with an observed NEGATIVE baseline
# IGRA (QFT-GIT "Negative") AND an observed negative baseline TST (< 5 mm).
#
# Two event definitions, both using observed (pre-imputation) results only:
#   same_visit  Event at the first follow-up visit where IGRA is Positive AND
#               TST >= 5 mm at that same visit.
#   cumulative  Event at the first follow-up visit by which each assay has
#               been positive at least once (the two conversions may occur at
#               different visits).
# Status at a visit is 1 (event), 0 (observed not an event: at least one
# relevant assay observed and negative) or NA (cannot be determined: censored).
# Follow-up stops at the first event or first undeterminable visit, matching
# build_person_period(). Indeterminate follow-up IGRA is treated as missing.
#
# Sanity check: the single-assay cohort sizes printed first should match the
# primary IGRA and TST cohorts (440 and 457 in the 2026-10-22 run).
#
# Output (aggregate only): tables_and_figures/dual-conversion-<date>/

library(tidyverse)
library(here)
source(here("mi_pooling_helpers.R"))   # load_latest()

out_dir <- file.path(here("tables_and_figures"), paste0("dual-conversion-", Sys.Date()))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

load(here("Data", "data.Rdata"))   # `complete`
load_latest("^list_imputed_dataframes_\\d{4}-\\d{2}-\\d{2}\\.Rdata$")
load_latest("^survival_person_period_qgit_\\d{4}-\\d{2}-\\d{2}\\.Rdata$")
analysis_ids <- unique(as.character(allImputations$imputations[[1]]$PID))

visits <- c("Entry", "Month4", "Month12", "Month24")
res <- complete %>%
  mutate(VISIT = str_replace_all(VISIT, " ", ""), PID = as.character(PID)) %>%
  filter(VISIT %in% visits, PID %in% analysis_ids) %>%
  distinct(PID, VISIT, .keep_all = TRUE) %>%
  transmute(PID, VISIT,
            igra = case_when(as.character(QGIT.Result) == "Positive" ~ 1,
                             as.character(QGIT.Result) == "Negative" ~ 0,
                             TRUE ~ NA_real_),   # missing or Indeterminate
            tst = case_when(is.na(TST) ~ NA_real_, TST >= 5 ~ 1, TRUE ~ 0)) %>%
  pivot_wider(id_cols = PID, names_from = VISIT, values_from = c(igra, tst))
for (v in c("igra", "tst")) for (vis in visits) {
  col <- paste0(v, "_", vis); if (!col %in% names(res)) res[[col]] <- NA_real_
}
res <- tibble(PID = analysis_ids) %>% left_join(res, by = "PID")

cat("\n== Single-assay cohorts (check against primary analysis) ==\n")
cat("Baseline IGRA negative:", sum(res$igra_Entry == 0, na.rm = TRUE), "\n")
cat("Baseline TST negative: ", sum(res$tst_Entry == 0, na.rm = TRUE), "\n")

dual <- res %>% filter(igra_Entry == 0, tst_Entry == 0)
cat("\n== Dual-negative cohort ==\nBaseline negative on both:", nrow(dual), "\n")

# Long format over follow-up intervals 1-3 (Month4, Month12, Month24)
long <- dual %>%
  pivot_longer(-PID, names_to = c(".value", "VISIT"), names_sep = "_") %>%
  filter(VISIT != "Entry") %>%
  mutate(interval = match(VISIT, visits) - 1) %>%
  arrange(PID, interval) %>%
  group_by(PID) %>%
  mutate(
    status_same = case_when(igra == 1 & tst == 1 ~ 1,
                            igra == 0 | tst == 0 ~ 0,
                            TRUE ~ NA_real_),
    # positive at some observed visit up to and including this one
    igra_ever = cummax(coalesce(igra, 0)) == 1,
    tst_ever = cummax(coalesce(tst, 0)) == 1,
    # 0 if an assay not yet positive is observed negative at this visit
    status_cum = case_when(igra_ever & tst_ever ~ 1,
                           (!igra_ever & igra %in% 0) | (!tst_ever & tst %in% 0) ~ 0,
                           TRUE ~ NA_real_)
  ) %>% ungroup()

truncate_pp <- function(df, status_col) {
  df %>% mutate(event = .data[[status_col]]) %>%
    group_by(PID) %>%
    mutate(first_event = suppressWarnings(min(interval[event %in% 1])),
           first_na = suppressWarnings(min(interval[is.na(event)])),
           keep_upto = case_when(is.finite(first_event) & (!is.finite(first_na) | first_event <= first_na) ~ first_event,
                                 is.finite(first_na) ~ first_na,
                                 TRUE ~ max(interval))) %>%
    filter(interval <= keep_upto) %>%
    ungroup() %>%
    dplyr::select(PID, interval, VISIT, igra, tst, event)
}

# PM2.5 availability per person-interval, from the primary IGRA person-period
# data (the dual cohort is a subset of the IGRA cohort; PM is assay-independent)
pm_avail <- allImputations_survival_qgit$imputations[[1]] %>%
  transmute(PID = as.character(PID), interval,
            pm_ok = !is.na(mean_pm_intervalavg) & !is.na(SE_pm_intervalavg)) %>%
  distinct(PID, interval, .keep_all = TRUE)

summarise_pp <- function(pp, def) {
  pp <- pp %>% left_join(pm_avail, by = c("PID", "interval"))
  # PM models drop person-intervals from the first one without PM onward
  pp <- pp %>% group_by(PID) %>%
    mutate(pm_kept = cumall(coalesce(pm_ok, FALSE))) %>% ungroup()
  pp %>% filter(interval <= 2) %>% group_by(interval) %>%
    summarise(at_risk = n_distinct(PID),
              events = sum(event %in% 1),
              censored = sum(is.na(event)),
              at_risk_pm_model = n_distinct(PID[pm_kept]),
              events_pm_model = sum(event %in% 1 & pm_kept),
              .groups = "drop") %>%
    mutate(definition = def, .before = 1)
}

pp_same <- truncate_pp(long, "status_same")
pp_cum <- truncate_pp(long, "status_cum")
counts <- bind_rows(summarise_pp(pp_same, "same_visit"), summarise_pp(pp_cum, "cumulative"))
cat("\n== Dual conversion: contacts at risk and events by interval ==\n",
    "(interval 1 = Entry to Month 4, interval 2 = Month 4 to Month 12)\n")
print(counts)

# Patterns at Month 4 among those with both assays observed
pattern_m4 <- dual %>%
  filter(!is.na(igra_Month4), !is.na(tst_Month4)) %>%
  count(IGRA_Month4 = ifelse(igra_Month4 == 1, "Pos", "Neg"),
        TST_Month4 = ifelse(tst_Month4 == 1, ">=5 mm", "<5 mm"))
cat("\n== Month 4 results among dual-negative contacts with both assays observed ==\n")
print(pattern_m4)
followup_any <- dual %>% summarise(
  no_followup_either = sum(is.na(igra_Month4) & is.na(tst_Month4) & is.na(igra_Month12) &
                             is.na(tst_Month12) & is.na(igra_Month24) & is.na(tst_Month24)))
print(followup_any)

write_csv(tibble(item = c("analysis_dataset", "baseline_igra_negative", "baseline_tst_negative",
                          "baseline_both_negative", "both_negative_no_followup"),
                 n = c(length(analysis_ids), sum(res$igra_Entry == 0, na.rm = TRUE),
                       sum(res$tst_Entry == 0, na.rm = TRUE), nrow(dual),
                       followup_any$no_followup_either)),
          file.path(out_dir, "dual_cohort_sizes.csv"))
write_csv(counts, file.path(out_dir, "dual_conversion_counts.csv"))
write_csv(pattern_m4, file.path(out_dir, "dual_month4_patterns.csv"))
message("Saved to ", out_dir)
