### Month 6 IGRA/TST results among baseline-negative contacts ##################
# Do some contacts have an IGRA or TST result at Month 6 but not at Month 4?
# If so, they could contribute to the conversion analysis. Writes aggregate
# counts only: tables_and_figures/month6-check-<date>/month6_tests.csv
library(tidyverse)
library(here)

load(here("Data", "data.Rdata"))
visits <- c("Entry", "Month4", "Month6", "Month12")
d <- complete %>%
  mutate(VISIT = gsub(" ", "", as.character(VISIT))) %>%
  filter(VISIT %in% visits) %>%
  mutate(VISIT = factor(VISIT, levels = visits)) %>%
  arrange(PID, VISIT) %>%
  distinct(PID, VISIT, .keep_all = TRUE) %>%
  transmute(PID, VISIT, igra = tolower(trimws(as.character(QGIT.Result))),
            tst = suppressWarnings(as.numeric(as.character(TST))))

cat("Recorded IGRA result values by visit:\n")
print(table(d$VISIT, d$igra, useNA = "ifany"))
cat("\nTST recorded (non-missing) by visit:\n")
print(table(d$VISIT, !is.na(d$tst)))

w <- d %>% pivot_wider(id_cols = PID, names_from = VISIT, values_from = c(igra, tst))
for (v in visits) for (a in c("igra", "tst")) {
  nm <- paste0(a, "_", v); if (!nm %in% names(w)) w[[nm]] <- if (a == "igra") NA_character_ else NA_real_
}
igra_obs <- function(x) grepl("^neg|^pos", x)
igra_pos <- function(x) grepl("^pos", x)
tst_obs  <- function(x) !is.na(x)
tst_pos  <- function(x) !is.na(x) & x >= 5

summ <- function(assay) {
  obs <- if (assay == "igra") igra_obs else tst_obs
  pos <- if (assay == "igra") igra_pos else tst_pos
  g <- function(v) w[[paste0(assay, "_", v)]]
  cohort <- obs(g("Entry")) & !pos(g("Entry"))
  m4 <- obs(g("Month4")); m6 <- obs(g("Month6")); m12 <- obs(g("Month12"))
  tibble(
    assay = toupper(assay),
    baseline_negative = sum(cohort),
    month4_observed = sum(cohort & m4),
    month6_observed = sum(cohort & m6),
    month12_observed = sum(cohort & m12),
    month6_and_month4 = sum(cohort & m6 & m4),
    month6_without_month4 = sum(cohort & m6 & !m4),
    month6_without_month4_positive = sum(cohort & m6 & !m4 & pos(g("Month6"))),
    month6_without_month4_or_month12 = sum(cohort & m6 & !m4 & !m12),
    no_month4_result = sum(cohort & !m4)
  )
}
out <- bind_rows(summ("igra"), summ("tst"))
cat("\nBaseline-negative cohorts (first record per visit):\n")
print(as.data.frame(out))

od <- here("tables_and_figures", paste0("month6-check-", Sys.Date()))
dir.create(od, recursive = TRUE, showWarnings = FALSE)
write_csv(out, file.path(od, "month6_tests.csv"))
message("Saved ", file.path(od, "month6_tests.csv"))
