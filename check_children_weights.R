### Diagnostics for the baseline PM2.5 weights in the children subgroup ######
# Checks whether unstable baseline PM2.5 estimates among children arise from
# the generalized propensity score (GPS) weights (extreme weights, poor
# balance, limited overlap) or are also present without weighting.
#
# Follows the continuous-PM2.5 model in propensity_scores.R (stabilized GPS
# from weightit(method = "glm"), 3-df natural spline, weighted logistic
# outcome model, observed outcomes only), using the observed city-month PM2.5
# mean instead of exposure draws. Point estimates are averaged across
# imputations; no confidence intervals are computed.
#
# Covariate lists are copied from propensity_scores.R and must be kept in sync.
#
# Outputs (aggregate only), tables_and_figures/children-weights-<date>/:
#   weight_summary.csv          - weight distribution and ESS by population
#   rd_by_method.csv            - RD and RR at 60 and 80 vs 40 ug/m3, by method
#   balance_children.csv        - covariate-exposure correlations before/after weighting
#   balance_overall.csv         - same, overall population
#   pm_by_site_children.csv     - PM2.5 distribution by site among children
#   prevalence_by_pm_band.csv   - observed positivity by PM2.5 band and site
#   top_weights_children.csv    - characteristics of the largest 5% of weights
#   weights_histogram.jpg       - weight distributions by population

library(tidyverse)
library(here)
library(mitools)
library(WeightIt)
library(cobalt)
library(splines)

source(here("mi_pooling_helpers.R"))   # load_latest(), collapse_covariates()
set.seed(123)

# PM2.5 models are adjusted for household fuel use unless
# adjust_pm_for_fuel <- FALSE is set beforehand.
if (!exists("adjust_pm_for_fuel")) adjust_pm_for_fuel <- TRUE
pm_fuel_adjustment <- if (adjust_pm_for_fuel) "solid_fuel_use" else character(0)

out_dir <- file.path(here("tables_and_figures"), paste0("children-weights-", if (adjust_pm_for_fuel) "" else "nofuel-", Sys.Date()))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

load_latest("^list_imputed_dataframes_\\d{4}-\\d{2}-\\d{2}\\.Rdata$")
imps <- lapply(allImputations$imputations, collapse_covariates)
m <- length(imps)

### Populations (copied from propensity_scores.R) ##############################
cov_overall <- c("Income_collapsed", "Education_collapsed", "Employment_collapsed",
                 "Living.Area", "B_Age", "B_Gender", "Site", "season",
                 "alcohol_disorder", "BMI", "B_Height",
                 "B_HbA1c", "B_BCG", "B_Sleepindex", "Type.of.House",
                 "smoking_status", "Household.size")
cov_children <- c("Income_collapsed",
                  "Living.Area_collapsed", "B_Age", "B_Gender", "Site", "season",
                  "BMI_for_age", "Household.size")

prepare_children <- function(df) {
  child_df <- df[df$B_Age < 15, ]
  zscorer::addWGSR(data = child_df, sex = "B_Gender_numeric",
                   firstPart = "B_Weight", secondPart = "B_Height",
                   thirdPart = "B_Age_months", index = "bfa",
                   output = "BMI_for_age")
}
pops <- list(
  overall  = list(prep = function(d) d,                 covs = c(cov_overall, pm_fuel_adjustment)),
  children = list(prep = prepare_children,              covs = c(cov_children, pm_fuel_adjustment)),
  adults   = list(prep = function(d) d[d$B_Age >= 15, ], covs = c(cov_overall, pm_fuel_adjustment))
)
outcomes <- c(IGRA = "qgit_positive_entry", TST = "tst_positive_entry")
ref_val <- 40
contrasts <- c(60, 80)

### Helpers ####################################################################
gps_weights <- function(d, covs) {
  f <- as.formula(paste("mean_pm_recentavg ~", paste(covs, collapse = " + ")))
  weightit(f, data = d, method = "glm", stabilize = TRUE)
}

weight_stats <- function(w) {
  n <- length(w); s <- sort(w, decreasing = TRUE)
  tibble(n = n, ess = sum(w)^2 / sum(w^2), ess_pct = 100 * sum(w)^2 / sum(w^2) / n,
         min = min(w), p01 = quantile(w, .01), p05 = quantile(w, .05), median = median(w),
         p95 = quantile(w, .95), p99 = quantile(w, .99), max = max(w),
         top5pct_share_of_weight = sum(s[seq_len(ceiling(.05 * n))]) / sum(w),
         n_weights_over_10 = sum(w > 10))
}

truncate_w <- function(w, lo = .01, hi = .99) {
  q <- quantile(w, c(lo, hi)); pmin(pmax(w, q[1]), q[2])
}

# RD and RR at each contrast vs ref_val. Spline basis is built on the whole
# population (as in propensity_scores.R); outcome rows with NA are dropped.
contrast_fit <- function(d, y, w = NULL, covs = NULL) {
  basis <- ns(d$mean_pm_recentavg, df = 3)
  B <- as.data.frame(predict(basis, d$mean_pm_recentavg)); names(B) <- c("ns1", "ns2", "ns3")
  dd <- bind_cols(d, B)
  keep <- !is.na(dd[[y]])
  dd <- dd[keep, ]; ww <- if (is.null(w)) NULL else w[keep]
  rhs <- paste(c("ns1 + ns2 + ns3", covs), collapse = " + ")
  fit <- suppressWarnings(glm(as.formula(paste(y, "~", rhs)), data = dd,
                              family = quasibinomial(), weights = ww))
  risk_at <- function(x) {
    b <- predict(basis, rep(x, nrow(dd)))
    nd <- dd; nd$ns1 <- b[, 1]; nd$ns2 <- b[, 2]; nd$ns3 <- b[, 3]
    # With covariates, standardize predictions over the sample (g-computation);
    # without, every row has the same prediction.
    mean(predict(fit, newdata = nd, type = "response"))
  }
  r0 <- risk_at(ref_val)
  map_dfr(contrasts, function(cv) { r1 <- risk_at(cv); tibble(contrast = cv, risk_ref = r0, risk_contrast = r1, rd = r1 - r0, rr = r1 / r0) })
}

pm_band <- function(x) cut(x, c(-Inf, 40, 60, 80, Inf), right = FALSE,
                           labels = c("<40", "40–59", "60–79", "≥80"))

### Run ########################################################################
wsum <- list(); rds <- list(); w_first <- list()
for (pn in names(pops)) {
  message("Population: ", pn)
  for (i in seq_len(m)) {
    d <- pops[[pn]]$prep(imps[[i]])
    wfit <- gps_weights(d, pops[[pn]]$covs)
    w <- wfit$weights
    if (i == 1) {
      w_first[[pn]] <- list(d = d, w = w, fit = wfit)
    }
    wsum[[length(wsum) + 1]] <- weight_stats(w) %>% mutate(population = pn, imputation = i)
    for (on in names(outcomes)) {
      y <- outcomes[[on]]
      methods <- list(
        `Unweighted` = function() contrast_fit(d, y),
        `GPS-weighted (as in paper)` = function() contrast_fit(d, y, w = w),
        `GPS-weighted, truncated 1st/99th pct` = function() contrast_fit(d, y, w = truncate_w(w)),
        `GPS-weighted, truncated 5th/95th pct` = function() contrast_fit(d, y, w = truncate_w(w, .05, .95)),
        `Outcome regression (covariate-adjusted)` = function() contrast_fit(d, y, covs = pops[[pn]]$covs),
        `GPS-weighted, Site-stratified average` = function() {
          bind_rows(lapply(split(seq_len(nrow(d)), d$Site), function(ix) {
            contrast_fit(d[ix, ], y, w = w[ix]) %>% mutate(n_site = length(ix))
          })) %>% group_by(contrast) %>%
            summarise(across(c(risk_ref, risk_contrast, rd, rr), ~ weighted.mean(.x, n_site)), .groups = "drop")
        }
      )
      for (mn in names(methods)) {
        res <- tryCatch(methods[[mn]](), error = function(e) NULL)
        if (!is.null(res)) rds[[length(rds) + 1]] <- res %>% mutate(population = pn, outcome = on, method = mn, imputation = i)
      }
    }
  }
}

weight_summary <- bind_rows(wsum) %>% group_by(population) %>%
  summarise(across(n:n_weights_over_10, median), .groups = "drop") %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))
write_csv(weight_summary, file.path(out_dir, "weight_summary.csv"))
print(weight_summary, width = Inf)

rd_by_method <- bind_rows(rds) %>%
  group_by(population, outcome, method, contrast) %>%
  summarise(n_imputations = n(), risk_ref = mean(risk_ref), risk_contrast = mean(risk_contrast),
            rd = mean(rd), rr = exp(mean(log(rr))), .groups = "drop") %>%
  mutate(across(c(risk_ref, risk_contrast, rd), ~ round(.x, 3)), rr = round(rr, 2)) %>%
  arrange(population, outcome, contrast, method)
write_csv(rd_by_method, file.path(out_dir, "rd_by_method.csv"))
print(rd_by_method, n = Inf, width = Inf)

### Balance (imputation 1) #####################################################
for (pn in c("children", "overall")) {
  bt <- bal.tab(w_first[[pn]]$fit, un = TRUE, stats = "correlations")
  b <- bt$Balance %>% rownames_to_column("covariate") %>%
    select(covariate, any_of(c("Corr.Un", "Corr.Adj"))) %>%
    mutate(across(where(is.numeric), ~ round(.x, 3)))
  write_csv(b, file.path(out_dir, paste0("balance_", pn, ".csv")))
  message("Balance, ", pn, ": max |corr| before = ", max(abs(b$Corr.Un), na.rm = TRUE),
          ", after = ", max(abs(b$Corr.Adj), na.rm = TRUE))
}

### Overlap and crude patterns (imputation 1) ##################################
ch <- w_first$children$d
pm_site <- ch %>% group_by(Site) %>%
  summarise(n = n(), min = min(mean_pm_recentavg), p25 = quantile(mean_pm_recentavg, .25),
            median = median(mean_pm_recentavg), p75 = quantile(mean_pm_recentavg, .75),
            max = max(mean_pm_recentavg), n_at_or_above_80 = sum(mean_pm_recentavg >= 80),
            n_below_40 = sum(mean_pm_recentavg < 40), .groups = "drop") %>%
  mutate(across(where(is.double), ~ round(.x, 1)))
write_csv(pm_site, file.path(out_dir, "pm_by_site_children.csv"))
print(pm_site)

prev <- bind_rows(
  w_first$children$d %>% mutate(group = "Children"),
  w_first$adults$d %>% mutate(group = "Adults")
) %>%
  mutate(pm_band = pm_band(mean_pm_recentavg)) %>%
  pivot_longer(all_of(unname(outcomes)), names_to = "outcome", values_to = "positive") %>%
  filter(!is.na(positive)) %>%
  mutate(outcome = names(outcomes)[match(outcome, outcomes)]) %>%
  group_by(group, outcome, Site, pm_band) %>%
  summarise(n = n(), n_positive = sum(positive), percent_positive = round(100 * mean(positive), 1), .groups = "drop")
write_csv(prev, file.path(out_dir, "prevalence_by_pm_band.csv"))
print(prev, n = Inf)

wc <- w_first$children$w
top <- ch %>% mutate(w = wc, top5 = w >= quantile(w, .95), pm_band = pm_band(mean_pm_recentavg)) %>%
  filter(top5) %>% count(Site, pm_band, igra_positive = qgit_positive_entry, tst_positive = tst_positive_entry)
write_csv(top, file.path(out_dir, "top_weights_children.csv"))
print(top)

hist_df <- bind_rows(lapply(names(w_first), function(pn) tibble(population = pn, w = w_first[[pn]]$w)))
jpeg(file.path(out_dir, "weights_histogram.jpg"), width = 9, height = 5, units = "in", res = 200)
print(ggplot(hist_df, aes(w)) + geom_histogram(bins = 50) + scale_x_log10() +
        facet_wrap(~ population, scales = "free_y") +
        labs(x = "Stabilized GPS weight (log scale), imputation 1", y = "Count") + theme_minimal())
dev.off()
message("Saved to ", out_dir)
