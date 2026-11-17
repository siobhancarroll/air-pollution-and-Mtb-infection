### Alternative PM2.5 exposure windows for the conversion models ##############
# pm_window_lookup(window) returns one row per PID x follow-up interval with
# mean_pm_intervalavg / SE_pm_intervalavg recomputed for an alternative window,
# from the same site x calendar-month PM2.5 series and visit dates used in
# data_processing.R (actual visit dates, otherwise Entry.Date + 4/12/24
# months). Calendar months are whole months; a window average is NA if any of
# its months has no PM2.5 value. The SE of a K-month mean is
# sqrt(mean(SE_i^2) / K), as for the primary window.
#
# Windows (interval k runs from the previous scheduled visit to visit k):
#   "pre_entry3"     the 3 calendar months before the Entry month, held fixed
#                    for every interval (exposure before enrolment, when much
#                    household transmission is likely to have occurred)
#   "early2"         the first 2 calendar months of each interval (exposure
#                    early enough for an infection acquired then to be
#                    detectable by the end-of-interval test)
#   "first_interval" the Entry-to-Month-4 average, held fixed for every
#                    interval (early follow-up exposure predicting later
#                    conversion)
library(dplyr); library(lubridate); library(here)

pm_windows <- c("pre_entry3", "early2", "first_interval")
# "interval" (the primary window) is also available here, so that the
# window-season covariate below can be computed for it as well.
# pm_window_<season>: share of the window's calendar months in each season
# (Winter Dec-Feb, Pre-Monsoon Mar-May, Monsoon Jun-Sep, Post-Monsoon Oct-Nov,
# as in data_processing.R), used when adjust_window_season is TRUE in
# discrete_time_survival.R to adjust for the season of the exposure months
# rather than only the season at the end of the interval.

pm_monthly_series <- function() {
  ch <- readxl::read_xlsx(here("Data", "Statistics_Chennai.xlsx"))
  pu <- readxl::read_xlsx(here("Data", "Statistics_pune.xlsx"))
  bind_rows(
    transmute(ch, Site = "NIRT",  month = floor_date(as.Date(Months), "month"), mean_pm = Chennai, SE_pm = `Standard error`),
    transmute(pu, Site = "BJGMC", month = floor_date(as.Date(Months), "month"), mean_pm = Pune,    SE_pm = SE)
  )
}

pm_visit_dates <- function() {
  env <- new.env(); load(here("Data", "data.Rdata"), envir = env)
  cmp <- env$complete %>% mutate(PID = as.character(PID), VISIT = gsub(" ", "", as.character(VISIT)))
  entry <- cmp %>% filter(VISIT == "Entry", !is.na(Entry.Date)) %>% distinct(PID, .keep_all = TRUE) %>%
    transmute(PID, Site = as.character(Site), entry = as.Date(Entry.Date))
  vis <- cmp %>% filter(VISIT %in% c("Month4", "Month12", "Month24")) %>% distinct(PID, VISIT, .keep_all = TRUE) %>%
    transmute(PID, VISIT, d = as.Date(visit_date)) %>%
    tidyr::pivot_wider(names_from = VISIT, values_from = d)
  for (v in c("Month4", "Month12", "Month24")) if (!v %in% names(vis)) vis[[v]] <- as.Date(NA)
  entry %>% left_join(vis, by = "PID") %>%
    mutate(m4 = coalesce(Month4, entry %m+% months(4)), m12 = coalesce(Month12, entry %m+% months(12)),
           m24 = coalesce(Month24, entry %m+% months(24))) %>%
    select(PID, Site, entry, m4, m12, m24)
}

# months: list column of calendar months (Dates) per PID x interval
pm_window_months <- function(dates, window) {
  em <- floor_date(dates$entry, "month")
  starts <- list(dates$entry, dates$m4, dates$m12)
  ends <- list(dates$m4, dates$m12, dates$m24)
  bind_rows(lapply(1:3, function(k) {
    months <- switch(window,
      interval = Map(function(a, b) if (is.na(a) || is.na(b) || b < a) as.Date(NA) else seq(a, b, by = "month"),
                     floor_date(starts[[k]], "month"), floor_date(ends[[k]], "month")),
      pre_entry3 = lapply(em, function(m) m %m-% months(1:3)),
      early2 = lapply(floor_date(starts[[k]], "month"), function(m) m %m+% months(0:1)),
      # Month 4 dated before Entry (recording irregularity, as excluded in
      # data_processing.R): no window, so the value is NA
      first_interval = Map(function(a, b) if (is.na(a) || is.na(b) || b < a) as.Date(NA) else seq(a, b, by = "month"),
                           em, floor_date(dates$m4, "month")),
      stop("Unknown pm_window: ", window))
    tibble(PID = dates$PID, Site = dates$Site, interval = k, month = months)
  }))
}

pm_window_lookup <- function(window, dates = pm_visit_dates(), series = pm_monthly_series()) {
  stopifnot(window %in% c("interval", pm_windows))
  long <- pm_window_months(dates, window) %>% tidyr::unnest(month) %>%
    left_join(series, by = c("Site", "month"))
  out <- long %>% group_by(PID, interval) %>%
    summarise(mean_pm_intervalavg = mean(mean_pm), SE_pm_intervalavg = sqrt(mean(SE_pm^2) / n()),
              n_months = n(), pm_window_monsoon = mean(month(month) %in% 6:9),
              pm_window_winter = mean(month(month) %in% c(12, 1, 2)), pm_window_premonsoon = mean(month(month) %in% 3:5),
              pm_window_postmonsoon = mean(month(month) %in% 10:11), .groups = "drop") %>%
    mutate(across(starts_with("pm_window_"), ~ if_else(is.na(mean_pm_intervalavg), NA_real_, .x)))
  message("pm_window ", window, ": ", sum(!is.na(out$mean_pm_intervalavg)), " of ", nrow(out),
          " PID x interval windows have complete PM2.5 coverage")
  out
}

# Replace the primary interval-average PM2.5 in a person-period data frame
# replace_pm = FALSE keeps the primary interval-average PM2.5 and only adds
# pm_window_monsoon. Where the window is undefined but PM2.5 is present
# (rare date irregularities), the monsoon share falls back to the
# interval-end season.
apply_pm_window <- function(pp, lookup, replace_pm = TRUE) {
  key <- as.character(pp$PID)
  m <- match(paste(key, pp$interval), paste(lookup$PID, lookup$interval))
  if (replace_pm) {
    pp$mean_pm_intervalavg <- lookup$mean_pm_intervalavg[m]
    pp$SE_pm_intervalavg <- lookup$SE_pm_intervalavg[m]
  }
  end_season <- as.character(pp$season)
  for (v in names(SEASON_MONTHS)) {
    col <- paste0("pm_window_", v)
    pp[[col]] <- lookup[[col]][m]
    fill <- is.na(pp[[col]]) & !is.na(pp$mean_pm_intervalavg)
    if (any(fill)) pp[[col]][fill] <- as.numeric(end_season[fill] == SEASON_LABELS[[v]])
  }
  pp
}

SEASON_MONTHS <- list(monsoon = 6:9, winter = c(12, 1, 2), premonsoon = 3:5, postmonsoon = 10:11)
SEASON_LABELS <- c(monsoon = "Monsoon", winter = "Winter", premonsoon = "Pre-Monsoon", postmonsoon = "Post-Monsoon")
# Window-season covariates: shares of the exposure window's months in each
# season. With four seasons, Monsoon is the omitted reference share.
# The four-season shares (3 parameters on top of the 4-level interval-end
# season) made the TST pre-enrolment model non-estimable (non-finite
# cumulative hazards in most imputations), so a single parameter is used:
# the share of the window's months in the monsoon, the low-PM2.5 season.
window_season_covariates <- function(season_detail) "pm_window_monsoon"
