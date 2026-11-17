# Descriptive audit only: no imputation, weights or outcome models are fitted.
# Run from the project root: source("scripts/audit_shared_exposure_structure.R")
# Only aggregate summaries are exported; participant/household IDs stay in memory.
local({
  stopifnot(dir.exists("Data"))
  if (!requireNamespace("lubridate", quietly = TRUE)) stop("lubridate is required for visit-date fallback")
  input_dir <- normalizePath("Data")
  out <- file.path("tables_and_figures", paste0("shared-exposure-audit-", Sys.Date()))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  manifest <- list()
  load_private <- function(path, object) {
    env <- new.env(parent = baseenv())
    load(path, envir = env)
    if (!exists(object, envir = env, inherits = FALSE)) stop("Missing object ", object, " in ", path)
    manifest[[length(manifest) + 1L]] <<- data.frame(file = normalizePath(path), md5 = unname(tools::md5sum(path)))
    get(object, envir = env, inherits = FALSE)
  }
  latest <- function(pattern) {
    files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
    if (!length(files)) stop("No input matching ", pattern)
    sort(files)[length(files)]
  }
  baseline_file <- latest("^list_imputed_dataframes_[0-9]{4}-[0-9]{2}-[0-9]{2}\\.Rdata$")
  input_date <- sub(".*_([0-9]{4}-[0-9]{2}-[0-9]{2})\\.Rdata$", "\\1", baseline_file)
  baseline <- load_private(baseline_file, "allImputations")$imputations[[1]]
  stopifnot(all(c("PID", "RELATED.INDEX.PID", "Site") %in% names(baseline)), !anyDuplicated(baseline$PID))
  raw <- load_private(file.path(input_dir, "data.Rdata"), "complete")
  stopifnot(all(c("PID", "VISIT", "Entry.Date", "visit_date") %in% names(raw)))
  raw$PID <- as.character(raw$PID)
  raw$visit <- gsub(" ", "", as.character(raw$VISIT))
  hh <- baseline[, c("PID", "RELATED.INDEX.PID", "Site")]
  hh$PID <- as.character(hh$PID)
  hh$household <- ifelse(is.na(hh$RELATED.INDEX.PID) | !nzchar(trimws(as.character(hh$RELATED.INDEX.PID))),
    NA_character_, paste(hh$Site, hh$RELATED.INDEX.PID, sep = ":"))
  entry <- raw[raw$visit == "Entry", c("PID", "Entry.Date", "visit_date")]
  if (anyDuplicated(entry$PID)) stop("Multiple Entry records per PID; reconcile before auditing")
  entry$entry_date <- as.Date(entry$Entry.Date)
  entry$entry_date[is.na(entry$entry_date)] <- as.Date(entry$visit_date[is.na(entry$entry_date)])
  hh <- merge(hh, entry[, c("PID", "entry_date")], by = "PID", all.x = TRUE, sort = FALSE)
  outcomes <- list(baseline_all = baseline)
  for (assay in c("qgit", "tst")) {
    file <- latest(paste0("^survival_person_period_", assay, "_[0-9]{4}-[0-9]{2}-[0-9]{2}\\.Rdata$"))
    if (!endsWith(file, paste0(input_date, ".Rdata"))) stop("Baseline and survival input dates differ; do not mix analysis datasets")
    value <- load_private(file, paste0("allImputations_survival_", assay))$imputations[[1]]
    outcomes[[paste0(assay, "_primary")]] <- value[value$interval <= 2, , drop = FALSE]
  }
  summaries <- list(); windows <- list(); months <- list()
  summarize <- function(df, label) {
    df$PID <- as.character(df$PID)
    matched <- match(df$PID, hh$PID)
    if (anyNA(matched)) stop("A cohort PID is absent from the baseline household lookup")
    household <- hh$household[matched]
    people <- !duplicated(df$PID)
    sizes <- table(household[people], useNA = "no")
    summaries[[length(summaries) + 1L]] <<- data.frame(
      analysis = label, n_rows = nrow(df), n_people = sum(people),
      n_people_missing_household = sum(is.na(household[people])), n_households = length(sizes),
      n_households_multiple_contacts = sum(sizes > 1),
      median_contacts_per_household = if (length(sizes)) median(as.numeric(sizes)) else NA_real_,
      max_contacts_per_household = if (length(sizes)) max(sizes) else NA_integer_)
    if (label == "baseline_all") {
      starts <- lubridate::floor_date(hh$entry_date[matched], "month") %m-% lubridate::months(3)
      ends <- lubridate::floor_date(hh$entry_date[matched], "month") %m-% lubridate::months(1)
    } else {
      # Preserve recorded dates; nominal dates are used only where absent.
      visit_date <- function(visit, offset) {
        r <- raw[raw$visit == visit, c("PID", "visit_date")]
        if (anyDuplicated(r$PID)) stop("Duplicate ", visit, " dates; reconcile before auditing")
        dates <- as.Date(r$visit_date[match(df$PID, r$PID)])
        nominal <- lubridate::add_with_rollback(hh$entry_date[matched], lubridate::months(offset))
        dates[is.na(dates)] <- nominal[is.na(dates)]
        dates
      }
      m4 <- visit_date("Month4", 4); m12 <- visit_date("Month12", 12)
      # ifelse() strips Date classes; reconstruct explicitly before floor_date.
      starts <- lubridate::floor_date(as.Date(ifelse(df$interval == 1, as.numeric(hh$entry_date[matched]), as.numeric(m4)), origin = "1970-01-01"), "month")
      ends <- lubridate::floor_date(as.Date(ifelse(df$interval == 1, as.numeric(m4), as.numeric(m12)), origin = "1970-01-01"), "month")
    }
    site <- as.character(hh$Site[matched])
    valid <- !is.na(starts) & !is.na(ends) & starts <= ends
    key <- paste(site, starts, ends, sep = "|")
    counts <- table(key[valid])
    windows[[length(windows) + 1L]] <<- data.frame(analysis = label,
      n_valid_window_rows = sum(valid), n_missing_or_reversed_windows = sum(!valid),
      n_distinct_site_windows = length(counts), n_windows_shared_by_multiple_rows = sum(counts > 1),
      max_rows_per_site_window = if (length(counts)) max(counts) else NA_integer_)
    pieces <- lapply(which(valid), function(i) data.frame(PID = df$PID[i],
      Site = site[i], month = as.character(seq(starts[i], ends[i], by = "month"))))
    if (length(pieces)) {
      expanded <- do.call(rbind, pieces)
      rows <- aggregate(PID ~ Site + month, expanded, length); names(rows)[3] <- "n_interval_rows"
      unique_people <- aggregate(PID ~ Site + month, unique(expanded), length); names(unique_people)[3] <- "n_people"
      value <- merge(rows, unique_people, by = c("Site", "month")); value$analysis <- label
      months[[length(months) + 1L]] <<- value
    }
  }
  # Access the month-subtraction operator without attaching packages globally.
  `%m-%` <- getExportedValue("lubridate", "%m-%")
  for (label in names(outcomes)) {
    df <- outcomes[[label]]; summarize(df, label)
    if (label != "baseline_all") {
      df <- df[!is.na(df$event) & !is.na(df$mean_pm_intervalavg), , drop = FALSE]
      if (nrow(df)) summarize(df, paste0(label, "_observed_pm_rows"))
    }
  }
  write.csv(do.call(rbind, summaries), file.path(out, "household_structure.csv"), row.names = FALSE)
  write.csv(do.call(rbind, windows), file.path(out, "shared_exposure_windows.csv"), row.names = FALSE)
  if (length(months)) write.csv(do.call(rbind, months), file.path(out, "site_month_exposure_use.csv"), row.names = FALSE)
  write.csv(do.call(rbind, manifest), file.path(out, "input_manifest.csv"), row.names = FALSE)
  writeLines(c("Descriptive audit only; no fitted models, weights, imputations or estimates were changed.",
    "The first imputation is used to describe IDs and cohort structure, not imputation variability.",
    "RELATED.INDEX.PID is an index-case household proxy; validate its meaning and completeness before clustering.",
    "Calendar-month overlap is structural, not a test of residual correlation or a reconstruction of daily-weighted PM.",
    "Observed PM rows are before weighting, not a reconstruction of every fitted comparison population.",
    "Temporal dependence between adjacent/overlapping windows is not resolved by counting site-month groups."),
    file.path(out, "interpretation_notes.txt"))
  message("Aggregate audit saved to ", out)
})
