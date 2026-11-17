# Questionnaire-aware fuel processing. Pure data transformations: no model fitting.
# Raw fields are preserved so cached Data/data.Rdata can be reprocessed safely.
fuel_secondary_columns <- paste0("secondary_fuel_", c(
  "kerosene", "gas", "coal", "electricity", "wood", "dung", "crop_waste", "other"))
fuel_polluting_secondary <- paste0("secondary_fuel_", c(
  "kerosene", "coal", "wood", "dung", "crop_waste", "other"))
fuel_source_columns <- c("primary_fuel", fuel_secondary_columns)
fuel_clean_primary <- c("Gas (cylinder/piped)", "Electricity")
fuel_polluting_primary <- c("Wood", "Paraffin/Kerosene", "Coal", "Crop Waste")

normalize_fuel_flag <- function(x) {
  value <- toupper(trimws(as.character(x)))
  value[is.na(x) | value == ""] <- NA_character_
  bad <- !is.na(value) & !value %in% c("YES", "NO")
  if (any(bad)) stop("Unexpected fuel response(s): ", paste(unique(value[bad]), collapse = ", "))
  ifelse(is.na(value), NA_character_, ifelse(value == "YES", "YES", "No"))
}

# Three-valued classification: confirmed positive overrides uncertainty;
# a negative classification requires the relevant components to be known.
derive_fuel_indicators <- function(df, suffix = "") {
  col <- function(name) as.character(df[[paste0(name, suffix)]])
  required <- paste0(fuel_source_columns, suffix)
  if (!all(required %in% names(df))) stop("Underlying fuel fields are missing: ",
    paste(setdiff(required, names(df)), collapse = ", "))
  primary <- col("primary_fuel")
  known_primary <- !is.na(primary)
  primary_clean <- primary %in% fuel_clean_primary
  primary_polluting <- primary %in% fuel_polluting_primary
  unsupported <- known_primary & !(primary_clean | primary_polluting)
  if (any(unsupported)) stop("Unrecognized primary-fuel category: ", paste(unique(primary[unsupported]), collapse = ", "))
  secondary <- as.data.frame(lapply(fuel_polluting_secondary, col), stringsAsFactors = FALSE)
  any_polluting_secondary <- rowSums(secondary == "YES", na.rm = TRUE) > 0
  all_polluting_secondary_no <- rowSums(!is.na(secondary) & secondary == "No") == ncol(secondary)
  exclusive <- rep(NA_character_, nrow(df))
  exclusive[primary_polluting | any_polluting_secondary] <- "No"
  exclusive[primary_clean & all_polluting_secondary_no] <- "YES"
  fuel_any <- function(fuel, category) {
    secondary_fuel <- col(paste0("secondary_fuel_", fuel))
    out <- rep(NA_character_, nrow(df))
    out[known_primary & secondary_fuel %in% "No"] <- "No"
    out[primary %in% category | secondary_fuel %in% "YES"] <- "YES"
    out
  }
  derived <- list(exclusive_clean = exclusive,
                  wood_use = fuel_any("wood", "Wood"),
                  kerosene_use = fuel_any("kerosene", "Paraffin/Kerosene"),
                  solid_fuel_primary = ifelse(!known_primary, NA_character_,
                                              ifelse(primary_polluting, "YES", "No")))
  for (name in names(derived)) df[[paste0(name, suffix)]] <- factor(derived[[name]], levels = c("No", "YES"))
  df
}

process_fuel_questionnaire <- function(df) {
  required <- c("PID", "RELATED.INDEX.PID", "VISIT", fuel_source_columns, "secondary_fuel_none")
  if (!all(required %in% names(df))) stop("Fuel audit requires: ", paste(setdiff(required, names(df)), collapse = ", "))
  for (name in c(fuel_source_columns, "secondary_fuel_none")) {
    raw <- paste0(name, "_raw")
    if (!raw %in% names(df)) df[[raw]] <- df[[name]]
    df[[name]] <- as.character(df[[raw]])
  }
  primary <- trimws(df$primary_fuel)
  primary[is.na(primary) | primary == ""] <- NA_character_
  df$primary_fuel <- primary
  for (name in c(fuel_secondary_columns, "secondary_fuel_none")) df[[name]] <- normalize_fuel_flag(df[[name]])
  secondary <- df[fuel_secondary_columns]
  all_no <- rowSums(!is.na(secondary) & secondary == "No") == length(fuel_secondary_columns)
  any_yes <- rowSums(secondary == "YES", na.rm = TRUE) > 0
  none_yes <- df$secondary_fuel_none %in% "YES"
  df$fuel_default_no_pattern <- all_no & df$secondary_fuel_none %in% "No"
  df$fuel_likely_questionnaire_nonresponse <- df$fuel_default_no_pattern & is.na(primary)
  df$fuel_none_conflict <- none_yes & any_yes
  df$fuel_primary_observed <- !is.na(primary)
  # An all-No secondary section is ambiguous even when primary fuel is known.
  for (name in fuel_secondary_columns) df[[name]][df$fuel_default_no_pattern] <- NA_character_
  # Explicit none supplies negative secondary responses only if not contradicted.
  for (name in fuel_secondary_columns) {
    fill <- none_yes & !any_yes & is.na(df[[name]])
    df[[name]][fill] <- "No"
  }
  df$fuel_secondary_missing_before_household <- rowSums(is.na(df[fuel_secondary_columns]))
  # Never combine different visits or unrelated participants with missing HH IDs.
  hh <- as.character(df$RELATED.INDEX.PID)
  hh[is.na(hh) | trimws(hh) == ""] <- paste0("PERSON:", df$PID[is.na(hh) | trimws(hh) == ""])
  groups <- split(seq_len(nrow(df)), interaction(hh, as.character(df$VISIT), drop = TRUE))
  df$fuel_household_conflict <- FALSE
  df$fuel_household_fields_filled <- 0L
  df$fuel_pre_imputation_conflict <- df$fuel_none_conflict
  for (idx in groups) {
    conflict <- any(vapply(fuel_source_columns, function(name)
      length(unique(stats::na.omit(df[[name]][idx]))) > 1L, logical(1)))
    df$fuel_household_conflict[idx] <- conflict
    # Preserve conflicting observed reports for audit; do not pick a respondent.
    if (conflict || any(df$fuel_none_conflict[idx])) next
    for (name in fuel_source_columns) {
      values <- unique(stats::na.omit(df[[name]][idx]))
      if (length(values) == 1L) {
        missing <- idx[is.na(df[[name]][idx])]
        df[[name]][missing] <- values
        df$fuel_household_fields_filled[missing] <- df$fuel_household_fields_filled[missing] + 1L
      }
    }
  }
  # Differing reports from household members are retained as individual
  # observations, including reports made on the same date. These differences
  # remain descriptive audit flags and prevent household borrowing, but do not
  # block imputation. A within-person none/positive contradiction still does.
  df$fuel_pre_imputation_conflict <- df$fuel_none_conflict
  for (name in fuel_source_columns) {
    lev <- if (name == "primary_fuel") c(fuel_clean_primary, fuel_polluting_primary) else c("No", "YES")
    if (name == "primary_fuel" && any(!is.na(df[[name]]) & !df[[name]] %in% lev))
      stop("Unrecognized primary fuel; review raw categories before imputation")
    df[[name]] <- factor(df[[name]], levels = lev)
  }
  df <- derive_fuel_indicators(df)
  attr(df, "fuel_processing_version") <- "questionnaire-aware-v1"
  df
}

write_fuel_questionnaire_audit <- function(df, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  flags <- c("fuel_primary_observed", "fuel_default_no_pattern",
             "fuel_likely_questionnaire_nonresponse", "fuel_none_conflict", "fuel_household_conflict")
  audit <- df %>% dplyr::group_by(Site, VISIT) %>%
    dplyr::summarise(n_records = dplyr::n(),
      # Count households before summarise() replaces record-level flags with
      # their sums; otherwise the sum becomes a numeric row subscript.
      households_default_no = dplyr::n_distinct(RELATED.INDEX.PID[fuel_default_no_pattern], na.rm = TRUE),
      households_with_conflicts = dplyr::n_distinct(RELATED.INDEX.PID[fuel_none_conflict | fuel_household_conflict], na.rm = TRUE),
      dplyr::across(dplyr::all_of(flags), ~sum(.x)),
      primary_missing_after_household = sum(is.na(primary_fuel)),
      records_with_secondary_missing = sum(rowSums(is.na(dplyr::pick(dplyr::all_of(fuel_secondary_columns)))) > 0),
      records_receiving_household_fields = sum(fuel_household_fields_filled > 0),
      exclusive_clean_observed = sum(exclusive_clean == "YES", na.rm = TRUE),
      any_polluting_observed = sum(exclusive_clean == "No", na.rm = TRUE),
      fuel_classification_unknown = sum(is.na(exclusive_clean)), .groups = "drop")
  readr::write_csv(audit, file.path(output_dir, "questionnaire_response_audit.csv"))
  patterns <- df %>% dplyr::count(Site, VISIT, fuel_primary_observed,
    fuel_default_no_pattern, fuel_none_conflict, fuel_household_conflict, name = "n_records")
  readr::write_csv(patterns, file.path(output_dir, "response_patterns.csv"))
  availability <- dplyr::bind_rows(lapply(fuel_source_columns, function(name) {
    df %>% dplyr::group_by(Site, VISIT) %>%
      dplyr::summarise(variable = name,
        n_observed = sum(!is.na(.data[[name]])), n_missing = sum(is.na(.data[[name]])),
        n_observed_categories = dplyr::n_distinct(.data[[name]], na.rm = TRUE),
        .groups = "drop")
  }))
  readr::write_csv(availability, file.path(output_dir, "source_response_availability.csv"))
  categories <- dplyr::bind_rows(lapply(fuel_source_columns, function(name) {
    df %>% dplyr::transmute(Site, VISIT, variable = name, response = as.character(.data[[name]])) %>%
      dplyr::count(Site, VISIT, variable, response, name = "n_records")
  }))
  readr::write_csv(categories, file.path(output_dir, "source_response_categories.csv"))
  # Identifiable tracing records stay in Data/, which is excluded from git.
  detail_dir <- here::here("Data", "fuel_questionnaire_audit")
  dir.create(detail_dir, recursive = TRUE, showWarnings = FALSE)
  flagged <- df[df$fuel_default_no_pattern | df$fuel_none_conflict | df$fuel_household_conflict, ]
  readr::write_csv(flagged %>% dplyr::select(PID, RELATED.INDEX.PID, VISIT, Site,
    dplyr::all_of(flags), fuel_household_fields_filled,
    dplyr::all_of(fuel_source_columns), secondary_fuel_none,
    dplyr::all_of(paste0(c(fuel_source_columns, "secondary_fuel_none"), "_raw"))),
    file.path(detail_dir, paste0("flagged_records_", Sys.Date(), ".csv")))
  baseline_conflicts <- df[gsub(" ", "", as.character(df$VISIT)) == "Entry" &
                            df$fuel_household_conflict, , drop = FALSE]
  baseline_conflicts <- baseline_conflicts[!is.na(baseline_conflicts$PID), , drop = FALSE]
  groups <- split(seq_len(nrow(baseline_conflicts)), interaction(
    baseline_conflicts$Site, baseline_conflicts$RELATED.INDEX.PID, drop = TRUE))
  conflict_review <- dplyr::bind_rows(lapply(groups, function(idx) {
    d <- baseline_conflicts[idx, , drop = FALSE]
    differing <- fuel_source_columns[vapply(fuel_source_columns,
      function(v) length(unique(stats::na.omit(d[[v]]))) > 1L, logical(1))]
    tibble::tibble(Site = as.character(d$Site[1]),
      household_id = as.character(d$RELATED.INDEX.PID[1]), n_contacts = nrow(d),
      conflict_type = if ("primary_fuel" %in% differing) "Primary fuel" else
        if (any(differing %in% fuel_polluting_secondary)) "Polluting secondary fuel" else "Clean secondary fuel only",
      conflicting_fields = paste(differing, collapse = "; "),
      reported_values = paste(vapply(differing, function(v)
        paste0(v, ": ", paste(sort(unique(stats::na.omit(as.character(d[[v]])))), collapse = " / ")),
        character(1)), collapse = "; "))
  }))
  readr::write_csv(conflict_review,
    file.path(detail_dir, paste0("baseline_household_conflicts_", Sys.Date(), ".csv")))
  print(audit, n = Inf, width = Inf)
  invisible(audit)
}

# Derived indicators are not separately imputed. Reconstruct them from sources
# whenever a completed dataset is extracted, including the .imp = 0 records.
complete_fuel_imputation <- function(data, action = 1L, include = FALSE) {
  result <- mice::complete(data, action = action, include = include)
  result <- derive_fuel_indicators(result, "_Entry")
  if (!".imp" %in% names(result) || any(result$.imp %in% seq_len(data$m))) {
    check <- if (".imp" %in% names(result)) result[result$.imp > 0, , drop = FALSE] else result
    if (anyNA(check[paste0(c("exclusive_clean", "wood_use", "kerosene_use", "solid_fuel_primary"), "_Entry")]))
      stop("Completed source fuel variables did not yield complete classifications")
  }
  result
}
