# Synthetic checks only: no participant data, imputation, or analysis fits.
source(here::here("fuel_processing_helpers.R"))
example <- function(n = 1L) {
  d <- data.frame(PID = paste0("P", seq_len(n)), RELATED.INDEX.PID = paste0("H", seq_len(n)),
                  VISIT = "Entry", Site = "TEST", primary_fuel = "Gas (cylinder/piped)",
                  secondary_fuel_none = "YES", stringsAsFactors = FALSE)
  for (name in fuel_secondary_columns) d[[name]] <- rep("No", n)
  d
}
x <- example(6)
x$secondary_fuel_none[c(2, 3, 4, 5)] <- "No"
x$primary_fuel[c(3, 4)] <- NA_character_
x$secondary_fuel_wood[4] <- "YES"
x$secondary_fuel_wood[5] <- "YES"
x$secondary_fuel_kerosene[5] <- "YES"
x$primary_fuel[6] <- "Wood"
a <- process_fuel_questionnaire(x)
stopifnot(as.character(a$exclusive_clean[1]) == "YES",
          all(a$fuel_default_no_pattern == c(FALSE, TRUE, TRUE, FALSE, FALSE, FALSE)),
          is.na(a$exclusive_clean[2]), !is.na(a$primary_fuel[2]),
          all(is.na(a[2, fuel_secondary_columns])),
          a$fuel_likely_questionnaire_nonresponse[3],
          as.character(a$wood_use[4]) == "YES", as.character(a$exclusive_clean[4]) == "No",
          is.na(a$kerosene_use[4]),
          as.character(a$wood_use[5]) == "YES", as.character(a$kerosene_use[5]) == "YES",
          as.character(a$exclusive_clean[5]) == "No",
          as.character(a$solid_fuel_primary[6]) == "YES")
# An unanswered household member receives compatible source responses.
x <- example(2);x$RELATED.INDEX.PID[] <- "SHARED"
x$primary_fuel[2] <- NA_character_;x$secondary_fuel_none[2] <- "No"
a <- process_fuel_questionnaire(x)
stopifnot(as.character(a$exclusive_clean[2]) == "YES", a$fuel_household_fields_filled[2] == 9L)
# Reprocessing uses raw responses, preserving provenance and the same results.
b <- process_fuel_questionnaire(a)
stopifnot(identical(a$fuel_default_no_pattern, b$fuel_default_no_pattern),
          identical(a$fuel_household_fields_filled, b$fuel_household_fields_filled))
# No cross-visit sharing, even in the same household.
x$VISIT[2] <- "Month4";a <- process_fuel_questionnaire(x)
stopifnot(is.na(a$primary_fuel[2]), is.na(a$exclusive_clean[2]))
# Different reports, including on the same date, retain each person's values
# and do not block imputation or copy a respondent's report to the other.
x <- example(2);x$RELATED.INDEX.PID[] <- "SHARED";x$primary_fuel[2] <- "Wood"
x$visit_date <- as.Date("2015-01-01")
a <- process_fuel_questionnaire(x)
stopifnot(all(a$fuel_household_conflict), all(a$fuel_household_fields_filled == 0L),
          !any(a$fuel_pre_imputation_conflict),
          as.character(a$primary_fuel[1]) == "Gas (cylinder/piped)",
          as.character(a$primary_fuel[2]) == "Wood",
          as.character(a$exclusive_clean[1]) == "YES",
          as.character(a$exclusive_clean[2]) == "No")
# Missing household IDs must not cause unrelated people to share responses.
x <- example(2);x$RELATED.INDEX.PID[] <- NA_character_
x$primary_fuel[2] <- NA_character_;x$secondary_fuel_none[2] <- "No"
a <- process_fuel_questionnaire(x);stopifnot(is.na(a$primary_fuel[2]))
# Explicit none plus a positive response is audited as a source conflict.
x <- example();x$secondary_fuel_wood <- "YES"
a <- process_fuel_questionnaire(x)
stopifnot(a$fuel_none_conflict, a$fuel_pre_imputation_conflict,
          as.character(a$exclusive_clean) == "No")
message("Fuel-processing synthetic checks passed")
