# Air pollution and *Mycobacterium tuberculosis* infection

Analysis code for the second manuscript of my PhD thesis. The manuscript estimates the effects of outdoor PM2.5 and household polluting fuel use on IGRA and TST conversion among household contacts of people with TB in Pune and Chennai, India (CTRIUMPH cohort).

The study data are not included in this repository.

## Requirements

R with the following packages: `here`, `tidyverse`, `mice`, `WeightIt`, `cobalt`, `survey`, `marginaleffects`, `mitools`, `splines`, `future.apply` and `zscorer`.

## Running the analyses

1. Prepare the analysis dataset with `data_processing.R`. This saves `Data/data.Rdata`.
2. From the project folder, run the full pipeline:

   ```r
   source("run_all.R")
   ```

   This runs multiple imputation, the conversion (survival) analyses, the baseline analyses and the sensitivity analyses, in that order. To reuse existing imputations, set `run_imputation <- FALSE` before sourcing.

Results are saved in date-stamped folders under `tables_and_figures/`, and a log for each step is saved under `logs/`.

## Main scripts

| Script | Purpose |
|---|---|
| `data_processing.R` | Builds the analysis dataset |
| `multiple_imputation.R` | Multiple imputation of missing covariates and fuel use |
| `discrete_time_survival.R` | IGRA and TST conversion analyses |
| `propensity_scores.R` | Baseline (prevalent positivity) analyses |
| `run_all.R` | Runs all of the above in order |
