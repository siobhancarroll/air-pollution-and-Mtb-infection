### Figure 2: household fuel use and IGRA/TST conversion (forest plot) ##########
# Risk differences (95% CI), by assay, at Month 4 (panel A) and Month 12
# (panel B). Any polluting fuel, wood and kerosene use are compared with
# exclusively clean fuels. Primary polluting fuel is compared with other
# primary fuels. Reads the
# pooled estimates written by discrete_time_survival.R
# (survival_hap_effects_summary.csv) from the most recent primary survival
# output folder, or from `survival_dir` if set before source()-ing.
#
# Output: tables_and_figures/manuscript-figures-<date>/
#   figure2_fuel_conversion.png (300 dpi) and figure2_fuel_conversion.pdf

library(tidyverse)
library(here)

fig_root <- here("tables_and_figures")
if (!exists("survival_dir")) {
  dirs <- list.files(fig_root, pattern = "^survival-\\d{4}-\\d{2}-\\d{2}$", full.names = TRUE)
  survival_dir <- dirs[order(basename(dirs), decreasing = TRUE)][1]
}
message("Reading estimates from ", survival_dir)
if (!exists("clean_survival_dir")) {
  clean_dirs <- list.files(fig_root, pattern = "^survival-cleanref-haponly-\\d{4}-\\d{2}-\\d{2}(-\\d{6})?$", full.names = TRUE)
  if (!length(clean_dirs)) stop("Run run_fuel_clean_reference.R with full adjustment, or set clean_survival_dir to its validated full-adjustment output")
  clean_survival_dir <- tail(sort(clean_dirs), 1)
}
for (folder in c(survival_dir, clean_survival_dir)) {
  settings <- file.path(folder, "run_settings.txt")
  if (!file.exists(settings) || !any(grepl("^outcome_adjustment_mode: full$", readLines(settings))))
    stop("Figure 2 requires validated full-adjustment results; set survival_dir and clean_survival_dir to the corresponding full-adjustment folders: ", folder)
}
out_dir <- file.path(fig_root, paste0("manuscript-figures-", Sys.Date()))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

exposure_labels <- c(solid_fuel_use = "Any polluting fuel use",
                     wood_use = "Wood use",
                     kerosene_use = "Kerosene use",
                     solid_fuel_primary = "Polluting fuel as\nprimary cooking fuel")

est <- bind_rows(
  read_csv(file.path(survival_dir, "survival_hap_effects_summary.csv"), show_col_types = FALSE) %>%
    filter(!exposure %in% c("wood_use", "kerosene_use")),
  read_csv(file.path(clean_survival_dir, "survival_hap_effects_summary.csv"), show_col_types = FALSE) %>%
    filter(exposure %in% c("wood_use", "kerosene_use"))) %>%
  filter(measure == "Risk difference", exposure %in% names(exposure_labels)) %>%
  mutate(assay = factor(recode(cohort, qgit = "IGRA", tst = "TST"), levels = c("IGRA", "TST")),
         exposure = factor(exposure_labels[exposure], levels = rev(exposure_labels)),
         y = as.numeric(exposure) + if_else(assay == "IGRA", 0.14, -0.14),
         horizon = factor(recode(horizon, month4 = "A  By Month 4", month12 = "B  By Month 12"),
                          levels = c("A  By Month 4", "B  By Month 12")))

# Exposures as numeric rows (first exposure at the top), IGRA just above TST
p <- ggplot(est, aes(x = estimate, y = y, colour = assay, shape = assay)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey35", linewidth = 0.3) +
  geom_segment(aes(x = ci_low, xend = ci_high, yend = y), linewidth = 0.7) +
  geom_point(size = 2.2) +
  facet_wrap(~ horizon, nrow = 1) +
  scale_colour_manual(values = c(IGRA = "grey5", TST = "grey55"), name = NULL) +
  scale_shape_manual(values = c(IGRA = 16, TST = 15), name = NULL) +
  scale_y_continuous(breaks = seq_along(exposure_labels), labels = rev(exposure_labels)) +
  coord_cartesian(xlim = c(-0.3, 0.45)) +
  labs(x = "Risk difference (95% CI)", y = NULL) +
  theme_classic(base_size = 10) +
  theme(legend.position = "top",
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0, face = "bold", size = 10),
        panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3),
        axis.text.y = element_text(colour = "black"),
        panel.spacing = unit(1.2, "lines"))

ggsave(file.path(out_dir, "figure2_fuel_conversion.png"), p, width = 8.5, height = 3.9, dpi = 300)
ggsave(file.path(out_dir, "figure2_fuel_conversion.pdf"), p, width = 8.5, height = 3.9)
message("Saved to ", out_dir)
