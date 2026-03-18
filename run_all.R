## -----------------------------------------------------------------------------
## Author: Francisco Garre-Frutos
## -----------------------------------------------------------------------------

# Edit to run specific analysis scripts
force_main <- FALSE
force_comp <- FALSE
force_reanalysis <- FALSE

dirs <- c("output/analysis", "output/figures", "output/models", "output/manuscript")
for (dir in dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
}


# Results from Experiment 1
if (force_main || !file.exists("output/analysis/analysis_e1.rds")) {
  message("Running Experiment 1 analysis...")
  source("scripts/analysis_e1.R", local = FALSE)
} else {
  message("Using cached Experiment 1 analysis: output/analysis/analysis_e1.rds")
}

# Results from Experiment 2a and 2b
if (force_main || !file.exists("output/analysis/analysis_e2.rds")) {
  message("Running Experiment 2 analysis...")
  source("scripts/analysis_e2.R", local = FALSE)
} else {
  message("Using cached Experiment 2 analysis: output/analysis/analysis_e2.rds")
}

# Re-analysis of Le et al. (2024)
if (force_reanalysis || !file.exists("output/figures/fig-7.png")) {
  message("Running Le et al. (2024) re-analysis...")
  source("scripts/re_analysis_Le2024.R", local = FALSE)
} else {
  message("Using cached Le et al. (2024) re-analysis: output/figures/fig-7.png")
}

# Computational modelling
if (force_comp || !file.exists("output/models/model_comparison_results.rds")) {
  message("Running Bayesian model comparison...")
  source("scripts/run_model_comparison.R", local = FALSE)
} else {
  message("Using cached Bayesian model comparison: output/models/model_comparison_results.rds")
}

# PPC plot.
message("Running posterior predictive plots...")
source("scripts/posterior_predictions.R", local = FALSE)

message("Pipeline completed.")

