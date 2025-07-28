# This is a script running targets pipeline of the project.

# Load packages required to define the pipeline:
library(targets)

# Set target options:
tar_option_set( packages = c(
  "here", # for path listing
  "openxlsx", # for data reading
  "tidyverse", # for data wrangling
  "psych", # for data description
  "emmeans", # for avg_comparisons in boxplots
  "performance", # for regression diagnostics
  "ggh4x", # for improved facet_grid
  "ggpubr", # for plotting boxplots with p values
  "brms", # for Bayesian regressions
  "bayesplot", # for posterior predictive checks
  "priorsense", # for checking sensitivity to powerscaling
  "patchwork" # for plots arranging
))

# Load all in-house functions:
tar_source()

# Use multiple cores for model fitting:
options( mc.cores = parallel::detectCores() )

# Replace the target list below with your own:
list(
  # Read data ----
  tar_target(files, unlist(data_paths()), format = "file"), # For tracking purposes
  tar_target(datafiles, data_paths()),
  tar_target(helpers, extract_helpers()),
  tar_target(raw_data, import_data(datafiles, helpers)),
  tar_target(rt_variables, extract_rt_variables(helpers)),
  tar_target(preprocessed_data, preprocess_data(raw_data, helpers, rt_variables, return = "df")),
  tar_target(scales, preprocess_data(raw_data, helpers, rt_variables, return = "scl")),

  # Data analysis ----
  ## ---- Description ----
  tar_target(descriptives, describe_data(raw_data)),
  ## ---- Regression ----
  tar_target(regressions, fit_regressions(preprocessed_data, helpers)),
  tar_target(summaries, summarise_regressions(regressions)),
  tar_target(formulas, set_formulas()),
  tar_target(bayesian_regressions, fit_bayesian(preprocessed_data, formulas)),
  tar_target(prior_sensitivity, model_check(bayesian_regressions, helpers, formulas, "prior_sense")),
  tar_target(posterior_predictive_checks, model_check(bayesian_regressions, helpers, formulas, "ppc")),
  tar_target(bayesian_coefficients, extract_coefficients(bayesian_regressions)),
  ## ---- Visualisation ----
  tar_target(boxplot_brains, boxplots(raw_data, preprocessed_data, regressions$subcortical, helpers, scales, rt_variables, which = "brains")),
  tar_target(boxplot_cognition_osa, boxplots(raw_data, preprocessed_data, regressions$cognition, helpers, scales, rt_variables, which = "cognition_1")),
  tar_target(boxplot_cognition_pd, boxplots(raw_data, preprocessed_data, regressions$cognition, helpers, scales, rt_variables, which = "cognition_2"))

)
