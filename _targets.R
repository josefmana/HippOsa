# This is a script running targets pipeline of the project.

# Load packages required to define the pipeline:
library(targets)

# Set target options:
tar_option_set( packages = c(
  
  "here", # for path listing
  "openxlsx", # for data reading
  "tidyverse", # for data wrangling
  "performance", # for regression diagnostics 
  "emmeans" # for models' marginal means estimation
  
) )

# Load all in-house functions:
tar_source()

# Replace the target list below with your own:
list(
  
  # read data ----
  tar_target(
    
    files = data.frame(
      lhippo = here("_raw", "Tabhipposubfields_lhx.xlsx"),
      rhippo = here("_raw", "Tabhipposubfields_rhx.xlsx"),
      subcor = here("_raw", "Tabsubcortex1_RBD.xlsx"),
      psych = here("_raw", "RBDBIOPDCON_DATA_2024-07-17_1146.csv"),
      motor = here("_raw", "BIOPD_MDSUPDRSIII.xlsx")
    ),

    helpers = data.frame(
      psychvar = here("helpers","psychs.csv"),
      calculat = here("_raw","calculator_final_v7_c_301116.xlsx"),
      calc_sheet = "equations"

    )
  ),
  
  # 
  
)
