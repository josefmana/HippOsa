#' Document path to data files
#'
#' This function takes no parameter but needs to be adjusted manually if desired.
#' It is used in the \pkg{targets} pipeline to check for changes in data.
#'
#' @exports A data.frame containing paths to data
data_paths <- function() {
  data.frame(
    lhippo = here("_raw", "Tabhipposubfields_lhx.xlsx"),
    rhippo = here("_raw", "Tabhipposubfields_rhx.xlsx"),
    subcor = here("_raw", "asegTab.xlsx"),
    psych = here("_raw", "RBDBIOPDCON_DATA_2024-07-17_1146.csv"),
    motor = here("_raw", "BIOPD_MDSUPDRSIII.xlsx")
  )
}
