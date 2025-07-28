#' Extracts variable names of response times indexes
#'
#' Takes in helpers files, finds neuropsychology helper and identifies tests in
#' domains that are based on response times. The function is very dataset specific!
#'
#' @param help A list containing data.frame named psych with a proper study-specific
#' structure.
#'
#' @exports A vector of variable names.
extract_rt_variables <- function(help) {
  with(help, {
    psych |>
      filter(domain %in% c("Attention", "Executive function", "Processing speed")) |>
      pull(variable)
  })
}
