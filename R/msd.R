#' Prints mean ± standard deviation
#'
#' Takes in a numeric vector, calculates mean and standard deviation (SD) and
#' returns these formatted.
#'
#' @param x A vector from which mean and SD ought to be calculated
#' @param d Decimals (defaults to 2)
#'
#' @exports A character with mean and SD
msd <- function(x, d = 2) {
  paste0(rprint(mean(x, na.rm = TRUE), d), " ± ", rprint(sd(x, na.rm = TRUE), d))
}
