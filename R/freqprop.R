#' Prints frequency (proportion)
#'
#' Takes in a character/factor vector, calculates frequency and proportion of its
#' groups and returns these formatted.
#'
#' @param x A vector from which frequencies and proportions ought to be calculated
#' @param d Decimals (defaults to 0)
#'
#' @exports A character with frequencies and proportions
freqprop <- function(x, d = 0) {
  paste0(table(x)[2], " (", rprint(100*prop.table(table(x))[2], d = d), "%)")
}
