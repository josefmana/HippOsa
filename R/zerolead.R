#' Get rid of the leading zero
#'
#' Takes in a number and prints it with specified decimal places without the
#' leading zero.
#'
#' @param x A numeric to be stripped of zeor
#' @param d Decimals (defaults to 3)
#'
#' @exports A character with the number
zerolead <- function(x, d = 3) {
  ifelse(x < .001, "< .001", sub("0.", ".", rprint(x, 3), fixed = TRUE))
}
