#' Print numeric as a character
#'
#' Takes in a number and prints it with specified decimal places.
#'
#' @param x A numeric to be printed
#' @param d Decimals (defaults to 2)
#'
#' @exports A character with the number
rprint <- function(x, d = 2) {
  sprintf(paste0("%.", d, "f"), round(x, d))
}
