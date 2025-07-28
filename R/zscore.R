#' Computes z-scores for neuropsychological tests
#'
#' Takes in linear calculator values, observed score, name of the scale,
#' age, gender and education and returns a vector of z-scores.
#'
#' @param calc A data.frame with calculator values
#' @param x A vector with test scores to be evaluated
#' @param nam A data.frame containing mapping of variable names
#' @param AGE A vector of observed age values in years
#' @param GEN A vector of observed gender values (female = 0, male = 1)
#' @param EDU A vector with observed education values in years
#'
#' @exports A vector of z-scores
zscore <- function(calc, x, nam, AGE, GEN, EDU) {
  with(calc, {
    # Prepare a matrix of data and vector of parameters:
    pars <- as.numeric(c( Constant[var == nam], age[var == nam], gender[var == nam], education[var == nam]))
    data <- as.matrix(cbind(rep(1, length(x)), AGE, GEN, EDU))
    # Compute linear predictions:
    x_bar <- data %*% pars
    # Return z-score:
    c(sign[var == nam] * (x - x_bar) / RMSE[var == nam])
  })
}
