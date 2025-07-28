#' Fit linear regression models
#'
#' Loops through outcomes to calculate a set of linear regressions.
#'
#' @param d The data
#' @param outcomes A character vector of the outcomes of interest
#' @param X The right side of the linear model
#' @param w Indicator whether regression weight shall be used (TRUE) or not(FALSE, default)
#'
#' @exports A list with linear models, one per outcome
fit_reg <- function(d, outcomes, X = "SUBJ * AHI.F + AGE + GENDER + SBTIV", w = FALSE) {
  lapply(set_names(outcomes), function(y) {
    if (w == TRUE) {
      lm(formula = as.formula(paste(y, X, sep = " ~ ")), data = d, weights = weights)
    } else {
      lm(formula = as.formula(paste(y, X, sep = " ~ ")), data = d, weights = NULL)
    }
  })
}
