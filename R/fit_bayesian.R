#' Fit Bayesian regression models
#'
#' Loops through outcomes to calculate a set of linear regressions.
#'
#' @param df The data
#' @param formulas A two-layered list with formulas for Bayesian regression.
#'
#' @exports A list with regressions, one per outcome
fit_bayesian <- function(df, formulas) {
  # Re-fit the models for TMT-B and GPT via brms with differing variance for the groups:
  priors <- list(
    # Priors for base models (assuming homoscedasticity):
    varequal = c(
      prior(normal(0, 1), class = Intercept),
      prior(normal(0, 1), class = b),
      prior(exponential(1), class = sigma)
    ),
    # Priors for variance adjusted models (allowing heteroscedasticity):
    heteroscedastic = c(
      prior(normal(0, 1), class = Intercept),
      prior(normal(0, 1), class = b),
      prior(normal(0, 1), class = Intercept, dpar = sigma),
      prior(normal(0, 1), class = b, dpar = sigma)
    )
  )
  # Fit the models:
  lapply(set_names(names(formulas)), function(i) {
    lapply(set_names(names(formulas[[i]])), function(j) {
      brm(
        formula = formulas[[i]][[j]],
        data = df,
        family = gaussian(link = "identity"),
        prior = priors[[i]],
        seed = 87542
      )
    })
  })
}
