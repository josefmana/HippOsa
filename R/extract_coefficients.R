#' Extracts linear model coefficients of Bayesian regressions
#'
#' From a Bayesian regression fit, extracts parameter estimate and with
#' it associated metrics.
#'
#' @param fit The regression model computed by \code{\link[brms]{brm}}
#'
#' @exports A table summarising models coefficient of interest
extract_coefficients <- function(fit) {
  map_dfr(names(fit), function(i) {
    do.call(rbind.data.frame, list(
      lm_coeff(fit[[i]], term = "SUBJ1", type = "Bayesian"),
      lm_coeff(fit[[i]], term = "AHI.F1", type = "Bayesian"),
      lm_coeff(fit[[i]], term = "SUBJ1:AHI.F1", type = "Bayesian")
    ))
  }) |>
    mutate(across(all_of(c("X","sigma","coefficient")), re_formulate))
}
