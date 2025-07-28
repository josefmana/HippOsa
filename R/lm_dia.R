#' Runs residual model diagnostics
#'
#' From a linear regression fit, extracts measures of the Breusch-Pagan test
#' for homoscedasticity, number of outliers according to the Cook distance,
#' and the Shapiro-Wilk test for normality.
#'
#' @param fit The linear regression model computed by \code{lm}
#'
#' @exports A table summarising model diagnostics
lm_dia <- function(fit) {
  sapply(names(fit), function(y) {
    data.frame(
      X = sub(paste0(y," ~ "), "", c(formula(fit[[y]]))),
      p_breusch_pagan = c(check_heteroscedasticity(fit[[y]])),
      n_cook = sum(check_outliers(fit[[y]]), na.rm = TRUE),
      p_shapiro_wilk = c(check_normality(fit[[y]]))
    )
  }) |>
    t() |>
    as.data.frame() |>
    mutate(across(everything(), ~unlist(.x, use.names = FALSE))) |>
    mutate(heteroscedasticity = ifelse(p_breusch_pagan < .05, "!", ""), .after = p_breusch_pagan) |>
    mutate(nonnormality = ifelse(p_shapiro_wilk < .05, "!", ""), .after = p_shapiro_wilk) |>
    rownames_to_column("y")
}
