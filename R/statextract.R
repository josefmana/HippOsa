#' Extracts and formats test statistics
#'
#' Takes in results of a regression analysis and extracts regression coefficients
#' with accompanying test statistics and p-values.
#'
#' @param coeffs Regression coefficients
#' @param y Outcome of interest
#' @param stat A character indicating whether t-values ("t", default) or
#' z-values ("z") should be looked for.
#'
#' @exports A tibble containing regression coefficients together with
#' test statistics and p-values associated with them.
statextract <- function(coeffs, y, stat = "t") {
  coeffs |>
    as.data.frame() |>
    mutate(y = y, .before = 1) |>
    rownames_to_column("coefficient") |>
    mutate(
      out = paste0(
        stat," = ", rprint(get(paste0(stat, " value")), 2),
        ", p ",
        if_else(
          get(paste0("Pr(>|", stat, "|)")) < .001,
          zerolead(get(paste0("Pr(>|", stat, "|)"))),
          paste0("= ", zerolead(get(paste0("Pr(>|", stat, "|)"))))
        )
      )
    ) |>
    select(-paste0(stat, " value"), -paste0("Pr(>|", stat, "|)") ) |>
    pivot_wider(names_from = coefficient, values_from = out) |>
    select(-`(Intercept)`)
}
