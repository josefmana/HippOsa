#' Summarise a set of linear regression models
#'
#' Loops through linear regression fits and summarise their parameters.
#'
#' @param fit A list of linear regression models
#'
#' @exports A list with summary of linear models, one per model
summarise_regressions <- function(fit) {
  lapply( set_names(names(fit)), function(y) {
    # regression coefficients with threshold-based decisions
    # do full analysis for subcortical structures and cognition, interaction only for hippocampal substructures
    if (y == "hippocampi") {
      lm_coeff(fit[[y]], term = "SUBJ1:AHI.F1") |>
        left_join(lm_dia(fit[[y]]), by = c("y","X")) |>
        mutate(across(all_of(c("X","coefficient")), re_formulate))
    } else {
      left_join(
        rbind.data.frame(
          lm_coeff(fit[[y]], term = "SUBJ1"),
          lm_coeff(fit[[y]], term = "AHI.F1"),
          lm_coeff(fit[[y]], term = "SUBJ1:AHI.F1")
        ),
        lm_dia(fit[[y]]),
        by = c("y","X")
      ) |>
        mutate(
          across(all_of( c("X","coefficient") ), re_formulate),
          sig_FDR = bh_adjust(`p value`) # re-calculate Benjamini-Hochberg adjusted significance statements
        )
    }
  })
}
