#' Prepare formulas for regressions
#'
#' Lists formulas for Bayesian models in the sensitivity-to-heteroscedasticity
#' analysis.
#'
#' @exports A list with formulas
set_formulas <- function() {
  list(
    # Formulas for base models (assuming homoscedasticity):
    varequal = lapply(set_names(c("tmt_b", "gpt_phk", "gpt_lhk")), function(i) {
      paste0(i, " ~ SUBJ * AHI.F + AGE + GENDER + EDU.Y") |>
        as.formula() |>
        bf()
    }),
    # Formulas for variance adjusted models (allowing heteroscedasticity)
    # (listing them one by one because as.formula does not want work with commas):
    heteroscedastic = list(
      tmt_b = bf(tmt_b ~ SUBJ * AHI.F + AGE + GENDER + EDU.Y, sigma ~ SUBJ * AHI.F + AGE + GENDER + EDU.Y),
      gpt_phk = bf(gpt_phk ~ SUBJ * AHI.F + AGE + GENDER + EDU.Y, sigma ~ SUBJ * AHI.F + AGE + GENDER + EDU.Y),
      gpt_lhk = bf(gpt_lhk ~ SUBJ * AHI.F + AGE + GENDER + EDU.Y, sigma ~ SUBJ * AHI.F + AGE + GENDER + EDU.Y)
    )
  )
}
