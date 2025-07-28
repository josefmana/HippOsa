#' Fit a set of linear regression models
#'
#' Loops through outcomes to calculate a set of linear regressions.
#'
#' @param df The data frame
#' @param help A list with helpers files prepared by \code{extract_helpers}
#'
#' @exports A list with linear models, one per outcome
fit_regressions <- function(df, help){
  # Extract helpers:
  for (i in names(help)) {
    assign(i, help[[i]])
  }
  # Set-up formulas:
  forms <- data.frame(
    object = c("subco","hippo","psych"),
    y = c("subcortical", "hippocampi","cognition"),
    X = c(rep("SUBJ * AHI.F + AGE + GENDER + SBTIV", 2), "SUBJ * AHI.F + AGE + GENDER + EDU.Y")
  )
  # Loop through types of regressions (base vs interaction), and structures
  lapply(set_names(seq_len(nrow(forms)), forms$y), function(r) {
    with(get(forms[r, "object"]), {
      y <- forms[r, "y"] # extract outcome
      if(y == "cognition") { # extract variables
        vars <- variable
      } else {
        vars <- unique(name)
      }
      fit_reg(df, vars, X = forms[r, "X"], w = FALSE) # fit the models
    })
  })
}
