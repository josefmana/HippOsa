#' Reformat model formula for presentation
#'
#' Takes in a model formula and changes the terms to presentation-ready
#' format.
#'
#' @param form A character containing model formula
#'
#' @exports Formated formula for presentation
re_formulate <- function(form) {
  # Prepare a data frame with changes to be made:
  map <- data.frame(
    old = c("SUBJ", "AGE", "GENDER", "AHI.F", "EDU.Y", "tmt_b", "gpt_phk", "gpt_lhk", "sigma"),
    new = c("Group", "Age", "Sex", "OSA", "Education", "TMT-B", "GPT right", "GPT left", "\u03C3")
  )
  # Change it
  for(i in seq_len(nrow(map))) {
    form <- gsub(map$old[i], map$new[i], form, fixed = TRUE)
  }
  form
}
