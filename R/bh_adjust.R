#' Adjust p-values for multiple comparisons
#'
#' Evaluate a vector of p-values for statistical significance on 5% False Discovery
#' Rate according to the Benjamini-Hochberg procedure.
#'
#' @param p A vector of p-values
#'
#' @exports A vector indicating statistical significance after adjustment.
bh_adjust <- function(p) {
  # Extract threshold:
  bh_thres <- data.frame(
    p = sort(p), # sort p values from smallest to largest
    thres = .05 * (seq_along(p)) / length(p) # prepare BH thresholds for each p value
  ) |>
    mutate(sig = if_else(p <= thres, TRUE, FALSE)) |>
    filter(sig == TRUE) |>
    select(thres) |>
    max()
  # Return stars based on this threshold:
  if_else(p < bh_thres, "*", "")
}
