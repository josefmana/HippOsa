#' Extracts a linear model coefficient
#'
#' From a linear regression fit, extracts parameter estimate and with
#' it associated metrics.
#'
#' @param fit The regression model computed by \code{lm} or \code{\link[brms]{brm}}
#' @param term Regression coefficient of interest
#' @param type A character indicating whether it is a "frequentist" (default)
#' or "Bayesian" model being evaluated
#'
#' @exports A table summarising models coefficient of interest
lm_coeff <- function(fit, term = "SUBJ1:AHI.F1", type = "frequentist") {
  # Coefficients from a frequentist model:
  if (type == "frequentist") {
    sapply(names(fit), function(y) {
      summary(fit[[y]])$coefficients[term, ] |>
        t() |>
        as.data.frame() |>
        rename("p value" = "Pr(>|t|)") |>
        cbind(t( confint(fit[[y]])[term, ])) |>
        relocate(`2.5 %`, .before = `t value`) |>
        relocate(`97.5 %`, .before = `t value`) |>
        mutate(X = sub(paste0(y," ~ "), "", c(formula(fit[[y]]))), .before = 1)
    }) |>
      t() |>
      as.data.frame() |>
      mutate(coefficient = term, .after = X) |>
      mutate_if(is.list, unlist) |>
      rownames_to_column("y") |>
      mutate(
        `q value` = p.adjust(`p value`, method = "BH"),
        `s value` = -log(`p value`, base = 2),
        sig_PCER = if_else(`p value` < .05, "*", ""),
        sig_FDR = bh_adjust(`p value`),
        sig_FWER = if_else(`p value` < .05/n(), "*", "")
      )
    # Coefficients from a Bayesian model:
  } else if (type == "Bayesian") {
    sapply(names(fit), function(y) {
      fixef(fit[[y]])[term, ] |>
        t() |>
        as.data.frame() |>
        mutate(
          X = sub(".* ~ ", "", as.character(formula(fit[[y]]))[1]),
          sigma = if_else(
            grepl("sigma", formula(fit[[y]])[2]),
            sub(")", "", sub( ".* ~ ", "", as.character(formula(fit[[y]]))[2])),
            "1"
          ),
          .before = 1
        )
    }) |>
      t() |>
      as.data.frame() |>
      mutate(coefficient = term, .after = sigma) |>
      mutate_if(is.list, unlist) |>
      rownames_to_column("y")
  }
}
