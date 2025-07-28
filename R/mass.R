#' Extracts effects from a linear model
#'
#' Using the \code{\link[emmeans]{emmeans}} function extracts different types
#' of effect estimates from fitted linear models.
#'
#' @param fit A list of regression models computed by \code{lm} or
#' \code{\link[brms]{brm}}
#' @param type A character indicating whether only a "moderation" (default)
#' effects or all comparisons ("full") should be estimated
#'
#' @exports A table summarising models coefficient of interest
mass <- function(fit, type = "moderation") {
  map_dfr(names(fit), function(y) {
    # Moderation comparisons only:
    if (type == "moderation") {
      full_join(
        emmeans(fit[[y]], specs = pairwise ~ AHI.F | SUBJ) |>
          contrast(interaction = "consec") |>
          as_tibble() |>
          mutate(X = sub( paste0(y," ~ "), "", c(formula(fit[[y]]))), term = "AHI.F", .before = 1),
        emmeans(fit[[y]], specs = pairwise ~ AHI.F * SUBJ) |>
          contrast(interaction = "consec") |>
          as_tibble() |>
          mutate(X = sub( paste0(y," ~ "), "", c(formula(fit[[y]]))), .before = 1) |>
          rename("term" = "SUBJ_consec")
      ) |>
        as.data.frame() |>
        rename("contrast" = "AHI.F_consec") |>
        mutate(contrast = if_else(term == "PD - CON", NA, contrast)) |> # for compatibility with previous versions of the script
        mutate(y = y, .before = 1)
      # All comparisons:
    } else if (type == "full") {
      reduce(list(
        # 'main effects'
        emmeans(fit[[y]], specs = pairwise ~ SUBJ) |> contrast("consec") |> as_tibble() |> mutate(X = sub(paste0(y," ~ "), "", c(formula(fit[[y]]))), term = "SUBJ", .before = 1),
        emmeans(fit[[y]], specs = pairwise ~ AHI.F) |> contrast("consec") |> as_tibble() |> mutate(X = sub(paste0(y," ~ "), "", c(formula(fit[[y]]))), term = "AHI.F", .before = 1),
        # 'simple main effects'
        emmeans(fit[[y]], specs = pairwise ~ AHI.F | SUBJ) |> contrast("consec") |> as_tibble() |> mutate( X = sub(paste0(y," ~ "), "", c(formula(fit[[y]]))), term = "AHI.F", .before = 1),
        emmeans(fit[[y]], specs = pairwise ~ SUBJ | AHI.F) |> contrast("consec") |> as_tibble() |> mutate( X = sub(paste0(y," ~ "), "", c(formula(fit[[y]]))), term = "SUBJ", .before = 1),
        # interactions
        emmeans(fit[[y]], specs = pairwise ~ AHI.F * SUBJ) |> contrast(interaction = "consec") |> as_tibble() |> mutate(X = sub(paste0(y," ~ "), "", c(formula(fit[[y]]))), .before = 1 ) |> rename("term" = "SUBJ_consec"),
        emmeans(fit[[y]], specs = pairwise ~ SUBJ * AHI.F) |> contrast(interaction = "consec") |> as_tibble() |> mutate(X = sub(paste0(y," ~ "), "", c(formula(fit[[y]]))), .before = 1 ) |> rename("term" = "AHI.F_consec")
      ),
      full_join
      ) |>
        mutate(y = y, .before = 1) |>
        select(-ends_with("consec"))
    }
  })
}
