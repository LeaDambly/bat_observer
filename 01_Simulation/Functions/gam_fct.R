indsp_func <- function(bats, dfvec = 6, q = 1) {
  # All rows with NA's (missing values) must be removed from the data before applying this function.
  if (length(bats$count[is.na(bats$count)]) > 0) {
    stop("no missing data allowed")
  }

  # fit the GAM using MGCV - get rid of any roosts that were only surveyed for one year
  bats <- bats %>%
    group_by(roost) %>%
    dplyr::filter(n() > 1) %>%
    arrange(year) %>%
    ungroup()

  gam_df <- gam(count ~ as.factor(roost) + s(year, fx = TRUE, k = dfvec + 1),
    family = poisson(link = log), data = bats
  )

  # Takes a fitted gam object produced by gam() and produces predictions given a new set of values for
  # the model covariates or the original values used for the model fit. Predictions can be accompanied
  # by standard errors, based on the posterior distribution of the model coefficients.
  fit <- as_tibble(predict(gam_df, type = "terms", se.fit = TRUE)$fit) %>%
    select("s(year)") %>%
    rename(smooth = "s(year)")

  se <- as_tibble(predict(gam_df, type = "terms", se.fit = TRUE)$se.fit) %>%
    select("s(year)") %>%
    rename(error = "s(year)")

  # calculate confidence intervals
  lcl <- as_tibble(fit - 1.96 * se) %>%
    rename(lower = smooth)

  ucl <- as_tibble(fit + 1.96 * se) %>%
    rename(upper = smooth) %>%
    mutate(row_name = row_number())

  # combine them to one tibble
  predict <- bind_cols(fit, lcl, ucl, se)

  # entr gives the entries corresponding to the first occurrence of each year in the tibble
  entr <- bats %>%
    mutate(row_name = row_number()) %>%
    dplyr::filter(!duplicated(year)) %>%
    select(year, row_name)

  predict <- predict %>%
    dplyr::filter(predict$row_name %in% entr$row_name) %>%
    add_column(year = entr$year) %>%
    mutate(
      index = exp(smooth) / exp(smooth[q]),
      upper = exp(upper) / exp(smooth[q]),
      lower = exp(lower) / exp(smooth[q])
    ) %>%
    select(year, index, upper, lower) %>%
    arrange(year)

  return(predict)
}
