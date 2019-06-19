# TO DO: describe function and parameters it takes

commitment <- function(obs, pr = 0.3) {
  mot <- obs %>%
    group_by(roost) %>%
    mutate(prob = rbinom(length(n()), 3, pr)) %>%
    ungroup()

  # no errors - it just gives a warning if there isn't anything
  # that fulfills the conditions,
  # which is absolutely fine but annoying
  suppressWarnings(
    zer <- mot %>%
      group_by(roost) %>%
      dplyr::filter(prob == 0) %>%
      dplyr::filter(row_number() <= min(which(count == 0 & dplyr::lag(count, n = 0) == 0))) %>%
      ungroup()
  )

  suppressWarnings(
    one <- mot %>%
      group_by(roost) %>%
      dplyr::filter(prob == 1) %>%
      dplyr::filter(row_number() <= min(which(count == 0 & dplyr::lag(count, n = 1) == 0))) %>%
      ungroup()
  )

  suppressWarnings(
    two <- mot %>%
      group_by(roost) %>%
      dplyr::filter(prob == 2) %>%
      dplyr::filter(row_number() <= min(which(count == 0 & dplyr::lag(count, n = 2) == 0))) %>%
      ungroup()
  )

  suppressWarnings(
    thr <- mot %>%
      group_by(roost) %>%
      dplyr::filter(prob == 3) %>%
      dplyr::filter(row_number() <= min(which(count == 0 & dplyr::lag(count, n = 3) == 0))) %>%
      ungroup()
  )

  mot <- rbind(zer, one, two, thr)


  return(mot)
}
