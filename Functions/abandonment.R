# abandon function: simulates the 'abandonment' of a roost by its population
# This is when the entire population leaves its monitored roost to go to a new, not yet monitored roost
#
# Takes the following parameters
# pop: a tibble that contains at least count, year, roost
# p: the fraction of roosts that are abandoned
#
# It returns pop: a tibble with columns count, roost, year, abyr, orgroost
# abyr records the year in which the population abandoned the original roost
# orgroost records a population's original roost that it left to move to a new roost

abandon <- function(pop, p = 0.1) {

  # create new tibble with populations that have abandoned their roost
  pop2 <- pop %>%
    dplyr::filter(roost <= z & count > 0) %>%
    sample_frac(p, replace = FALSE) %>%
    arrange(roost)

  # create new tibble without the populations that have abandoned their roost
  pop3 <- anti_join(pop, pop2, by = c("count", "roost", "year", "abyr")) %>%
    arrange(roost)

  # create new tibble where count is zero (because the original roosts are now empty)
  pop4 <- pop2 %>%
    select(roost, year, abyr) %>%
    add_column(count = 0)

  # add pop4 back to pop3 - now the original roosts are empty
  pop3 <- pop3 %>%
    bind_rows(pop4) %>%
    arrange(roost)

  # rename roost to orgroost (original roost) - new roost sites will have a 'memory' of where the population came from
  pop2 <- pop2 %>%
    select(count, roost, abyr, year) %>%
    rename(orgroost = roost)

  # add the newly establised roosts back to the entire tibble
  pop <- pop3 %>%
    bind_rows(pop2)

  # number of new populations without respective roost number
  num <- sum(is.na(pop$roost))

  # the new roosts have NAs, we need to make them 0
  pop$roost[is.na(pop$roost)] <- 0

  # now we index the new roosts
  pop <- pop %>%
    mutate(roost = ifelse(roost == 0, seq(
      max(pop$roost) + 1, max(pop$roost) + num, 1
    ), roost)) %>%
    arrange(roost) %>%
    # and give the original z roosts an NA for the orgroost variable
    mutate(orgroost = ifelse(roost <= z, NA, orgroost)) %>%
    # finally add another year to abyr
    mutate(abyr = ifelse((roost > z & count > 0), abyr + 1, 0))

  return(pop)
}
