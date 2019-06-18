# reoccupy function: simulates the reoccupation of a roost by its former population
#
# Takes the following parameters
  # pop: a tibble with count, roost, year, abyr, orgroost
  # p2: the fraction of populations that have previously abandoned their original roost and are now moving back
#
# It returns pop: a tibble with columns count, roost, year, abyr, orgroost

reoccupy <- function(pop, p2 = 0.7) {
  # create new tibble with the populations that are coming back - it has to have been 1 yr since the abandonment (this is based on BCT report on abandonments for pipistrelles)
  pop2 <- pop %>%
    dplyr::filter(!is.na(orgroost) & abyr >= 1) %>%
    sample_frac(p2)
  
  # if there aren't any to move back, the function will end
  if (dim(pop2)[1] == 0) {
    return(pop)
  } else{
    # create new tibble where populations that will move back are removed
    pop3 <- anti_join(pop, pop2, by = c("roost", "count", "year", "orgroost", "abyr"))
    
    # make the orgroosts their inhabited roosts again
    pop2 <- pop2 %>%
      mutate(roost = orgroost) %>%
      select(-orgroost)
    
    # and put everything back together
    pop <- left_join(pop3, pop2, by = "roost") %>%
      select(-c(year.y, abyr.y)) %>%
      mutate(count = ifelse(is.na(count.y), count.x, count.y)) %>%
      rename(year = year.x) %>%
      rename(abyr = abyr.x) %>%
      select(count, roost, year, abyr, orgroost)
    
    return(pop)
  }
}