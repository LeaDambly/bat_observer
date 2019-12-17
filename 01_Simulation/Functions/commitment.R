commitment <- function(obs, pr = 0.2) {
  
  # calculate after how many years that roost would be dropped
  mot <- obs %>%
    group_by(roost) %>%
    mutate(prob = rbinom(length(n()), 2, pr) + 1) %>%
    ungroup()
  
  # and now drop any once their abyr corresponds with prob
  mot2 <- mot %>%
    mutate(drop = ifelse(lag(prob) <= lag(abyr), year, NA)) %>%
    group_by(roost) %>%
    tidyr::fill(drop) %>%
    ungroup() %>%
    dplyr::filter(is.na(drop))
  
  return(mot2)
}
