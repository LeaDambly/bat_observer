# TO DO: describe function and parameters it takes

adss <- function(obs, a = 0, b = 1, c = 1) {
  # Only take the first year that roosts are monitored in
  # This is because, if we use this function on all years, any kind of 0 would be removed
  # It also makes more sense that once a roost is known, it'll stay known
  det <- obs %>% 
    group_by(roost) %>% 
    dplyr::filter(count > 0) %>%
    top_n(year, n = -1) %>%
    ungroup()
  
  # create minimum (a) and maximum (b) values for accessibility rating
  seqq <- seq(a, b, 1)
  
  # calculate slope: m = (y2 - y1) / (x2 - x1)
  slope <- (max(seqq) - min(seqq)) / (max(det$count) - min(det$count))
  
  # calculate intercept: b = (m * x1 - y1) * (-1)
  int <- (slope * min(det$count) - min(seqq)) * (-1)
  
  # Create detectability 'score'
  det <- det %>% 
    rowwise(.) %>% 
    mutate(score = ((slope * count) + int) ^ c) %>%
    arrange(score)
  
  nroost <- (length(det$roost)/100) * 75
  
  det <- det[sample(nrow(det), nroost, prob = det$score),]
  
  # now only keep the roosts that were sampled continuing from the year in det
  final <- obs %>%
    dplyr::filter(roost %in% det$roost) %>%
    left_join(., det, by = "roost") %>%
    dplyr::filter(year.x >= year.y) %>%
    select(roost, year.x, count.x) %>%
    rename(year = year.x, count = count.x)
  
  return(final)
}  