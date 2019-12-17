adss <- function(obs, frac = 0.5, yr = 20) {

  # kick out any roost that wasn't part of the original set
  obs <- obs %>%
    dplyr::filter(roost <= z)
  
  # get roosts which are selected based on their max abundance  
    maxab <- obs %>%
      group_by(roost) %>%
      arrange(roost, year) %>%
      dplyr::filter(count == max(count)) %>% # potential for duplicated values so take earliest year only
      slice(1) %>%
      rename(fstyear = year) %>%
      ungroup() %>%
      sample_frac(., size = frac) %>%
      arrange(roost)
    
    # get the entire time series for roosts that were selected in maxab
    maxall <- left_join(maxab, obs, by = "roost") %>%
      group_by(roost) %>%
      dplyr::filter(year >= fstyear) %>%
      select(roost, year, count.y, abyr.y) %>%
      rename(count = count.y, abyr = abyr.y) %>%
      ungroup()
    
    # get roosts which aren't selected based on their max abundance
    lftovr <- obs %>%
      group_by(roost) %>%
      dplyr::filter(!roost %in% maxab$roost) %>%
      ungroup()
    
    # get leftover roosts and assign them random years at which their sampling starts
    vec <- lftovr %>%
      dplyr::filter(count != 0) %>%
      group_by(roost) %>%
      sample_n(1) %>%
      select(year, roost) %>%
      rename(fstyr = year) %>%
      ungroup()

    lftovr <- left_join(vec, obs, by = "roost") %>%
      dplyr::filter(year >= fstyr) %>%
      select(roost, year, count, abyr)
      
    adss <- bind_rows(lftovr, maxall) %>%
      arrange(roost, year)
    
    
    return(adss)
}
