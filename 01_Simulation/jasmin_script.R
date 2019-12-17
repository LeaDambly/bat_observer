library(foreach)
library(doParallel)

numCores <- detectCores()
cl <- makeCluster(numCores)
registerDoParallel(cl)

ab_high <- foreach(i = 1:500, .packages = c("tibble", "dplyr", "mgcv")) %dopar% {
  srcf <- "/home/users/leadam/bat_sim/scripts/functions"
  R.utils::sourceDirectory(srcf)
  z <- 400
  x <- 72
  y <- 1
  
  N0 <- rpois(lam = x, n = z)
  
  bats_1 <- tibble(count = N0) %>%
    add_column(roost = 1:z, year = y, abyr = NA, reocc = NA)
  
  y <- y + 1
  
  pogr <- bats_1
  
  while (y <= 20) {
    pogr <- pogr %>%
      mutate(year = year + 1)
    
    pogr <- pogr %>%
      group_by(roost) %>%
      mutate(count = growth(Nt = count)) %>%
      ungroup()
    
    pogr <- reoccupy(pogr, p2 = 0.7)
    pogr <- abandon(pogr, p = 0.22)
    
    assign(paste("bats", y, sep = "_"), pogr)
    
    y <- y + 1
    gc()
  }
  
  rm(pogr)
  
  all <- mget(ls(pattern = "bats_*")) %>%
    bind_rows() %>%
    as_tibble()
  
  rm(list = ls(pattern = "bats_*"))
  
  # instead of acast and melt - to 'fill in' zeroes for all roosts
  temp <- enframe(rep(seq(1, length(unique(all$roost)), 1), times = 20), name = NULL) %>%
    rename(roost = value) %>%
    arrange(roost) %>%
    add_column(year = rep(seq(1, y-1, 1), times = length(unique(all$roost)))) %>%
    add_column(count = 0)
  
  act <- right_join(all, temp, by = c("roost", "year")) %>%
    mutate(count = ifelse(is.na(count.x), count.y, count.x)) %>%
    select(roost, count, year, abyr, reocc, orgroost)
  
  rm(all)
  rm(temp)
  
  adss_low <- adss2(act, frac = 0.2)
  adss_med <- adss2(act, frac = 0.5)
  adss_high <- adss2(act, frac = 0.8)
  
  adss_low_com_high <- commitment2(adss_low, pr = 0.8)
  adss_med_com_high <- commitment2(adss_med, pr = 0.8)
  adss_high_com_high <- commitment2(adss_high, pr = 0.8)
  
  adss_low_com_med <- commitment2(adss_low, pr = 0.5)
  adss_med_com_med <- commitment2(adss_med, pr = 0.5)
  adss_high_com_med <- commitment2(adss_high, pr = 0.5)
  
  adss_low_com_low <- commitment2(adss_low, pr = 0.2)
  adss_med_com_low <- commitment2(adss_med, pr = 0.2)
  adss_high_com_low <- commitment2(adss_high, pr = 0.2)
  
  inda <- indsp_func(act) %>%
    mutate(type = "Actual")
  
  inda <- inda %>%
    mutate(rmse = sqrt(mean((inda$index - index)^2)))
  
  indb <- indsp_func(adss_low_com_high) %>%
    mutate(type = "ADSS low commitment high") %>%
    mutate(rmse = sqrt(mean((inda$index - index)^2)))
  
  indc <- indsp_func(adss_med_com_high) %>%
    mutate(type = "ADSS med commitment high") %>%
    mutate(rmse = sqrt(mean((inda$index - index)^2)))
  
  indd <- indsp_func(adss_high_com_high) %>%
    mutate(type = "ADSS high commitment high") %>%
    mutate(rmse = sqrt(mean((inda$index - index)^2)))
  
  inde <- indsp_func(adss_low_com_med) %>%
    mutate(type = "ADSS low commitment med") %>%
    mutate(rmse = sqrt(mean((inda$index - index)^2)))
  
  indf <- indsp_func(adss_med_com_med) %>%
    mutate(type = "ADSS med commitment med") %>%
    mutate(rmse = sqrt(mean((inda$index - index)^2)))
  
  indg <- indsp_func(adss_high_com_med) %>%
    mutate(type = "ADSS high commitment med") %>%
    mutate(rmse = sqrt(mean((inda$index - index)^2)))
  
  indh <- indsp_func(adss_low_com_low) %>%
    mutate(type = "ADSS low commitment low") %>%
    mutate(rmse = sqrt(mean((inda$index - index)^2)))
  
  indi <- indsp_func(adss_med_com_low) %>%
    mutate(type = "ADSS med commitment low") %>%
    mutate(rmse = sqrt(mean((inda$index - index)^2)))
  
  indj <- indsp_func(adss_high_com_low) %>%
    mutate(type = "ADSS high commitment low") %>%
    mutate(rmse = sqrt(mean((inda$index - index)^2)))
  
  allind <- bind_rows(inda, indb, indc, indd, inde, indf, indg, indh, indi, indj)
  allind$type <- as.factor(allind$type)
  allind$rep <- i
  
  return(allind)
  gc()
}

save(ab_high, file = "/home/users/leadam/bat_sim/ab_high.RData")
stopCluster(cl)
