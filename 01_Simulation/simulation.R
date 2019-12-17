# a script to simulate virtual observers monitoring virtual bat populations
# here we have four different 'types' of bat populations being monitored
# first, a bat population is simulated with a set inter-annual variation in abundance and
# a set level site fidelity (this is the 'actual' system state)
## note: Site fidelity is termed 'abandonment' and 'reoccupation' in this script which may lead
## to confusion as high site fidelity == low rate of abandonment
# secondly, a monitoring scheme is simulated with different levels of biased site selection
## (here termed ADSS - Abundance Dependent Site Selection)
# and different levels of observer retention (this is the 'observed' system state)
# thirdly, population trends for all datasets are modelled
# and finally, root mean square errors are calculated to measure the difference between
# actual and observed population trends


# run in parallel
library(foreach)
library(doParallel)

numCores <- detectCores()
cl <- makeCluster(numCores)
registerDoParallel(cl)

# Bat population: Inter-annual variation low, Abandonment low
var_low_ab_low <- foreach(i = 1:1000, .packages = c("tibble", "dplyr", "mgcv", "lme4")) %dopar% {
  srcf <- "N:/RProjects/Bat_Observer/01_Simulation/Functions"
  R.utils::sourceDirectory(srcf)
  
  z <- 400 # number of roosts
  x <- 72 # mean abundance
  y <- 1 # starting year
  
  N0 <- rpois(lam = x, n = z) # creates z number of roosts with mean abundance x
  
  bats_1 <- tibble(count = N0) %>%
    add_column(roost = 1:z, year = y, abyr = NA, reocc = NA)
  
  y <- y + 1
  
  pogr <- bats_1
  
  # loop for 20 years simulating population growth and abandonment/reoccupation for each year
  while (y <= 20) {
    pogr <- pogr %>%
      mutate(year = year + 1)
    
    pogr <- pogr %>%
      group_by(roost) %>%
      mutate(count = growth(Nt = count, se.r = 0.01)) %>%
      ungroup()
    
    pogr <- reoccupy(pogr, p2 = 0.7)
    pogr <- abandon(pogr, p = 0.02)
    
    assign(paste("bats", y, sep = "_"), pogr)
    
    y <- y + 1
    gc()
  }
  
  rm(pogr)
  
  all <- mget(ls(pattern = "bats_*")) %>%
    bind_rows() %>%
    as_tibble()
  
  rm(list = ls(pattern = "bats_*"))
  
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
  
  adss_low <- adss(act, frac = 0.2)
  adss_med <- adss(act, frac = 0.5)
  adss_high <- adss(act, frac = 0.8)
  
  adss_low_com_high <- commitment(adss_low, pr = 0.8)
  adss_med_com_high <- commitment(adss_med, pr = 0.8)
  adss_high_com_high <- commitment(adss_high, pr = 0.8)
  
  adss_low_com_med <- commitment(adss_low, pr = 0.5)
  adss_med_com_med <- commitment(adss_med, pr = 0.5)
  adss_high_com_med <- commitment(adss_high, pr = 0.5)
  
  adss_low_com_low <- commitment(adss_low, pr = 0.2)
  adss_med_com_low <- commitment(adss_med, pr = 0.2)
  adss_high_com_low <- commitment(adss_high, pr = 0.2)
  

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
  
  all <- allind
  
  
  return(all)
  gc()
}
save(var_low_ab_low, file = "N:/RProjects/Bat_Observer/02_Results/1000_runs/var-low_ab-low.RData")
gc()

var_low_ab_high <- foreach(i = 1:1000, .packages = c("tibble", "dplyr", "mgcv", "lme4")) %dopar% {
  srcf <- "N:/RProjects/Bat_Observer/01_Simulation/Functions"
  R.utils::sourceDirectory(srcf)
  
  z <- 400 # number of roosts
  x <- 72 # mean abundance
  y <- 1 # starting year
  
  N0 <- rpois(lam = x, n = z) # creates z number of roosts with mean abundance x
  
  bats_1 <- tibble(count = N0) %>%
    add_column(roost = 1:z, year = y, abyr = NA, reocc = NA)
  
  y <- y + 1
  
  pogr <- bats_1
  
  # loop for 20 years simulating population growth and abandonment/reoccupation for each year
  while (y <= 20) {
    pogr <- pogr %>%
      mutate(year = year + 1)
    
    pogr <- pogr %>%
      group_by(roost) %>%
      mutate(count = growth(Nt = count, se.r = 0.01)) %>%
      ungroup()
    
    pogr <- reoccupy(pogr, p2 = 0.7)
    pogr <- abandon(pogr, p = 0.1)
    
    assign(paste("bats", y, sep = "_"), pogr)
    
    y <- y + 1
    gc()
  }
  
  rm(pogr)
  
  all <- mget(ls(pattern = "bats_*")) %>%
    bind_rows() %>%
    as_tibble()
  
  rm(list = ls(pattern = "bats_*"))
  
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
  
  adss_low <- adss(act, frac = 0.2)
  adss_med <- adss(act, frac = 0.5)
  adss_high <- adss(act, frac = 0.8)
  
  adss_low_com_high <- commitment(adss_low, pr = 0.8)
  adss_med_com_high <- commitment(adss_med, pr = 0.8)
  adss_high_com_high <- commitment(adss_high, pr = 0.8)
  
  adss_low_com_med <- commitment(adss_low, pr = 0.5)
  adss_med_com_med <- commitment(adss_med, pr = 0.5)
  adss_high_com_med <- commitment(adss_high, pr = 0.5)
  
  adss_low_com_low <- commitment(adss_low, pr = 0.2)
  adss_med_com_low <- commitment(adss_med, pr = 0.2)
  adss_high_com_low <- commitment(adss_high, pr = 0.2)
  
  
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
  
  all <- allind
  
  
  return(all)
  gc()
}
save(var_low_ab_high, file = "N:/RProjects/Bat_Observer/02_Results/1000_runs/var-low_ab-high.RData")
gc()

var_high_ab_low <- foreach(i = 1:1000, .packages = c("tibble", "dplyr", "mgcv", "lme4")) %dopar% {
  srcf <- "N:/RProjects/Bat_Observer/01_Simulation/Functions"
  R.utils::sourceDirectory(srcf)
  
  z <- 400 # number of roosts
  x <- 72 # mean abundance
  y <- 1 # starting year
  
  N0 <- rpois(lam = x, n = z) # creates z number of roosts with mean abundance x
  
  bats_1 <- tibble(count = N0) %>%
    add_column(roost = 1:z, year = y, abyr = NA, reocc = NA)
  
  y <- y + 1
  
  pogr <- bats_1
  
  # loop for 20 years simulating population growth and abandonment/reoccupation for each year
  while (y <= 20) {
    pogr <- pogr %>%
      mutate(year = year + 1)
    
    pogr <- pogr %>%
      group_by(roost) %>%
      mutate(count = growth(Nt = count, se.r = 0.1)) %>%
      ungroup()
    
    pogr <- reoccupy(pogr, p2 = 0.7)
    pogr <- abandon(pogr, p = 0.02)
    
    assign(paste("bats", y, sep = "_"), pogr)
    
    y <- y + 1
    gc()
  }
  
  rm(pogr)
  
  all <- mget(ls(pattern = "bats_*")) %>%
    bind_rows() %>%
    as_tibble()
  
  rm(list = ls(pattern = "bats_*"))
  
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
  
  adss_low <- adss(act, frac = 0.2)
  adss_med <- adss(act, frac = 0.5)
  adss_high <- adss(act, frac = 0.8)
  
  adss_low_com_high <- commitment(adss_low, pr = 0.8)
  adss_med_com_high <- commitment(adss_med, pr = 0.8)
  adss_high_com_high <- commitment(adss_high, pr = 0.8)
  
  adss_low_com_med <- commitment(adss_low, pr = 0.5)
  adss_med_com_med <- commitment(adss_med, pr = 0.5)
  adss_high_com_med <- commitment(adss_high, pr = 0.5)
  
  adss_low_com_low <- commitment(adss_low, pr = 0.2)
  adss_med_com_low <- commitment(adss_med, pr = 0.2)
  adss_high_com_low <- commitment(adss_high, pr = 0.2)
  
  
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
  
  all <- allind
  
  
  return(all)
  gc()
}
save(var_high_ab_low, file = "N:/RProjects/Bat_Observer/02_Results/1000_runs/var-high_ab-low.RData")
gc()

var_high_ab_high <- foreach(i = 1:1000, .packages = c("tibble", "dplyr", "mgcv", "lme4")) %dopar% {
  srcf <- "N:/RProjects/Bat_Observer/01_Simulation/Functions"
  R.utils::sourceDirectory(srcf)
  
  z <- 400 # number of roosts
  x <- 72 # mean abundance
  y <- 1 # starting year
  
  N0 <- rpois(lam = x, n = z) # creates z number of roosts with mean abundance x
  
  bats_1 <- tibble(count = N0) %>%
    add_column(roost = 1:z, year = y, abyr = NA, reocc = NA)
  
  y <- y + 1
  
  pogr <- bats_1
  
  # loop for 20 years simulating population growth and abandonment/reoccupation for each year
  while (y <= 20) {
    pogr <- pogr %>%
      mutate(year = year + 1)
    
    pogr <- pogr %>%
      group_by(roost) %>%
      mutate(count = growth(Nt = count, se.r = 0.1)) %>%
      ungroup()
    
    pogr <- reoccupy(pogr, p2 = 0.7)
    pogr <- abandon(pogr, p = 0.1)
    
    assign(paste("bats", y, sep = "_"), pogr)
    
    y <- y + 1
    gc()
  }
  
  rm(pogr)
  
  all <- mget(ls(pattern = "bats_*")) %>%
    bind_rows() %>%
    as_tibble()
  
  rm(list = ls(pattern = "bats_*"))
  
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
  
  adss_low <- adss(act, frac = 0.2)
  adss_med <- adss(act, frac = 0.5)
  adss_high <- adss(act, frac = 0.8)
  
  adss_low_com_high <- commitment(adss_low, pr = 0.8)
  adss_med_com_high <- commitment(adss_med, pr = 0.8)
  adss_high_com_high <- commitment(adss_high, pr = 0.8)
  
  adss_low_com_med <- commitment(adss_low, pr = 0.5)
  adss_med_com_med <- commitment(adss_med, pr = 0.5)
  adss_high_com_med <- commitment(adss_high, pr = 0.5)
  
  adss_low_com_low <- commitment(adss_low, pr = 0.2)
  adss_med_com_low <- commitment(adss_med, pr = 0.2)
  adss_high_com_low <- commitment(adss_high, pr = 0.2)
  
  
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
  
  all <- allind
  
  
  return(all)
  gc()
}
save(var_high_ab_high, file = "N:/RProjects/Bat_Observer/02_Results/1000_runs/var-high_ab-high.RData")

stopCluster(cl)