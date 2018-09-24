# This code runs bat roost simulations for debugging

# libraries----
library(tidyr)
library(dplyr)
library(mgcv)
library(reshape2)

growth <- function(Nt = N0, mu.r = 0, se.r = 0.0) {
  lambda <- exp(rnorm(n = 1, mean = mu.r, sd = se.r)) #growth rate lambda
  Ntplus1 <- round(Nt * lambda, 0) #growth for that year
  return(Ntplus1)
}

switch1 <- function(pogr){
  perc <-  runif(1, 0.01, 0.1)
  
  w <- pogr # the population for that year per site
  colnames(w) <- c("count", "roost", "year")
  
  w2 <- w %>% sample_frac(perc) #roosts that switch
  
  w3 <- suppressMessages(anti_join(w, w2)) #remove switched roosts from original df
  
  colnames(w2) <- c("count", "orgroost", "year")
  w4 <- w2 #copy switched roosts
  w4$count = 0  #set count to zero (this is the count of their origin roost)
  
  w3 <- w3 %>% bind_rows(w4) #add 0 counts to orig
  
  w3 <- w3 %>% mutate(r = ifelse(is.na(orgroost), roost, orgroost)) %>% arrange(r)

  w3 <- w3[c(1,3)]
  
  w <-  w3 %>% 
    
    bind_rows(w2) %>% #top newly switched roosts on top of orig roosts
    
    tibble::rownames_to_column(., var = "roost")
  w$roost <- as.numeric(w$roost)
  
  return(w)
}

switch2 <- function(ltsw){
  perc <-  runif(1, 0.8, 0.9)
  
  w <- ltsw
  
  w2 <- w %>% filter(!is.na(orgroost)) %>% sample_frac(perc)
  
  w <- suppressMessages(anti_join(w, w2))
  
  w2$roost <- w2$orgroost
  
  w2$orgroost <- NA
  
  w$count2 <- w2[match(w$roost, w2$roost),2]
  
  w <- w %>% mutate(a = ifelse(is.na(count2), count, count2))
  
  w <- w[c(1,3,4,6)]
  colnames(w) <- c("roost", "year", "orgroost", "count")
  
  return(w)
  
  }


# sim----
  z = 300 #number of monitored roosts
  x = 66 #roost size avg
  
  # first year populations
  N0 <- rpois(lam = x, n = z)
  ltsw_1 <- data.frame(N0)
  ltsw_1$roost <- 1:z
  ltsw_1$year <- 1
  colnames(ltsw_1) <- c("count", "roost", "year")
  
  nyr = 2
  
  while(nyr <= 20){
    pogr <- lapply(N0, growth) # population growth for that year
    pogr <- melt(pogr)
    pogr$year <- nyr
    
    ltsw <- switch1(pogr)
    ltsw <- switch2(ltsw)
    
    N0 = ltsw$count
    assign(paste("ltsw", nyr, sep = "_"), ltsw)
    
    nyr = nyr+1
    gc()
  }
  
  rm(ltsw)
  rm(pogr)
  
  all <- mget(ls(pattern="ltsw_*")) %>%
    bind_rows()
  
  act <- all[c(1:3)]
  act2 <- acast(act, roost~year, value.var = "count", fill = 0)
  act2 <- melt(act2, value.name = "count")
  colnames(act2) <- c("roost", "year", "count")
  
  org <- act %>% filter(roost <= z)
  rm(all)
  rm(list = ls(pattern="ltsw_*"))
  
  obs <- accessability(org)
  obs <- obs[c(1:3)]
  obs <- start_stop_monitoring(obs)
  
  
  gama <- gam_func(act)
  gamb <- gam_func(act2)
  gamc <- gam_func(org)
  gamd <- gam_func(obs)
  
  inda <- index_func(act, gama)
  indb <- index_func(act2, gamb)
  indc <- index_func(org, gamc)
  indd <- index_func(obs, gamd)
  
  return(list(inda, indb, indc, indd))
  
  gc()

library(ggplot2)
ggplot(data = inda, aes(x = years_ord, y = index_df)) +
  geom_line()
ggplot(data = indb, aes(x = years_ord, y = index_df)) +
  geom_line()
ggplot(data = indc, aes(x = years_ord, y = index_df)) +
  geom_line()
