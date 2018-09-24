# This code runs bat roost simulations in parallel

# libraries----
library(foreach)
library(doParallel)

# make parallel----
numCores <- detectCores()
cl <- makeCluster(numCores)
registerDoParallel(cl)

clusterCall(cl, function() { source("01_functions_180614.R") })

# run sim----
splcom <- foreach(i = 1:100, .packages = c("tidyr", "dplyr", "mgcv", "reshape2")) %do%{
  z = 300 #number of monitored roosts
  x = 66 #roost size avg
  
  # first year populations
  N0 <- rpois(lam = x, n = z)
  pogr_1 <- data.frame(N0)
  pogr_1$roost <- 1:z
  pogr_1$year <- 1
  colnames(pogr_1) <- c("count", "roost", "year")

  nyr = 2
  
  while(nyr <= 20){
    pogr <- lapply(N0, growth) # population growth for that year
    pogr <- melt(spl)
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
  
  act$rep <- i
  act2$rep <- i
  org$rep <- i
  obs$rep <- i
  
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
}

stopCluster(cl)

a <- sapply(splcom,function(x) x[1])
a <- do.call("rbind", a)

b <- sapply(splcom,function(x) x[2])
b <- do.call("rbind", b)

c <- sapply(splcom,function(x) x[3])
c <- do.call("rbind", c)

d <- sapply(splcom,function(x) x[4])
d <- do.call("rbind", d)

write.csv(a, file = "a.csv")
write.csv(b, file = "b.csv")
write.csv(c, file = "c.csv")
write.csv(d, file = "d.csv")
