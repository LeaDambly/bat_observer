library(foreach)
library(doParallel)
numCores <- detectCores()
cl <- makeCluster(numCores)
registerDoParallel(cl)

clusterCall(cl, function() { source("01_functions_180614.R") })

spl120 <- foreach(i = 1:100, .packages = c("tidyr", "dplyr", "mgcv", "reshape2")) %dopar%{
  z = 300 #number of monitored roosts
  x = 66 #roost size avg
  
  N0 <- rpois(lam = x, n = z)
  fis_1 <- data.frame(N0)
  fis_1$roost <- 1:z
  fis_1$year <- 1
  colnames(fis_1) <- c("count", "roost", "year")
  
  nyr = 2
  
  while(nyr <= 20){
    spl <- lapply(N0, growth) # populations for that year
    spl <- melt(spl)
    spl$year <- nyr
    colnames(spl) <- c("count", "roost", "year")
    
    fis <- split(spl, k = 120, n = 3, yr = nyr)
    
    N0 = fis$count
    assign(paste("fis", nyr, sep = "_"), fis)
    
    nyr = nyr+1
    gc()
  }
  
  rm(spl)
  rm(fis)
  
  all <- mget(ls(pattern="fis_*")) %>%
    bind_rows()
  
  act <- all[c(1:3)]
  obs <- act %>% filter(roost <= z)
  rm(all)
  rm(list = ls(pattern="fis_*"))
  
  obs <- accessability(obs)
  obs <- obs[c(1:3)]
  obs <- start_stop_monitoring(obs)
  
  gama <- gam_func(act)
  gamb <- gam_func(obs)
  
  inda <- index_func(act, gama)
  indb <- index_func(obs, gamb)
  
  return(list(act, obs, inda, indb))
  
  gc()
}

stopCluster(cl)

a <- sapply(spl120,function(x) x[1])
a <- do.call("rbind", a)

b <- sapply(spl120,function(x) x[2])
b <- do.call("rbind", b)

c <- sapply(spl120,function(x) x[3])
c <- do.call("rbind", c)

d <- sapply(spl120,function(x) x[4])
d <- do.call("rbind", d)

write.csv(a, file = "a.csv")
write.csv(b, file = "b.csv")
write.csv(c, file = "c.csv")
write.csv(d, file = "d.csv")

