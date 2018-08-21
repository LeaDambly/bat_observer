library(reshape2)
library(dplyr)
library(tidyr)
library(mgcv)
library(foreach)
library(doParallel)

source("01_functions_180614.R")

numCores <- detectCores()
cl <- makeCluster(numCores)
registerDoParallel(cl)

gro <- foreach(i = 1:1000) %dopar%{
  library(reshape2)
  library(dplyr)
  library(tidyr)
  library(mgcv)
  source("01_functions_180614.R")
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
    
    #fis <- split(spl, k = 110, n = 3) # split pops
    
    N0 = spl$count
    #N0 = fis$count
    
    
    #assign(paste("fis", nyr, sep = "_"), fis)
    assign(paste("fis", nyr, sep = "_"), spl)
    
    nyr = nyr+1
    gc()
  }
  
  rm(spl)
  
  all <- mget(ls(pattern="fis_*")) %>%
    bind_rows()
  
  act <- all[c(1:3)]
  obs <- act %>% filter(roost <= z)
  rm(all)
  rm(list = ls(pattern="fis_*"))
  
  obs <- accessability(obs)
  obs <- obs[c(1:3)]
  obs2 <- start_stop_monitoring(obs)
  
  gama <- gam_func(act)
  gamb <- gam_func(obs)
  gamc <- gam_func(obs2)
  
  inda <- index_func(act, gama)
  indb <- index_func(obs, gamb)
  indc <- index_func(obs2, gamc)
  
  return(list(act, obs, obs2, gama, gamb, gamc, inda, indb, indc))
  
  gc()
}

a <- sapply(gro,function(x) x[1])
a <- do.call("rbind", a)

b <- sapply(gro,function(x) x[2])
b <- do.call("rbind", b)

c <- sapply(gro,function(x) x[3])
c <- do.call("rbind", c)

d <- sapply(gro,function(x) x[4])
d <- do.call("rbind", d)

e <- sapply(gro,function(x) x[5])
e <- do.call("rbind", e)

f <- sapply(gro,function(x) x[6])
f <- do.call("rbind", f)

g <- sapply(gro,function(x) x[7])
g <- do.call("rbind", g)

h <- sapply(gro,function(x) x[8])
h <- do.call("rbind", h)

i <- sapply(gro,function(x) x[9])
i <- do.call("rbind", i)


write.csv(a, file="a.csv") 
write.csv(b, file="b.csv") 
write.csv(c, file="c.csv") 
write.csv(d, file="d.csv") 
write.csv(e, file="e.csv") 
write.csv(f, file="f.csv") 
write.csv(g, file="g.csv") 
write.csv(h, file="h.csv")
write.csv(i, file="i.csv")