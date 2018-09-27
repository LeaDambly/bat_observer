# This code runs bat roost simulations in parallel

# libraries----
library(foreach)
library(doParallel)

# make parallel----
numCores <- detectCores()-2
cl <- makeCluster(numCores)
registerDoParallel(cl)

# run sim----
swicom <- foreach(i = 1:2, .packages = c("tibble", "tidyr", "dplyr", "mgcv", "reshape2")) %do%{
  z = 300 #number of monitored roosts
  x = 66 #roost size avg
  
  # first year populations
  N0 <- rpois(lam = x, n = z)
  ltsw_1 <- tibble(count = N0) %>% add_column(., roost = 1:z) %>% add_column(year = 1)
 
  nyr = 2
  
  while(nyr <= 20){
    pogr <- lapply(N0, growth) # population growth for that year
    pogr <- melt(pogr) %>% rename(., count = value, roost = L1) %>% add_column(year = nyr)
    
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
  
  act <- all[c("roost", "year", "count")]
  act2 <- acast(act, roost~year, value.var = "count", fill = 0)
  act2 <- melt(act2, value.name = "count") %>% rename(roost = Var1, year = Var2)
  
  org <- act %>% filter(roost <= z)
  
  rm(all)
  rm(list = ls(pattern="ltsw_*"))
  
  obs <- accessability(org)
  obs <- obs[c("roost", "year", "count")]
  obs <- start_stop_monitoring(obs)
  
  gama <- gam_func(act)
  gamb <- gam_func(act2)
  gamc <- gam_func(org)
  gamd <- gam_func(obs)
  
  rmsea <- rmse_func(gama)
  rmseb <- rmse_func(gamb)
  rmsec <- rmse_func(gamc)
  rmsed <- rmse_func(gamd)
  
  inda <- index_func(act, gama)
  indb <- index_func(act2, gamb)
  indc <- index_func(org, gamc)
  indd <- index_func(obs, gamd)
  
  return(list(inda, indb, indc, indd, rmsea, rmseb, rmsec, rmsed))
  
  gc()
}

stopCluster(cl)

a <- sapply(swicom,function(x) x[1])
a <- do.call("rbind", a)

b <- sapply(swicom,function(x) x[2])
b <- do.call("rbind", b)

c <- sapply(swicom,function(x) x[3])
c <- do.call("rbind", c)

d <- sapply(swicom,function(x) x[4])
d <- do.call("rbind", d)

e <- sapply(swicom,function(x) x[5])
e <- do.call("rbind", e)

f <- sapply(swicom,function(x) x[6])
f <- do.call("rbind", f)

g <- sapply(swicom,function(x) x[7])
g <- do.call("rbind", g)

h <- sapply(swicom,function(x) x[8])
h <- do.call("rbind", h)

write.csv(a, file = "a.csv")
write.csv(b, file = "b.csv")
write.csv(c, file = "c.csv")
write.csv(d, file = "d.csv")
write.csv(e, file = "e.csv")
write.csv(f, file = "f.csv")
write.csv(g, file = "g.csv")
write.csv(h, file = "h.csv")