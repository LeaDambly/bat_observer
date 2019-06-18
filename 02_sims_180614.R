# This code runs bat roost simulations in parallel

# libraries----
library(foreach)
library(doParallel)

# make parallel----
numCores <- detectCores() - 1
cl <- makeCluster(numCores)
registerDoParallel(cl)

# run sim----
swicom <-
  foreach(
    i = 1:1,
    .packages = c("tibble", "tidyr", "dplyr", "mgcv", "reshape2")
  ) %do% {
    z = 300 #number of monitored roosts
    x = 66 #roost size avg
    
    # first year populations
    N0 <- rpois(lam = x, n = z)
    ltsw_1 <- tibble(count = N0) %>%
      add_column(., roost = 1:z) %>%
      add_column(year = 1) %>%
      add_column(swyr = 0)
    
    pogr <- ltsw_1 %>%
      mutate(year = year + 1)
    nyr = 2
    
    while (nyr <= 20) {
      pogr <- pogr %>%
        group_by(1:n()) %>% #apply function row by row
        mutate(count = growth(Nt = count, mu.r = 0.79, se.r = 0.099)) %>%
        ungroup() %>%
        select(-"1:n()")
      
      #ltsw <- switch1(pogr)
      #ltsw <- switch2(ltsw)
      
      assign(paste("ltsw", nyr, sep = "_"), pogr)
      
      pogr <- ltsw %>%
        mutate(year = year + 1)
      nyr <- nyr + 1
      
      gc()
    }
    
    #rm(ltsw)
    rm(pogr)
    
    all <- mget(ls(pattern = "ltsw_*")) %>%
      bind_rows()
    rm(list = ls(pattern = "ltsw_*"))
    
    act <- acast(all, roost ~ year, value.var = "count", fill = 0)
    act <- melt(act, value.name = "count") %>%
      rename(roost = Var1, year = Var2)
    
    
    #mot <- motivation(act)
    acc <- accessibility(act, c = 0.8)
    #mot_acc <- accessibility(act, c = 0.8)
    #obs <- obs[c("roost", "year", "count")]
    #obs <- start_stop_monitoring(obs)

    gama <- gam_func(act)
    gamb <- gam_func(mot)
    gamc <- gam_func(acc)
    gamd <- gam_func(mot_acc)
    
    rmsea <- rmse_func(gama)
    rmseb <- rmse_func(gamb)
    
    inda <- index_func(act, gama)
    indb <- index_func(mot, gamb)
    indc <- index_func(acc, gamc)
    indd <- index_func(mot_acc, gamd)
    
    ggplot(NULL) +
      geom_line(data = inda,
                aes(years_ord, index_df),
                colour = "red") +
      geom_line(data = indb,
                aes(years_ord, index_df),
                colour = "blue") +
      geom_line(data = indc,
                aes(years_ord, index_df),
                colour = "green") +
      geom_line(data = indd,
                aes(years_ord, index_df),
                colour = "black")
    
    return(list(inda, indb, indc, indd, rmsea, rmseb, rmsec, rmsed))
    
    gc()
  }

stopCluster(cl)

a <- sapply(swicom, function(x)
  x[1])
a <- do.call("rbind", a)

b <- sapply(swicom, function(x)
  x[2])
b <- do.call("rbind", b)

c <- sapply(swicom, function(x)
  x[3])
c <- do.call("rbind", c)

d <- sapply(swicom, function(x)
  x[4])
d <- do.call("rbind", d)

e <- sapply(swicom, function(x)
  x[5])
e <- do.call("rbind", e)

f <- sapply(swicom, function(x)
  x[6])
f <- do.call("rbind", f)

g <- sapply(swicom, function(x)
  x[7])
g <- do.call("rbind", g)

h <- sapply(swicom, function(x)
  x[8])
h <- do.call("rbind", h)

write.csv(a, file = "a.csv")
write.csv(b, file = "b.csv")
write.csv(c, file = "c.csv")
write.csv(d, file = "d.csv")
write.csv(e, file = "e.csv")
write.csv(f, file = "f.csv")
write.csv(g, file = "g.csv")
write.csv(h, file = "h.csv")