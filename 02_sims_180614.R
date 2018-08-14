# simulations
# 14 06 18 - Lea I Dambly
# Libraries----
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(mgcv)
library(gridExtra)
source("01_functions_180614.R")


### SIM----
i <- 1
while(i <= 1000){

z = 100 #number of monitored roosts
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

  fis <- split(spl, k = 80, n = 3) # split pops

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
obs <- act %>% group_by(year) %>% filter(roost <= z) %>% ungroup()
rm(all)

obs <- accessability(obs)

obs <- start_stop_monitoring(obs)

gama <- gam_func(act, c(6)) #actual
gamb <- gam_func(obs, c(6)) #bias
inda <- index_func(act, gama, c(6))
indb <- index_func(obs, gamb, c(6))

assign(paste("act", i, sep = "_"), act)
assign(paste("obs", i, sep = "_"), obs)
assign(paste("inda", i, sep = "_"), inda)
assign(paste("indb", i, sep = "_"), indb)
rm(act)
rm(obs)
rm(gama)
rm(gamb)
rm(inda)
rm(indb)
gc()


i <- i+1
}

all_act <- mget(ls(pattern="inda_*")) %>%
  bind_rows()
all_obs <- mget(ls(pattern="indb_*")) %>%
  bind_rows()

all_act$reps <- c(0, rep(1:(nrow(all_act)-1) %/% 20))
all_act$reps <- as.factor(all_act$reps)
all_obs$reps <- c(0, rep(1:(nrow(all_obs)-1) %/% 20))
all_obs$reps <- as.factor(all_obs$reps)

pl1 <- ggplot(NULL) +
  geom_line(data = all_act, aes(years_ord, index_df, group = reps), colour = "red", alpha = 0.5) +
  geom_hline(yintercept = 1, size = 0.5)

pl2 <- ggplot(NULL) +
  geom_line(data = all_obs, aes(years_ord, index_df, group = reps), colour = "blue",  alpha = 0.5) +
  geom_hline(yintercept = 1, size = 0.5)

grid.arrange(pl1, pl2, ncol = 2)


### Roost switch
### High abundance: likelihood of selecting site depends on size of roost
### Wrong counts