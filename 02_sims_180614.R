# simulations
# 14 06 18 - Lea I Dambly
# Libraries----
library(ggplot2)
library(reshape2)
library(dplyr)
library(mgcv)
library(gridExtra)
library(grid)
library(rmarkdown)

source("01_functions_180614.R")


### Pop growth AND roost split----
z = 600 #number of roosts
x = 66 #roost size avg

N0 <- rpois(lam = x, n = z)
fis_1 <- data.frame(N0)
fis_1$roost <- 1:z
fis_1$year <- 1
colnames(fis_1) <- c("count", "roost", "year")

nyr = 2

while(nyr <= 20){
  obs <- lapply(N0, growth)
  fis <- lapply(obs, split)
  
  fis <- melt(fis)
  fis$year <- nyr
  colnames(fis) <- c("count", "L1", "year")
  
  obs <- melt(obs)
  obs$year <- nyr
  colnames(obs) <- c("count", "L1", "year")
  obs <- obs%>%
    distinct(.) %>%
    arrange(L1)
  obs <- obs[1:z,]
  
  dup <- fis %>% 
    group_by(L1) %>% 
    filter(n()>1) %>%
    distinct(.)
    
  fis <- bind_rows(obs, dup, id = NULL)
  
  fis <- fis %>% 
   mutate(roost = 1:nrow(fis))
  
  N0 = fis$count
  
  assign(paste("fis", nyr, sep = "_"), fis)
  
  nyr = nyr+1
}
rm(obs)
rm(fis)
rm(dup)
rm(all)

all <- mget(ls(pattern="fis_*")) %>%
  bind_rows()

#100----
all_100_70 <- all[c(1:3)]
obs_100_70 <- all_100_70 %>% group_by(year) %>% filter(roost <= z) %>% ungroup()

all_100_90 <- all[c(1:3)]
obs_100_90 <- all_100_90 %>% group_by(year) %>% filter(roost <= z) %>% ungroup()

all_100_110 <- all[c(1:3)]
obs_100_110 <- all_100_110 %>% group_by(year) %>% filter(roost <= z) %>% ungroup()

#500----
all_500_70 <- all[c(1:3)]
obs_500_70 <- all_500_70 %>% group_by(year) %>% filter(roost <= z) %>% ungroup()

all_500_90 <- all[c(1:3)]
obs_500_90 <- all_500_90 %>% group_by(year) %>% filter(roost <= z) %>% ungroup()

all_500_110 <- all[c(1:3)]
obs_500_110 <- all_500_110 %>% group_by(year) %>% filter(roost <= z) %>% ungroup()

#900----
all_900_70 <- all[c(1:3)]
obs_900_70 <- all_900_70 %>% group_by(year) %>% filter(roost <= z) %>% ungroup()

all_900_90 <- all[c(1:3)]
obs_900_90 <- all_900_90 %>% group_by(year) %>% filter(roost <= z) %>% ungroup()

all_900_110 <- all[c(1:3)]
obs_900_110 <- all_900_110 %>% group_by(year) %>% filter(roost <= z) %>% ungroup()

# fit gam----
gam_act_100_70 <- gam_func(all_100_70, c(6))
gam_obs_100_70 <- gam_func(obs_100_70, c(6))

gam_act_100_90 <- gam_func(all_100_90, c(6))
gam_obs_100_90 <- gam_func(obs_100_90, c(6))

gam_act_100_110 <- gam_func(all_100_110, c(6))
gam_obs_100_110 <- gam_func(obs_100_110, c(6))

gam_act_500_70 <- gam_func(all_500_70, c(6))
gam_obs_500_70 <- gam_func(obs_500_70, c(6))

gam_act_500_90 <- gam_func(all_500_90, c(6))
gam_obs_500_90 <- gam_func(obs_500_90, c(6))

gam_act_500_110 <- gam_func(all_500_110, c(6))
gam_obs_500_110 <- gam_func(obs_500_110, c(6))

gam_act_900_70 <- gam_func(all_900_70, c(6))
gam_obs_900_70 <- gam_func(obs_900_70, c(6))

gam_act_900_90 <- gam_func(all_900_90, c(6))
gam_obs_900_90 <- gam_func(obs_900_90, c(6))

gam_act_900_110 <- gam_func(all_900_110, c(6))
gam_obs_900_110 <- gam_func(obs_900_110, c(6))

# produce pop index----
ind_act_100_70 <- index_func(all_100_70, gam_act_100_70, c(6))
ind_obs_100_70 <- index_func(obs_100_70, gam_obs_100_70, c(6))

ind_act_100_90 <- index_func(all_100_90, gam_act_100_90, c(6))
ind_obs_100_90 <- index_func(obs_100_90, gam_obs_100_90, c(6))

ind_act_100_110 <- index_func(all_100_110, gam_act_100_110, c(6))
ind_obs_100_110 <- index_func(obs_100_110, gam_obs_100_110, c(6))

ind_act_500_70 <- index_func(all_500_70, gam_act_500_70, c(6))
ind_obs_500_70 <- index_func(obs_500_70, gam_obs_500_70, c(6))

ind_act_500_90 <- index_func(all_500_90, gam_act_500_90, c(6))
ind_obs_500_90 <- index_func(obs_500_90, gam_obs_500_90, c(6))

ind_act_500_110 <- index_func(all_500_110, gam_act_500_110, c(6))
ind_obs_500_110 <- index_func(obs_500_110, gam_obs_500_110, c(6))

ind_act_900_70 <- index_func(all_900_70, gam_act_900_70, c(6))
ind_obs_900_70 <- index_func(obs_900_70, gam_obs_900_70, c(6))

ind_act_900_90 <- index_func(all_900_90, gam_act_900_90, c(6))
ind_obs_900_90 <- index_func(obs_900_90, gam_obs_900_90, c(6))

ind_act_900_110 <- index_func(all_900_110, gam_act_900_110, c(6))
ind_obs_900_110 <- index_func(obs_900_110, gam_obs_900_110, c(6))

# bootstrap (should be higher when time)----
boot_act_100_70 <- outer_boot_func(all_100_70, 6, 39)
boot_obs_100_70 <- outer_boot_func(obs_100_70, 6, 39)

boot_act_100_90 <- outer_boot_func(all_100_90, 6, 39)
boot_obs_100_90 <- outer_boot_func(obs_100_90, 6, 39)

boot_act_100_110 <- outer_boot_func(all_100_110, 6, 39)
boot_obs_100_110 <- outer_boot_func(obs_100_110, 6, 39)

boot_act_500_70 <- outer_boot_func(all_500_70, 6, 39)
boot_obs_500_70 <- outer_boot_func(obs_500_70, 6, 39)

boot_act_500_90 <- outer_boot_func(all_500_90, 6, 39)
boot_obs_500_90 <- outer_boot_func(obs_500_90, 6, 39)

boot_act_500_110 <- outer_boot_func(all_500_110, 6, 39)
boot_obs_500_110 <- outer_boot_func(obs_500_110, 6, 39)

boot_act_900_70 <- outer_boot_func(all_900_70, 6, 39)
boot_obs_900_70 <- outer_boot_func(obs_900_70, 6, 39)

boot_act_900_90 <- outer_boot_func(all_900_90, 6, 39)
boot_obs_900_90 <- outer_boot_func(obs_900_90, 6, 39)

boot_act_900_110 <- outer_boot_func(all_900_110, 6, 39)
boot_obs_900_110 <- outer_boot_func(obs_900_110, 6, 39)

# confidence intervals----
ci_act_100_70 <- ci_func(boot_act_100_70, conf = 0.95)
ci_obs_100_70 <- ci_func(boot_obs_100_70, conf = 0.95)

ci_act_100_90 <- ci_func(boot_act_100_90, conf = 0.95)
ci_obs_100_90 <- ci_func(boot_obs_100_90, conf = 0.95)

ci_act_100_110 <- ci_func(boot_act_100_110, conf = 0.95)
ci_obs_100_110 <- ci_func(boot_obs_100_110, conf = 0.95)

ci_act_500_70 <- ci_func(boot_act_500_70, conf = 0.95)
ci_obs_500_70 <- ci_func(boot_obs_500_70, conf = 0.95)

ci_act_500_90 <- ci_func(boot_act_500_90, conf = 0.95)
ci_obs_500_90 <- ci_func(boot_obs_500_90, conf = 0.95)

ci_act_500_110 <- ci_func(boot_act_500_110, conf = 0.95)
ci_obs_500_110 <- ci_func(boot_obs_500_110, conf = 0.95)

ci_act_900_70 <- ci_func(boot_act_900_70, conf = 0.95)
ci_obs_900_70 <- ci_func(boot_obs_900_70, conf = 0.95)

ci_act_900_90 <- ci_func(boot_act_900_90, conf = 0.95)
ci_obs_900_90 <- ci_func(boot_obs_900_90, conf = 0.95)

ci_act_900_110 <- ci_func(boot_act_900_110, conf = 0.95)
ci_obs_900_110 <- ci_func(boot_obs_900_110, conf = 0.95)


# RMSE----
rmse_act_100_70 <- rmse_func(gam_act_100_70)
rmse_obs_100_70 <- rmse_func(gam_obs_100_70)

rmse_act_100_90 <- rmse_func(gam_act_100_90)
rmse_obs_100_90 <- rmse_func(gam_obs_100_90)

rmse_act_100_110 <- rmse_func(gam_act_100_110)
rmse_obs_100_110 <- rmse_func(gam_obs_100_110)

rmse_act_500_70 <- rmse_func(gam_act_500_70)
rmse_obs_500_70 <- rmse_func(gam_obs_500_70)

rmse_act_500_90 <- rmse_func(gam_act_500_90)
rmse_obs_500_90 <- rmse_func(gam_obs_500_90)

rmse_act_500_110 <- rmse_func(gam_act_500_110)
rmse_obs_500_110 <- rmse_func(gam_obs_500_110)

rmse_act_900_70 <- rmse_func(gam_act_900_70)
rmse_obs_900_70 <- rmse_func(gam_obs_900_70)

rmse_act_900_90 <- rmse_func(gam_act_900_90)
rmse_obs_900_90 <- rmse_func(gam_obs_900_90)

rmse_act_900_110 <- rmse_func(gam_act_900_110)
rmse_obs_900_110 <- rmse_func(gam_obs_900_110)

# combine----
act_100 <- data.frame(year = ind_act_100_70[,1],
                     pop_index70 = ind_act_100_70[,2]*100,
                     pop_index90 = ind_act_100_90[,2]*100,
                     pop_index110 = ind_act_100_110[,2]*100)
                     
obs_100 <- data.frame(year = ind_obs_100_70[,1],
                      pop_index70 = ind_obs_100_70[,2]*100,
                      pop_index90 = ind_obs_100_90[,2]*100,
                      pop_index110 = ind_obs_100_110[,2]*100)

act_500 <- data.frame(year = ind_act_500_70[,1],
                      pop_index70 = ind_act_500_70[,2]*100,
                      pop_index90 = ind_act_500_90[,2]*100,
                      pop_index110 = ind_act_500_110[,2]*100)

obs_500 <- data.frame(year = ind_obs_500_70[,1],
                      pop_index70 = ind_obs_500_70[,2]*100,
                      pop_index90 = ind_obs_500_90[,2]*100,
                      pop_index110 = ind_obs_500_110[,2]*100)

act_900 <- data.frame(year = ind_act_900_70[,1],
                      pop_index70 = ind_act_900_70[,2]*100,
                      pop_index90 = ind_act_900_90[,2]*100,
                      pop_index110 = ind_act_900_110[,2]*100)

obs_900 <- data.frame(year = ind_obs_900_70[,1],
                      pop_index70 = ind_obs_900_70[,2]*100,
                      pop_index90 = ind_obs_900_90[,2]*100,
                      pop_index110 = ind_obs_900_110[,2]*100)


# plotting----
plot_act_100 <- ggplot(data = act_100) +
  geom_line(aes(year, pop_index70, color = "70")) +
  geom_line(aes(year, pop_index90, color = "90")) +
  geom_line(aes(year, pop_index110, color = "110")) +
  scale_x_continuous(expand = c(0, 0), limits = c(1, 20.5), breaks = seq(0, 20, 2)) +
  scale_y_continuous(expand = c(0, 0), limits = c(95, 115)) +
  ylab(" ") +
  geom_hline(yintercept = 100, colour = "#9FA3A7") +
  theme_bw() + 
  theme(axis.text.x = element_text(face = "bold", color = "black", size = 10),
        axis.text.y = element_text(face = "bold", color = "black", size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 13, margin = margin(r = 15)),
        plot.title = element_text(color = "black", hjust = 0.5),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none") +
  scale_colour_manual("Roost split threshold", 
                      labels = c("70", "90", "110"),
                      values = c("#1f78b4", "#33a02c", "#d01c8b"), 
                      breaks = c("70", "90", "110"))
plot_act_100

plot_obs_100 <- ggplot(data = obs_100) +
  geom_line(aes(year, pop_index70, color = "70")) +
  geom_line(aes(year, pop_index90, color = "90")) +
  geom_line(aes(year, pop_index110, color = "110")) +
  scale_x_continuous(expand = c(0, 0), limits = c(1, 20.5), breaks = seq(0, 20, 2)) +
  scale_y_continuous(expand = c(0, 0), limits = c(95, 115)) +
  geom_hline(yintercept = 100, colour = "#9FA3A7") +
  theme_bw() + 
  theme(axis.text.x = element_text(face = "bold", color = "black", size = 10),
        axis.text.y = element_text(face = "bold", color = "black", size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(color = "black", hjust = 0.5),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none") +
  scale_colour_manual("Roost split threshold", 
                      labels = c("70", "90", "110"),
                      values = c("#1f78b4", "#33a02c", "#d01c8b"), 
                      breaks = c("70", "90", "110"))
plot_obs_100

plot_act_500 <- ggplot(data = act_500) +
  geom_line(aes(year, pop_index70, color = "70")) +
  geom_line(aes(year, pop_index90, color = "90")) +
  geom_line(aes(year, pop_index110, color = "110")) +
  ylab(" ") +
  scale_x_continuous(expand = c(0, 0), limits = c(1, 20.5), breaks = seq(0, 20, 2)) +
  scale_y_continuous(expand = c(0, 0), limits = c(95, 115)) +
  geom_hline(yintercept = 100, colour = "#9FA3A7") +
  theme_bw() + 
  theme(axis.text.x = element_text(face = "bold", color = "black", size = 10),
        axis.text.y = element_text(face = "bold", color = "black", size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 13, margin = margin(r = 15)),
        plot.title = element_text(color = "black", hjust = 0.5),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none") +
  scale_colour_manual("Roost split threshold", 
                      labels = c("70", "90", "110"),
                      values = c("#1f78b4", "#33a02c", "#d01c8b"), 
                      breaks = c("70", "90", "110"))
plot_act_500

plot_obs_500 <- ggplot(data = obs_500) +
  geom_line(aes(year, pop_index70, color = "70")) +
  geom_line(aes(year, pop_index90, color = "90")) +
  geom_line(aes(year, pop_index110, color = "110")) +
  scale_x_continuous(expand = c(0, 0), limits = c(1, 20.5), breaks = seq(0, 20, 2)) +
  scale_y_continuous(expand = c(0, 0), limits = c(95, 115)) +
  geom_hline(yintercept = 100, colour = "#9FA3A7") +
  theme_bw() + 
  theme(axis.text.x = element_text(face = "bold", color = "black", size = 10),
        axis.text.y = element_text(face = "bold", color = "black", size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(color = "black", hjust = 0.5),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none") +
  scale_colour_manual("Roost split threshold", 
                      labels = c("70", "90", "110"),
                      values = c("#1f78b4", "#33a02c", "#d01c8b"), 
                      breaks = c("70", "90", "110"))
plot_obs_500

plot_act_900 <- ggplot(data = act_900) +
  geom_line(aes(year, pop_index70, color = "70")) +
  geom_line(aes(year, pop_index90, color = "90")) +
  geom_line(aes(year, pop_index110, color = "110")) +
  ylab(" ") +
  scale_x_continuous(expand = c(0, 0), limits = c(1, 20.5), breaks = seq(0, 20, 2)) +
  scale_y_continuous(expand = c(0, 0), limits = c(95, 115)) +
  geom_hline(yintercept = 100, colour = "#9FA3A7") +
  theme_bw() + 
  theme(axis.text.x = element_text(face = "bold", color = "black", size = 10),
        axis.text.y = element_text(face = "bold", color = "black", size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 13, margin = margin(r = 15)),
        plot.title = element_text(color = "black", hjust = 0.5),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none") +
  scale_colour_manual("Roost split threshold", 
                      labels = c("70", "90", "110"),
                      values = c("#1f78b4", "#33a02c", "#d01c8b"), 
                      breaks = c("70", "90", "110"))
plot_act_900

plot_obs_900 <- ggplot(data = obs_900) +
  geom_line(aes(year, pop_index70, color = "70")) +
  geom_line(aes(year, pop_index90, color = "90")) +
  geom_line(aes(year, pop_index110, color = "110")) +
  scale_x_continuous(expand = c(0, 0), limits = c(1, 20.5), breaks = seq(0, 20, 2)) +
  scale_y_continuous(expand = c(0, 0), limits = c(95, 115)) +
  geom_hline(yintercept = 100, colour = "#9FA3A7") +
  theme_bw() + 
  theme(axis.text.x = element_text(face = "bold", color = "black", size = 10),
        axis.text.y = element_text(face = "bold", color = "black", size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(color = "black", hjust = 0.5),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none") +
  scale_colour_manual("Roost split threshold", 
                      labels = c("70", "90", "110"),
                      values = c("#1f78b4", "#33a02c", "#d01c8b"), 
                      breaks = c("70", "90", "110"))
plot_obs_900



# all together plot----
all <- grid.arrange(arrangeGrob(plot_act_100,plot_act_500,plot_act_900,
            ncol=1), 
            arrangeGrob(plot_obs_100,plot_obs_500,plot_obs_900,
            ncol=1),   
            ncol = 2)
all
