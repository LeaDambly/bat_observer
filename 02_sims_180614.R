# simulations
# 14 06 18 - Lea I Dambly

# Libraries----
library(ggplot2)
library(reshape2)
library(dplyr)
library(mgcv)

source("01_functions_180614.R")

# geometric population growth model----
z = 600 #number of roosts
i = 0 #iterations
# each time a roost splits, a new roost is created. As long as there are roosts, the loop is running
# tdm_1 are the original 600 roosts, therefore the ones that are being observed
while(z != 0){
  i <- i+1
  N0 <- rpois(lam = 66, n = z)
  z <- 0
  tdm <- melt(sapply(N0, sim, nyr = 20))
  assign(paste("tdm", i, sep = "_"), tdm)
  rm(tdm)
}

tdm_all <- as.list(mget(paste("tdm_", 1:i, sep="")))
tdm_all <- bind_rows(tdm_all)

names(tdm_1) <- c('year', 'site', 'count')
names(tdm_all) <- c('year', 'site', 'count')
qplot(data = tdm_1, x = year, y = count, group = site, geom = 'line')
qplot(data = tdm_all, x = year, y = count, group = site, geom = 'line')



# fit gam----
gam1 <- gam_func(tdm_1, c(6))
gam2 <- gam_func(tdm_all, c(6))


# produce pop index----
ind1 <- index_func(tdm_1, gam1, c(6))
ind2 <- index_func(tdm_all, gam2, c(6))

plot(ind1, type='o', pch=20, cex=1, col='black')
plot(ind2, type='o', pch=20, cex=1, col='black')


# bootstrap (should be 399 when time)----
boot1 <- outer_boot_func(tdm, 6, 39)
boot2 <- outer_boot_func(mon, 6, 39)


# confidence intervals----
ci1 <- ci_func(boot1, conf = 0.95)
ci2 <- ci_func(boot2, conf = 0.95)


#combine----
actual <- data.frame(year = indmon1[,1],
                     pop_index = indmon1[,2]*100,
                     upper_ci = ci1[2,]*100, 
                     lower_ci = ci1[1,]*100)

observed <- data.frame(year = indmon2[,1],
                       pop_index = indmon2[,2]*100,
                       upper_ci = ci2[2,]*100, 
                       lower_ci = ci2[1,]*100)

#plotting----
plot.actual <- ggplot(data = actual) +
  geom_line(aes(year, pop_index)) +
  geom_line(aes(year, lower_ci), linetype = "dashed") +
  geom_line(aes(year, upper_ci), linetype = "dashed") +
  xlab("Years") +
  ylab("Index of Abundance (base year = 3)") +
  ggtitle("Actual state of the system") +
  scale_x_continuous(expand = c(0, 0), limits = c(1, 20.5), breaks = seq(0, 20, 2)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 200)) +
  geom_hline(yintercept = 100, colour = "#9FA3A7") +
  theme_bw() + 
  theme(axis.text.x = element_text(face = "bold", color = "black", size = 10),
        axis.text.y = element_text(face = "bold", color = "black", size = 10),
        axis.title.x = element_text(color = "black", size = 13, margin = margin(t = 15)),
        axis.title.y = element_text(color = "black", size = 13, margin = margin(r = 15)),
        plot.title = element_text(color = "black", hjust = 0.5),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank())
plot.actual

plot.observed <- ggplot(data = observed) +
  geom_line(aes(year, pop_index)) +
  geom_line(aes(year, lower_ci), linetype = "dashed") +
  geom_line(aes(year, upper_ci), linetype = "dashed") +
  xlab("Years") +
  ylab("Index of Abundance (base year = 3)") +
  ggtitle("Observed state of the system") +
  scale_x_continuous(expand = c(0, 0), limits = c(1, 20.5), breaks = seq(0, 20, 2)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 200)) +
  geom_hline(yintercept = 100, colour = "#9FA3A7") +
  theme_bw() + 
  theme(axis.text.x = element_text(face = "bold", color = "black", size = 10),
        axis.text.y = element_text(face = "bold", color = "black", size = 10),
        axis.title.x = element_text(color = "black", size = 13, margin = margin(t = 15)),
        axis.title.y = element_text(color = "black", size = 13, margin = margin(r = 15)),
        plot.title = element_text(color = "black", hjust = 0.5),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank())
plot.observed

comparison <- ggplot(NULL) +
  geom_line(data = actual, aes(year, pop_index, color = "act"), size = 1) +
  geom_line(data = observed, aes(year, pop_index, color = "obs"), size = 1) +
    xlab("Years") +
  ylab("Index of Abundance (base year = 3)") +
  ggtitle("Comparison") +
  scale_x_continuous(expand = c(0, 0), limits = c(1, 20.5), breaks = seq(0, 20, 2)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 200)) +
  geom_hline(yintercept = 100, colour = "#9FA3A7") +
  theme_bw() +
  theme(axis.text.x = element_text(face = "bold", color = "black", size = 10),
        axis.text.y = element_text(face = "bold", color = "black", size = 10),
        axis.title.x = element_text(color = "black", size = 13, margin = margin(t = 15)),
        axis.title.y = element_text(color = "black", size = 13, margin = margin(r = 15)),
        plot.title = element_text(color = "black", hjust = 0.5),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_colour_manual("Population Indeces", 
                      labels = c("Actual", "Observed"),
                      values = c("#FF567D", "#5699FF"), 
                      breaks = c("act", "obs"))
comparison
