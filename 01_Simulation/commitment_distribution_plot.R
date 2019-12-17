library(ggplot2)
library(jpeg)

low <- as.data.frame(rbinom(400, 2, 0.2)+1) %>%
  rename(years = "rbinom(400, 2, 0.2) + 1") %>%
  mutate(scenario = "low")

med <- as.data.frame(rbinom(400, 2, 0.5)+1) %>%
  rename(years = "rbinom(400, 2, 0.5) + 1") %>%
  mutate(scenario = "med")

high <- as.data.frame(rbinom(400, 2, 0.8)+1) %>%
  rename(years = "rbinom(400, 2, 0.8) + 1") %>%
  mutate(scenario = "high")

all <- bind_rows(low, med, high)

all$scenario <- factor(all$scenario, levels = unique(all$scenario))
all$scenario <- factor(all$scenario, levels = c("low", "med", "high"))

tiff("com_dist.tiff", units = "in", width = 10, height = 5, res = 300)
ggplot(all) +
  geom_bar(aes(years)) +
  facet_wrap(vars(scenario)) +
  theme_bw()
dev.off()
