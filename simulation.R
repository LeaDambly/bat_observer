# This code runs bat simulations in parallel

# libraries----
library(foreach)
library(doParallel)

# make parallel----
# choose number of cores to run on
numCores <- detectCores() - 1
# make local cluster
cl <- makeCluster(numCores)
registerDoParallel(cl)

# run sim----
# i is number of repetitions
# packages don't need to be loaded manually but do have to be installed
simulation <- foreach(i = 1:10, .packages = c("tibble", "dplyr", "mgcv", "reshape2")) %dopar% {
  # define path to folder where functions are
  srcf <- "N:/RProjects/Bat_Observer/01_Simulation/Functions"
  R.utils::sourceDirectory(srcf)

  ## set up parameters
  # z is the number of roosts
  z <- 200

  # x is the average population size per roost
  x <- 72

  # y is the current year of the simulation
  y <- 1

  # N0: a vector of random numbers that are Poisson distributed and based on our starting parameters x and z
  N0 <- rpois(lam = x, n = z)

  # bats_1: a tibble with 4 columns
  # count: the population size
  # roost: the respective roost number
  # year: the respective year
  # abyr ("abandonment year"): a variable that tallies how long it has been since a population has abandoned a roost
  bats_1 <- tibble(count = N0) %>%
    add_column(roost = 1:z, year = y, abyr = 0)

  # move ahead one year
  y <- y + 1

  # 'pogr' stands for population growth
  pogr <- bats_1

  # 20 year loop of population growth and movement
  while (y <= 20) {
    pogr <- pogr %>%
      mutate(year = year + 1)

    pogr <- pogr %>%
      group_by(roost) %>%
      mutate(count = growth(Nt = count)) %>%
      ungroup()

    pogr <- abandon(pogr, p = 0.01)
    pogr <- reoccupy(pogr)

    assign(paste("bats", y, sep = "_"), pogr)

    y <- y + 1
    gc()
  }

  rm(pogr)

  all <- mget(ls(pattern = "bats_*")) %>%
    bind_rows() %>%
    as.tibble()

  rm(list = ls(pattern = "bats_*"))

  act <- acast(all, roost ~ year, value.var = "count", fill = 0)
  act <- melt(act, value.name = "count") %>%
    rename(roost = Var1, year = Var2) %>%
    as.tibble()

  dl <- adss(act, c = 0.5)
  dl_ml <- commitment(dl, pr = 0.1)
  dl_mh <- commitment(dl, pr = 0.9)

  dh <- adss(act, c = 1.5)
  dh_ml <- commitment(dh, pr = 0.1)
  dh_mh <- commitment(dh, pr = 0.9)

  inda <- indsp_func(act) %>%
    mutate(type = "Actual")

  inda <- inda %>%
    mutate(rmse = sqrt(mean((inda$index - index)^2)))

  indb <- indsp_func(dl_ml) %>%
    mutate(type = "ADSSS low + commitment low") %>%
    mutate(rmse = sqrt(mean((inda$index - index)^2)))

  indc <- indsp_func(dl_mh) %>%
    mutate(type = "ADSS low + commitment high") %>%
    mutate(rmse = sqrt(mean((inda$index - index)^2)))

  indd <- indsp_func(dh_ml) %>%
    mutate(type = "ADSS high + commitment low") %>%
    mutate(rmse = sqrt(mean((inda$index - index)^2)))

  inde <- indsp_func(dh_mh) %>%
    mutate(type = "ADSS high + commitment high") %>%
    mutate(rmse = sqrt(mean((inda$index - index)^2)))


  allind <- bind_rows(inda, indb, indc, indd, inde)
  allind$type <- as.factor(allind$type)
  allind$rep <- i

  return(allind)
  gc()
}

save(simulation, file = "N:/RProjects/Bat_Observer/02_Results/test1.RData")
stopCluster(cl)
