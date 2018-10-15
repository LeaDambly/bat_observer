# Simulation functions
# 14 06 18 - Lea I Dambly
growth <- function(Nt = N0, mu.r = 1, se.r = 0.1) {
  # set the growth rate Lambda
  lambda <- rnorm(n = 1, mean = mu.r, sd = se.r)
  
  # calculate new population size
  Ntplus1 <- round(Nt * lambda, 0)
  return(Ntplus1)
}

switch1 <- function(swi){
  # percentage of roosts that will switch
  perc <-  runif(1, 0.1, 0.3)

  # create new tibble with roosts that have switched
  swi2 <- swi %>% 
    dplyr::filter(roost <= z & count > 0) %>% 
    sample_frac(perc, replace = FALSE) %>% 
    arrange(roost)
  
  # create new tibble without the roosts that have switched
  swi3 <- anti_join(swi, swi2, by = c("count", "roost", "year", "swyr")) %>% 
    arrange(roost)
  
  # create new tibble where count is zero (because the original roosts are now empty)
  swi4 <- swi2 %>% 
    select(roost, year, swyr) %>% 
    add_column(count = 0)
  
  # add swi4 back to swi3 - now the original roosts are empty
  swi3 <- swi3 %>% 
    bind_rows(swi4) %>% 
    arrange(roost) 
  
  # rename roost to orgroost for the newly switched roosts
  swi2 <- swi2 %>%
    select(count, roost, swyr, year) %>% 
    rename(orgroost = roost)
  
  # add the newly switched roosts back to the entire tibble
  swi <- swi3 %>% bind_rows(swi2)
  
  # number of new roosts without roost number
  num <- sum(is.na(swi$roost))
  
  # the new roosts have NAs, we need to make them 0
  swi$roost[is.na(swi$roost)] <- 0
  
  # now we index the new roosts
  swi <- swi %>%
    mutate(roost = ifelse(roost == 0, seq(max(swi$roost)+1, max(swi$roost)+num, 1), roost)) %>%
    arrange(roost) %>%
  # and give the original z roosts an NA for the orgroost variable
    mutate(orgroost = ifelse(roost <= z, NA, orgroost)) %>%
  # finally add another year to swyr
  mutate(swyr = ifelse((roost > z & count > 0), swyr+1, 0))
  
  return(swi)
}

switch2 <- function(swi){
  # percentage of roosts that will switch back to their original
  perc <-  runif(1, 0.5, 0.8)
  
  # create new tibble with the roosts that are coming back - has to have been 2 yrs since the switch
  swi2 <- swi %>%
    dplyr::filter(!is.na(orgroost) & swyr >= 2) %>% 
    sample_frac(perc)
  
  # if there aren't any to switch back, the function will end
  if (dim(swi2)[1] == 0){
    return(swi)
  }else{

    # create new tibble where switched roosts that will switch back now are removed
    swi3 <- anti_join(swi, swi2, by = c("roost", "count", "year", "orgroost", "swyr"))

    # make the orgroosts their actual roosts again
    swi2 <- swi2 %>% 
      mutate(roost = orgroost) %>% 
      select(-orgroost)
  
    # and put everything back together
    swi <- left_join(swi3, swi2, by = "roost") %>%
      select(-c(year.y, swyr.y)) %>%
      mutate(count = ifelse(is.na(count.y), count.x, count.y)) %>%
      rename(year = year.x) %>%
      rename(swyr = swyr.x) %>%
      select(count, roost, year, swyr, orgroost)
      
    return(swi)
  }
}

split <- function(spl, k = 120, n = 3, yr = nyr){
  perc <-  runif(1, 0.6, 0.9)
  
  w <- spl # the population for that year per site
  
  w2 <- w %>%
    
    sample_frac(perc) %>%
    
    mutate(a = ifelse(count >= k, round(count/n), NA)) %>% 
    
    mutate(b = ifelse(is.na(a), count, (a*n)-a)) %>% 
    
    select(-count)
  
  w <- full_join(w, w2, by = c("roost", "year"))
  
  w <- w %>%
    
    mutate(c = ifelse(is.na(b), count, b))
  
  w2 <- w %>% select(a) %>% na.omit(a) %>% rename(., "c" = "a")
  
  w <- w[c(3,6)]
  
  w <-  w %>% 
    
    bind_rows(w2) %>%
    
    replace_na(list(year = yr)) %>%
    
    tibble::rownames_to_column(., var = "roost")
  
  w$roost <- as.numeric(w$roost)
  
  colnames(w)[3] <- c("count")
  
  return(w)
}

accessibility <- function(obs, a = 0, b = 1, c = 1, d = 4){
  # people will likely only monitor sites that are accessible
  # TO DO: Comment on code

  seqq <- seq(a,b,1)
  
  slope <- (max(seqq) - min(seqq))/(max(obs$count)-min(obs$count))
  
  int <- (slope*min(obs$count)-min(seqq))*(-1)
  
  obs <- obs %>% rowwise(.) %>% mutate(acc = ((slope*count)+int)^c)
  
  lth <- length(obs$count)/d
  
  perc <- runif(1, lth*2, lth*3)
  
  acc <- obs[sample(obs$count, perc, replace = FALSE, prob = obs$acc),]
  
  return(acc)
}

# TO DO sample with increasing probability? (ie the longer they sample, the more likely they are to go on)
start_stop_monitoring <- function(mon) {
  # start monitoring at random and stop at random
  sam <- mon %>% group_by(roost) %>% sample_frac(., size = runif(1, 0.65, 0.85)) %>% arrange(roost, year)
  sam <- as.data.frame(sam)
  return(sam)
  }

rmse_func <- function(gamsp) {
  error <- residuals.gam(gamsp)
  rmse <- sqrt(mean((error)^2))
  rmse
}

# GAM functions----
# These functions have been taken and adapted from Rachel Fewsters's GAM code 
# which is available online at https://www.stat.auckland.ac.nz/~fewster/gams/R/
# Fewster et al. 2000. Ecology

# function to fit GAM for a range of degrees of freedom to assess
gam_func <- function(sp, dfvec = c(6)){
  gam <- gam(count~as.factor(roost) + s(year, fx=TRUE, k=dfvec+1), 
             family = poisson(link = log), data = sp)
  gam
}


index_func <- function(sp, gam, dfvec = c(6)) {
  if(length(sp$count[is.na(sp$count)]) > 0) stop(
    "no missing data allowed")
  
  fit_func <- function(dfval)
  {
    pred <- predict.gam(gam,type="terms")
    
    srow <- length(pred[1,])   
    Nentries <- length(sp$count)
    
    # setting index year
    # gives the position corresponding to the first occurrence of each year
    yrs_pos <- (1:Nentries)[!duplicated(sp$year)]
    
    # gives you the value of 'year' of yrs_pos
    yrs_entries <- sp$year[!duplicated(sp$year)]
    
    
    yrs_pred <- pred[yrs_pos, srow]
    
    
    # Both are then ordered  
    years_ord <- sort(yrs_entries)
    pred_ord <- yrs_pred[order(yrs_entries)]
    
    
    # calculate index based on year 3 (observed vs expected)
    index_df <- exp(pred_ord)/exp(pred_ord[3])
    
    index_df <- unname(index_df)
    
    index_df <- data.frame(years_ord, index_df)    
    
    index_df
    
    
  }
  
  index_all <- fit_func(dfvec)
  
  index_all
}

boot_func <- function(sp, defree = 6) {
  uniq_roost <- unique(sp$roost)
  Nentries <- length(sp$count)
  Nyears <- length(unique(sp$year))
  Nroosts <- length(uniq_roost)
  
  # take sample of roosts to be included in the resample 
  sam <- sample(uniq_roost, replace = T)
  
  # lists the rows of data frame sp for inclusion in the resample
  # ragged list of rows to be included in resample (with repetitions listed separately) 
  elements_func <- function(x) {
    (1:Nentries)[sp$roost == x]
  }
  elements <- lapply(sam, elements_func)
  
  # vector of roost levels for the resample (should go from 1:nroosts because want to pick nroosts w/ repetition)
  dr_roost <- rep(1:Nroosts, as.vector(unlist(lapply(elements, length))))
  
  # extract year and count data for resample 
  elements <- as.vector(unlist(elements))
  data_resample <- data.frame(roost = dr_roost, 
                              year = sp$year[elements],
                              count = sp$count[elements])
  
  # fit GAM on resample
  dr_gam <- gam(count~ s(year, fx=TRUE, k=defree+1) + as.factor(roost), family = 
                  poisson(link = log), data = data_resample)	#
  
  # Data frame for prediction w/ all years, roosts and covariates
  year_base <- min(sp$year)-1
  dr_new <- expand.grid(year=year_base+1:Nyears, roost=sort(unique(dr_roost)))
  
  # predict from model for all years
  result <- predict.gam(dr_gam, newdata=dr_new, type = "terms")                
  srow <- length(result[1,])   
  pred_allyears <-  result[1:Nyears, srow]
  
  indices <- exp(pred_allyears)/exp(pred_allyears[3])     
  indices
}

outer_boot_func <- function(sp, defree = 6, Nrep = 39) {
  #if(!is.character(indexfile)) stop("filenames have to be character strings")	
  if(length(sp$count[is.na(sp$count)]) > 0) stop("no missing data allowed")
  
  replicate_func <- function(i) {
    print(i)
    rep_indices <- boot_func(sp, defree = defree)
    #write(rep.indices, file = indexfile, ncolumns = 4, append = T)
    rep_indices
  }
  index_matrix <- sapply(1:Nrep, replicate_func)
  matrix(unlist(index_matrix), byrow = T, nrow = Nrep)
}

ci_func <- function(index_matrix, conf = 0.95) {
  # index.ci.func: 
  # Given a matrix of bootstrapped index values of the form 
  # index.matrix = sp.bootind.Nrep.df, this function extracts 
  # the lower and upper (100*conf)% confidence limits for the abundance
  # indices. 
  
  Nyears <- ncol(index_matrix)
  Nrep <- nrow(index_matrix)
  
  alp <- (1 - conf)/2 
  
  # Interpolation is not ideal, so discard some rows of the bootstrap
  # matrix if (Nrep+1)*alp is not an integer:
  if(abs((Nrep + 1) * alp - round((Nrep + 1) * alp)) > 1e-05)
    stop("need to discard rows of index.matrix, or change conf, 
         so that (Nrep+1)*(1-conf)/2 is an integer.")
  
  lowerpt <- (Nrep + 1) * alp
  upperpt <- (Nrep + 1) * (1 - alp)
  
  
  inner_func <- function(yr) {
    sort(index_matrix[, yr])[c(lowerpt, upperpt)]
  }
  
  index_ci <- sapply(seq(1, Nyears), inner_func)
  dimnames(index_ci) <- list(c("lower", "upper"), NULL)
  index_ci
}
