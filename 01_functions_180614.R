# Simulation functions
# 14 06 18 - Lea I Dambly
growth <- function(Nt = N0, mu.r = 0, se.r = 0.1) {
  lambda <- exp(rnorm(n = 1, mean = mu.r, sd = se.r)) #growth rate lambda
  Ntplus1 <- round(Nt * lambda, 0) #growth for that year
  return(Ntplus1)
}

split <- function(spl, k = 80, n = 3){
  perc <-  runif(1, 0.7, 1.0)
  
  x <- spl # the population for that year per site
  
  x2 <- x %>%
    
    sample_frac(perc) %>%
    
    rowwise() %>%
    
    mutate(a = ifelse(count >= k, round(count/n), NA)) %>% 
    
    mutate(b = ifelse(is.na(a), count, (a*n)-a)) %>% 
    
    select(-count)
 
  x <- full_join(x, x2, by = c("roost", "year"))
  
  x <- x %>%
    
    mutate(c = ifelse(is.na(b), count, b))
 
  x2 <- x %>% select(a) %>% na.omit(a) %>% rename(., "c" = "a")
   
  x <- x[c(3,6)]
  
  x <-  x %>% 
    
    bind_rows(x2) %>%
    
    replace_na(list(year = nyr)) %>%
    
    tibble::rownames_to_column(., var = "roost")
    
    x$roost <- as.numeric(x$roost)
  
  return(x)
}

# TO DO attractiveness may change with time
accessability <- function(obs){
  # people will likely only monitor sites that are easily accessible/more attractive
  roost <- obs %>% filter(year == 1) %>% select(roost)
  
  # 1 = not attractive, 3 = very attractive
  seqq <- c(1,1,2,2,2,3,3,3,3,3)
  r <- length(roost$roost)/length(seqq)
  acc <- sample(rep.int(seqq, r))
  roost$acc <- acc
  roost <- as.data.frame(roost)
  obs$acc <- roost[match(obs$roost, roost$roost), 2]
  
  # weighted random sized (between 0.7-0.9%) random sample
  perc <- runif(1, 0.7, 0.9)
  obs2 <- obs %>% filter(year == 1) %>% sample_frac(., size = perc, weight = acc)
  obs <- semi_join(obs, obs2, by = "roost")
  mon <- as.data.frame(obs[c(1:3)])
  return(mon)
  }

# TO DO sample with increasing probability? (ie the longer they sample, the more likely they are to go on)
# something with weights?
start_stop_monitoring <- function(mon) {
  # start monitoring at random and stop at random
  sam <- mon %>% group_by(roost) %>% sample_frac(., size = runif(1, 0.6, 0.9)) %>% arrange(roost, year)
  sam <- as.data.frame(sam)
  return(sam)
  }

rmse_func <- function(gamsp) {
  error <- residuals.gam(gamsp)
  rmse <- sqrt(mean((error)^2))
  rmse
}

# not modular version----
sim1 <- function(N0 = 100, z = 600) {
  growth <- function(Nt = N0, mu.r = 0, se.r = 0.1, k = 150) {
    lambda <- exp(rnorm(n = 1, mean = mu.r, sd = se.r)) #growth rate lambda
    Ntplus1 <- round(Nt * lambda, 0) #growth for that year
    if(Ntplus1 >= k){ #if N is larger than k, it'll split and creates new roost z
    Ntplus1 <- round(Ntplus1/2); 
    # the new roost needs to start at the size of ntplus1...
    # plus, it can't just run for 20 yrs again. it needs to run 20-t years
    z <<- z + 1
    }
    return(Ntplus1)
  }
  N <- list()
  N <- growth(N0)
  for(t in 2:nyr){
    N[t] <- growth(N[t-1])
  }
  return(N)
}


# Geomtric growth V2----
# creates the whole time series in one go (ONE GROWTH!)
gro <- function(N0 = 100, nyr = 20, mu.r = 0, se.r = 0.1){
  yrs <- seq(1, nyr, 1)
  lambda <- exp(rnorm(n = 1, mean = mu.r, sd = se.r))
  Nt <- round(N0 * lambda^yrs)
  return(Nt)
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
