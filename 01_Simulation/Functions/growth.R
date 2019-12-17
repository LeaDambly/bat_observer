# growth function: calculates the population change from one year to the next
#
# Takes the following parameters
# Nt: the previous year's population size of a roost
# mu.r: the mean growth rate
# se.r: the standard deviation of the growth rate
#
# It returns Ntplus1: The current population size of a roost

growth <- function(Nt, mu.r = 1, se.r = 0.1) {

  # set the growth rate Lambda
  lambda <- rnorm(n = 1, mean = mu.r, sd = se.r)

  # calculate new population size
  Ntplus1 <- round(Nt * lambda, 0)
  return(Ntplus1)
}
