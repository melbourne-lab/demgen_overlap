##### Source code for model 2 (single-stage model, viability selection on offspring)

##### Load packages
library(dplyr)
library(tidyr)

### Suppress output from summarise()
options(dplyr.summarise.inform = FALSE)

##### Aux functions

### Wrapper for adding data to data frame
# (a pre-allocated data frame to improve performance)
dim.add = function(df, rows, addition) {
  if (sum(is.na(df$i)) < nrow(addition)) {
    df = df %>%
      rbind(df %>%
              sample_n(size = rows) %>%
              mutate_all(function(x) NA))
  }
  if (nrow(addition)) {
    df[which.max(is.na(df$i)) + 0:(nrow(addition)-1),] = addition
  }
  return(df)
}

### Wrapper for initializing sim

init.sim2 = function(params, theta0) {
  
  # initial population size
  size0 = params$n.pop0
  # optimal phenotypic value
  s     = params$s 
  # mean fecundity per individual
  r.max = params$r.max
  # selection strength (this seems like a good value to use)
  wfitn = params$wfitn
  # genetic variance (should be approx. static)
  sig.a = params$sig.a
  # standard dev. of environmental phenoytpic noise
  sig.e = params$sig.e
  # initial genotype
  gbar0 = ifelse(any(grepl('gbar0', names(params))), params$gbar0, 0)
  # density dependence strength
  alpha = ifelse(any(grepl('alpha', names(params))), params$alpha, 0)
  # ceiling-like carrying capacity term
  kceil = ifelse(any(grepl('ceil' , names(params))), params$kceil, Inf)
  # equilibrium population growth rate 
  lstar = ifelse(any(grepl('lstar', names(params))), params$lstar, 1)
  
  popn = data.frame(
    # Unique identifier
    i   = 1:size0,
    # Generation
    gen = 0,
    # Age
    age = rgeom(n = size0, prob = 1 - (s / lstar)),
    # Sex (TRUE = female)
    fem = as.logical(sample(0:1, size0, replace = TRUE))
  ) %>%
    mutate(
      b_i = rnorm(size0, mean = gbar0, sd = sig.a),
      # Phenotype
      z_i = rnorm(size0, mean = b_i, sd = sig.e),
      # Offspring (0 for males, Poisson draw for females)
      r_i = rpois(size0, lambda = ifelse(fem, 2 * r.max, 0)),
      # Phenotypic optimum in this time step
      theta_t = theta0
    )
  
  # If initial size is above the specified ceiling, truncate
  if (nrow(popn) > kceil) { 
    popn = popn %>% 
      slice_sample(kceil, replace = FALSE) %>%
      # Sort rows to order individuals by index
      arrange(i)
  }
  
  return(popn)
  
}

### Function for porpagating one generation of simulated populations forward
# NOTE: here, survival (and selection) are first
# under this model, replacement rate of individual i should be s_i (1 + r_i)
propagate.sim2 = function(popn, params, theta) {
  
  # Pull out parameters of interest
  
  # Maximum survival
  s     = params$s
  # mean fecundity per individual
  r.max = params$r.max
  # selection strength (this seems like a good value to use)
  wfitn = params$wfitn
  # genetic variance (should be approx. static)
  sig.a = params$sig.a
  # standard dev. of environmental phenoytpic noise
  sig.e = params$sig.e
  # mutation rate per capita
  mu    = ifelse(any(grepl('mu', names(params))), params$mu, 0)
  # standard dev. of mutations
  sig.m = ifelse(any(grepl('sig.m', names(params))), params$sig.m, 0)
  # Strength of density dependence
  alpha = ifelse(any(grepl('alpha', names(params))), params$alpha, 0) 
  # ceiling-like carrying capacity term
  kceil = ifelse(any(grepl('ceil' , names(params))), params$kceil, Inf)
  
  # Current time step (useful for adding time information to finished data frame)
  cur.gen = max(popn$gen)
  
  # Survival (of all individuals)
  popn.surv = popn %>% 
    # simulate survival step using Bernoulli draw for each individual
    filter(as.logical(rbinom(n = nrow(.), size = 1, prob = s * exp(-alpha * nrow(.))))) %>%
    # Handle parents - iterating forward age, current generation
    mutate(age = age + 1, gen = cur.gen + 1)  %>%
    # get new number of offspring per parent
    mutate(r_i = rpois(nrow(.), lambda = ifelse(fem, 2 * r.max, 0)))
  
  
  # Perform mating to get offspring (if there are fecund females)
  if (sum(popn.surv$r_i) & sum(popn.surv$fem) & sum(!popn.surv$fem)) {
    
    offspring = cbind(
      # Mom info
      popn.surv %>%
        filter(fem) %>%
        select(-c(i, gen, age, fem, z_i, theta_t)) %>%
        rename(mom_b_i = b_i),
      # Dad info
      popn.surv %>%
        sample_n(size = sum(fem),
                 weight = as.numeric(!fem) / sum(as.numeric(!fem)),
                 replace = TRUE) %>%
        select(b_i) %>%
        rename(dad_b_i = b_i)
    ) %>%
      # Get midparent value for each mating pair
      mutate(midp_b_i = (mom_b_i + dad_b_i) / 2) %>%
      # Duplicate midparent value according to r_i, number of offspring per pair
      uncount(r_i) %>%
      mutate(i = max(popn$i) + 1:nrow(.),           # assign ID number
             gen = cur.gen + 1,                     # add generation number
             age = 0,                               # all offspring are age 0
             fem = as.logical(sample(0:1, nrow(.), replace = TRUE)), # assign sex
             b_i = rnorm(nrow(.), midp_b_i, sqrt(var(popn$b_i) / 2)) +
               rbinom(nrow(.), 1, mu) * rnorm(nrow(.), 0, sig.m),
             # assign breeding value (mean is midparent, sd is sqrt(additive var))
             z_i = rnorm(nrow(.), b_i, sqrt(sig.e^2 * (1 + var(b_i))/(1+var(b_i) - sig.e^2))),      # assign phenotype (mean is breeding value, sd is sd of env. variance)
             r_i = rpois(nrow(.), lambda = ifelse(fem, 2 * r.max, 0)), # re-draw number of offspring per mother
             theta_t = theta) %>%                   # phenotypic optimum in this time step
      select(i, gen, age, fem, b_i, z_i, r_i, theta_t) %>%
      # Perform viability selection on offspring
      filter(as.logical(rbinom(n = nrow(.), size = 1, prob = exp(-(z_i - theta_t)^2 / (2*wfitn^2)))))
    
  } else {
    # If sum(r_i) is zero, there are zero offspring from mating this generation
    offspring = popn.surv %>% sample_n(size = 0)
  }
  
  # Combine offspring and parents
  popn.out = rbind(popn.surv, offspring) 
  
  # If offspring + parents is above the ceiling carrying capacity, truncate
  if (nrow(popn.out) > kceil) {
    popn.out = popn.out %>% 
      # If population is above the ceiling, sample just the ceiling
      slice_sample(n = kceil, replace = FALSE) %>%
      # Sort rows to order individuals by index
      arrange(i)
  }
  
  # Return
  return(popn.out)
  
}

### Wrapper function for simulating a population

sim2 = function(params, theta.t, init.rows, init.popn = NULL) {
  
  ### Declaring variables
  # (must this be done here?)
  
  # length of simulations
  timesteps = params$timesteps
  # initial population size
  size0 = params$n.pop0
  # optimal phenotypic value
  s     = params$s 
  # mean fecundity per individual
  r.max = params$r.max
  # selection strength (this seems like a good value to use)
  wfitn = params$wfitn
  # genetic variance (should be approx. static)
  sig.a = params$sig.a
  # standard dev. of environmental phenoytpic noise
  sig.e = params$sig.e
  # initial genotype
  gbar0 = ifelse(any(grepl('gbar0', names(params))), params$gbar0, 0)
  # density dependence strength
  alpha = ifelse(any(grepl('alpha', names(params))), params$alpha, 0)
  # ceiling-like carrying capacity term
  kceil = ifelse(any(grepl('ceil' , names(params))), params$kceil, Inf)
  
  ### Initialize data frame
  all.data = data.frame(
    i   = rep(NA, init.rows),
    gen = rep(NA, init.rows),
    age = rep(NA, init.rows),
    fem = rep(NA, init.rows),
    b_i = rep(NA, init.rows),
    z_i = rep(NA, init.rows),
    r_i = rep(NA, init.rows),
    theta_t = rep(NA, init.rows)
  )
  
  ### Check lengths of theta parameter
  if (!length(theta.t) %in% c(timesteps+1, 1)) {
    stop("Length of theta vector does not match other inputs")
  } else if (length(theta.t) == 1) {
    theta.t = rep(theta.t, timesteps+1)
  }
  
  ### Initialize population
  if (!is.null(init.popn)) {
    init.popn = init.popn %>% mutate(theta_t = theta.t[1])
  } else {
    init.popn = init.sim2(params = params, theta0 = theta.t[1])
  }
  
  ### Add initial population to data frame
  all.data = dim.add(df = all.data, rows = init.rows, addition = init.popn)
  
  ### Iterate simulation
  prev.gen = init.popn
  
  for (tt in 1:timesteps) {
    if (nrow(prev.gen)) {
      # Propagate sim (one time step)
      popn = propagate.sim2(popn = prev.gen, params = params, theta = theta.t[tt+1])
      # Add to data frame
      all.data = dim.add(df = all.data, rows = init.rows, addition = popn)
      # Update "previous generation"
      prev.gen = popn
    }
  }
  
  ### Remove empty rows at the end
  all.data = all.data %>% filter(!is.na(i))
  
  return(all.data)
  
}

