library(dplyr)
library(tidyr)

rm(list = ls())

# initial population size
n.pop0 = 100
# optimal phenotypic value
s.max = 0.9
# mean fecundity per individual
r     = (1.2 / (s.max)) - 1
# selection strength (this seems like a good value to use)
wfitn = 3
# genetic variance (should be approx. static)
sig.a = 0.5
# standard deviation of the selection pressure
sig.e = 0
# standard dev. of environmental phenoytpic noise
p.pos = 0.5
# Initial frequency of the positive allele
alpha = 0
# optimal genotype
theta_t = 2.5
# # maximum age
# lspan = 2

###

pop0 = data.frame(
    i   = 1:n.pop0,
    gen = 1,
    age = 1 + rgeom(n.pop0, (1-s.max)),
    fem = as.logical(sample(0:1, n.pop0, replace = TRUE)),
    b_i = rnorm(n.pop0, 0, sig.a)
) %>%
  mutate(
    z_i = rnorm(n.pop0, b_i, sig.e),
    s_i = s.max * exp(-(z_i - theta_t)^2 / (2*wfitn^2)),
    r_i = rpois(n.pop0, lambda = ifelse(fem, 2 * r, 0))
  )

# Mate (n.b., mating occurs after selection - does this matter?)

# survival
# pop = pop %>%
#     mutate(surv = as.logical(rbinom(n = nrow(.), size = 1, prob = s_i))) %>%
#     filter(surv) %>%
#     select(-surv)

offspring = cbind(
  # Mom info
  pop0 %>%
    filter(fem) %>%
    select(-c(i, gen, age, fem, z_i, s_i)) %>%
    rename(mom_b_i = b_i),
  # Dad info
  pop0 %>%
    sample_n(size = sum(fem),
             weight = as.numeric(!fem) / sum(as.numeric(!fem)),
             replace = TRUE) %>%
    select(b_i) %>%
    rename(dad_b_i = b_i)
) %>%
  mutate(midp_b_i = (mom_b_i + dad_b_i) / 2) %>%
  uncount(r_i) %>%
  mutate(i = max(pop0$i) + 1:nrow(.),
         gen = max(pop0$gen) + 1,
         age = 1,
         fem = as.logical(sample(0:1, nrow(.), replace = TRUE)),
         b_i = rnorm(nrow(.), midp_b_i, sig.a),
         z_i = rnorm(nrow(.), b_i, sig.e),
         s_i = s.max * exp(-(z_i - theta_t)^2 / (2*wfitn^2)),
         r_i = rpois(nrow(.), lambda = ifelse(fem, 2 * r, 0))) %>%
  select(i, gen, age, fem, b_i, z_i, s_i, r_i)

# Survival (of all? or just of those alive previously?)

pop1 = rbind(pop0 %>% mutate(age = age + 1), offspring) %>%
      mutate(surv = as.logical(rbinom(n = nrow(.), size = 1, prob = s_i))) %>%
      filter(surv) %>%
      select(-surv)

nrow(pop0)
nrow(pop1) # hmm... this is big.

sum(pop0$r_i)

# Function for propagation

propagate.sim = function(popn, theta) {
  
  cur.gen = max(popn$gen)
  
  if (sum(popn$r_i)) {
  
  offspring = cbind(
    # Mom info
    popn %>%
      filter(fem) %>%
      select(-c(i, gen, age, fem, z_i, s_i, theta_t)) %>%
      rename(mom_b_i = b_i),
    # Dad info
    popn %>%
      sample_n(size = sum(fem),
               weight = as.numeric(!fem) / sum(as.numeric(!fem)),
               replace = TRUE) %>%
      select(b_i) %>%
      rename(dad_b_i = b_i)
  ) %>%
    mutate(midp_b_i = (mom_b_i + dad_b_i) / 2) %>%
    uncount(r_i) %>%
    mutate(i = max(popn$i) + 1:nrow(.),
           gen = cur.gen + 1,
           age = 1,
           fem = as.logical(sample(0:1, nrow(.), replace = TRUE)),
           b_i = rnorm(nrow(.), midp_b_i, sig.a),
           z_i = rnorm(nrow(.), b_i, sig.e),
           s_i = s.max * exp(-(z_i - theta_t)^2 / (2*wfitn^2)),
           r_i = rpois(nrow(.), lambda = ifelse(fem, 2 * r, 0)),
           theta_t = theta) %>%
    select(i, gen, age, fem, b_i, z_i, s_i, r_i, theta_t)
  
  } else {
    offspring = popn %>% sample_n(size = 0)
  }
  
  # Survival (of all? or just of those alive previously?)
  
  popn.out = rbind(
      popn %>% mutate(age = age + 1, gen = cur.gen + 1), 
      offspring
    ) %>%
    # recalculate s_i based on current theta value
    mutate(theta_t = theta) %>%
    mutate(s_i = s.max * exp(-(z_i - theta_t)^2 / (2*wfitn^2))) %>%
    mutate(surv = as.logical(rbinom(n = nrow(.), size = 1, prob = s_i))) %>%
    filter(surv) %>%
    select(-surv)
  
  return(popn.out)
  
}

propagate.sim(pop0, 2.5)

# Look at variation in size from trial to trial
# (with the same initial population size)

list.one = vector('list', length = 400)

# Feel like this is still kinda slow...
for (trial in 1:400) list.one[[trial]] = propagate.sim(pop0, 2.5) %>% mutate(k = trial)

sapply(list.one, function(x) x %>% summarise(n = nrow(.), b = mean(b_i))) %>%
  t() %>%
  plot()

sapply(list.one, nrow) %>% mean()
sapply(list.one, nrow) %>% (function(x) sd(x) / sqrt(length(x)))
sapply(list.one, function(x) mean(x$b_i)) %>% mean()
sapply(list.one, function(x) mean(x$b_i)) %>% (function(x) sd(x) / sqrt(length(x)))


# hmm... interesting
# what exactly is that growth rate though? not what I expected.
# also genotypic change is very slow here... high longevity I guess

# Try repeated propagation 

pop1 = propagate.sim(pop0, theta = 2.5)
pop2 = propagate.sim(pop1, theta = 2.5)
pop3 = propagate.sim(pop2, theta = 2.5)

pop3 %>% nrow()

rbind(pop0, pop1, pop2, pop3) %>% group_by(gen) %>% summarise(n = n())

# advancing phenotype
pop1 = propagate.sim(pop0, theta = 2.5)
pop2 = propagate.sim(pop1, theta = 2.55)
pop3 = propagate.sim(pop2, theta = 2.6)

pop3 %>% nrow()

rbind(pop0, pop1, pop2, pop3) %>% group_by(gen) %>% summarise(n = n())
# hmm...

### Set up parameters

pars = data.frame(
  # number of time steps
  timesteps = 20,
  # initial population size
  n.pop0 = 100,
  # optimal phenotypic value
  s.max = 0.9,
  # mean fecundity per individual
  r     = (1.2 / (s.max)) - 1,
  # selection strength (this seems like a good value to use)
  wfitn = 3,
  # genetic variance (should be approx. static)
  sig.a = 0.5,
  # standard deviation of the selection pressure
  sig.e = 0,
  # standard dev. of environmental phenoytpic noise
  p.pos = 0.5,
  # Initial frequency of the positive allele
  alpha = 0,
  # optimal genotype
)

### Wrapper for adding data to data frame

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

init.sim = function(params, theta0) {
  
  # initial population size
  size0 = params$n.pop0
  # optimal phenotypic value
  s.max = params$s.max 
  # mean fecundity per individual
  r     = params$r
  # selection strength (this seems like a good value to use)
  wfitn = params$wfitn
  # genetic variance (should be approx. static)
  sig.a = params$sig.a
  # standard dev. of environmental phenoytpic noise
  sig.e = params$sig.e
  # initial genotype
  gbar0 = ifelse(grepl('gbar0', names(params)), params$gbar0, 0)
  # density dependence strength
  # alpha = ifelse(grepl('alpha', names(params)), params$alpha, 0)
  
  popn = data.frame(
    # Unique identifier
    i   = 1:size0,
    # Generation
    gen = 0,
    # Age (age distibution for bernoulli survival is geometric)
    age = 1 + rgeom(size0, (1-s.max)),
    # Sex (TRUE = female)
    fem = as.logical(sample(0:1, size0, replace = TRUE)),
    # Breeding value (sig.a is additive genetic variance, according to Lynch & Walsh)
    b_i = rnorm(size0, 0, sig.a)
  ) %>%
  mutate(
    # Phenotype
    z_i = rnorm(size0, b_i, sig.e),
    # Survival ()
    s_i = s.max * exp(-(z_i - theta0)^2 / (2*wfitn^2)),
    # Offspring (0 for males, Poisson draw for females)
    r_i = rpois(size0, lambda = ifelse(fem, 2 * r, 0)),
    # Phenotypic optimum in this time step
    theta_t = theta0
  )
  
  return(popn)
  
}

pop0 = init.sim(params = pars, theta0 = 2.5)

head(pop0)
mean(pop0$r_i)
mean(pop0$b_i)

pop0 = init.sim(params = pars %>% mutate(sig.e = 0.5), theta0 = 2.5)

head(pop0)
mean(pop0$r_i)
mean(pop0$b_i)
mean(pop0$z_i)
var(pop0$b_i)
var(pop0$z_i)

### Write a function pulling it all together

sim = function(params, theta_t, init.rows, init.popn = NULL) {
  
  ### Declaring variables
  # (must this be done here?)
  
  # length of simulations
  timesteps = params$timesteps
  # initial population size
  size0 = params$n.pop0
  # optimal phenotypic value
  s.max = params$s.max 
  # mean fecundity per individual
  r     = params$r
  # selection strength (this seems like a good value to use)
  wfitn = params$wfitn
  # genetic variance (should be approx. static)
  sig.a = params$sig.a
  # standard dev. of environmental phenoytpic noise
  sig.e = params$sig.e
  # initial genotype
  gbar0 = ifelse(grepl('gbar0', names(params)), params$gbar0, 0)
  # density dependence strength
  # alpha = ifelse(grepl('alpha', names(params)), params$alpha, 0)
  
  ### Initialize data frame
  all.data = data.frame(
    i   = rep(NA, init.rows),
    gen = rep(NA, init.rows),
    fem = rep(NA, init.rows),
    b_i = rep(NA, init.rows),
    z_i = rep(NA, init.rows),
    s_i = rep(NA, init.rows),
    r_i = rep(NA, init.rows),
    theta_t = rep(NA, init.rows)
  )
  
  ### Check lengths of theta parameter
  if (!length(theta_t) %in% c(timesteps + 1, 1)) {
    stop("Length of theta vector does not match other inputs")
  } else if (length(theta_t) == 1) {
    theta_t = rep(theta_t, timesteps)
  }
    
  ### Initialize population
  if (!is.null(init.popn)) {
    init.popn = init.popn %>%
      mutate(theta_t = theta_t[1]) %>%
      mutate(s_i = s.max * exp(-(z_i - theta_t)^2 / (2*wfitn^2)))
  } else {
    init.popn = init.sim(params = pars, theta0 = theta_t[1])
  }
  
  ### Add initial population to data frame
  all.data = dim.add(df = all.data, rows = init.rows, addition = init.popn)
  
  ### Iterate simulation
  
  prev.gen = init.popn
  
  for (tt in 1:timesteps) {
    if (nrow(prev.gen)) {
      # Propagate sim (one time step)
      popn = propagate.sim(popn = prev.gen, theta = theta[tt + 1])
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

### HAS NOT YET BEEN TESTED