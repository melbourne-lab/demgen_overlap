### Script comparing sims in a slowly moving environmnt with no variation
### SN - 14 Apr 2021

library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyselect)

# Clear namespace
rm(list = ls())

# Source files allowing sims
source('model_source/sim_functions_perennial.R')
source('model_source/sim_functions_annual.R')
source('model_source/sim_aux_handling_functions.R')

# Define parameters
pars = data.frame(
  end.time = 30,
  init.row = 1e5,
  max.age = 5,
  n.pop0 = 100,
  n.loci = 25,
  w.max = 2,
  wfitn = sqrt(1 / 0.14 / 2),
  sig.e = sqrt(0.5),
  pos.p = 0.5,
  alpha = 0.005
)

# Define number of trials
n.trials = 1000

# Define variance of environmental fluctuations
sig.theta = sqrt(0.5)
# Define mean temporal trend in environment
delta.tht = 0.1

# Define source populations.
pop0 = init.simp(params = pars %>% mutate(n.pop0 = pars$n.pop0 * n.trials),
                 theta0 = 0) %>%
  mutate(trial = ((0:(nrow(.)-1)) %/% pars$n.pop0) + 1)

# Initialize objects for storing sim results
list.age1 = vector('list', n.trials)
list.age3 = vector('list', n.trials)
list.age5 = vector('list', n.trials)

# Run simulations
set.seed(440089)

for (tr in 1:n.trials) {
  
  # Note: for each trio of life histories,
  # the initial population (sans age distribution)
  # and environment are the same
  
  # Setup
  pop.init = pop0 %>% filter(trial %in% tr)
  theta_t  = (1:pars$end.time)*delta.tht + rnorm(pars$end.time, 0, sig.theta)
  
  ### Run annual trial
  sim.out = sim1(params = pars, 
                 theta_t = theta_t,
                 init.popn = pop.init %>% select(-trial, age))
  
  list.age1[[tr]] = cbind(
    summarise.demo(sim.out),
    summarise.geno(sim.out, pars) %>% select(-gen)
  )
  
  ### Run age three trials
  
  sim.out = simp(params = pars %>% mutate(max.age = 3),
                 theta_t = theta_t,
                 init.popn = pop.init %>% 
                   select(-trial) %>% 
                   mutate(age = sample(3, size = nrow(.), replace = TRUE)))
  
  list.age3[[tr]] = cbind(
    summarise.demo(sim.out),
    summarise.geno(sim.out, pars) %>% select(-gen)
  )
  
  ### Run age five trials
  
  sim.out = simp(params = pars,
                 theta_t = theta_t,
                 init.popn = pop.init %>% select(-trial))
  
  list.age5[[tr]] = cbind(
    summarise.demo(sim.out),
    summarise.geno(sim.out, pars) %>% select(-gen)
  )
  
  print(tr)
  
}

# Cat outputs together and write csv
rbind(
  unroll.sums(list.age1) %>% mutate(age = 1),
  unroll.sums(list.age3) %>% mutate(age = 3),
  unroll.sums(list.age5) %>% mutate(age = 5)
) %>%
  write.csv('run_sims/dynamic/dynamic_lovar_out.csv',
            row.names = FALSE)
