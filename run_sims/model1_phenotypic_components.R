##########
# Simulated experiment
# Here: recording population-level phenotypic components
# SN - init 9 Mar 2023
##########

# Packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(parallel)

# Clear namespace
rm(list = ls())

# Load source code
source('model_source/sim_model1_functions.R')

### Load in parameters

# Trials per parameter combo
trys.per = 200

# Parameters
pars = expand.grid(
  # (Equilibrium) size of initial cohort
  p0    = c(.5, .75, .95),
  # Equilibrium growth rate relative to maximum
  l.rat = c(.8),
  # Heritability of fitness
  h2    = c(.25, .5, 1)
) %>%
  # Demographic rates
  mutate(
    # Maximum growth rate
    l.max = 2,
    # Maximum survival
    s.max = l.max * (1 - p0),
    # Equilibrium growth rate
    lstar = l.rat * l.max,
    # Mean fecundity
    r     = p0 / (1-p0),
    # Initial population size
    n.pop0 = 20000,
    # Strength of density dependence
    alpha = log(lstar) / n.pop0,
    # Ceciling-type carrying capacity just in case
    kceil = 30000
  ) %>%
  filter(s.max <= 1) %>%
  # Genetic info
  group_by(p0, lstar, s.max, h2) %>%
  mutate(
    # Gamma-parameterization
    # wfitn = 1 in gamma parameterization
    wfitn = 1,
    # Phenotypic standard deviation in new cohorts
    sig.0 = sqrt(newt.method.g1(.1, 1e-8, s.max / lstar, r)),
    # Breeding value standard deviation in new cohorts
    sig.a = sqrt(h2 * sig.0^2),
    # Non-inherited standard dxeviation in new cohorts
    sig.e = sqrt((1-h2) * sig.0^2),
    # Population-wide breeding value standard deviation
    sig.p = sqrt(gamma.a.calc(sig.a^2, s.max / lstar, r, sig.e^2)),
    mu    = 1,
    sig.m = sqrt(wfitn^2 * (sig.a^2 - (sig.p^2 - p0*sig.a^2)/(1-p0))),
    gbar0 = 2
  ) %>%
  ungroup() %>%
  # Other junk
  mutate(
    timesteps = 50
  )

# Run simulations
set.seed(4523)

sim.out2 = mclapply(
  pars %>% uncount(trys.per) %>% mutate(try.no = 1:(nrow(.))) %>% split(.$try.no),
  function(pars) {
    sim(pars, theta.t = 0, init.rows = 50 * 30000) %>%
      mutate(
        e_i = z_i - b_i
      ) %>%
      group_by(gen) %>%
      summarise(
        n = n(),
        bbar = mean(b_i),
        bvar = ifelse(any(age > 0), var(b_i[age > 0]), NA),
        ebar = mean(e_i),
        evar = ifelse(any(age > 0), var(e_i[age > 0]), NA),
        zbar = mean(z_i),
        zvar = ifelse(any(age > 0), var(z_i[age > 0]), NA)
      )  %>%
      mutate(
        trial = pars$try.no, 
        p0    = pars$p0,
        h2    = pars$h2
      )
  },
  mc.cores = 12
) %>%
  do.call(rbind, .)

# nrow(sim.out2)

write.csv(
  sim.out2,
  file = 'run_sims/out/sim_results_m1_phtype.csv',
  row.names = FALSE
)
