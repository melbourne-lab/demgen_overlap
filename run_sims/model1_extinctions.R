##########
# Simulated experiment
# Here: trying to capture age structure for populations
# SN - init 6 Mar 2023
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
trys.per = 1000

# Parameters
pars = expand.grid(
  # (Equilibrium) size of initial cohort
  p0    = c(.5, .75, .95),
  # Equilibrium growth rate relative to maximum
  l.rat = c(.8, .9, .95),
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

sim.out = mclapply(
  pars %>% uncount(trys.per) %>% mutate(try.no = 1:(nrow(.))) %>% split(.$try.no),
  function(pars) {
    sim(pars, theta.t = 0, init.rows = 50 * 30000) %>%
      group_by(gen) %>%
      summarise(
        n = n()
      ) %>%
      mutate(
        trial = pars$try.no, 
        p0    = pars$p0,
        lstar = pars$lstar,
        h2    = pars$h2
      )
  },
  mc.cores = 16
) %>%
  do.call(rbind, .)
    
sim.r = sim.out %>%
  group_by(trial) %>%
  mutate(log.lam = c(diff(log(n)), NA)) %>%
  filter(!is.na(log.lam)) %>%
  group_by(gen, p0, lstar, h2) %>%
  summarise(
    llbar = mean(log.lam),
    llvar = var(log.lam),
    n     = mean(n),
    nn = n()
  )

sim.n = sim.out %>%
  merge(
    expand.grid(
      # (t=0 not needed here - everyone initialized with non-zero)
      gen   = 1:pars$timesteps[1],
      p0    = c(.5, .75, .95),
      lstar = c(.8, .9, .95) * pars$l.max[1],
      h2    = c(.25, .5, 1)
    ),
    all.x = TRUE, all.y = TRUE
  ) %>%
  mutate(n = ifelse(is.na(n), 0, n)) %>%
  group_by(p0, lstar, h2, gen, trial) %>%
  summarise(n = sum(n)) %>%
  group_by(p0, lstar, h2, gen) %>%
  summarise(
    nbar = mean(n),
    nvar = var(n),
    psrv = mean(n > 0),
    nn   = n()
  )

# sim.n %>%
#   mutate(p0 = factor(p0)) %>%
#   ggplot(aes(x = gen)) +
#   geom_line(
#     aes(
#       y = nbar,
#       group = p0,
#       colour = p0
#     )
#   ) +
#   geom_ribbon(
#     aes(
#       ymin = nbar - 2 * sqrt(nvar / trys.per),
#       ymax = nbar + 2 * sqrt(nvar / trys.per),
#       group = p0,
#       fill = p0
#     ),
#     alpha = 0.2
#   ) +
#   facet_wrap(h2 ~ lstar)

write.csv(
  sim.r,
  row.names = FALSE,
  file = 'run_sims/out/sim_results_m1_sizes_r.csv'
)

write.csv(
  sim.n,
  row.names = FALSE,
  file = 'run_sims/out/sim_results_m1_sizes_n.csv'
)
