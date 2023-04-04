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
  s.max = c(0.1, 0.5, 0.9),
  # Heritability of fitness
  h2    = c(.25, .5, 1),
  # Gamma squared (pheno variance / sel pressure)
  sig.z = sqrt(c(.1, .2, .4))
) %>%
  # Demographic rates
  mutate(
    # Maximum expected lifetime fitness
    w.max = 3,
    # Equilibrium lifetime fitness
    wstar = w.max * (1 - s.max) / (sqrt(1 + sig.z^2) - s.max),
    # Mean fecundity
    r     = w.max * (1 - s.max) / s.max,
    # Equilibrium population growth rate
    lstar = (s.max + w.max * (1 - s.max)) / (s.max + (w.max/wstar) * (1 - s.max)),
    # Initial population size
    n.pop0 = 20000,
    # Strength of density dependence
    alpha = log(lstar) / n.pop0,
    # Ceciling-type carrying capacity just in case
    kceil = 30000,
    p0    = (w.max * (1 - s.max)) / (w.max * (1 - s.max) + s.max)
  ) %>%
  # Genetic info
  group_by(lstar, s.max, h2, p0) %>%
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
  mutate(timesteps = 100)

# Run simulations
set.seed(4523)

sim.out = mclapply(
  pars %>% uncount(trys.per) %>% mutate(try.no = 1:(nrow(.))) %>% split(.$try.no),
  function(pars) {
    sim(pars, theta.t = 0, init.rows = 100 * 30000) %>%
      group_by(gen) %>%
      summarise(n = n())  %>%
      mutate(
        trial = pars$try.no, 
        p0    = pars$p0,
        h2    = pars$h2,
        var.z = pars$sig.z^2
      )
  },
  mc.cores = 12
) %>%
  do.call(rbind, .)
    
sim.r = sim.out %>%
  group_by(trial) %>%
  mutate(log.lam = c(diff(log(n)), NA)) %>%
  filter(!is.na(log.lam)) %>%
  group_by(gen, p0, var.z, h2) %>%
  summarise(
    llbar = mean(log.lam),
    llvar = var(log.lam),
    n     = mean(n),
    nn = n()
  )

sim.n.all = sim.out %>%
  merge(
    expand.grid(
      gen   = 0:pars$timesteps[1],
      trial = 1:(nrow(pars) * trys.per)
    ),
    all.y = TRUE
  ) %>%
  mutate(n = ifelse(is.na(n), 0, n)) %>%
  group_by(trial) %>%
  arrange(gen) %>%
  mutate(
    p0 =    p0[1],
    h2 =    h2[1],
    var.z = var.z[1]
  ) %>%
  group_by(p0, var.z, h2, gen, trial) %>%
  summarise(n = sum(n)) %>%
  group_by(p0, var.z, h2, gen) %>%
  summarise(
    nbar = mean(n),
    nvar = var(n),
    psrv = mean(n > 0),
    nn   = n()
  )

sim.n.surv = sim.out %>%
  group_by(trial) %>%
  filter(max(gen) == pars$timesteps[1]) %>%
  group_by(h2, var.z, p0, gen) %>%
  summarise(
    nbar = mean(n),
    nvar = var(n),
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
  sim.n.all,
  row.names = FALSE,
  file = 'run_sims/out/sim_results_m1_allsizes_n.csv'
)

write.csv(
  sim.n.surv,
  row.names = FALSE,
  file = 'run_sims/out/sim_results_m1_survsizes_n.csv'
)
