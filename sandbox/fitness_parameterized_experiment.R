# script looking at an experimental design where lifetime fitness is held constant
# rather than maximum intrinsic growth rate

library(ggplot2)
library(dplyr)
library(tidyr)

rm(list = ls())

source('model_source/sim_model1_functions.R')

trys.per = 10

pars = expand.grid(
  # (Equilibrium) size of initial cohort
  s.max = c(0.1, 0.5, 0.9),
  # Heritability of fitness
  h2    = c(.25, .5, 1)
) %>%
  # Demographic rates
  mutate(
    # Gamma squared (pheno variance / sel pressure)
    sig.z = sqrt(.25),
    # Maximum expected lifetime fitness
    w.max = 2.5,
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
    wfitn = 2.5,
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
    timesteps = 30
  )

# Run simulations
set.seed(4523)

sim.out = mclapply(
  pars %>% uncount(trys.per) %>% mutate(try.no = 1:(nrow(.))) %>% split(.$try.no),
  function(pars) {
    sim(pars, theta.t = 0, init.rows = 30 * 30000) %>%
      group_by(gen) %>%
      mutate(e_i = z_i - b_i) %>%
      summarise(
        n = n(),
        bbar = mean(b_i),
        bvar = var(b_i),
        ebar = mean(e_i),
        evar = var(e_i),
        zbar = mean(z_i),
        zvar = var(z_i)
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

sim.out %>%
  mutate(p0 = factor(round(p0, 2))) %>%
  ggplot(aes(x = gen, y = n, group = trial)) +
  geom_line(aes(colour = p0)) +
  scale_y_log10() +
  facet_wrap(~ h2)

sim.out %>%
  mutate(p0 = factor(round(p0, 2))) %>%
  ggplot(aes(x = gen, y = zbar, group = trial)) +
  geom_line(aes(colour = p0)) +
  facet_wrap(~ h2)

sim.out %>%
  mutate(p0 = factor(round(p0, 2))) %>%
  ggplot(aes(x = gen, y = bbar, group = trial)) +
  geom_line(aes(colour = p0)) +
  facet_wrap(~ h2)

sim.out %>%
  mutate(p0 = factor(round(p0, 2))) %>%
  ggplot(aes(x = gen, y = zvar, group = trial)) +
  geom_line(aes(colour = p0)) +
  facet_wrap(~ h2)
