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
trys.per = 25

# Parameters
pars = expand.grid(
  # (Equilibrium) size of initial cohort
  p0    = c(.55, .75, .95),
  # Equilibrium growth rate relative to maximum
  l.rat = c(.9, .95, .975),
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
    n.pop0 = 5000,
    # PStrength of density dependence
    alpha = log(lstar) / n.pop0,
    # Ceciling-type carrying capacity just in case
    kceil = 5500
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
    sig.m = sqrt(wfitn^2 * (sig.a^2 - (sig.p^2 - p0*sig.a^2)/(1-p0)))
  ) %>%
  ungroup() %>%
  # Other junk
  mutate(
    timesteps = 10
  )

# Run simulations
set.seed(4523)

sim.out2 = mclapply(
  pars %>% uncount(trys.per) %>% mutate(try.no = 1:(nrow(.))) %>% split(.$try.no),
  function(pars) {
    sim(pars, theta.t = 1.75, init.rows = 5 * 1e5) %>%
      mutate(
        e_i = z_i - b_i,
        b_i = theta_t - b_i,
        z_i = theta_t - z_i
      ) %>%
      group_by(gen) %>%
      summarise(
        n = n(),
        kbar = mean(age),
        bbar = mean(b_i),
        bvar = var(b_i),
        zbar = mean(e_i),
        zvar = var(e_i),
        zbar = mean(z_i),
        zvar = var(z_i)
      )  %>%
      mutate(
        trial = pars$try.no, 
        p0    = pars$p0,
        lstar = pars$lstar,
        h2    = pars$h2
      )
  },
  mc.cores = 4
) %>%
  do.call(rbind, .)

# Population size
sim.out2 %>%
  ggplot(aes(x = gen, y = n, group = interaction(trial, lstar, h2))) +
  geom_line(aes(colour = h2), linewidth = 0.25) +
  scale_colour_gradient2(low = 'red', mid = 'yellow', midpoint = .5, high = 'blue') +
  facet_wrap(paste('p0', p0) ~ paste('lstar', lstar), nrow = 3) +
  scale_y_log10()

sim.out2 %>%
  ggplot(aes(x = gen, y = zbar, group = interaction(trial, lstar, h2))) +
  geom_line(aes(colour = h2), linewidth = 0.25) +
  scale_colour_gradient2(low = 'red', mid = 'yellow', midpoint = .5, high = 'blue') +
  facet_wrap(paste('p0', p0) ~ paste('lstar', lstar), nrow = 3) 

sim.out2 %>%
  ggplot(aes(x = gen, y = zvar, group = interaction(trial, lstar, h2))) +
  geom_line(aes(colour = h2), linewidth = 0.15) +
  scale_colour_gradient2(low = 'red', mid = 'yellow', midpoint = .5, high = 'blue') +
  facet_wrap(paste('p0', p0) ~ paste('lstar', lstar), nrow = 3) 

sim.out2 %>%
  ggplot(aes(x = gen, y = n, group = interaction(trial, lstar, h2))) +
  geom_line(aes(colour = p0), linewidth = 0.25) +
  scale_colour_gradient2(low = 'red', mid = 'yellow', midpoint = .5, high = 'blue') +
  facet_wrap(paste('h2', h2) ~ paste('lstar', lstar), nrow = 3) +
  scale_y_log10()

sim.out2 %>%
  ggplot(aes(x = gen, y = bbar, group = interaction(trial, lstar, h2))) +
  geom_line(aes(colour = p0), linewidth = 0.25) +
  scale_colour_gradient2(low = 'red', mid = 'yellow', midpoint = .5, high = 'blue') +
  facet_wrap(paste('h2', h2) ~ paste('lstar', lstar), nrow = 3) 

sim.out2 %>%
  ggplot(aes(x = gen, y = zbar, group = interaction(trial, lstar, h2))) +
  geom_line(aes(colour = p0), linewidth = 0.25) +
  scale_colour_gradient2(low = 'red', mid = 'yellow', midpoint = .5, high = 'blue') +
  facet_wrap(paste('h2', h2) ~ paste('lstar', lstar), nrow = 3) 

# Look at some individual output.

sim.out3 = mclapply(
  pars %>% uncount(5) %>% mutate(try.no = 1:(nrow(.))) %>% split(.$try.no),
  function(pars) {
    sim(pars, theta.t = 1.75, init.rows = 5 * 1e6) %>%
      mutate(
        trial = pars$try.no, 
        p0    = pars$p0,
        lstar = pars$lstar,
        h2    = pars$h2
      )
  },
  mc.cores = 4
) %>%
  do.call(rbind, .)

head(sim.out3)
nrow(sim.out3)

t1 = sim.out3 %>% filter(trial %in% 1)

s1 = t1 %>%
  filter(gen < 10) %>%
  mutate(e_i = z_i - b_i, b_i = theta_t - b_i, z_i = theta_t - z_i) %>%
  group_by(gen, age) %>%
  summarise(
    n = n(),
    bbar = mean(b_i),
    bvar = var(b_i),
    zbar = mean(e_i),
    zvar = var(e_i),
    zbar = mean(z_i),
    zvar = var(z_i)
  )

s1 %>%
  ggplot(aes(x = age, y = bbar)) +
  geom_point(aes(colour = gen)) +
  geom_line(aes(colour = gen, group = gen)) +
  scale_colour_viridis_c()

s1 %>%
  ggplot(aes(x = gen, y = bbar)) +
  geom_point(aes(colour = age)) +
  geom_line(aes(colour = age, group = age)) +
  scale_colour_viridis_c()

c1 = t1 %>%
  filter(gen < 10) %>%
  mutate(e_i = z_i - b_i, b_i = theta_t - b_i, z_i = theta_t - z_i) %>%
  group_by(gen, cohort = gen - age) %>%
  summarise(
    n = n(),
    bbar = mean(b_i),
    bvar = var(b_i),
    ebar = mean(e_i),
    evar = var(e_i),
    zbar = mean(z_i),
    zvar = var(z_i)
  )

c1 %>%
  mutate(age = gen - cohort) %>%
  ggplot(aes(x = gen, y = bbar)) +
  geom_line(aes(colour = cohort, group = cohort)) +
  geom_point(aes(fill = age), shape = 21, size = 3) +
  scale_colour_viridis_c() +
  scale_fill_gradient(low = 'black', high = 'white')

c1long = c1 %>%
  select(-contains("var")) %>%
  select(-n) %>%
  pivot_longer(cols = contains("bar"), names_to = "vartype", values_to = "varval") %>%
  mutate(age = gen - cohort)

c1long %>%
  ggplot(aes(x = gen, y = varval)) +
  geom_line(aes(colour = cohort, group = cohort)) +
  geom_point(aes(fill = age), shape = 21, size = 3.5) +
  scale_colour_viridis_c() +
  scale_fill_gradient(low = 'black', high = 'white') +
  facet_wrap(~ vartype, nrow = 1)

##### Try a big batch now, get some averages

# Trials per parameter combo
trys.per = 400

set.seed(407)

sim.out4 = mclapply(
  pars %>% uncount(trys.per) %>% mutate(try.no = 1:(nrow(.))) %>% split(.$try.no),
  function(pars) {
    sim(pars, theta.t = 1.875, init.rows = 5 * 1e5) %>%
      mutate(
        e_i = z_i - b_i,
        b_i = theta_t - b_i,
        z_i = theta_t - z_i
      ) %>%
      group_by(gen) %>%
      summarise(
        n = n(),
        kbar = mean(age),
        bbar = mean(b_i),
        bvar = var(b_i),
        zbar = mean(e_i),
        zvar = var(e_i),
        zbar = mean(z_i),
        zvar = var(z_i)
      )  %>%
      mutate(
        trial = pars$try.no, 
        p0    = pars$p0,
        lstar = pars$lstar,
        h2    = pars$h2
      )
  },
  mc.cores = 4
) %>%
  do.call(rbind, .)
