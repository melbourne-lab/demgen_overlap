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
    # Strength of density dependence
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
  scale_colour_gradient2(low = 'red', mid = 'yellow', midpoint = .75, high = 'blue') +
  facet_wrap(paste('h2', h2) ~ paste('lstar', lstar), nrow = 3) +
  scale_y_log10()

sim.out2 %>%
  ggplot(aes(x = gen, y = bbar, group = interaction(trial, lstar, h2))) +
  geom_line(aes(colour = p0), linewidth = 0.25) +
  scale_colour_gradient2(low = 'red', mid = 'yellow', midpoint = .75, high = 'blue') +
  facet_wrap(paste('h2', h2) ~ paste('lstar', lstar), nrow = 3) 

sim.out2 %>%
  ggplot(aes(x = gen, y = zbar, group = interaction(trial, lstar, h2))) +
  geom_line(aes(colour = p0), linewidth = 0.25) +
  scale_colour_gradient2(low = 'red', mid = 'yellow', midpoint = .75, high = 'blue') +
  facet_wrap(paste('h2', h2) ~ paste('lstar', lstar), nrow = 3) 

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
        ebar = mean(e_i),
        evar = var(e_i),
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

sim.out4 %>%
  ggplot(aes(x = gen, y = n, group = interaction(trial, lstar, h2))) +
  geom_line(aes(colour = h2), linewidth = 0.05) +
  scale_colour_gradient2(low = 'red', mid = 'yellow', midpoint = .5, high = 'blue') +
  facet_wrap(paste('p0', p0) ~ paste('lstar', lstar), nrow = 3) +
  scale_y_log10()

sim.sum = 
  rbind(
    expand.grid(gen = 1:pars$timesteps[1], trial = unique(sim.out4$trial)) %>%
      merge(sim.out4 %>% distinct(trial, p0, lstar, h2)) %>%
      mutate(n = 0, kbar = NA, bbar = NA, bvar = NA, zbar = NA, zvar = NA) %>%
      select(names(sim.out4)),
    sim.out4
  ) %>%
  group_by(gen, trial, h2, p0, lstar) %>%
  summarise(across(c(n, kbar, bbar, bvar, zbar, zvar), function(x) sum(x, na.rm = TRUE))) %>%
  group_by(gen, h2, p0, lstar) %>%
  summarise(
    across(
      c(n, kbar, bbar, bvar, zbar, zvar), 
      list(bar = function(x) mean(x, na.rm = TRUE), var = function(x) var(x, na.rm = TRUE)),
      .names = "{.col}{.fn}"
    ),
    nn = sum(n > 0)
  )

head(sim.sum)
tail(sim.sum)

sim.sum %>%
  mutate(p0 = factor(p0)) %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = nbar, 
      group = p0, 
      colour = p0
    )
  ) +
  geom_ribbon(
    aes(
      ymin = nbar - sqrt(nvar / trys.per),
      ymax = nbar + sqrt(nvar / trys.per),
      fill = p0,
      group = p0
    ),
    alpha = 0.25
  ) +
  scale_y_log10() +
  facet_wrap(h2 ~ lstar)

sim.sum %>%
  mutate(p0 = factor(p0), ext = factor(nn < trys.per)) %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = bbarbar, 
      group = interaction(p0, ext), 
      colour = p0,
      linetype = ext
    )
  ) +
  geom_ribbon(
    aes(
      ymin = bbarbar - sqrt(bbarvar / nn),
      ymax = bbarbar + sqrt(bbarvar / nn),
      fill = p0,
      group = interaction(p0, ext),
      alpha = ext
    ),
  ) +
  scale_alpha_manual(values = c(.1, .5)) +
  facet_wrap(h2 ~ lstar)

sim.sum %>%
  filter(nn %in% trys.per) %>%
  mutate(p0 = factor(p0)) %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = bbarbar, 
      group = p0, 
      colour = p0
    )
  ) +
  geom_ribbon(
    aes(
      ymin = bbarbar - sqrt(bbarvar / nn),
      ymax = bbarbar + sqrt(bbarvar / nn),
      fill = p0,
      group = p0,
    ),
    alpha = 0.5
  ) +
  facet_wrap(h2 ~ lstar)

sim.sum %>%
  filter(nn %in% trys.per) %>%
  mutate(p0 = factor(p0)) %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = zbarbar, 
      group = p0, 
      colour = p0
    )
  ) +
  geom_ribbon(
    aes(
      ymin = zbarbar - sqrt(zbarvar / nn),
      ymax = zbarbar + sqrt(zbarvar / nn),
      fill = p0,
      group = p0,
    ),
    alpha = 0.5
  ) +
  facet_wrap(h2 ~ lstar)

sim.sum %>%
  filter(nn %in% trys.per) %>%
  mutate(p0 = factor(p0)) %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = zvarbar, 
      group = p0, 
      colour = p0
    )
  ) +
  geom_ribbon(
    aes(
      ymin = zvarbar - sqrt(zvarvar / nn),
      ymax = zvarbar + sqrt(zvarvar / nn),
      fill = p0,
      group = p0,
    ),
    alpha = 0.5
  ) +
  facet_wrap(h2 ~ lstar)

sim.sum %>%
  filter(nn %in% trys.per) %>%
  
  mutate(p0 = factor(p0)) %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = zvarbar, 
      group = p0, 
      colour = p0
    )
  ) +
  geom_ribbon(
    aes(
      ymin = zvarbar - sqrt(zvarvar / nn),
      ymax = zvarbar + sqrt(zvarvar / nn),
      fill = p0,
      group = p0,
    ),
    alpha = 0.5
  ) +
  facet_wrap(h2 ~ lstar)
