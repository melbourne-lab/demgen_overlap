# Here: is there a difference between setting theta (gbar = 0) and setting gbar = theta?
# My work with the math suggests that there should be.
# (I think setting theta w/gbar = 0 produces non-normal phenotypes...)
# Here, just simulating that out (w = 1) to see if that is true
# Looking at the results, it looks like actually there isn't really a difference.
# hmm... might be because of w choice...
# SN - Feb 1 2023

# Packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(parallel)

# Clear namespace
rm(list = ls())

# Load source code
source('model_source/sim_model1_functions.R')

# Define parameters

# Trials per parameter combo
trys.per = 10
# 
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
    n.pop0 = 10000,
    # Strength of density dependence
    alpha = log(lstar) / n.pop0,
    # Ceciling-type carrying capacity just in case
    kceil = 10500
  ) %>%
  filter(s.max <= 1) %>%
  # Genetic info
  group_by(p0, lstar, s.max, h2) %>%
  mutate(
    # Gamma-parameterization
    # Initial genotype
    gbar0 = 2,
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
    timesteps = 15
  )

# 
set.seed(4523)

sim.out1 = mclapply(
  pars %>% uncount(trys.per) %>% mutate(try.no = 1:(nrow(.))) %>% split(.$try.no),
  function(pars) {
    sim(pars, theta.t = 0, init.rows = 5 * 1e5) %>%
      mutate(e_i = z_i - b_i) %>%
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

sim.sum1 = sim.out1 %>%
  group_by(gen, h2, lstar, p0) %>%
  summarise(
    across(
      everything(),
      list(bar = function(x) mean(x, na.rm = TRUE), var = function(x) var(x, na.rm = TRUE)),
      .names = "{.col}{.fn}"
    ),
    nn = n()
  )

sim.sum1 %>%
  mutate(p0 = factor(p0)) %>%
  filter(nn %in% trys.per) %>%
  ggplot(aes(x = gen, y = bbarbar)) +
  geom_line(
    aes(
      colour = p0,
      group = p0
    )
  ) +
  geom_ribbon(
    aes(
      ymin = bbarbar - 2 * sqrt(bbarvar / nn),
      ymax = bbarbar + 2 * sqrt(bbarvar / nn),
      fill = p0,
      group = p0
    ),
    alpha = 0.1
  ) +
  facet_wrap(~ paste('h2', h2, 'lstar', lstar))

sim.sum1 %>%
  mutate(p0 = factor(p0)) %>%
  filter(nn %in% trys.per) %>%
  ggplot(aes(x = gen, y = zbarbar)) +
  geom_line(
    aes(
      colour = p0,
      group = p0
    )
  ) +
  geom_ribbon(
    aes(
      ymin = zbarbar - 2 * sqrt(zbarvar / nn),
      ymax = zbarbar + 2 * sqrt(zbarvar / nn),
      fill = p0,
      group = p0
    ),
    alpha = 0.1
  ) +
  facet_wrap(~ paste('h2', h2, 'lstar', lstar))

sim.sum1 %>%
  mutate(p0 = factor(p0)) %>%
  filter(nn %in% trys.per) %>%
  ggplot(aes(x = gen, y = nbar)) +
  geom_line(
    aes(
      colour = p0,
      group = p0
    )
  ) +
  geom_ribbon(
    aes(
      ymin = nbar - 2 * sqrt(nvar / nn),
      ymax = nbar + 2 * sqrt(nvar / nn),
      fill = p0,
      group = p0
    ),
    alpha = 0.1
  ) +
  facet_wrap(~ paste('h2', h2, 'lstar', lstar))

sim.sum1 %>%
  mutate(p0 = factor(p0)) %>%
  filter(nn %in% trys.per) %>%
  ggplot(aes(x = gen, y = zvarbar)) +
  geom_line(
    aes(
      colour = p0,
      group = p0
    )
  ) +
  geom_ribbon(
    aes(
      ymin = zvarbar - 2 * sqrt(zvarvar / nn),
      ymax = zvarbar + 2 * sqrt(zvarvar / nn),
      fill = p0,
      group = p0
    ),
    alpha = 0.1
  ) +
  facet_wrap(~ paste('h2', h2, 'lstar', lstar))

###

set.seed(45230)

sim.out2 = mclapply(
  pars %>% uncount(trys.per) %>% mutate(try.no = 1:(nrow(.))) %>% split(.$try.no),
  function(pars) {
    sim(pars %>% mutate(gbar0 = 0), theta.t = 2, init.rows = 5 * 1e5) %>%
      mutate(e_i = z_i - b_i, z_i = theta_t - z_i, b_i = theta_t - b_i) %>%
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

sim.sum2 = sim.out2 %>%
  group_by(gen, h2, lstar, p0) %>%
  summarise(
    across(
      everything(),
      list(bar = function(x) mean(x, na.rm = TRUE), var = function(x) var(x, na.rm = TRUE)),
      .names = "{.col}{.fn}"
    ),
    nn = n()
  )

sim.sum2 %>%
  mutate(p0 = factor(p0)) %>%
  filter(nn %in% trys.per) %>%
  ggplot(aes(x = gen, y = bbarbar)) +
  geom_line(
    aes(
      colour = p0,
      group = p0
    )
  ) +
  geom_ribbon(
    aes(
      ymin = bbarbar - 2 * sqrt(bbarvar / nn),
      ymax = bbarbar + 2 * sqrt(bbarvar / nn),
      fill = p0,
      group = p0
    ),
    alpha = 0.1
  ) +
  facet_wrap(~ paste('h2', h2, 'lstar', lstar))

sim.sum2 %>%
  mutate(p0 = factor(p0)) %>%
  filter(nn %in% trys.per) %>%
  ggplot(aes(x = gen, y = zbarbar)) +
  geom_line(
    aes(
      colour = p0,
      group = p0
    )
  ) +
  geom_ribbon(
    aes(
      ymin = zbarbar - 2 * sqrt(zbarvar / nn),
      ymax = zbarbar + 2 * sqrt(zbarvar / nn),
      fill = p0,
      group = p0
    ),
    alpha = 0.1
  ) +
  facet_wrap(~ paste('h2', h2, 'lstar', lstar))

sim.sum2 %>%
  mutate(p0 = factor(p0)) %>%
  filter(nn %in% trys.per) %>%
  ggplot(aes(x = gen, y = nbar)) +
  geom_line(
    aes(
      colour = p0,
      group = p0
    )
  ) +
  geom_ribbon(
    aes(
      ymin = nbar - 2 * sqrt(nvar / nn),
      ymax = nbar + 2 * sqrt(nvar / nn),
      fill = p0,
      group = p0
    ),
    alpha = 0.1
  ) +
  facet_wrap(~ paste('h2', h2, 'lstar', lstar))

sim.sum1 %>%
  mutate(p0 = factor(p0)) %>%
  filter(nn %in% trys.per) %>%
  ggplot(aes(x = gen, y = zvarbar)) +
  geom_line(
    aes(
      colour = p0,
      group = p0
    )
  ) +
  geom_ribbon(
    aes(
      ymin = zvarbar - 2 * sqrt(zvarvar / nn),
      ymax = zvarbar + 2 * sqrt(zvarvar / nn),
      fill = p0,
      group = p0
    ),
    alpha = 0.1
  ) +
  facet_wrap(~ paste('h2', h2, 'lstar', lstar))

# neat.

###
# Genotypic distributions - are they the same?

set.seed(1997)

sim.out3 = mclapply(
  pars[1:3,] %>% uncount(2*trys.per) %>% mutate(try.no = 1:(nrow(.))) %>% split(.$try.no),
  function(pars) {
    sim(
      pars %>% mutate(gbar0 = ifelse(pars$try.no %% 2, 2, 0), timesteps = 4), 
      theta.t = ifelse(pars$try.no %% 2, 0, 2), 
      init.rows = 2 * 1e5
    ) %>%
      mutate(trial = pars$try.no, p0 = pars$p0)
  },
  mc.cores = 4
) %>%
  do.call(rbind, .)

sim.out3 %>%
  mutate(theta_t = factor(theta_t)) %>%
  group_by(trial, gen) %>%
  mutate(zct = z_i - mean(z_i)) %>%
  ungroup() %>%
  ggplot(aes(x = zct, group = interaction(trial, gen))) +
  geom_density(aes(colour = theta_t)) +
  facet_grid(rows = vars(p0), cols = vars(gen))

# departure from normality, but actually only for slow LH
# which actually may just be the buffering from e_i

sim.out3 %>%
  mutate(theta_t = factor(theta_t)) %>%
  # group_by(trial, gen) %>%
  # mutate(bct = b_i - mean(b_i)) %>%
  # ungroup() %>%
  ggplot(aes(x = b_i, group = interaction(trial, gen))) +
  geom_density(aes(colour = theta_t)) +
  facet_grid(rows = vars(p0), cols = vars(gen))
# actually considerable variation in breeding value var
# as time goes on... (wonder if this is just drift/small sample size?)
# conclusion... hmm might not be so bad  

# Not seeing much difference in the distribution here!
# Speed of adaptation?

sim.out3 %>%
  mutate(zdist = (z_i - theta_t) * ifelse(theta_t > 0, 1, -1)) %>%
  mutate(theta_t = factor(theta_t)) %>%
  ggplot(aes(x = zdist, group = interaction(trial, gen))) +
  geom_density(aes(colour = theta_t)) +
  facet_grid(rows = vars(p0), cols = vars(gen))

# yeah honestly these are pretty much the same lmao
# cool!
