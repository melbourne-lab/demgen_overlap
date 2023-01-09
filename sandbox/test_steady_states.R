# Checking steady states:
# Last week, I figured out an age distribution with an equilibrium pop growth rate of lambda*
# Using this, we have populations at equilibrium population size using alpha.
# Here, I'm just testing to see if these populations remain at:
#   (1) stable population age distribution
#   (2) stable population-wide pheno-variance
#   (3) stable pheno-variance x age distribution
# SN, 9 Jan 2023

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
trys.per = 50
# 
pars = expand.grid(
  p0 = c(0.25, 0.5)
) %>%
  # Demographic rates
  mutate(
    # Max growth rate for perfectly-adapted opulation
    l.max = 1.20,
    # Growth rate with phenotypic variance
    lstar = 1.08,
    # Maximum attainable survival
    s.max = l.max * (1 - p0),
    # Mean fecundity
    r     = 1/(1-p0) - 1,
    # Initial population size
    n.pop0 = 1000,
    # PStrength of density dependence
    alpha = log(lstar) / n.pop0,
    # Ceciling-type carrying capacity just in case
    kceil = 1500
  ) %>%
  # Genetic info
  group_by(p0) %>%
  mutate(
    # Gamma-parameterization
    wfitn = 1,
    sig.e = 0,
    sig.a = sqrt(newt.method.g1(.1, 1e-8, s.max / lstar, r)),
    sig.p = sqrt(gamma.calc(sig.a^2, s.max / lstar, r)),
    mu    = 1,
    sig.m = sqrt(wfitn^2 * (sig.a^2 - (sig.p^2 - p0*sig.a^2)/(1-p0)))
  ) %>%
  ungroup() %>%
  # Other junk
  mutate(
    timesteps = 15
  )

# Run sims!

set.seed(29029)

sim.out1 = mclapply(
  pars %>% uncount(trys.per) %>% mutate(try.no = 1:(nrow(.))) %>% split(.$try.no),
  function(pars) {
    sim(pars, theta.t = 0, init.rows = 5 * 1e5) %>%
      group_by(gen) %>%
      summarise(
        n = n(),
        bvar = var(b_i)
      ) %>%
      mutate(trial = pars$try.no, s.max = pars$p0)
  },
  mc.cores = 8
) %>%
  do.call(rbind, .)
# why did this take so long...?

sim.out1 %>% head()

sim.out1 %>%
  mutate(p0 = factor(s.max)) %>%
  ggplot(aes(x = gen, y = n)) +
  geom_line(aes(colour = p0, group = trial), linewidth = 0.1)
# bingo
# okay lol so this was just because of sig.a in the params list lmao

sim.out1 %>%
  mutate(p0 = factor(s.max)) %>%
  ggplot(aes(x = gen, y = bvar)) +
  geom_line(aes(colour = p0, group = trial), linewidth = 0.1)
# bingo

sim.sum1 = sim.out1 %>%
  mutate(bvar = ifelse(is.na(bvar), 0, bvar)) %>%
  group_by(gen, s.max) %>%
  summarise(
    nbar = mean(n),
    nvar = var(n),
    bvarbar = mean(bvar),
    bvarvar = var(bvar),
    nn   = n()
  ) %>%
  filter(nn %in% max(nn))

sim.sum1 %>%
  mutate(p0 = factor(s.max)) %>%
  ggplot(aes(x = gen)) +
  geom_line(aes(y = nbar, colour = p0, group = p0)) +
  geom_ribbon(
    aes(
      ymin = nbar - 2 * sqrt(nvar / nn),
      ymax = nbar + 2 * sqrt(nvar / nn),
      fill = p0, group = p0
    ),
    alpha = 0.2
  )
# yeah some movement but this looks fine to me

sim.sum1 %>%
  mutate(p0 = factor(s.max)) %>%
  ggplot(aes(x = gen)) +
  geom_line(aes(y = bvarbar, colour = p0, group = p0)) +
  geom_ribbon(
    aes(
      ymin = bvarbar - 2 * sqrt(bvarvar / nn),
      ymax = bvarbar + 2 * sqrt(bvarvar / nn),
      fill = p0, group = p0
    ),
    alpha = 0.2
  )
# hmm... not great but not horrendous...

set.seed(44)

test.sim = sim(pars[1,], theta.t = 0, init.rows = 5 * 1e5)

test.sim %>% group_by(gen) %>% summarise(n = n()) %>%
  ggplot(aes(x = gen, y = n)) +
  geom_line()

test.sim %>% group_by(gen) %>% summarise(bvar = var(b_i)) %>%
  ggplot(aes(x = gen, y = bvar)) +
  geom_line()
# but this is def wrong...

test.sim %>%
  group_by(gen, age) %>%
  summarise(n = n()) %>%
  group_by(gen) %>%
  mutate(p = n / sum(n)) %>%
  ggplot(aes(x = age, y = p)) +
  geom_point(aes(colour = gen > 0)) +
  geom_line(aes(colour = gen > 0, group = gen), linewidth = 0.5) + theme()
  # scale_colour_viridis_c()
# hmm this looks fine

test.sim %>%
  ggplot(aes(x = b_i)) +
  geom_density(aes(group = gen, colour = gen)) +
  scale_colour_viridis_c()

### Look at 

set.seed(29029)

sim.out1 = mclapply(
  pars %>% uncount(trys.per) %>% mutate(try.no = 1:(nrow(.))) %>% split(.$try.no),
  function(pars) {
    sim(pars, theta.t = 0, init.rows = 5 * 1e5) %>%
      group_by(gen, age) %>%
      summarise(
        n = n(),
        bvar = var(b_i)
      ) %>%
      mutate(trial = pars$try.no, s.max = pars$p0)
  },
  mc.cores = 8
) %>%
  do.call(rbind, .)

sims.p.age = sim.out1 %>%
  group_by(gen, trial) %>%
  mutate(p = n / sum(n)) %>%
  group_by(gen, age, s.max) %>%
  summarise(
    pbar = mean(p),
    pvar = var(p),
    nn = n()
  ) %>%
  ungroup()

# DOES the age distribution stay constant over time?
sims.p.age %>%
  filter(age < 20, nn %in% trys.per) %>%
  ggplot(aes(x = age)) +
  geom_point(aes(y = pbar, colour = gen)) +
  geom_line(aes(y = pbar, colour = gen, group = gen)) +
  scale_colour_viridis_c() +
  facet_wrap(~ s.max)
# hell yeah it does
# wow that is beautiful.

sims.p.var = sim.out1 %>%
  mutate(bvar = ifelse(is.na(bvar), 0, bvar)) %>%
  group_by(gen, age, s.max) %>%
  summarise(
    bvarbar = mean(bvar),
    bvarvar = var(bvar),
    nn = n()
  ) %>%
  ungroup()

sims.p.var %>%
  filter(age < 20, nn %in% trys.per) %>%
  ggplot(aes(x = age)) +
  geom_point(aes(y = bvarbar, colour = gen)) +
  geom_line(aes(y = bvarbar, colour = gen, group = gen)) +
  scale_colour_viridis_c() +
  facet_wrap(~ s.max)
# it does look like it is falling but otherwise yeah this looks fine to me!

### What have we learned?

# In an adaptive environment (theta = 0), these expressions produce constant
# age distributions, variance, and age-variance distributions!
# Very cool!

### Now: do these hold under a shifting environment?

# using theta = 1.25

set.seed(14194)

tim.out1 = mclapply(
  pars %>% uncount(trys.per) %>% mutate(try.no = 1:(nrow(.))) %>% split(.$try.no),
  function(pars) {
    sim(pars, theta.t = 1.15, init.rows = 5 * 1e5) %>%
      group_by(gen) %>%
      summarise(
        n = n(),
        bvar = var(b_i)
      ) %>%
      mutate(trial = pars$try.no, s.max = pars$p0)
  },
  mc.cores = 8
) %>%
  do.call(rbind, .)

tim.sum1 = tim.out1 %>%
  mutate(bvar = ifelse(is.na(bvar), 0, bvar)) %>%
  group_by(gen, s.max) %>%
  summarise(
    nbar = mean(n),
    nvar = var(n),
    bvarbar = mean(bvar),
    bvarvar = var(bvar),
    nn   = n()
  ) %>%
  filter(nn %in% max(nn))

tim.sum1 %>% head()
tim.sum1 %>% tail()
# no extinctions? very nice.

tim.sum1 %>%
  mutate(p0 = factor(s.max)) %>%
  ggplot(aes(x = gen)) +
  geom_line(aes(y = nbar, colour = p0, group = p0)) +
  geom_ribbon(
    aes(
      ymin = nbar - 2 * sqrt(nvar / nn),
      ymax = nbar + 2 * sqrt(nvar / nn),
      fill = p0, group = p0
    ),
    alpha = 0.2
  )
# whoah... surprisingly similar
tim.sum1 %>% filter(nn < trys.per)
# and no extinctions!

tim.sum1 %>%
  mutate(p0 = factor(s.max)) %>%
  ggplot(aes(x = gen)) +
  geom_line(aes(y = bvarbar, colour = p0, group = p0)) +
  geom_ribbon(
    aes(
      ymin = bvarbar - 2 * sqrt(bvarvar / nn),
      ymax = bvarbar + 2 * sqrt(bvarvar / nn),
      fill = p0, group = p0
    ),
    alpha = 0.2
  )
# huh... interesting
# what is happening in the first few timesteps? 

set.seed(9922861)

test.sim2 = sim(pars[1,], theta.t = 1.15, init.rows = 5 * 1e5)

test.sim2 %>% group_by(gen) %>% summarise(n = n()) %>%
  ggplot(aes(x = gen, y = n)) +
  geom_line()

test.sim2 %>% group_by(gen) %>% summarise(bvar = var(b_i)) %>%
  ggplot(aes(x = gen, y = bvar)) +
  geom_line()
# but this is def wrong...

test.sim2 %>%
  group_by(gen, age) %>%
  summarise(n = n()) %>%
  group_by(gen) %>%
  mutate(p = n / sum(n)) %>%
  ggplot(aes(x = age, y = p)) +
  geom_point(aes(colour = gen > 0)) +
  geom_line(aes(colour = gen > 0, group = gen), linewidth = 0.5) + theme()
# this age distribution actually looks good to me!

test.sim2 %>%
  ggplot(aes(x = b_i)) +
  geom_density(aes(group = gen, colour = gen)) +
  scale_colour_viridis_c()
# ahhhh look at that holy shit!

test.sim2 %>%
  ggplot(aes(x = age, y = z_i)) +
  geom_segment(aes(x = 0, xend = 20, y = 1.17, yend = 1.17),
               colour = 'blue', linetype = 2) +
  geom_point(position = position_jitter(width = 0.25), size = 0.25) +
  facet_wrap(~ gen)

test.sim2 %>% group_by(gen) %>% summarise(bbar = mean(b_i)) %>%
  ggplot(aes(x = gen, y = bbar)) +
  geom_segment(aes(x = 0, xend = 15, y = 1.17, yend = 1.17),
               colour = 'blue', linetype = 2) +
  geom_line()

test.sim2 %>% 
  group_by(cohort = gen - age, gen) %>%
  summarise(bbar = mean(b_i), n = n()) %>%
  filter(n() > 3) %>%
  ggplot(aes(x = gen, y = bbar, colour = cohort)) +
  geom_segment(aes(x = 0, xend = 15, y = 1.17, yend = 1.17),
               colour = 'blue', linetype = 2) +
  geom_point(aes(size = n)) +
  geom_line(aes(group = cohort)) +
  scale_color_viridis_c()

# hmm...

# Try looking now at age breakdowns

set.seed(14194)

tim.out2 = mclapply(
  pars %>% uncount(trys.per) %>% mutate(try.no = 1:(nrow(.))) %>% split(.$try.no),
  function(pars) {
    sim(pars, theta.t = 1.15, init.rows = 5 * 1e5) %>%
      group_by(gen, age) %>%
      summarise(
        n = n(),
        bvar = var(b_i)
      ) %>%
      mutate(trial = pars$try.no, s.max = pars$p0)
  },
  mc.cores = 8
) %>%
  do.call(rbind, .)

tims.p.age = tim.out2 %>%
  group_by(gen, trial) %>%
  mutate(p = n / sum(n)) %>%
  group_by(gen, age, s.max) %>%
  summarise(
    pbar = mean(p),
    pvar = var(p),
    nn = n()
  ) %>%
  ungroup()

# DOES the age distribution stay constant over time?
tims.p.age %>%
  filter(age < 20, nn %in% trys.per) %>%
  ggplot(aes(x = age)) +
  geom_point(aes(y = pbar, colour = gen)) +
  scale_colour_viridis_c() +
  facet_wrap(~ s.max)
# interesting -
# there is some change in the "slow" group but it might be temporary?

tims.p.age %>%
  filter(nn %in% trys.per) %>%
  ggplot(aes(x = gen)) +
  geom_line(aes(y = pbar, colour = age, group = age)) +
  geom_point(aes(y = pbar, colour = age)) +
  scale_colour_viridis_c() +
  facet_wrap(~ s.max)
# whoooah shit this is cool!
# travelling wave - there is a cohort that seems very large
# (who are they?) but otherwise the structure looks constant

tims.p.var = tim.out2 %>%
  mutate(bvar = ifelse(is.na(bvar), 0, bvar)) %>%
  group_by(gen, age, s.max) %>%
  summarise(
    bvarbar = mean(bvar),
    bvarvar = var(bvar),
    nn = n()
  ) %>%
  ungroup()

tims.p.var %>%
  filter(age < 20, nn %in% trys.per) %>%
  ggplot(aes(x = age)) +
  geom_point(aes(y = bvarbar, colour = gen)) +
  geom_line(aes(y = bvarbar, colour = gen, group = gen)) +
  scale_colour_viridis_c() +
  facet_wrap(~ s.max)
# variance structure *looks* similar over time

tims.p.var %>%
  filter(nn %in% trys.per) %>%
  ggplot(aes(x = gen)) +
  geom_line(aes(y = bvarbar, colour = age, group = age)) +
  geom_point(aes(y = bvarbar, colour = age)) +
  scale_colour_viridis_c() +
  facet_wrap(~ s.max)

tims.p.var = tim.out2 %>%
  mutate(bvar = ifelse(is.na(bvar), 0, bvar)) %>%
  group_by(gen, cohort = gen - age, s.max) %>%
  summarise(
    bvarbar = mean(bvar),
    bvarvar = var(bvar),
    nn = n()
  ) %>%
  ungroup()

tims.p.var %>%
  filter(nn %in% trys.per) %>%
  ggplot(aes(x = gen)) +
  geom_line(aes(y = bvarbar, colour = cohort, group = cohort)) +
  geom_point(aes(y = bvarbar, colour = cohort)) +
  scale_colour_viridis_c() +
  facet_wrap(~ s.max)
# meh not so useful

### LESSONS:
# - Age structure, variance, variance structure ARE ~constant
#   in a stable environment
#   (so long as you're careful about the parameters you define in the initial
#    pars list)
# - Under a single environmental change, age structure and variance *DO*
#   seem to change but only slightly (and only for slow LH)
#   - age structure looks like there is a cohort that makes up a "travelling wave"
#   - for some reason I haven't figured out yet, there is also a temporary var
#     increase
# 
# Neat shit!
