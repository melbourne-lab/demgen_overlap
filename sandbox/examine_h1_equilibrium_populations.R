###
# Okay - I have accepted (although I still don't like) the fact that genetic
# variance will vary across treatments. But with that said, do the other things
# I've worked out (population age structure and genetic variance) look like they
# produce equilibrium conditions?
#   (looks like the answer to this question is - yes)
# And how do they respond differently to a single environmental shift?
#   (obviously only looked at this super crudely, but, there are differences
#    but another weird thing I'm seeing is that after they reach the new steady state
#    they don't return to their original size - think the mutations might be wrong?
#    might be because it is all determined by genetic variance?
#    maybe mutations need to be different... but how?)

##### Setup

# Load packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(parallel)

# Clear namespace 
rm(list = ls())

# Load in functions
source('model_source/sim_model1_functions.R')

### Define parameters

l.max = 1.2
ntrys = 5000

pars = expand.grid(p.0 = c(0.25, 0.5)) %>%
  # Define survival and reproduction
  mutate(
    s.max = (1 - p.0) * l.max,
    r     = (l.max / s.max) - 1,
  ) %>%
  # Get gamma^2_0
  group_by(p.0) %>%
  mutate(
    gsqr0 = newt.method.g1(1, 1e-8, s.max, r)
  ) %>%
  ungroup() %>%
  # Get genetic and phenotypic variances, selection strength
  mutate(
    wfitn = sqrt(10),
    sig.a = sqrt(wfitn^2 * gsqr0),
    sig.e = 0
  ) %>%
  # Get other control parameters
  mutate(
    timesteps = 5,
    n.pop0 = 1000,
    alpha = 0,
    gsqrd = gamma.calc(gsqr0, s.max, r),
    kceil = 1200,
    mu = 1, 
    sig.m = sqrt(wfitn^2 * (gsqr0 - (gsqrd - p.0*gsqr0)/(1-p.0)))
  )

# Okay let's see what happens lol

#######-------------------------------------------------

set.seed(55)

sim.out = mclapply(
  pars %>% uncount(ntrys) %>% mutate(try.no = 1:(nrow(.))) %>% split(.$try.no),
  function(pars) {
    sim(pars, theta.t = 0, init.rows = 7 * 1e3) %>%
      group_by(gen) %>%
      summarise(
        n = n(),
        bvar = var(b_i)
      ) %>%
      mutate(trial = pars$try.no, s.max = pars$s.max)
  },
  mc.cores = 6
) %>%
  do.call(rbind, .)

head(sim.out)

sim.sum = sim.out %>% 
  group_by(s.max, gen) %>%
  summarise(
    nbar = mean(n),
    nvar = var(n),
    bvarbar = mean(bvar),
    bvarvar = var(bvar),
    nn = n()
  )

sim.sum %>%
  mutate(s.max = factor(s.max)) %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = nbar,
      colour = s.max
    )
  ) +
  geom_ribbon(
    aes(
      ymin = nbar - 2 * sqrt(nvar / nn),
      ymax = nbar + 2 * sqrt(nvar / nn),
      fill = s.max
    ),
    alpha = 0.1
  )

sim.sum %>%
  mutate(s.max = factor(s.max)) %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = bvarbar,
      colour = s.max
    )
  ) +
  geom_ribbon(
    aes(
      ymin = bvarbar - 2 * sqrt(bvarvar / nn),
      ymax = bvarbar + 2 * sqrt(bvarvar / nn),
      fill = s.max
    ),
    alpha = 0.1
  )

# Okay I guess this is good enough!
# Let's see what happens with an environmental change...
# (eyes emoji)

#######-------------------------------------------------

ntrys = 1000

set.seed(9201)

sim.out = mclapply(
  pars %>% uncount(ntrys) %>% mutate(try.no = 1:(nrow(.))) %>% split(.$try.no),
  function(pars) {
    sim(pars %>% mutate(timesteps = 9), theta.t = 3, init.rows = 7 * 1e3) %>%
      group_by(gen) %>%
      summarise(
        n = n(),
        bvar = var(b_i)
      ) %>%
      mutate(trial = pars$try.no, s.max = pars$s.max)
  },
  mc.cores = 6
) %>%
  do.call(rbind, .)

head(sim.out)

sim.sum = sim.out %>% 
  group_by(s.max, gen) %>%
  summarise(
    nbar = mean(n),
    nvar = var(n),
    bvarbar = mean(bvar),
    bvarvar = var(bvar),
    nn = n()
  )

sim.sum %>%
  mutate(s.max = factor(s.max)) %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = nbar,
      colour = s.max
    )
  ) +
  geom_ribbon(
    aes(
      ymin = nbar - 2 * sqrt(nvar / nn),
      ymax = nbar + 2 * sqrt(nvar / nn),
      fill = s.max
    ),
    alpha = 0.1
  )

sim.sum %>%
  mutate(s.max = factor(s.max)) %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = bvarbar,
      colour = s.max
    )
  ) +
  geom_ribbon(
    aes(
      ymin = bvarbar - 2 * sqrt(bvarvar / nn),
      ymax = bvarbar + 2 * sqrt(bvarvar / nn),
      fill = s.max
    ),
    alpha = 0.1
  )

# What happens if we are looking/thinking by age class?

ntrys = 100

set.seed(9201)

sim.out = mclapply(
  pars %>% uncount(ntrys) %>% mutate(try.no = 1:(nrow(.))) %>% split(.$try.no),
  function(pars) {
    sim(pars %>% mutate(timesteps = 9), theta.t = 3, init.rows = 7 * 1e3) %>%
      group_by(gen, age) %>%
      summarise(
        n = n(),
        bvar = var(b_i)
      ) %>%
      mutate(trial = pars$try.no, s.max = pars$s.max)
  },
  mc.cores = 6
) %>%
  do.call(rbind, .)

age.p = rbind(
  expand.grid(
    gen = 0:9,
    age = 0:max(sim.out$age),
    n = 0,
    trial = 1:max(sim.out$trial)
  ) %>%
    merge(y = sim.out %>% ungroup() %>% distinct(trial, s.max)),
  sim.out %>% select(-bvar)
) %>%
  group_by(age, gen, trial, s.max) %>%
  summarise(n = sum(n)) %>%
  group_by(gen, trial, s.max) %>%
  mutate(p = n / sum(n))

age.p %>%
  filter(trial %/% ntrys < 10, age < 30) %>%
  ggplot(aes(x = age, y = p, colour = gen)) +
  # geom_point(alpha = 0.1) +
  geom_line(aes(group = interaction(trial, gen)), linewidth = 0.1) +
  scale_colour_viridis_c() +
  theme(panel.background = element_rect(fill = 'black')) +
  facet_wrap(~ s.max, nrow = 2)

age.p %>%
  filter(trial %/% ntrys < 20, age < 30) %>%
  ggplot(aes(x = age, y = n, colour = gen)) +
  # geom_point(alpha = 0.1) +
  geom_line(aes(group = interaction(trial, gen)), linewidth = 0.1) +
  scale_colour_viridis_c() +
  theme(panel.background = element_rect(fill = 'black')) +
  facet_wrap(~ s.max, nrow = 2)

age.p %>%
  filter(age < 20) %>%
  group_by(age, gen, s.max) %>%
  filter(any(n > 0)) %>%
  summarise(pbar = mean(p)) %>%
  ggplot(aes(x = age, y = pbar, group = gen, colour = gen)) +
  #geom_point() +
  geom_line(linewidth = 2) +
  scale_colour_viridis_c() +
  theme(panel.background = element_rect(fill = 'black')) +
  facet_wrap(~ s.max, nrow = 2)
  
  # ggplot(aes(x = pbar, y = bvarbar, colour = gen)) +
  # geom_point(size = 3) +
  # geom_line(aes(group = gen), linewidth = 0.5) +
  # theme(panel.background = element_rect(fill = 'black')) +
  # facet_wrap(~ s.max)

set.seed(8008)

test.sim = sim(pars %>% filter(p.0 %in% 0.25) %>% mutate(timesteps = 9), 
               theta.t = 3, 
               init.rows = 7 * 1e3)

test.sim %>%
  filter(gen < 9, age < 20) %>%
  ggplot(aes(x = age, y = z_i)) +
  geom_segment(aes(x = 0, xend = 19, y = 3, yend = 3), linetype = 2, colour = 'blue') +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.5) +
  facet_wrap(~ gen)

test.sim %>%
  filter(age < 16) %>%
  ggplot(aes(x = gen, y = z_i)) +
  geom_segment(aes(x = 0, xend = 15, y = 3, yend = 3), linetype = 2, colour = 'blue') +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.5) +
  facet_wrap(~ age)

test.sim %>%
  filter(gen < 9) %>%
  group_by(i) %>% mutate(first.gen = min(gen)) %>% ungroup() %>%
  ggplot(aes(x = i, y = z_i)) +
  geom_segment(aes(x = 0, xend = max(i), y = 3, yend = 3), linetype = 2, colour = 'blue') +
  geom_point(aes(colour = first.gen)) +
  scale_colour_viridis_c() +
  facet_wrap(~ gen)


set.seed(8008)
test.sim = sim(pars %>% filter(p.0 %in% 0.25) %>% mutate(timesteps = 25), 
               theta.t = 3, 
               init.rows = 7 * 1e3)

test.sim %>% group_by(gen) %>% summarise(n = n()) %>% 
  ggplot(aes(x = gen, y = n)) +
  geom_line()

test.sim %>% group_by(gen) %>% summarise(zbar = mean(z_i)) %>% 
  ggplot(aes(x = gen, y = zbar)) +
  geom_line()

test.sim %>% group_by(gen) %>% summarise(zvar = var(z_i)) %>% 
  ggplot(aes(x = gen, y = zvar)) +
  geom_line()

test.sim %>%
  mutate(p.0 = 0.25) %>%
  merge(y = pars %>% filter(p.0 %in% 0.25) %>% select(p.0, s.max, r, gsqrd, wfitn)) %>%
  group_by(gen) %>%
  summarise(
    sbar = mean(s_i),
    expl = (1+r)*sbar,
    exps = s.max * sqrt(wfitn^2 / (wfitn^2 + var(z_i)))
  ) %>% 
  distinct(gen, .keep_all = TRUE) %>%
  as.data.frame()

#----

set.seed(21)

sim.out = mclapply(
  pars %>% uncount(10) %>% mutate(try.no = 1:(nrow(.))) %>% split(.$try.no),
  function(pars) {
    sim(pars %>% mutate(timesteps = 24), theta.t = 3, init.rows = 5 * 1e5) %>%
      group_by(gen) %>%
      summarise(
        n = n(),
        bbar = mean(b_i),
        bvar = var(b_i)
      ) %>%
      mutate(trial = pars$try.no, s.max = pars$s.max)
  },
  mc.cores = 6
) %>%
  do.call(rbind, .)

sim.out %>%
  mutate(s.max = factor(s.max)) %>%
  ggplot(aes(x = gen, y = n, colour = s.max, group = trial)) +
  geom_line() +
  scale_y_log10() +
  scale_colour_manual(values = c('purple', 'orange'))

sim.out %>%
  mutate(s.max = factor(s.max)) %>%
  ggplot(aes(x = gen, y = bbar, colour = s.max, group = trial)) +
  geom_line() +
  scale_colour_manual(values = c('purple', 'orange'))

sim.out %>%
  mutate(s.max = factor(s.max)) %>%
  ggplot(aes(x = gen, y = bvar, colour = s.max, group = trial)) +
  geom_line() +
  scale_colour_manual(values = c('purple', 'orange'))

sim.sum = sim.out %>%
  group_by(gen, s.max) %>%
  summarise(
    nbar = mean(n),
    nvar = var(n),
    bbarbar = mean(bbar),
    bbarvar = var(bbar),
    bvarbar = mean(bvar),
    bvarvar = var(bvar),
    nn = n()
  )

sim.sum %>%
  mutate(s.max = factor(s.max)) %>%
  ggplot(aes(x = gen, y = nbar)) +
  geom_line(aes(colour = s.max)) +
  geom_ribbon(
    aes(
      ymin = nbar - 2*sqrt(nvar / nn),
      ymax = nbar + 2*sqrt(nvar / nn),
      fill = s.max
    ),
    alpha = 0.1
  )

sim.sum %>%
  mutate(s.max = factor(s.max)) %>%
  ggplot(aes(x = gen, y = bbarbar)) +
  geom_line(aes(colour = s.max)) +
  geom_ribbon(
    aes(
      ymin = bbarbar - 2*sqrt(bbarvar / nn),
      ymax = bbarbar + 2*sqrt(bbarvar / nn),
      fill = s.max
    ),
    alpha = 0.1
  )

sim.sum %>%
  mutate(s.max = factor(s.max)) %>%
  ggplot(aes(x = gen, y = bvarbar)) +
  geom_line(aes(colour = s.max)) +
  geom_ribbon(
    aes(
      ymin = bvarbar - 2*sqrt(bvarvar / nn),
      ymax = bvarbar + 2*sqrt(bvarvar / nn),
      fill = s.max
    ),
    alpha = 0.1
  )

# genetic variance doesn't seem like it's changing... so why is there a different equilibrium now??
# hmm... fuck... wonder if it has something to do with those mutation rates being fucked up?
# mutation rates being independent of the population size... hmm... maybe should be more probabalistic?
