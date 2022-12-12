# script for probing age structure including proposed new pmf for stable age
# dist'n
# SN Dec 2022

library(ggplot2)
library(dplyr)
library(tidyr)
library(parallel)

rm(list = ls())

######---------------------------------------
### Load in functions

source('model_source/sim_model1_functions.R')

### Other aux functions

# f we are solving in Newton's method
# (loss-function )
f_fun = function(x, w, s, r) {
  # NOTE: s here is the *target* population variance (or, here, s.d.),
  # x that is returned is the new cohort variance
  sumo = 0
  for (k in 0:1000) {
    sumo = sumo + 
      (1/(1+r))^k * (w^2 * x^2) / (w^2 + k * x^2)
  }
  sumo = sumo * (r / (1 + r))
  return(sumo - s^2)
  # NOTE sumo aka x is sigma, NOT sigma^2
}

# derivative of f wrt x
f_prm = function(x, w, r) {
  sumo = 0
  for (k in 0:1000) {
    nasty.quot = (2*x*w^2 * (w^2 + k*x^2) - (x^2 * w^2) * (2*k*x)) / (w^2 + k*x^2)^2
    sumo = sumo + ((1/(1+r))^k * nasty.quot)
  }
  sumo = sumo * (r / (1 + r))
  return(sumo)
}

# Newton's method wrapper for solving for initial additive genetic variance
newt.sigma.a0 = function(x0, tol, w, s, r) {
  
  xold = x0
  
  while(abs(f_fun(xold, w, s, r)) > tol) xold = xold - (f_fun(xold, w, s, r) / f_prm(xold, w, r))
  
  return(xold)
  
}

######---------------------------------------
### Defining parameters for init sims

ntrials = 10000
l.max = 1.1

pars.list = expand.grid(s.max = c(0.5, 0.9)) %>%
  mutate(
    n.pop0 = 1000,
    sig.a.pop = 1,
    sig.e = 0, alpha = 0.0000,
    kceil = 2000,
    timesteps = 1,
    mu = 0, sig.m = 0
  ) %>%
  # Get selection pressure and fecundity
  mutate(
    wfitn = sqrt((sig.a.pop^2 + sig.e^2) / (l.max^2 - 1)),
    r = (l.max / s.max) - 1
  ) %>%
  # Need this to for vectorizing functions listed above
  group_by(s.max) %>%
  mutate(
    sig.a = newt.sigma.a0(1, 1e-6, wfitn, sig.a.pop, r)
  ) %>%
  ungroup()

######---------------------------------------
### Run initial batch of sims

set.seed(5418)

sim.out = mclapply(
  pars.list %>% uncount(ntrials) %>% mutate(try.no = 1:(nrow(.))) %>% split(.$try.no),
  function(pars) {
    sim(pars, theta.t = 0, init.rows = 3 * 1e3) %>%
      group_by(gen) %>%
      summarise(
        n = n(),
        bvar = var(b_i)
      ) %>%
      mutate(trial = pars$try.no, s.max = pars$s.max)
  },
  mc.cores = 8
) %>%
  do.call(rbind, .)

sim.sum = sim.out %>%
  group_by(s.max, gen) %>%
  summarise(
    nbar = mean(n),
    nvar = var(n),
    bvarbar = mean(bvar),
    bvarvar = var(bvar),
    nn = n()
  )

head(sim.sum)

sim.sum %>% filter(gen > 0) %>% select(s.max, bvarbar) %>% as.data.frame()

######---------------------------------------
### RNow run with values above used for mutation

ntrials = 5000
l.max = 1.1

parlist = expand.grid(s.max = c(0.5, 0.9)) %>%
  mutate(
    n.pop0 = 1000,
    sig.a.pop = 1,
    sig.e = 0, alpha = 0.0000,
    kceil = 2000,
    timesteps = 2
  ) %>%
  # Get selection pressure and fecundity
  mutate(
    wfitn = sqrt((sig.a.pop^2 + sig.e^2) / (l.max^2 - 1)),
    r = (l.max / s.max) - 1
  ) %>%
  # Need this to for vectorizing functions listed above
  group_by(s.max) %>%
  mutate(
    sig.a = newt.sigma.a0(1, 1e-6, wfitn, sig.a.pop, r)
  ) %>%
  ungroup() %>%
  # Get variance lost
  merge(sim.sum %>% filter(gen > 0) %>% select(s.max, bvarbar)) %>%
  mutate(
    mu = 1, 
    sig.m = sqrt( ((1 + r) / r) * (sig.a.pop^2 - bvarbar) )
  )

with(parlist, (bvarbar / (1 + r)) + (r / (1+r)) * (bvarbar + sig.m^2))

set.seed(260)

jim.out = mclapply(
  parlist %>% uncount(ntrials) %>% mutate(try.no = 1:(nrow(.))) %>% split(.$try.no),
  function(pars) {
    sim(pars, theta.t = 0, init.rows = 3 * 1e3) %>%
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

jim.sum = jim.out %>%
  group_by(s.max, gen) %>%
  summarise(
    nbar = mean(n),
    nvar = var(n),
    bvarbar = mean(bvar),
    bvarvar = var(bvar),
    nn = n()
  )

jim.sum %>%
  mutate(s.max = factor(s.max)) %>%
  ggplot(aes(x = gen, group = s.max)) +
  geom_line(aes(y = bvarbar, colour = s.max)) +
  geom_ribbon(
    aes(
      ymin = bvarbar - 2 * sqrt(bvarvar / nn),
      ymax = bvarbar + 2 * sqrt(bvarvar / nn),
      fill = s.max
    ),
    alpha = 0.25
  )
# yep...  okay so something is still not quite right with the longer-gen sims

### So what's going on...?

# Get some individual trials

set.seed(260)

tim.out = mclapply(
  parlist %>% uncount(200) %>% mutate(try.no = 1:(nrow(.))) %>% split(.$try.no),
  function(pars) {
    sim(pars %>% mutate(timesteps = 10), theta.t = 0, init.rows = 4 * 1e3) %>%
      mutate(
        trial = pars$try.no, 
        s.max = pars$s.max
      )
  },
  mc.cores = 6
) %>%
  do.call(rbind, .)

tim.sum = tim.out %>%
  group_by(trial, s.max, gen, age) %>%
  mutate(
    bbar = mean(b_i),
    bvar = var(b_i),
    n = n()
  ) %>%
  group_by(trial, s.max, gen) %>%
  mutate(
    bvar.all = var(b_i),
    p = n / sum(n)
  ) %>%
  ungroup() %>%
  select(
    trial, gen, age, s.max,
    bbar, bvar, n, p, bvar.all
  ) %>%
  distinct(trial, gen, age, .keep_all = TRUE)

tim.sum %>%
  filter((trial %% 100) < 10) %>%
  mutate(cohort = gen - age) %>%
  ggplot(aes(x = gen, y = bbar, group = cohort)) +
  geom_line() +
  geom_point(aes(colour = age)) +
  scale_colour_viridis_c() +
  facet_wrap(~ trial)
# not a super helpful plot...

tim.sum %>%
  group_by(gen, s.max) %>%
  summarise(
    bbarbar = sum(p * bbar),
    bbarvar = var(bbar),
    wbarvar = sum(p * (bbar - bbarbar)^2),
    nn = n()
  ) %>%
  ggplot(aes(x = gen, y = wbarvar, group = s.max)) +
  geom_line(aes(colour = s.max))
# well... they both seem to increase!!! makes sense...
# shit.

tim.sum %>%
  group_by(trial, gen, s.max) %>%
  summarise(
    bbarbar = sum(p * bbar),
    bbarvar = var(bbar),
    wbarvar = sum(p * (bbar - bbarbar)^2),
    n = n()
  ) %>%
  group_by(gen, s.max) %>%
  summarise(
    wbvbar = mean(wbarvar),
    wbvvar = var(wbarvar),
    nn = n()
  ) %>%
  ggplot(aes(x = gen, y = wbvbar, group = s.max)) +
  geom_line(aes(colour = s.max)) +
  geom_ribbon(
    aes(
      ymin = wbvbar - 2 * sqrt(wbvvar / nn),
      ymax = wbvbar + 2 * sqrt(wbvvar / nn),
      fill = s.max
    ),
    alpha = 0.1
  )

tim.out %>%
  group_by(trial, s.max, gen) %>%
  summarise(bvar = var(b_i)) %>%
  group_by(s.max = factor(s.max), gen) %>%
  summarise(
    bvarbar = mean(bvar),
    bvarvar = var(bvar),
    n = n()
  ) %>%
  ggplot(aes(x = gen, group = s.max)) +
  geom_line(aes(y = bvarbar, colour = s.max)) +
  geom_ribbon(
    aes(
      ymin = bvarbar - 2 * sqrt(bvarvar / n),
      ymax = bvarbar + 2 * sqrt(bvarvar / n),
      fill = s.max
    ),
    alpha = 0.2
  )
# yeah... seeing decline for the longer-lived ones...
# what the heck man

# so there is variance among population means
# but more variance for shorter-lived populations... (sampling variance?)

# how does this fit in with more variance loss in longer-lived populations?
# do these means even matter?

######---------------------------------------
### Okay well now try to brute force?

ntrials = 1000

parlist = expand.grid(s.max = c(0.5, 0.9), badj = ((0:5)*2)/1000) %>%
  mutate(
    n.pop0 = 1000,
    sig.a.pop = 1,
    sig.e = 0, alpha = 0.0000,
    kceil = 1200,
    timesteps = 3
  ) %>%
  # Get selection pressure and fecundity
  mutate(
    wfitn = sqrt((sig.a.pop^2 + sig.e^2) / (l.max^2 - 1)),
    r = (l.max / s.max) - 1
  ) %>%
  # Get variance lost
  merge(sim.sum %>% filter(gen > 0) %>% select(s.max, bvarbar)) %>%
  group_by(s.max, badj) %>%
  mutate(
    sig.a = newt.sigma.a0(1, 1e-6, wfitn, sig.a.pop, r)
  ) %>%
  mutate(
    mu = 1, 
    sig.m = sqrt( ((1 + r) / r) * (sig.a.pop^2 - bvarbar) + badj)
  ) %>%
  ungroup()

with(parlist, (bvarbar / (1 + r)) + (r / (1+r)) * (bvarbar + sig.m^2))

set.seed(260)

rim.out = mclapply(
  parlist %>% uncount(ntrials) %>% mutate(try.no = 1:(nrow(.))) %>% split(.$try.no),
  function(pars) {
    sim(pars, theta.t = 0, init.rows = 5 * 1e3) %>%
      group_by(gen) %>%
      summarise(
        n = n(),
        bvar = var(b_i)
      ) %>%
      mutate(trial = pars$try.no, s.max = pars$s.max, 
             sig.m = pars$sig.m,  badj  = pars$badj)
  },
  mc.cores = 6
) %>%
  do.call(rbind, .)

rim.sum = rim.out %>%
  group_by(gen, s.max, sig.m, badj) %>%
  summarise(
    bvarbar = mean(bvar),
    bvarvar = var(bvar),
    nn = n()
  )

rim.sum %>%
  ggplot(aes(x = gen, group = badj, fill = badj)) +
  geom_line(aes(y = bvarbar, colour = badj), size = 1.1) +
  geom_ribbon(
    aes(
      ymin = bvarbar - 2 * sqrt(bvarvar / nn),
      ymax = bvarbar + 2 * sqrt(bvarvar / nn)
    ),
    alpha = 0.1
  ) +
  scale_colour_viridis_c() +
  scale_fill_viridis_c() +
  facet_wrap(~ s.max)

ntrials = 2000

parlist = expand.grid(badj = ((0:5)*4)/1000) %>%
  mutate(
    n.pop0 = 1000,
    sig.a.pop = 1,
    sig.e = 0, alpha = 0.0000,
    kceil = 1200,
    timesteps = 3,
    s.max = 0.9
  ) %>%
  # Get selection pressure and fecundity
  mutate(
    wfitn = sqrt((sig.a.pop^2 + sig.e^2) / (l.max^2 - 1)),
    r = (l.max / s.max) - 1
  ) %>%
  # Get variance lost
  merge(sim.sum %>% filter(gen > 0) %>% select(s.max, bvarbar)) %>%
  group_by(s.max, badj) %>%
  mutate(
    sig.a = newt.sigma.a0(1, 1e-6, wfitn, sig.a.pop, r)
  ) %>%
  mutate(
    mu = 1, 
    sig.m = sqrt( ((1 + r) / r) * (sig.a.pop^2 - bvarbar) + badj)
  ) %>%
  ungroup()

set.seed(260)

rim.out = mclapply(
  parlist %>% uncount(ntrials) %>% mutate(try.no = 1:(nrow(.))) %>% split(.$try.no),
  function(pars) {
    sim(pars, theta.t = 0, init.rows = 5 * 1e3) %>%
      group_by(gen) %>%
      summarise(
        n = n(),
        bvar = var(b_i)
      ) %>%
      mutate(trial = pars$try.no, sig.m = pars$sig.m,  badj  = pars$badj)
  },
  mc.cores = 6
) %>%
  do.call(rbind, .)

rim.sum = rim.out %>%
  group_by(gen, sig.m, badj) %>%
  summarise(
    bvarbar = mean(bvar),
    bvarvar = var(bvar),
    nn = n()
  )

rim.sum %>%
  ggplot(aes(x = gen, group = badj, fill = badj)) +
  geom_line(aes(y = bvarbar, colour = badj), size = 1.1) +
  geom_ribbon(
    aes(
      ymin = bvarbar - 2 * sqrt(bvarvar / nn),
      ymax = bvarbar + 2 * sqrt(bvarvar / nn)
    ),
    alpha = 0.1
  ) +
  scale_colour_viridis_c() +
  scale_fill_viridis_c() 

parlist = expand.grid(badj = ((0:5)*4)/100) %>%
  mutate(
    n.pop0 = 1000,
    sig.a.pop = 1,
    sig.e = 0, alpha = 0.0000,
    kceil = 2000,
    timesteps = 10,
    s.max = 0.9
  ) %>%
  # Get selection pressure and fecundity
  mutate(
    wfitn = sqrt((sig.a.pop^2 + sig.e^2) / (l.max^2 - 1)),
    r = (l.max / s.max) - 1
  ) %>%
  # Get variance lost
  merge(sim.sum %>% filter(gen > 0) %>% select(s.max, bvarbar)) %>%
  group_by(s.max, badj) %>%
  mutate(
    sig.a = newt.sigma.a0(1, 1e-6, wfitn, sig.a.pop, r)
  ) %>%
  mutate(
    mu = 1, 
    sig.m = sqrt( ((1 + r) / r) * (sig.a.pop^2 - bvarbar) + badj)
  ) %>%
  ungroup()

set.seed(260)

rim.out = mclapply(
  parlist %>% uncount(ntrials) %>% mutate(try.no = 1:(nrow(.))) %>% split(.$try.no),
  function(pars) {
    sim(pars, theta.t = 0, init.rows = 5 * 1e4) %>%
      group_by(gen) %>%
      summarise(
        n = n(),
        bvar = var(b_i)
      ) %>%
      mutate(trial = pars$try.no, sig.m = pars$sig.m,  badj  = pars$badj)
  },
  mc.cores = 6
) %>%
  do.call(rbind, .)

rim.sum = rim.out %>%
  group_by(gen, sig.m, badj) %>%
  summarise(
    bvarbar = mean(bvar),
    bvarvar = var(bvar),
    nn = n()
  )

rim.sum %>%
  ggplot(aes(x = gen, group = badj, fill = badj)) +
  geom_line(aes(y = bvarbar, colour = badj), size = 1.1) +
  geom_ribbon(
    aes(
      ymin = bvarbar - 2 * sqrt(bvarvar / nn),
      ymax = bvarbar + 2 * sqrt(bvarvar / nn)
    ),
    alpha = 0.1
  ) +
  scale_colour_viridis_c() +
  scale_fill_viridis_c() 

# interesting it *does* look like pop sizestabilizes...?
# but what the fuck is up with this initial change? something transient is happening
# fucking christ

set.seed(2090)

rim.out = mclapply(
  parlist %>% slice(2:3) %>% uncount(ntrials) %>% mutate(try.no = 1:(nrow(.))) %>% split(.$try.no),
  function(pars) {
    sim(pars %>% mutate(timesteps = 24), theta.t = 0, init.rows = 5 * 1e4) %>%
      group_by(gen) %>%
      summarise(
        n = n(),
        bvar = var(b_i)
      ) %>%
      mutate(trial = pars$try.no, sig.m = pars$sig.m,  badj  = pars$badj)
  },
  mc.cores = 6
) %>%
  do.call(rbind, .)

rim.sum = rim.out %>%
  group_by(gen, sig.m, badj) %>%
  summarise(
    bvarbar = mean(bvar),
    bvarvar = var(bvar),
    nn = n()
  )

rim.sum %>%
  ggplot(aes(x = gen, group = badj, fill = badj)) +
  geom_line(aes(y = bvarbar, colour = badj), size = 1.1) +
  geom_ribbon(
    aes(
      ymin = bvarbar - 2 * sqrt(bvarvar / nn),
      ymax = bvarbar + 2 * sqrt(bvarvar / nn)
    ),
    alpha = 0.1
  ) +
  scale_colour_viridis_c() +
  scale_fill_viridis_c() 
# ugh... fucking jesus
# what the fuck is going on here, fr

# look at a few of these with the big adj to see what happens

set.seed(44)

pim.out = mclapply(
  parlist %>% slice(3) %>% uncount(50) %>% mutate(try.no = 1:(nrow(.))) %>% split(.$try.no),
  function(pars) {
    sim(pars %>% mutate(timesteps = 15), theta.t = 0, init.rows = 5 * 1e4) %>%
      mutate(trial = pars$try.no, sig.m = pars$sig.m,  badj  = pars$badj)
  },
  mc.cores = 6
) %>%
  do.call(rbind, .)

nrow(pim.out)

pim.out %>%
  group_by(trial, gen, sig.m) %>%
  summarise(bvar = var(b_i)) %>%
  group_by(gen, sig.m) %>%
  summarise(
    bvarbar = mean(bvar),
    bvarvar = var(bvar),
    nn = n()
  ) %>%
  ggplot(aes(x = gen, y = bvarbar)) +
  geom_line()

# 

pim.out %>%
  group_by(trial, gen, sig.m) %>%
  summarise(bvar = var(b_i)) %>%
  ggplot(aes(x = gen, y = bvar, group = trial)) +
  geom_line(size = 0.1)

pim.out %>%
  group_by(trial, gen, sig.m) %>%
  summarise(bvar = var(b_i)) %>%
  ungroup() %>%
  filter(gen %in% 1) %>%
  slice_max(bvar)

pim28 = pim.out %>% filter(trial %in% 28)

pim28 %>%
  ggplot(aes(x = b_i, y = age)) +
  geom_point(position = position_jitter(height = 0.2), alpha = 0.5) +
  facet_wrap(~ gen)

pim28 %>%
  ggplot(aes(x = b_i, group = gen)) +
  geom_freqpoly(aes(colour = gen)) +
  scale_colour_viridis_c()

pim28 %>%
  ggplot(aes(x = b_i, group = interaction(gen, age))) +
  geom_freqpoly(aes(colour = age)) +
  scale_colour_viridis_c() +
  facet_wrap(~ gen)

pim28 %>%
  filter(gen < 4, age < 25) %>%
  ggplot(aes(x = b_i, group = interaction(gen, age))) +
  geom_freqpoly(aes(colour = gen), binwidth = 0.5) +
  scale_colour_viridis_c() +
  facet_wrap(~ age)

pim28 %>%
  filter(gen < 4) %>%
  mutate(cohort = gen - age) %>%
  ggplot(aes(x = b_i, group = interaction(gen, cohort))) +
  geom_freqpoly(aes(colour = gen), binwidth = 0.5) +
  scale_colour_viridis_c() +
  facet_wrap(~ cohort)

pim.out %>%
  filter(gen < 2) %>%
  mutate(isnew = gen %in% 1 & !age) %>%
  group_by(trial, gen, isnew) %>%
  summarise(bvar = var(b_i), n = n()) %>%
  ggplot(aes(x = gen, y = bvar, group = interaction(trial, isnew), colour = isnew)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = c('black', 'blue'))

muto.fun = function(w2, s2, r) {
  sumo = 0
  for (k in 0:1000) {
    sumo = sumo +
      (1 / (1+r))^k * ( (w2 + k*s2)^(-1) - (w2 + (k+1)*s2)^(-1) )
  }
  sumo = sumo * s2 * w2
  return(sumo)
}

p1out = pim.out %>%
  filter(trial < 2) %>%
  group_by(gen, age) %>%
  summarise(
    bvar = var(b_i),
    n = n()
  ) %>%
  group_by(gen) %>%
  mutate(p = n / sum(n)) %>%
  ungroup()

pout = pim.out %>%
  group_by(trial, gen, age) %>%
  summarise(
    bvar = var(b_i),
    n = n()
  ) %>%
  mutate(bvar = ifelse(is.na(bvar), 0, bvar)) %>%
  group_by(trial, gen) %>%
  mutate(p = n / sum(n)) %>%
  ungroup() %>%
  group_by(gen, age) %>%
  summarise(
    bvar = mean(bvar),
    p = mean(p)
  )

pout %>%
  #filter(gen < 3) %>%
  ggplot(aes(x = p, y = bvar, group = age)) +
  geom_path() +
  geom_point(aes(fill = gen), shape = 21, size = 3) +
  scale_colour_viridis_c(option = 'C') +
  scale_fill_viridis_c(option = 'C')
# wait... here, p looks like it shifts a lot more than v! v maybe shifts a wee bit...

pout %>%
  ggplot(aes(x = age, y = p, group = gen, colour = gen)) +
  geom_line() +
  geom_point() +
  geom_line(
    data = data.frame(
      age = 0:60,
      p = with(parlist[1,], (r / (1+r)) * (1 / (1+r))^(0:60))
    ),
    aes(x = age, y = p),
    inherit.aes = FALSE,
    linetype = 2, colour = 'red'
  )
# oh that must be why... I'm not actually at SSD shoot.

pout %>%
  ggplot(aes(x = age, y = p, group = gen, colour = gen)) +
  geom_line() +
  geom_point() +
  geom_line(
    data = data.frame(
      age = 0:60,
      p = with(parlist[1,], (r / (1+r)) * (1 / (1+r))^(0:60))
    ),
    aes(x = age, y = p),
    inherit.aes = FALSE,
    linetype = 2, colour = 'red'
  ) +
  scale_y_log10()

##############-----------------------------------------------
# am I off of SSD because of the mutations

set.seed(5418)

sim.out = mclapply(
  parlist %>% slice(1) %>% uncount(50) %>% mutate(try.no = 1:(nrow(.))) %>% split(.$try.no),
  function(pars) {
    sim(pars %>% mutate(timesteps = 24), theta.t = 0, init.rows = 5 * 1e4) %>%
      mutate(trial = pars$try.no)
  },
  mc.cores = 8
) %>%
  do.call(rbind, .)

sout = sim.out %>%
  group_by(trial, gen, age) %>%
  summarise(
    bvar = var(b_i),
    n = n()
  ) %>%
  mutate(bvar = ifelse(is.na(bvar), 0, bvar)) %>%
  group_by(trial, gen) %>%
  mutate(p = n / sum(n)) %>%
  ungroup() %>%
  group_by(gen, age) %>%
  summarise(
    bvar = mean(bvar),
    p = mean(p)
  )

sout %>% 
  ggplot(aes(x = age, y = p, group = gen, colour = gen)) +
  # geom_line() +
  geom_point() +
  geom_line(
    data = data.frame(
      age = 0:60,
      p = with(parlist[1,], (r / (1+r)) * (1 / (1+r))^(0:60))
    ),
    aes(x = age, y = p),
    inherit.aes = FALSE,
    linetype = 2, colour = 'red'
  ) 

# oh... I guess we're *not* at SSD.. fuck

set.seed(5418)

mim.out = mclapply(
  parlist %>% slice(1) %>% uncount(50) %>% mutate(try.no = 1:(nrow(.))) %>% split(.$try.no),
  function(pars) {
    sim(pars %>% mutate(timesteps = 19, sig.m = 0), theta.t = 0, init.rows = 5 * 1e4) %>%
      mutate(trial = pars$try.no)
  },
  mc.cores = 8
) %>%
  do.call(rbind, .)

mout = mim.out %>%
  group_by(trial, gen, age) %>%
  summarise(
    bvar = var(b_i),
    n = n()
  ) %>%
  mutate(bvar = ifelse(is.na(bvar), 0, bvar)) %>%
  group_by(trial, gen) %>%
  mutate(p = n / sum(n)) %>%
  ungroup() %>%
  group_by(gen, age) %>%
  summarise(
    bvar = mean(bvar),
    p = mean(p)
  )

mout %>% 
  ggplot(aes(x = age, y = p, group = gen, colour = gen > 5)) +
  # geom_line() +
  geom_point() +
  geom_line(
    data = data.frame(
      age = 0:60,
      p = with(parlist[1,], (r / (1+r)) * (1 / (1+r))^(0:60))
    ),
    aes(x = age, y = p),
    inherit.aes = FALSE,
    linetype = 2, colour = 'red'
  )

# so... I'm not at SSD...?'
# I guess the mutation rate fs with this... shit shit fuck fuck

with(list(k = 0:30, s = .9, w = 2.18, sig.a.0 = 2.5, sig.a = 1, r = 1.1/.9 - 1),
     data.frame(
       age = k,
       pk1 = s^k * sqrt(w^2 / (w^2 + (k+1)*sig.a.0^2)),
       pk2 = (s/1.1)^k * (1 - s/1.1),
       pk3 = sqrt(w^2 / (w^2 + (k+1)*sig.a.0^2)) / ((1+r)*sqrt(w^2 / (w^2 + sig.a^2)))^k
       )
     ) %>%
  mutate(pk1 = pk1 / sum(pk1)) %>%
  mutate(pk3 = pk3 / sum(pk3)) %>%
  ggplot(aes(x = age)) +
  geom_line(aes(y = pk1), colour = 'blue') +
  geom_line(aes(y = pk2), colour = 'red') +
  geom_line(aes(y = pk3), colour = 'green')

with(list(k = 0:100, r = parlist$r[1], s = parlist$s.max[1], 
          sig.a = parlist$sig.a[1], sig.a.pop = parlist$sig.a.pop[1],
          w = parlist$wfitn[1]),
     data.frame(
       age = k,
       cumsurv = sqrt(w^2 / (w^2 + k*sig.a^2)) / ((1+r)*sqrt(w^2 / (w^2 + sig.a.pop^2)))^k,
       psurvk  = sqrt((w^2 + k*sig.a^2) / (w^2 + (k+1)*sig.a^2)) / ((1+r)*sqrt(w^2 / (w^2 + sig.a.pop^2)))
     ) %>%
       mutate(
         pk = cumsurv * (1 - psurvk)
       )
    ) %>%
  summarise(haha = sum(pk))


try.probs = with(
  list(k = 0:100, r = parlist$r[1], s = parlist$s.max[1], 
       sig.a = parlist$sig.a[1], sig.a.pop = parlist$sig.a.pop[1],
       w = parlist$wfitn[1]),
     data.frame(
       age = k,
       cumsurv = sqrt(w^2 / (w^2 + k*sig.a^2)) / ((1+r)*sqrt(w^2 / (w^2 + sig.a.pop^2)))^k,
       psurvk  = sqrt((w^2 + k*sig.a^2) / (w^2 + (k+1)*sig.a^2)) / ((1+r)*sqrt(w^2 / (w^2 + sig.a.pop^2)))
     ) %>%
       mutate(
         pk = cumsurv * (1 - psurvk)
       )
)

mout %>% 
  #filter(gen > 0) %>%
  ggplot(aes(x = age, y = p, group = gen, colour = gen)) +
  #geom_line() +
  geom_point() +
  geom_line(
    data = data.frame(
      age = 0:60,
      p = with(parlist[1,], (r / (1+r)) * (1 / (1+r))^(0:60))
    ),
    aes(x = age, y = p),
    inherit.aes = FALSE,
    linetype = 2, colour = 'red'
  ) +
  # geom_line(
  #   data = try.probs %>% filter(age < 40),
  #   aes(x = age, y = pk),
  #   inherit.aes = FALSE,
  #   linetype = 2, colour = 'forestgreen'
  # ) +
  # scale_colour_viridis_c() +
  theme(legend.position = 'bottom')

# okay so my expression looks like it's quite wrong
# (hard to know because sig.a is calculated assuming geometric though)
# when mutations are zero it looks like SSD actually is geometric
# when there are mutations, there's a divergence from this structure
# (although my expression is still not right)

cohort0.timeseries = mim.out %>%
  group_by(trial, gen) %>%
  mutate(n.overall = n()) %>%
  ungroup() %>%
  filter(gen == age) %>%
  group_by(trial, gen) %>%
  summarise(n = n(), n.overall = n.overall[1])

cohort0.timeseries %>%
  ggplot(aes(x = gen, y = n, group = trial)) +
  geom_line() +
  geom_line(
    data = try.probs %>% mutate(expec = cumsurv * 175) %>% filter(age < 19),
    aes(x = age, y = expec),
    inherit.aes = FALSE,
    colour = 'blue'
  )
# that looks pretty good to me...

cohort0.timeseries %>%
  group_by(trial) %>%
  mutate(n.norm = n / max(n)) %>%
  ungroup() %>%
  ggplot(aes(x = gen, y = n.norm, group = trial)) +
  geom_line() +
  geom_line(
    data = try.probs %>% filter(age < 20),
    aes(x = age, y = cumsurv),
    inherit.aes = FALSE,
    colour = 'blue'
  )
# honestly looks great up until time step ~18
# (although this could be due to the ceiling!)

cohort0.timeseries %>%
  group_by(trial) %>%
  mutate(n.dead = c(diff(n), NA),
         p.dead = -n.dead / n) %>%
  ungroup() %>%
  ggplot(aes(x = gen, y = p.dead, group = trial)) +
  geom_point(position = position_jitter(width = 0.1, height = 0.1)) +
  geom_line(
    data = try.probs %>% filter(age < 20),
    aes(x = age, y = 1 - psurvk),
    inherit.aes = FALSE,
    colour = 'blue'
  )
# (oh yeah... could also be the ceiling, lol)

# I think there's an issue with the assumption that you can just multiply by p(mortality)
# that gives us an initial size of ~0.25 for cohort 1
# whereas I'm seeing, consistently, proportions close to ~.18...

cohort0.timeseries %>%
  mutate(p = n / n.overall) %>%
  ggplot(aes(x = gen, y = p, group = trial)) +
  geom_line()

