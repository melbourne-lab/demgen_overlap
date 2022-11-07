# Implement scripts that might (might!) get proper equilibrium variance
# in a stable environment

library(parallel)
library(ggplot2)
library(dplyr)
library(tidyr)

rm(list = ls())

#####----------------------------------------
# Load in functions

# Load in simulation functions
source('model_source/sim_model1_functions.R')

# Wrapper for determining the mutation rate (given s, w, r)
muto.fun = function(w, s, r) {
  # NOTE: s is the *new cohort variance*, not the population-level variance (s.d.)
  sumo = 0
  for (k in 0:1000) {
    sumo = sumo +
      (1 / (1+r))^k * ( (w^2 + k*s^2)^(-1) - (w^2 + (k+1)*s^2)^(-1) )
  }
  sumo = sumo * s^2 * w^2
  return(sumo)
  # sumo is sigma_m^2, not sigma_m
}

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

#####----------------------------------------
# Setup simulations

ntrials = 100
l.max = 1.1

pars.list = data.frame(
  n.pop0 = 1000, s.max = 0.9,
  wfitn = sqrt(10),
  sig.a.pop = 1,
  sig.e = 0, alpha = 0.0000,
  kceil = 2000,
  timesteps = 40
) %>%
  # Get selection pressure and fecundity
  mutate(
    # wfitn = sqrt((sig.a.pop^2 + sig.e^2) / (l.max^2 - 1)),
    r = (l.max / s.max) - 1
  ) %>%
  mutate(
    sig.a = newt.sigma.a0(1, 1e-6, wfitn, sig.a.pop, r)
  ) %>%
  mutate(
    mu = 1, 
    sig.m = sqrt(muto.fun(wfitn, sig.a, r))
  )

set.seed(5418)

out1 = mclapply(
  pars.list %>% uncount(ntrials) %>% mutate(try.no = 1:ntrials) %>% split(.$try.no),
  function(pars) {
    sim(params = pars, theta.t = 0, init.row = 5 * 1000 * 10) %>%
      group_by(gen) %>%
      summarise(n = n(), bvar = var(b_i)) %>%
      mutate(trial = pars$try.no)
  },
  mc.cores = 4
) %>%
  do.call(rbind, .)

sum1 = out1 %>%
  group_by(gen) %>%
  summarise(
    nbar = mean(n),
    nvar = var(n),
    bvarbar = mean(bvar),
    bvarvar = var(bvar),
    nn = n()
  )

sum1 %>%
  ggplot(aes(x = gen, y = nbar)) +
  geom_line() +
  scale_y_log10()

sum1 %>%
  ggplot(aes(x = gen, y = bvarbar)) +
  geom_line() +
  geom_ribbon(
    aes(
      ymin = bvarbar - 2 * sqrt(bvarvar / nn),
      ymax = bvarbar - 2 * sqrt(bvarvar / nn)
    ),
    alpha = 0.5
  )

# argghhhhhhhhhh so close and yet still not there... what's missing? jeez

# needs work but much closer than ever before...
# (variance is still the problem here... declining a wee bit still, somehow)

######--------------------------------

pars.list = data.frame(
  n.pop0 = 1000, s.max = 0.9,
  # wfitn = sqrt(10),
  sig.a.pop = 1,
  sig.e = 0, alpha = 0.0000,
  kceil = 2000,
  timesteps = 40
) %>%
  # Get selection pressure and fecundity
  mutate(
    wfitn = sqrt((sig.a.pop^2 + sig.e^2) / (l.max^2 - 1)),
    r = (l.max / s.max) - 1
  ) %>%
  mutate(
    sig.a = newt.sigma.a0(1, 1e-6, wfitn, sig.a.pop, r)
  ) %>%
  mutate(
    mu = 1, 
    sig.m = sqrt(muto.fun(wfitn, sig.a, r))
  )

set.seed(5418)

out2 = mclapply(
  pars.list %>% uncount(ntrials) %>% mutate(try.no = 1:ntrials) %>% split(.$try.no),
  function(pars) {
    sim(params = pars, theta.t = 0, init.row = 5 * 1000 * 10) %>%
      group_by(gen) %>%
      summarise(n = n(), bvar = var(b_i)) %>%
      mutate(trial = pars$try.no)
  },
  mc.cores = 4
) %>%
  do.call(rbind, .)

sum2 = out2 %>%
  group_by(gen) %>%
  summarise(
    nbar = mean(n),
    nvar = var(n),
    bvarbar = mean(bvar),
    bvarvar = var(bvar),
    nn = n()
  )

sum2 %>%
  ggplot(aes(x = gen, y = nbar)) +
  geom_line() +
  scale_y_log10()

sum2 %>%
  ggplot(aes(x = gen, y = bvarbar)) +
  geom_line() +
  geom_ribbon(
    aes(
      ymin = bvarbar - 2 * sqrt(bvarvar / nn),
      ymax = bvarbar - 2 * sqrt(bvarvar / nn)
    ),
    alpha = 0.5
  )

exp(diff(log(sum2$nbar)))

with(pars.list, s.max * (1 + r) * sqrt(wfitn^2 / (wfitn^2 + sum2$bvarbar)))
# man...

out2 %>%
  ggplot(aes(x = gen, y = bvar, group = trial)) +
  geom_line(size = 0.2)
# yep.

########_-----------------------
#
#

set.seed(4012277)

test1 = sim(params = pars.list, theta.t = 0, init.row = 5 * 1000 * 10)

test1 %>%
  group_by(gen) %>%
  summarise(bvar = var(b_i) * (1 - 1/n())) %>%
  ggplot(aes(x = gen, y = bvar)) +
  geom_line()

test1 %>%
  group_by(gen) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = gen, y = n)) +
  geom_line()

test1 %>%
  filter(gen < 36) %>%
  group_by(i) %>%
  mutate(coh1 = ifelse(!min(gen), 0, gen - age)) %>%
  ggplot(aes(x = i, y = b_i)) +
  geom_point(aes(colour = coh1), size = 0.1) +
  scale_color_viridis_c() +
  # scale_color_manual(values = c('green', 'blue')) +
  facet_wrap(~ gen)

test1 = test1 %>%  
  group_by(i) %>%
  mutate(coh1 = ifelse(!min(gen), 0, gen - age)) %>%
  ungroup()

test1 %>%
  distinct(i, coh1, .keep_all = TRUE) %>%
  ggplot(aes(i, y = b_i)) +
  geom_point(aes(colour = coh1), size = 0.2) +
  scale_colour_viridis_c()

test1 %>%
  distinct(i, coh1, .keep_all = TRUE) %>%
  group_by(coh1) %>%
  summarise(bvar = var(b_i) * (1 - 1/n()),
            bvar.uncor = var(b_i)) %>%
  ggplot(aes(x = coh1, y = bvar.uncor)) +
  geom_line()

test1 %>%
  group_by(gen, new.old = coh1 > 0) %>%
  summarise(bvar = var(b_i)) %>%
  ggplot(aes(x = gen, y = bvar, group = new.old)) +
  geom_line(aes(colour = new.old))

test1 %>%
  filter(!gen) %>%
  ggplot(aes(x = age)) + 
  geom_histogram(binwidth = 1)

test1 %>%
  filter(!coh1) %>%
  ggplot(aes(x = age, group = gen, fill = gen)) + 
  geom_histogram(position = "identity", alpha = 0.2, binwidth = 1) +
  scale_fill_viridis_c()

test1 %>%
  filter(!coh1) %>%
  ggplot(aes(x = age, group = gen, colour = gen)) + 
  geom_freqpoly(position = "identity", binwidth = 1) +
  scale_colour_viridis_c()

test1 %>%
  filter(!coh1) %>%
  ggplot(aes(x = age, after_stat(density), group = gen, colour = gen)) + 
  geom_freqpoly(position = "identity", binwidth = 1) +
  scale_colour_viridis_c()

#############---------------------------------------

set.seed(5418)

out3 = mclapply(
  pars.list %>% uncount(ntrials) %>% mutate(try.no = 1:ntrials) %>% split(.$try.no),
  function(pars) {
    sim(params = pars, theta.t = 0, init.row = 5 * 1000 * 10) %>%
      group_by(i) %>%
      mutate(cohort = ifelse(!min(gen), 0, gen - age)) %>%
      group_by(gen, cohort) %>%
      summarise(n = n(), bvar = var(b_i)) %>%
      mutate(trial = pars$try.no)
  },
  mc.cores = 4
) %>%
  do.call(rbind, .)

out3 %>%
  group_by(gen, cohort) %>%
  summarise(
    bvarbar = mean(bvar),
    bvarvar = var(bvar),
    nn = n()
  ) %>%
  ggplot(aes(x = gen, y = bvarbar, group = cohort)) +
  geom_line(size = 0.1)

################----------------------------------

set.seed(29014)

out4 = mclapply(
  pars.list %>% uncount(ntrials) %>% mutate(try.no = 1:ntrials) %>% split(.$try.no),
  function(pars) {
    sim(params = pars, theta.t = 0, init.row = 5 * 1000 * 10) %>%
      group_by(i) %>%
      mutate(cohort = ifelse(!min(gen), 0, gen - age)) %>%
      group_by(gen, new.old = cohort > 0) %>%
      summarise(n = n(), bvar = var(b_i)) %>%
      mutate(trial = pars$try.no)
  },
  mc.cores = 4
) %>%
  do.call(rbind, .)

out4 %>%
  group_by(gen, new.old) %>%
  summarise(
    bvarbar = mean(bvar),
    bvarvar = var(bvar),
    nn = n()
  ) %>%
  ggplot(aes(x = gen, y = bvarbar, group = new.old, colour = new.old)) +
  geom_line()

# yeah... something's not right, even the new/old is going below original

out4 %>%
  filter(!new.old) %>%
  group_by(gen) %>%
  summarise(
    bvarbar = mean(bvar * (1 - 1/n())),
    bvarvar = var(bvar),
    nn = n()
  ) %>%
  mutate(wfitn = pars.list$wfitn) %>%
  mutate(exp.bvar = bvarbar[1] * (wfitn^2 / (wfitn^2 + gen * bvarbar[1]))) %>% 
  ggplot(aes(x = exp.bvar, y = bvarbar)) + 
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), linetype = 2) +
  geom_line() + 
  geom_point(aes(colour = gen), size = 4) + 
  scale_colour_viridis_c()
# so... my expression for variance over time is not quite right for first cohort
# (age structure maybe?)
# is it right for subsequent ones?

out3 %>%
  filter(cohort %in% 1) %>%
  group_by(gen) %>%
  summarise(
    bvarbar = mean(bvar * (1 - 1/n())),
    bvarvar = var(bvar),
    nn = n()
  ) %>%
  mutate(wfitn = pars.list$wfitn) %>%
  mutate(exp.bvar = bvarbar[1] * (wfitn^2 / (wfitn^2 + (gen-1) * bvarbar[1]))) %>% 
  ggplot(aes(x = exp.bvar, y = bvarbar)) + 
  geom_segment(aes(x = 0, xend = 2, y = 0, yend = 2), linetype = 2) +
  geom_line() + 
  geom_point(aes(colour = gen), size = 4) + 
  scale_colour_viridis_c()
# that actually looks about right...

# Look at all of the cohorts
out3 %>%
  group_by(gen, cohort) %>%
  summarise(
    bvarbar  = mean(bvar * (1 - 1/n())),
    bvarvar  = var(bvar),
    nn = n()
  ) %>%
  # Screen out cohorts with extinctions
  filter(nn %in% 100) %>%
  group_by(cohort) %>%
  # Get initial for each cohort
  mutate(bvarbar0 = bvarbar[1]) %>%
  mutate(wfitn = pars.list$wfitn) %>%
  mutate(exp.bvar = bvarbar0 * (wfitn^2 / (wfitn^2 + (gen-cohort)*bvarbar0))) %>% 
  ggplot(aes(x = exp.bvar, y = bvarbar)) + 
  geom_line(aes(group = cohort, colour = cohort)) + 
  geom_segment(aes(x = 0, xend = 2, y = 0, yend = 2), linetype = 2) +
  scale_colour_viridis_c()
# looks actually dead on, with the exception of cohort 1 I think
# wait... no, the expectations for cohort 1 are just wrong, that's all that's happening

out3 %>%
  group_by(gen, cohort) %>%
  summarise(
    bvarbar  = mean(bvar * (1 - 1/n())),
    bvarvar  = var(bvar),
    nn = n()
  ) %>%
  # Screen out cohorts with extinctions
  filter(nn %in% 100) %>%
  group_by(cohort) %>%
  # Get initial for each cohort
  mutate(bvarbar0 = bvarbar[1]) %>%
  mutate(wfitn = pars.list$wfitn) %>%
  mutate(exp.bvar = bvarbar0 * (wfitn^2 / (wfitn^2 + (gen-cohort)*bvarbar0))) %>% 
  ggplot(aes(x = exp.bvar, y = bvarbar / exp.bvar)) + 
  geom_line(aes(group = cohort, colour = cohort > 0)) +
  geom_point(aes(colour = cohort > 0)) #+ 
#  geom_segment(aes(x = 0, xend = 2, y = 0, yend = 2), linetype = 2) +
#   scale_colour_viridis_c()

# yep - it's something with the first cohort
# coming in way below expectation...

# look at a single sim...

test1 %>%
  filter(!coh1) %>%
  group_by(i) %>%
  mutate(age.coh = min(age)) %>%
  group_by(gen, age.coh) %>%
  summarise(
    n = n(),
    bvar = var(b_i) * (1 - 1/n),
  ) %>%
  filter(n > 1) %>%
  ggplot(aes(x = gen, y = bvar, group = age.coh)) +
  geom_line(aes(colour = age.coh)) +
  scale_colour_viridis_c()
# hmm okay not super helpful

test1 %>%
  filter(gen < 5) %>%
  group_by(gen, age) %>%
  summarise(n = n()) %>%
  group_by(gen) %>%
  mutate(p.age = n / sum(n)) %>%
  ggplot(aes(x = age, y = p.age)) +
  geom_line(aes(colour = gen, group = gen))

test1 %>%
  group_by(gen, age) %>%
  summarise(n = n()) %>%
  group_by(gen) %>%
  mutate(p.age = n / sum(n)) %>%
  ggplot(aes(x = age, y = p.age)) +
  geom_line(data = . %>% filter(gen > 0), aes(group = gen)) +
  geom_line(data = . %>% filter(!gen), colour = 'skyblue', size = 3)
# that looks fine to me???

######----------------------------

set.seed(5418)

out5 = mclapply(
  pars.list %>% uncount(ntrials) %>% mutate(try.no = 1:ntrials) %>% split(.$try.no),
  function(pars) {
    sim(params = pars, theta.t = 0, init.row = 5 * 1000 * 10) %>%
      group_by(i) %>%
      filter(!min(gen)) %>%
      mutate(age.cohort = min(age)) %>%
      group_by(gen, age.cohort) %>%
      summarise(n = n(), bvar = var(b_i) * (1 - 1/n)) %>%
      mutate(trial = pars$try.no)
  },
  mc.cores = 4
) %>%
  do.call(rbind, .)

sum5 = out5 %>%
  group_by(gen, age.cohort) %>%
  summarise(
    bvarbar = mean(bvar),
    bvarvar = var(bvar),
    nn = n()
  ) %>%
  filter(nn %in% 100)

sum5 %>%
  ggplot(aes(x = gen, y = bvarbar, group = age.cohort)) +
  geom_line(aes(colour = age.cohort)) +
  scale_colour_viridis_c()

sum5 %>%
  mutate(age = age.cohort + gen) %>%
  mutate(w = pars.list$wfitn, sig.a = pars.list$sig.a) %>%
  mutate(exp.bvar = w^2 * sig.a^2 / (w^2 + age * sig.a^2)) %>%
  ggplot(aes(x = age, bvarbar / exp.bvar)) +
  geom_line(aes(group = age.cohort, colour = age.cohort)) +
  geom_point(aes(colour = age.cohort), size = 3) +
  scale_colour_viridis_c('')
# yeah... this is bad isn't it?

sum5 %>%
  mutate(age = age.cohort + gen) %>%
  mutate(w = pars.list$wfitn, sig.a = pars.list$sig.a) %>%
  mutate(exp.bvar = w^2 * sig.a^2 / (w^2 + age * sig.a^2)) %>%
  ggplot(aes(x = age, bvarbar / exp.bvar)) +
  geom_line(aes(group = age.cohort, colour = age.cohort)) +
  geom_point(aes(fill = gen), size = 3, shape = 21) +
  scale_colour_viridis_c() +
  scale_fill_viridis_c(option = 'B')
# why is this happening...?

sum5 %>%
  mutate(age = age.cohort + gen) %>%
  mutate(w = pars.list$wfitn, sig.a = pars.list$sig.a) %>%
  mutate(exp.bvar = w^2 * sig.a^2 / (w^2 + age * sig.a^2)) %>%
  pivot_longer(cols = c(bvarbar, exp.bvar), names_to = 'vartype', values_to = 'value') %>%
  filter(age.cohort < 16) %>%
  ggplot(aes(x = gen, y = value, group = interaction(vartype, age.cohort))) +
  geom_line(aes(colour = age.cohort, linetype = vartype)) +
  scale_colour_viridis_c() +
  facet_wrap(~ age.cohort)
# this does look like there's a little undershooting...

# any chance this is due to the correction?
sum5.5 = out5 %>%
  group_by(gen, age.cohort) %>%
  summarise(
    bvarbar = mean(bvar * (1 + 1/n)),
    bvarvar = var(bvar * (1 + 1/n)),
    nn = n()
  ) %>%
  filter(nn %in% 100)

sum5.5 %>%
  mutate(age = age.cohort + gen) %>%
  mutate(w = pars.list$wfitn, sig.a = pars.list$sig.a) %>%
  mutate(exp.bvar = w^2 * sig.a^2 / (w^2 + age * sig.a^2)) %>%
  pivot_longer(cols = c(bvarbar, exp.bvar), names_to = 'vartype', values_to = 'value') %>%
  filter(age.cohort < 16) %>%
  ggplot(aes(x = gen, y = value, group = interaction(vartype, age.cohort))) +
  geom_line(aes(colour = age.cohort, linetype = vartype)) +
  scale_colour_viridis_c() +
  facet_wrap(~ age.cohort)
# yean lol these look spot on
# hnnnnngh

sum5.5 %>%
  mutate(age = age.cohort + gen) %>%
  mutate(w = pars.list$wfitn, sig.a = pars.list$sig.a) %>%
  mutate(exp.bvar = w^2 * sig.a^2 / (w^2 + age * sig.a^2)) %>%
  ggplot(aes(x = age, bvarbar / exp.bvar)) +
  geom_line(aes(group = age.cohort, colour = age.cohort)) +
  geom_point(aes(colour = age.cohort), size = 3) +
  scale_colour_viridis_c('')
# yeah it's the correction... shit!

############---------------------------------

set.seed(5418)

out6 = mclapply(
  pars.list %>% uncount(ntrials) %>% mutate(try.no = 1:ntrials) %>% split(.$try.no),
  function(pars) {
    sim(params = pars, theta.t = 0, init.row = 5 * 1000 * 10) %>%
      group_by(gen, age) %>%
      summarise(n = n(), bvar = var(b_i) * (1 - 1/n)) %>%
      mutate(trial = pars$try.no)
  },
  mc.cores = 4
) %>%
  do.call(rbind, .)

sum6 = out6 %>%
  filter(n > 1) %>%
  group_by(gen, age) %>%
  summarise(
    bvarbar = mean(bvar),
    bvarvar = var(bvar),
    nn = n()
  )

# Look at variance of each cohort, where now cohorts represent true age classes
# (rather than lumping in all initialized populations together)
sum6 %>%
  filter(nn > 99) %>%
  mutate(cohort = gen - age) %>%
  ggplot(aes(x = gen, y = bvarbar, group = cohort, colour = cohort)) +
  geom_line(size = 3) +
  scale_color_gradient2(low = 'goldenrod', high = 'royalblue', mid = 'white', midpoint = 0) +
  # scale_y_log10() +
  theme(panel.background = element_rect(fill = 'black'))

# look to see if mean variance at a given age differs between initialized and non-initialized cohorts

sum6 %>%
  filter(nn > 99) %>%
  mutate(cohort = gen - age) %>%
  ggplot(aes(x = gen, y = bvarbar, group = cohort, colour = cohort < 0)) +
  geom_line(size = 3)

sum6 %>%
  filter(nn > 99) %>%
  mutate(cohort.a = gen - age) %>%
  group_by(age, initialized = cohort.a < 0) %>%
  summarise(bvarbar = mean(bvarbar)) %>% 
  ggplot(aes(x = age, y = bvarbar, group = initialized)) +
  geom_line(aes(colour = initialized))
# this looks right too...
# what the heck is going on man????

sum6 %>%
  filter(nn > 99) %>%
  mutate(cohort = gen - age) %>%
  ggplot(aes(x = age, y = bvarbar, group = cohort, colour = abs(cohort))) +
  geom_line() +
  scale_colour_viridis_c() +
  facet_wrap(~ cohort < 0)

# how does this compare to expectations?
sum6 %>%
  filter(nn > 99) %>%
  mutate(cohort = gen - age) %>%
  mutate(w = pars.list$wfitn, sig.a = pars.list$sig.a) %>%
  mutate(exp.bvar = w^2 * sig.a^2 / (w^2 + age * sig.a^2)) %>%
  ggplot(aes(x = bvarbar, y = exp.bvar)) +
  geom_segment(aes(x = 0, xend = sig.a^2, y = 0, yend = sig.a^2),
               linetype = 2) +
  geom_point(aes(colour = cohort), size = 3) +
  scale_colour_viridis_c('') +
  facet_wrap(~ cohort < 0)
# I mean... it does kinda look like the expected is bigger than observed
# but could that be because of weirdness...
  
sum6 %>%
  filter(nn > 99, gen < 10) %>%
  mutate(cohort = gen - age) %>%
  mutate(w = pars.list$wfitn, sig.a = pars.list$sig.a) %>%
  mutate(exp.bvar = w^2 * sig.a^2 / (w^2 + age * sig.a^2)) %>%
  ggplot(aes(x = bvarbar, y = exp.bvar)) +
  geom_segment(aes(x = 0, xend = sig.a^2, y = 0, yend = sig.a^2),
               linetype = 2) +
  geom_point(aes(colour = abs(cohort)), size = 3) +
  scale_colour_viridis_c('') +
  facet_wrap(~ cohort < 0)
# more deviation for higher cohorts... 

sum6 %>%
  filter(nn > 99, gen < 10) %>%
  mutate(cohort = gen - age) %>%
  mutate(w = pars.list$wfitn, sig.a = pars.list$sig.a) %>%
  mutate(exp.bvar = w^2 * sig.a^2 / (w^2 + age * sig.a^2)) %>%
  ggplot(aes(x = age, bvarbar / exp.bvar)) +
  geom_line(aes(group = cohort, colour = abs(cohort))) +
  geom_point(aes(colour = abs(cohort)), size = 3) +
  scale_colour_viridis_c('') +
  facet_wrap(~ cohort < 0)
# okay, so even in the initial cohort... some bad stuff happening yes?
# why are these way undershooting *their* expectations?

sum6.5 = out6 %>%
  filter(n > 1) %>%
  group_by(gen, age) %>%
  summarise(
    bvarbar = mean(bvar * (1 + 1/n)),
    bvarvar = var(bvar * (1 + 1/n)),
    nn = n()
  )

sum6.5 %>%
  filter(nn > 99, gen < 10) %>%
  mutate(cohort = gen - age) %>%
  mutate(w = pars.list$wfitn, sig.a = pars.list$sig.a) %>%
  mutate(exp.bvar = w^2 * sig.a^2 / (w^2 + age * sig.a^2)) %>%
  ggplot(aes(x = age, bvarbar / exp.bvar)) +
  geom_line(aes(group = cohort, colour = abs(cohort))) +
  geom_point(aes(colour = abs(cohort)), size = 3) +
  scale_colour_viridis_c('') +
  facet_wrap(~ cohort < 0)
# ehhhh yeah this looks fine I think.
# so... calculating the varianees, we should leave off that correction
# (which seems kinda fucked to me but w/e)

########===============================================

# 11/7
# I tried a bunch of shit to try to figure out what's wrong and none of it worked...
# - age structure looks fine and consistent across generations
# - I think my expression for expected variance in a cohort of a certain age is right
# - mutation rate... seems right? idk, I can't tell where the math would be wrong
# - bessel correction warps observed variances... what does this mean for sims??