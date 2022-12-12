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

pvar = function(x) mean(x^2) - mean(x)^2

#####----------------------------------------
# Setup simulations

ntrials = 500
l.max = 1.1

pars.list = data.frame(
  n.pop0 = 1000, s.max = 0.9,
  wfitn = sqrt(10),
  sig.a.pop = 1,
  sig.e = 0, alpha = 0.0000,
  kceil = 2000,
  timesteps = 20
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
  geom_ribbon(
    aes(
      ymin = nbar - 2 * sqrt(nvar / nn),
      ymax = nbar + 2 * sqrt(nvar / nn)
    ),
    alpha = 0.1
  ) +
  scale_y_log10()

sum1 %>%
  ggplot(aes(x = gen, y = bvarbar)) +
  geom_line() +
  geom_ribbon(
    aes(
      ymin = bvarbar - 2 * sqrt(bvarvar / nn),
      ymax = bvarbar + 2 * sqrt(bvarvar / nn)
    ),
    fill = 'blue',
    alpha = 0.2
  )

# argghhhhhhhhhh so close and yet still not there... what's missing? jeez

# needs work but much closer than ever before...
# (variance is still the problem here... declining a wee bit still, somehow)

# Population growth against genetic diversity
sum1 %>% 
  mutate(lambda = c(diff(log(nbar)), NA)) %>% 
  filter(!is.na(lambda)) %>% 
  ggplot(aes(x = bvarbar, y = lambda)) + 
  geom_path() + 
  geom_point(aes(colour = gen)) + 
  scale_colour_viridis_c()
# well... the declining growth at the end is due to the ceiling
# but also, bvarbar doesn't look like it actually influences growth that much??>

out1 %>%
  group_by(trial) %>%
  mutate(log_lam = c(diff(log(n)), NA)) %>%
  filter(!is.na(log_lam)) %>%
  group_by(gen) %>%
  summarise(
    bvarbar = mean(bvar),
    bvarvar = var(bvar),
    llbar = mean(log_lam),
    llvar = var(log_lam),
    nn = n()
  ) %>%
  ggplot(aes(x = bvarbar, y = llbar)) +
  geom_path(aes(colour = gen), size = 3) +
  scale_color_viridis_c()

out1 %>%
  group_by(trial) %>%
  mutate(log_lam = c(diff(log(n)), NA)) %>%
  filter(!is.na(log_lam)) %>%
  ggplot(aes(x = bvar, y = log_lam)) +
  geom_point(alpha = 0.5) +
  scale_color_viridis_c() +
  facet_wrap(~ gen)
# ??

out1 %>%
  mutate(wfitn = pars.list$wfitn,
         s.max = pars.list$s.max,
         r     = pars.list$r) %>%
  group_by(trial) %>%
  mutate(lamda = exp(c(diff(log(n)), NA))) %>%
  filter(!is.na(lamda)) %>%
  mutate(elbar = sqrt(wfitn^2 / (wfitn^2 + bvar)) * s.max * (1 + r)) %>%
  ggplot(aes(x = elbar, y = lamda / elbar)) +
  geom_point(aes(colour = lamda > elbar), alpha = 0.5) +
  facet_wrap(~ gen)
# yeah okay this looks fine
# (dang wtf was I seeing above then?)

######--------------------------------

pars.list = data.frame(
  n.pop0 = 1000, s.max = 0.9,
  # wfitn = sqrt(10),
  sig.a.pop = 1,
  sig.e = 0, alpha = 0.0000,
  kceil = 2000,
  timesteps = 15
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
      ymax = bvarbar + 2 * sqrt(bvarvar / nn)
    ),
    alpha = 0.2
  )

plot(
  x = with(pars.list, s.max * (1 + r) * sqrt(wfitn^2 / (wfitn^2 + sum2$bvarbar)))[-pars.list$timesteps],
  y = exp(diff(log(sum2$nbar))),
  xlab = 'predicted growth rate',
  ylab = 'observed growth rate',
  asp = 1
)
abline(a = 0, b = 1)
# oh... shit

lam2 = out2 %>%
  mutate(wfitn = pars.list$wfitn,
         s.max = pars.list$s.max,
         r     = pars.list$r) %>%
  group_by(trial) %>%
  mutate(lamda = exp(c(diff(log(n)), NA))) %>%
  filter(!is.na(lamda)) %>%
  mutate(elbar = sqrt(wfitn^2 / (wfitn^2 + bvar)) * s.max * (1 + r))

lam2 %>%
  ggplot(aes(x = elbar, y = lamda)) +
  geom_point(aes(colour = gen), alpha = 0.1) +
  scale_colour_viridis_c()

lam2 %>%
  ggplot(aes(x = elbar, y = lamda)) +
  geom_point() +
  facet_wrap(~ gen)

lam2 %>%
  ggplot(aes(x = elbar, y = lamda / elbar)) +
  geom_point(aes(colour = lamda > elbar)) +
  facet_wrap(~ gen)
# looks fine...
# I think based on this my expression for expected growth is right
# (but why did it look wrong before...?)

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
  summarise(bvar = var(b_i)) %>%
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

test1 %>%
  filter(gen < 20) %>%
  mutate(r = pars.list$r,
         wfitn = pars.list$wfitn,
         s.max = pars.list$s.max,
         l.max = s.max * (1 + r),
         bthrs = sqrt(2 * wfitn^2 * log(s.max * (1 + r)))) %>%
  ggplot(aes(x = i, y = b_i)) +
  geom_segment(aes(x = 0, xend = max(i), y = bthrs, yend = bthrs), col = 'white', linetype = 2) +
  geom_segment(aes(x = 0, xend = max(i), y = -bthrs, yend = -bthrs), col = 'white', linetype = 2) +
  geom_point(aes(colour = s_i * (1 + r)), size = 0.1) +
  scale_color_gradient2(low = 'red', high = 'lightblue', mid = 'white', midpoint = 1, '') +
  facet_wrap(~ gen) +
  theme(panel.background = element_rect(fill = 'black'),
        panel.grid = element_blank()) 
# whoa... this is interesting

test1 %>%
  filter(gen %in% c(0, 9)) %>%
  mutate(r = pars.list$r,
         wfitn = pars.list$wfitn,
         s.max = pars.list$s.max,
         l.max = s.max * (1 + r),
         bthrs = sqrt(2 * wfitn^2 * log(s.max * (1 + r)))) %>%
  ggplot(aes(x = age, y = b_i)) +
  geom_point(position = position_jitter(width = 0.25), size = 0.5, colour = 'white') +
  geom_segment(aes(x = 0, xend = max(age), y = bthrs, yend = bthrs), col = 'white', linetype = 2) +
  geom_segment(aes(x = 0, xend = max(age), y = -bthrs, yend = -bthrs), col = 'white', linetype = 2) +
  theme(panel.background = element_rect(fill = 'black'),
        panel.grid = element_blank()) +
  facet_wrap(~ gen, nrow = 2)

# turn that pyramid on its side lmao
test1 %>%
  filter(gen %in% c(0, 20), age < 10) %>%
  mutate(r = pars.list$r,
         wfitn = pars.list$wfitn,
         s.max = pars.list$s.max,
         l.max = s.max * (1 + r),
         bthrs = sqrt(2 * wfitn^2 * log(s.max * (1 + r)))) %>%
  ggplot(aes(x = b_i, y = age)) +
  geom_point(position = position_jitter(height = 0.2), size = 0.5, colour = 'white') +
  geom_segment(aes(y = 0, yend = max(age), x = bthrs, xend = bthrs), col = 'white', linetype = 2) +
  geom_segment(aes(y = 0, yend = max(age), x = -bthrs, xend = -bthrs), col = 'white', linetype = 2) +
  theme(panel.background = element_rect(fill = 'black'),
        panel.grid = element_blank()) +
  facet_wrap(~ gen, nrow = 1)

test1 %>%
  filter(gen %in% c(0, 1, 2, 5, 10, 15, 20)) %>%
  mutate(r = pars.list$r,
         wfitn = pars.list$wfitn,
         s.max = pars.list$s.max,
         l.max = s.max * (1 + r),
         bthrs = sqrt(2 * wfitn^2 * log(s.max * (1 + r)))) %>%
  ggplot(aes(x = b_i)) +
  geom_density(aes(group = gen, color = factor(gen))) +
  geom_segment(aes(y = 0, yend = 0.5, x = bthrs, xend = bthrs), col = 'white', linetype = 3) +
  geom_segment(aes(y = 0, yend = 0.5, x = -bthrs, xend = -bthrs), col = 'white', linetype = 3) +
  geom_line(
    data = data.frame(b_i = (-26:26)/10, p_b = dnorm((-26:26)/10)),
    aes(x = b_i, y = p_b),
    colour = 'white', linetype = 2
  ) +
  theme(panel.background = element_rect(fill = 'black'),
        panel.grid = element_blank())

test0 = test1 %>%
  mutate(
    r = pars.list$r,
    wfitn = pars.list$wfitn,
    s.max = pars.list$s.max,
    l.max = s.max * (1 + r),
    bthrs = sqrt(2 * wfitn^2 * log(s.max * (1 + r)))
  )

test0 %>%
  group_by(gen) %>%
  summarise(p.out = mean(abs(z_i) < bthrs)) %>%
  ggplot(aes(x = gen, y = p.out)) +
  geom_line()

weird = test0 %>%
  filter(fem) %>%
  ungroup() %>%
  group_by(i) %>%
  summarise(
    age.cohort = gen[1] - age[1],
    in.ztar = sum(abs(z_i) < bthrs),
    n.offsp = sum(r_i)/2
  )

weird %>%
  ggplot(aes(x = factor(in.ztar > 0), y = n.offsp)) +
  geom_point(aes(colour = age.cohort < 0),
             position = position_jitterdodge(
               jitter.height = 0.125,
               dodge.width = 1
             ))

test0 %>%
  distinct(i, .keep_all = TRUE) %>%
  group_by(age.cohort = gen - age) %>%
  summarise(p.out = mean(abs(z_i) > bthrs),
            n = n()) %>%
  ggplot(aes(x = age.cohort, y = p.out)) +
  geom_point(aes(size = n))

test0 %>%
  group_by(gen, age.cohort = gen - age) %>%
  summarise(p.out = mean(abs(z_i) > bthrs),
            n = n()) %>%
  ggplot(aes(x = age.cohort, y = p.out)) +
  geom_point(aes(size = n)) +
  facet_wrap(~ gen)

test0 %>%
  group_by(gen, age.cohort = gen - age) %>%
  summarise(bvar = pvar(b_i),
            n = n()) %>%
  ggplot(aes(x = age.cohort, y = bvar)) +
  geom_point(aes(size = n)) +
  facet_wrap(~ gen)

test0 %>%
  group_by(gen, age.cohort = gen - age) %>%
  summarise(bvar = pvar(b_i),
            n = n()) %>%
  group_by(gen) %>%
  summarise(bvarbar = weighted.mean(bvar, w = n))

test0 %>%
  group_by(gen) %>%
  summarise(bvar = var(b_i))
# oh huh... is there some rounding discrepancies between these two??
# (but test0 is observed, which is still falling below 1...)

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
      summarise(n = n(), bvar = var(b_i), bvar.adj = mean(b_i^2) - mean(b_i)^2) %>%     
      mutate(trial = pars$try.no)
  },
  mc.cores = 4
) %>%
  do.call(rbind, .)

sum6 = out6 %>%
  group_by(gen, age) %>%
  summarise(
    bvarbar = mean(bvar.adj),
    bvarvar = var(bvar.adj),
    nn = n()
  )

# Look at variance of each cohort, where now cohorts represent true age classes
# (rather than lumping in all initialized populations together)
sum6 %>%
  filter(nn == ntrials) %>%
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

sum6 %>%
  filter(nn %in% ntrials, age >= gen) %>%
  mutate(expect.var = with(pars.list, sig.a^2 * wfitn^2 / (wfitn^2 + sig.a^2 * age))) %>%
  ggplot(aes(x = age, y = bvarbar / expect.var)) +
  geom_line(aes(group = gen, colour = gen)) +
  geom_point(aes(group = gen, colour = gen), size = 3) +
  scale_colour_viridis_c()

sum6.5 = out6 %>%
  filter(n > 1) %>%
  group_by(gen, age) %>%
  summarise(
    bvarbar = mean(bvar.adj),
    bvarvar = var(bvar.adj),
    nn = n()
  )

sum6.5 %>%
  filter(nn %in% ntrials, age >= gen) %>%
  mutate(expect.var = with(pars.list, sig.a^2 * wfitn^2 / (wfitn^2 + sig.a^2 * age))) %>%
  ggplot(aes(x = age, y = bvarbar / expect.var)) +
  geom_line(aes(group = gen, colour = gen)) +
  geom_point(aes(group = gen, colour = gen), size = 3) +
  scale_colour_viridis_c()

sum6.5 %>%
  filter(nn %in% ntrials, age >= gen) %>%
  mutate(expect.var = with(pars.list, sig.a^2 * wfitn^2 / (wfitn^2 + sig.a^2 * age))) %>%
  ggplot(aes(x = age, y = bvarbar)) +
  geom_line(aes(group = gen, colour = gen)) +
  geom_point(aes(group = gen, colour = gen), size = 3) +
  geom_line(aes(y = expect.var), linetype = 2) +
  scale_colour_viridis_c()

sum6.6 = out6 %>%
  group_by(gen, age) %>%
  summarise(
    bvarbar = mean(bvar),
    bvarvar = var(bvar),
    nn = n()
  )

sum6.6 %>%
  filter(nn %in% ntrials, age >= gen) %>%
  mutate(expect.var = with(pars.list, sig.a^2 * wfitn^2 / (wfitn^2 + sig.a^2 * age))) %>%
  ggplot(aes(x = age, y = bvarbar)) +
  geom_line(aes(group = gen, colour = gen)) +
  geom_point(aes(group = gen, colour = gen), size = 3) +
  geom_line(aes(y = expect.var), linetype = 2) +
  scale_colour_viridis_c()
# lmao gotta stop with this stupid fucking bessel correction shit goddamn

sum6.6 %>%
  filter(nn %in% ntrials, age >= gen) %>%
  mutate(expect.var = with(pars.list, sig.a^2 * wfitn^2 / (wfitn^2 + sig.a^2 * age))) %>%
  ggplot(aes(x = age, y = bvarbar / expect.var)) +
  geom_segment(aes(x = 0, xend = 20, y = 1, yend = 1), linetype = 2) +
  geom_line(aes(group = gen, colour = gen)) +
  geom_point(aes(group = gen, colour = gen), size = 3) +
  scale_colour_viridis_c()



#####---------------------

# kinda digging this threshold stuff...

set.seed(5418)

out7 = mclapply(
  pars.list %>% uncount(ntrials) %>% mutate(try.no = 1:ntrials) %>% split(.$try.no),
  function(pars) {
    sim(params = pars, theta.t = 0, init.row = 5 * 1000 * 10) %>%
      mutate(p.out = abs(z_i) > with(pars.list, sqrt(2*wfitn^2*log(s.max*(1+r))))) %>%
      group_by(gen) %>%
      summarise(n = n(), bvar = var(b_i) * (1 - 1/n), p.out = mean(p.out)) %>%
      mutate(trial = pars$try.no)
  },
  mc.cores = 4
) %>%
  do.call(rbind, .)

out7 %>%
  ggplot(aes(x = gen, y = p.out)) +
  geom_point(position = position_jitter(width = 0.2))

out7 %>%
  group_by(trial) %>%
  mutate(log_lam = c(diff(log(n)), NA)) %>%
  filter(!is.na(log_lam)) %>%
  ggplot(aes(x = gen, y = p.out, colour = log_lam)) +
  geom_point(position = position_jitter(width = 0.4)) +
  scale_colour_gradient2(low = 'red', high = 'royalblue', mid = 'white', midpoint = 0) +
  theme(panel.background = element_rect(fill = 'black'),
        panel.grid = element_blank())
# doesn't look like a strong association to me...

out7 %>%
  group_by(trial) %>%
  mutate(log_lam = c(diff(log(n)), NA)) %>%
  filter(!is.na(log_lam)) %>%
  ggplot(aes(x = gen, y = p.out, colour = log_lam > 0)) +
  geom_point(position = position_jitter(width = 0.4)) +
  theme(panel.background = element_rect(fill = 'black'),
        panel.grid = element_blank())
# hmm... I guess... maybe...?

########===============================================

set.seed(5418)

out8 = mclapply(
  pars.list %>% uncount(10000) %>% mutate(try.no = 1:n()) %>% split(.$try.no),
  function(pars) {
    sim(params = pars %>% mutate(timesteps = 3), theta.t = 0, init.row = 5 * 1000 * 10) %>%
      group_by(gen, age) %>%
      summarise(n = n(), bvar = var(b_i), bvar.adj = mean(b_i^2) - mean(b_i)^2) %>%
      mutate(trial = pars$try.no)
  },
  mc.cores = 4
) %>%
  do.call(rbind, .)

sum8 = out8 %>%
  group_by(gen, age) %>%
  summarise(
    bvarbar = mean(bvar.adj),
    bvarvar = var(bvar.adj),
    n = n()
  )

sum8 %>%
  ggplot(aes(x = age, y = bvarbar, group = gen)) +
  geom_line(aes(colour = gen))

# Relative to mean (of gens for timestep)
sum8 %>% 
  group_by(age) %>% 
  filter(all(n %in% max(out8$trial))) %>% 
  mutate(bvarbar.resid = bvarbar - mean(bvarbar)) %>% 
  ggplot(aes(x = age, y = bvarbar.resid)) + 
  geom_ribbon(
    aes(
      group = gen,
      fill = gen,
      ymin = bvarbar.resid - 2*sqrt(bvarvar / n),
      ymax = bvarbar.resid + 2*sqrt(bvarvar / n)
    ),
    alpha = 0.1
  ) +
  geom_segment(aes(x = 0, xend = 14, y = 0, yend = 0), linetype = 2, colour = 'gray44') +
  geom_point(aes(colour = gen), size = 3) + 
  geom_line(aes(colour = gen, group = gen)) + 
  scale_colour_viridis_c() +
  scale_fill_viridis_c()
  
# relative to my analytical expectation
sum8 %>% 
  group_by(age) %>% 
  filter(all(n %in% max(out8$trial))) %>%
  ungroup() %>%
  mutate(
    pred.var = with(pars.list, sig.a^2 * wfitn^2 / (wfitn^2 + sig.a^2 * age))
  ) %>%
  ggplot(aes(x = age, y = bvarbar / pred.var)) + 
  geom_segment(aes(x = 0, xend = max(age), y = 1, yend = 1), linetype = 2, colour = 'gray44') +
  geom_point(aes(colour = gen), size = 3) + 
  geom_line(aes(colour = gen, group = gen)) + 
  scale_colour_viridis_c() +
  scale_fill_viridis_c()

# getting dragged way below... wonder if that's due to all of the zeros though???

# relative to my analytical expectation
# but looking at cohorts that are pretty unlikely to have zeros...
sum8 %>% 
  filter(age < 11) %>%
  mutate(
    pred.var = with(pars.list, sig.a^2 * wfitn^2 / (wfitn^2 + sig.a^2 * age))
  ) %>%
  ggplot(aes(x = age, y = bvarbar / pred.var)) + 
  geom_segment(aes(x = 0, xend = 10, y = 1, yend = 1), linetype = 2, colour = 'gray44') +
  geom_point(aes(colour = gen), size = 3) + 
  geom_line(aes(colour = gen, group = gen)) + 
  scale_colour_viridis_c() +
  scale_fill_viridis_c()
# uh... looks inconsistent?

sum8 %>% 
  filter(age < 11) %>%
  mutate(pred.var = with(pars.list, sig.a^2 * wfitn^2 / (wfitn^2 + sig.a^2 * age))) %>%
  group_by(age) %>%
  mutate(pred.var.resid = (bvarbar / pred.var) - mean(bvarbar / pred.var)) %>%
  ggplot(aes(x = age, y = pred.var.resid)) + 
  geom_segment(aes(x = 0, xend = 10, y = 0, yend = 0), linetype = 2, colour = 'gray44') +
  geom_point(aes(colour = gen), size = 3) + 
  geom_line(aes(colour = gen, group = gen)) + 
  scale_colour_viridis_c() +
  scale_fill_viridis_c()
# it does kinda just look like... random fluctuation??
# at least for the later generations...
# but there is an obvious pattern for the first cohort!

# I think it's possi ble the mutation rate is just fucked up?

sum8.5 = out8 %>%
  group_by(gen, age) %>%
  summarise(
    bvarbar = mean(bvar),
    bvarvar = var(bvar),
    n = n()
  )

sum8.5 %>%
  filter(age < 11) %>%
  mutate(pred.var = with(pars.list, sig.a^2 * wfitn^2 / (wfitn^2 + sig.a^2 * age))) %>%
  group_by(age) %>%
  mutate(pred.var.resid = (bvarbar / pred.var) - mean(bvarbar / pred.var)) %>%
  ggplot(aes(x = age, y = pred.var.resid)) + 
  geom_segment(aes(x = 0, xend = 10, y = 0, yend = 0), linetype = 2, colour = 'gray44') +
  geom_point(aes(colour = gen), size = 3) + 
  geom_line(aes(colour = gen, group = gen)) + 
  scale_colour_viridis_c() +
  scale_fill_viridis_c()

sum8.5 %>%
  filter(n %in% 10000) %>%
  mutate(
    r = pars.list$r,
    sig.a = pars.list$sig.a,
    wfitn = pars.list$wfitn
  ) %>%
  mutate(
    prob.age = (1 / (1+r))^age * (r / (1+r)),
    exp.varn = sig.a^2 * wfitn^2 / (wfitn^2 + sig.a^2 * age)
  ) %>%
  # group_by(gen) %>%
  # mutate(gen.mean = (r/(1+r)) * (1/(1+r))^k * exp.varn) %>%
  # ungroup() %>%
  ggplot(aes(x = age, y = bvarbar / exp.varn)) +
  geom_segment(aes(x = 0, xend = 15, y = 1, yend = 1), linetype = 2, colour = 'gray44') +
  geom_point(aes(colour = gen), size = 3) + 
  geom_line(aes(colour = gen, group = gen)) + 
  scale_colour_viridis_c() #+
  # scale_fill_viridis_c()

sum8.5 %>%
  mutate(
    r = pars.list$r,
    sig.a = pars.list$sig.a,
    wfitn = pars.list$wfitn
  ) %>%
  mutate(
    prob.age = (1 / (1+r))^age * (r / (1+r)),
    exp.varn = sig.a^2 * wfitn^2 / (wfitn^2 + sig.a^2 * age)
  ) %>%
  # group_by(gen) %>%
  # mutate(gen.mean = (r/(1+r)) * (1/(1+r))^k * exp.varn) %>%
  # ungroup() %>%
  ggplot(aes(x = age, y = bvarbar / exp.varn)) +
  geom_segment(aes(x = 0, xend = 15, y = 1, yend = 1), linetype = 2, colour = 'gray44') +
  geom_point(aes(colour = gen), size = 3) + 
  geom_line(aes(colour = gen, group = gen)) + 
  scale_colour_viridis_c() #+
# scale_fill_viridis_c()

sum8.5 %>%
  mutate(
    r = pars.list$r,
    sig.a = pars.list$sig.a,
    wfitn = pars.list$wfitn
  ) %>%
  mutate(
    prob.age = (1 / (1+r))^age * (r / (1+r)),
    exp.varn = sig.a^2 * wfitn^2 / (wfitn^2 + sig.a^2 * age)
  ) %>%
  # group_by(gen) %>%
  # mutate(gen.mean = (r/(1+r)) * (1/(1+r))^k * exp.varn) %>%
  # ungroup() %>%
  ggplot(aes(x = age, y = bvarbar / exp.varn)) +
  geom_segment(aes(x = 0, xend = 15, y = 1, yend = 1), linetype = 2, colour = 'gray44') +
  geom_point(aes(colour = gen), size = 3) + 
  geom_line(aes(colour = gen, group = gen)) + 
  scale_colour_viridis_c()

out8 %>% 
  filter(age >= gen) %>% 
  # group_by(trial, cohort = gen - age) %>% 
  # mutate(var.loss = c(diff(bvar), NA)) %>% 
  # filter(!is.na(var.loss)) %>% 
  group_by(gen, age) %>% 
  summarise(bvar = mean(bvar), nn = n()) %>% 
  filter(nn %in% max(out8$trial)) %>%
  mutate(
    r = pars.list$r,
    sig.a = pars.list$sig.a,
    wfitn = pars.list$wfitn
  ) %>%
  mutate(
    prob.age = (1 / (1+r))^age * (r / (1+r)),
    exp.varn = sig.a^2 * wfitn^2 / (wfitn^2 + sig.a^2 * age)
  ) %>%
  ggplot(aes(x = age, y = bvar / exp.varn, group = gen)) +
  geom_point(aes(colour = gen), size = 3) +
  geom_line(aes(colour = gen)) + 
  scale_colour_viridis_c()

#######-------------------

# some shit is wrong with the mutation rate
# (best guess is... var(pop.surv$b_i)/2 is not right for the sig.m formula?)
set.seed(5418)

out9 = mclapply(
  pars.list %>% uncount(1000) %>% mutate(try.no = 1:n()) %>% split(.$try.no),
  function(pars) {
    sim(params = pars %>% mutate(timesteps = 2), theta.t = 0, init.row = 5 * 1000 * 10) %>%
      group_by(gen, new.coho = ifelse(age > gen, -1, gen - age)) %>%
      summarise(bvar = var(b_i)) %>%
      mutate(trial = pars$try.no)
  },
  mc.cores = 4
) %>%
  do.call(rbind, .)

out9 %>%
  group_by(gen, new.coho) %>%
  summarise(
    bvarbar = mean(bvar),
    bvarvar = var(bvar),
    nn = n()
  ) %>%
  mutate(try.varb = bvarbar + pars.list$sig.m^2) %>%
  as.data.frame()
# goddamn do I just have zero fucking idea

sum9 = out9 %>%
  group_by(gen, new.coho) %>%
  summarise(
    bvarbar = mean(bvar),
    bvarvar = var(bvar),
    nn = n()
  )

sum9 %>%
  filter(new.coho >= 0) %>%
  mutate(new.coho = factor(new.coho)) %>%
  ggplot(aes(x = gen, y = bvarbar)) + 
  geom_line(aes(group = new.coho, colour = new.coho)) +
  geom_ribbon(
    aes(
      ymin = bvarbar - 2 * sqrt(bvarvar / nn),
      ymax = bvarbar + 2 * sqrt(bvarvar / nn),
      fill = new.coho,
      group = new.coho
    ),
    alpha = 0.1
  ) +
  geom_line(aes(group = as.numeric(new.coho) - gen)) +
  geom_ribbon(
    aes(
      ymin = bvarbar - 2 * sqrt(bvarvar / nn),
      ymax = bvarbar + 2 * sqrt(bvarvar / nn),
      group = as.numeric(new.coho) - gen
    ),
    alpha = 0.1
  ) +
  geom_point(aes(group = new.coho, colour = new.coho), size = 3) 

# is that 0 cohort (gen0, age0) losing variance as expected? how much is it adding via mutation?

#########-----------------------------------------

set.seed(5418)

out10 = mclapply(
  pars.list %>% uncount(1000) %>% mutate(try.no = 1:n()) %>% split(.$try.no),
  function(pars) {
    sim(params = pars %>% mutate(n.pop0 = 10000, kceil = 12000, timesteps = 3), theta.t = 0, init.row = 5 * 1000 * 10) %>%
      group_by(gen) %>%
      summarise(n = n(), bvar = var(b_i)) %>%
      mutate(trial = pars$try.no)
  },
  mc.cores = 4
) %>%
  do.call(rbind, .)

sum10 = out10 %>%
  group_by(gen) %>%
  summarise(
    bvarbar = mean(bvar),
    bvarvar = var(bvar),
    nn = n()
  )

sum10 %>%
  ggplot(aes(x = gen, y = bvarbar)) + 
  geom_line() +
  geom_ribbon(
    aes(
      ymin = bvarbar - 2 * sqrt(bvarvar / nn),
      ymax = bvarbar + 2 * sqrt(bvarvar / nn)
    ),
    alpha = 0.1
  ) +
  geom_point(size = 3) 
  
out10.5 = mclapply(
  pars.list %>% uncount(1000) %>% mutate(try.no = 1:n()) %>% split(.$try.no),
  function(pars) {
    sim(params = pars %>% mutate(n.pop0 = 1000, kceil = 1200, timesteps = 3), theta.t = 0, init.row = 5 * 1000 * 10) %>%
      group_by(gen) %>%
      summarise(n = n(), bvar = var(b_i)) %>%
      mutate(trial = pars$try.no)
  },
  mc.cores = 4
) %>%
  do.call(rbind, .)

sum10.5 = out10.5 %>%
  group_by(gen) %>%
  summarise(
    bvarbar = mean(bvar),
    bvarvar = var(bvar),
    nn = n()
  )

rbind(
  sum10 %>% mutate(popsize = 'large'),
  sum10.5 %>% mutate(popsize = 'small')
) %>%
  ggplot(aes(x = gen, y = bvarbar, group = popsize, colour = popsize, fill = popsize)) + 
  geom_line() +
  geom_ribbon(
    aes(
      ymin = bvarbar - 2 * sqrt(bvarvar / nn),
      ymax = bvarbar + 2 * sqrt(bvarvar / nn),
    ),
    alpha = 0.1
  ) +
  geom_point(size = 3) 


######------------------------------------------

# if I wantd to brute-force it... what would I do??

set.seed(5418)

out11 = mclapply(
  pars.list %>% uncount(100000) %>% mutate(try.no = 1:n()) %>% split(.$try.no),
  function(pars) {
    sim(params = pars %>% mutate(n.pop0 = 1000, mu = 0, timesteps = 1), theta.t = 0, init.row = 5 * 1000 * 10) %>%
      group_by(gen) %>%
      summarise(n = n(), bvar = var(b_i)) %>%
      mutate(trial = pars$try.no)
  },
  mc.cores = 4
) %>%
  do.call(rbind, .)

sum11 = out11 %>%
  group_by(gen) %>%
  summarise(
    nbar = mean(n),
    nvar = var(n),
    bvarbar = mean(bvar),
    bvarvar = var(bvar),
    nn = n()
  )

sum11 %>% as.data.frame()
sum11$bvarbar %>% diff() %>% abs() %>% 
  (function(x) with(pars.list, (1+r) / r * x)) %>%
  sqrt()
# 1.115644, compared with my estimate of 1.07...

out11 %>%
  group_by(trial) %>%
  mutate(vardiff = c(diff))


####3#------------------------------------------

# Wrapper for determining the mutation rate (given s, w, r)
mutofun2 = function(w, s, r, ps) {
  # NOTE: s is the *new cohort variance*, not the population-level variance (s.d.)
  sumo = 0
  for (k in 0:1000) {
    sumo = sumo +
      (1 / (1+r))^k * ( (w^2 + k*s^2)^(-1) - (w^2 + (k+1)*s^2)^(-1) )
  }
  sumo = sumo * s^2 * w^2 * (r / (1+r))
  return((ps - sumo) / r)
  # sumo is sigma_m^2, not sigma_m
}

# 11/7
# I tried a bunch of shit to try to figure out what's wrong and none of it worked...
# - age structure looks fine and consistent across generations
# - I think my expression for expected variance in a cohort of a certain age is right
# - mutation rate... seems right? idk, I can't tell where the math would be wrong
# - bessel correction warps observed variances... what does this mean for sims??

# 11/10
# - it's gotta be the mutation rate, man, not sure what else it could be
# - but what about the mutation rate is wrong...?

# 11/11
# - Brett suggests that it's a second-order term missing from variance
# - (after some thought: my expr for the variance assumes all means are 0
# - but I suppose there is sampling variance in these means?)
# - (and/or variance in the size of each cohort... multinomial...)


