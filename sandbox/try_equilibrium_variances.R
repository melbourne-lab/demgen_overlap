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
  }
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
    alpha = 0.1
  )

# argghhhhhhhhhh so close and yet still not there... what's missing? jeez

# needs work but much closer than ever before...
# (variance is still the problem here... declining a wee bit still, somehow)
