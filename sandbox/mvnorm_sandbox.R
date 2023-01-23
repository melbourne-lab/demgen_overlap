# Trying out options for mvnorm function
# (Multi-variate normal) - needed for phenos and genos
# using this to look for proper non-heritable phenotypic variance components
# (and phenotypic/genotypic change)

library(ggplot2)
library(tidyr)
library(dplyr)

rm(list = ls())

g2a0 = .1
g2e0 = .3

vcvt = vector(mode = "list", length = 100)

for (k in 1:length(vcvt)) {
  vcvk = c(g2a0 * (1 + (k-1)*g2e0) / (1 + (k-1)*(g2a0 + g2e0)),
           -(k-1) * g2a0 * g2e0 / (1 + (k-1)*(g2a0 + g2e0)),
           -(k-1) * g2a0 * g2e0 / (1 + (k-1)*(g2a0 + g2e0)),
           g2e0 * (1 + (k-1)*g2a0) / (1 + (k-1)*(g2a0 + g2e0))
  )
  vcvt[[k]] = matrix(vcvk, nrow = 1)
}

vcvt

mnvt = vector(mode = "list", length = 100)

library(mvtnorm)

rmvnorm(100, mean = c(0, 0), sigma = vcvt[[1]]) %>%
  as.data.frame() %>%
  ggplot(aes(x = V1, y = V2)) +
  geom_point()

rmvnorm(100, mean = c(0, 0), sigma = vcvt[[1:100]]) %>%
  as.data.frame() %>%
  ggplot(aes(x = V1, y = V2)) +
  geom_point() +
  geom_point(aes(x = 0, y = V1 + V2), colour = 'red')
# oh wow
#

# but... not vectorized...

library(mc2d)

test.mv = rmultinormal(100, mean = rep(0, 2), sigma = do.call(rbind, vcvt))

test.mv %>%
  as.data.frame() %>%
  mutate(age = 1:100) %>%
  ggplot(aes(x = V1, y = V2, age = age)) +
  geom_point(aes(colour = age))

test.mv2 = rmultinormal(200, mean = rep(0, 2), sigma = do.call(rbind, vcvt[c(rep(1, 100), 1:100)]))

test.mv2 %>%
  as.data.frame() %>%
  mutate(age = c(rep(1, 100), 1:100)) %>%
  ggplot(aes(x = V1, y = V2, age = age)) +
  geom_point(aes(colour = age))

sig.a = sqrt(g2a0)
sig.e = sqrt(g2e0)
wfitn = 1
s.max = .9
r     = 2 / .9 - 1
lstar = 1.9
size0 = 10000
theta0 = 0
alpha = 0
  
p.age = data.frame(age = 0:500) %>%
  mutate(p.age = (r / (1 + r)) * (s.max / lstar)^age * sqrt(wfitn^2 / (wfitn^2 + age*(sig.a^2 + sig.e^2))))

# Build variance-covariance matrix
vcv.a = vector('list', length = nrow(p.age))

for (k in 1:length(vcv.a)) {
  vcvk = c(sig.a^2 * (1 + (k-1)*sig.e^2) / (1 + (k-1)*(sig.a^2 + sig.e^2)),
           -(k-1) * sig.a^2 * sig.e^2    / (1 + (k-1)*(sig.a^2 + sig.e^2)),
           -(k-1) * sig.a^2 * sig.e^2    / (1 + (k-1)*(sig.a^2 + sig.e^2)),
           sig.e^2 * (1 + (k-1)*sig.a^2) / (1 + (k-1)*(sig.a^2 + sig.e^2))
  )
  vcv.a[[k]] = matrix(vcvk, nrow = 1)
}
# might be smarter as data frame...?

popn = data.frame(
  # Unique identifier
  i   = 1:size0,
  # Generation
  gen = 0,
  # Age
  age = sample(size = size0, x = p.age$age, prob = p.age$p.age, replace = TRUE),
  # Sex (TRUE = female)
  fem = as.logical(sample(0:1, size0, replace = TRUE))
) %>%
  cbind(
    rmultinormal(size0, mean = c(0, 0), sigma = do.call(rbind, vcv.a[.$age + 1]))
  ) %>%
  rename(b_i = `1`, e_i = `2`)

ggplot(popn, aes(x = b_i, y = e_i)) +
  geom_point(aes(colour = age)) +
  scale_colour_viridis_c()

ggplot(popn, aes(x = b_i, y = e_i)) +
  geom_point(alpha = 0.1) +
  facet_wrap(~ age)

popn %>% group_by(age) %>% summarise(bbar = mean(b_i), ebar = mean(e_i), n = n())

# it's slow.

popn %>%
  mutate(
    # Phenotype
    z_i = b_i + e_i,
    # Survival
    s_i = s.max * exp(-(z_i - theta0)^2 / (2*wfitn^2)) * exp(-alpha * size0),
    # Offspring (0 for males, Poisson draw for females)
    r_i = rpois(size0, lambda = ifelse(fem, 2 * r, 0)),
    # Phenotypic optimum in this time step
    theta_t = theta0
) %>%
  ggplot(aes(x = age, y = z_i, colour = z_i - b_i)) +
  geom_point(position = position_jitter(width = 0.5))

# seems good...

init.sim = function(params, theta0) {
  
  # initial population size
  size0 = params$n.pop0
  # optimal phenotypic value
  s.max = params$s.max 
  # mean fecundity per individual
  r     = params$r
  # selection strength (this seems like a good value to use)
  wfitn = params$wfitn
  # genetic variance (should be approx. static)
  sig.a = params$sig.a
  # standard dev. of environmental phenoytpic noise
  sig.e = params$sig.e
  # initial genotype
  gbar0 = ifelse(any(grepl('gbar0', names(params))), params$gbar0, 0)
  # density dependence strength
  alpha = ifelse(any(grepl('alpha', names(params))), params$alpha, 0)
  # ceiling-like carrying capacity term
  kceil = ifelse(any(grepl('ceil' , names(params))), params$kceil, Inf)
  # equilibrium population growth rate 
  lstar = ifelse(any(grepl('lstar', names(params))), params$lstar, 1)
  
  # Get age distribution
  p.age = data.frame(age = 0:500) %>%
    mutate(p.age = (r / (1 + r)) * (s.max / lstar)^age * sqrt(wfitn^2 / (wfitn^2 + age*(sig.a^2 + sig.e^2))))
  
  # Get phenotpyic age variance/covariance matrices
  vcv.a = vector('list', length = nrow(p.age))
  
  for (k in 1:length(vcv.a)) {
    vcvk = c(sig.a^2 * (1 + (k-1)*sig.e^2) / (1 + (k-1)*(sig.a^2 + sig.e^2)),
             -(k-1) * sig.a^2 * sig.e^2    / (1 + (k-1)*(sig.a^2 + sig.e^2)),
             -(k-1) * sig.a^2 * sig.e^2    / (1 + (k-1)*(sig.a^2 + sig.e^2)),
             sig.e^2 * (1 + (k-1)*sig.a^2) / (1 + (k-1)*(sig.a^2 + sig.e^2))
    )
    vcv.a[[k]] = matrix(vcvk, nrow = 1)
  }
  
  popn = data.frame(
    # Unique identifier
    i   = 1:size0,
    # Generation
    gen = 0,
    # Age
    age = sample(size = size0, x = p.age$age, prob = p.age$p.age, replace = TRUE),
    # Sex (TRUE = female)
    fem = as.logical(sample(0:1, size0, replace = TRUE))
  ) %>%
  cbind(
    rmultinormal(size0, mean = rep(gbar0, 2), sigma = do.call(rbind, vcv.a[.$age + 1]))
  ) %>%
  rename(b_i = `1`, e_i = `2`) %>%
  mutate(
      # Phenotype
      z_i = b_i + e_i,
      # Survival
      s_i = s.max * exp(-(z_i - theta0)^2 / (2*wfitn^2)) * exp(-alpha * size0),
      # Offspring (0 for males, Poisson draw for females)
      r_i = rpois(size0, lambda = ifelse(fem, 2 * r, 0)),
      # Phenotypic optimum in this time step
      theta_t = theta0
    ) %>%
  select(-e_i)
  
  # If initial size is above the specified ceiling, truncate
  if (nrow(popn) > kceil) { 
    popn = popn %>% 
      slice_sample(kceil, replace = FALSE) 
    # Sort rows to order individuals by index
    arrange(i)
  }
  
  return(popn)
  
}

# ooh,... need to figure out mutational variance...

# Trials per parameter combo
trys.per = 5
# 
pars = expand.grid(
  p0    = c(.55, .75),
  l.rat = c(.9, .95, .975),
  h2    = c(.25, .5, 1)
) %>%
  # Demographic rates
  mutate(
    l.max = 2,
    # Steady state size of newborn cohort
    s.max = l.max * (1 - p0),
    # Equilibrium growth rate
    lstar = l.rat * l.max,
    # Mean fecundity
    r     = p0 / (1-p0),
    # Initial population size
    n.pop0 = 1000,
    # PStrength of density dependence
    alpha = log(lstar) / n.pop0,
    # Ceciling-type carrying capacity just in case
    kceil = 1200
  ) %>%
  filter(s.max <= 1) %>%
  # Genetic info
  group_by(p0, lstar, s.max, h2) %>%
  mutate(
    # Gamma-parameterization
    wfitn = 1,
    sig.0 = sqrt(newt.method.g1(.1, 1e-8, s.max / lstar, r)),
    sig.a = sqrt(h2 * sig.0),
    sig.e = sqrt((1-h2) * sig.0),
    sig.p = sqrt(gamma.calc(sig.a^2, s.max / lstar, r)),
    mu    = 1,
    sig.m = sqrt(wfitn^2 * (sig.a^2 - (sig.p^2 - p0*sig.a^2)/(1-p0)))
  ) %>%
  ungroup() %>%
  # Other junk
  mutate(
    timesteps = 10
  )

# Run sims

set.seed(4523)

sim.out2 = mclapply(
  pars %>% uncount(40) %>% mutate(try.no = 1:(nrow(.))) %>% split(.$try.no),
  function(pars) {
    sim(pars, theta.t = 0, init.rows = 5 * 1e5) %>%
      group_by(gen) %>%
      summarise(
        n = n(),
        kbar = mean(age),
        bbar = mean(b_i),
        bvar = var(b_i),
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

sim.out2 %>%
  ggplot(aes(x = gen, y = n, group = interaction(trial, lstar, h2))) +
  geom_line(aes(colour = h2), linewidth = 0.5) +
  facet_wrap(paste('p0', p0) ~ paste('lstar', lstar), nrow = 2) +
  scale_y_log10()
# nope... something's wrong

sim.out2 %>%
  ggplot(aes(x = gen, y = zbar, group = interaction(trial, lstar, h2))) +
  geom_line(aes(colour = h2), linewidth = 0.5) +
  facet_wrap(paste('p0', p0) ~ paste('lstar', lstar), nrow = 2) 
# ooks okay actually

sim.out2 %>%
  ggplot(aes(x = gen, y = zvar, group = interaction(trial, lstar, h2))) +
  geom_line(aes(colour = h2), linewidth = 0.5) +
  facet_wrap(paste('p0', p0) ~ paste('lstar', lstar), nrow = 2) 
# looks fine actually!

sim.out2 %>%
  group_by(gen, h2, lstar, p0) %>%
  summarise(
    bvarbar = mean(bvar),
    bvarvar = var(bvar),
    zvarbar = mean(zvar),
    zvarvar = var(zvar),
    n = n()
  ) %>%
  ggplot(aes(x = gen,  group = interaction(lstar, h2))) +
  geom_line(
    aes(
      y = bvarbar,
      colour = h2
      )
  ) +
  geom_ribbon(
    aes(
      ymin = bvarbar - 2 * sqrt(bvarvar / n),
      ymax = bvarbar + 2 * sqrt(bvarvar / n),
      group = h2,
      fill = h2
    ),
    alpha = 0.2
  ) +
  facet_wrap(paste('p0', p0) ~ paste('lstar', lstar), nrow = 2) 
# hmm... some very slight changes...

sim.out2 %>%
  group_by(gen, h2, lstar, p0) %>%
  summarise(
    bvarbar = mean(bvar),
    bvarvar = var(bvar),
    zvarbar = mean(zvar),
    zvarvar = var(zvar),
    n = n()
  ) %>%
  ggplot(aes(x = gen,  group = interaction(lstar, h2))) +
  geom_line(
    aes(
      y = zvarbar,
      colour = h2
    )
  ) +
  geom_ribbon(
    aes(
      ymin = zvarbar - 2 * sqrt(zvarvar / n),
      ymax = zvarbar + 2 * sqrt(zvarvar / n),
      group = h2,
      fill = h2
    ),
    alpha = 0.2
  ) +
  facet_wrap(paste('p0', p0) ~ paste('lstar', lstar), nrow = 2) 
# hmm... phenotypic variance does change


sim.out2 %>%
  ggplot(aes(x = gen, y = bvar, group = interaction(trial, lstar, h2))) +
  geom_line(aes(colour = h2), linewidth = 0.5) +
  facet_wrap(paste('p0', p0) ~ paste('lstar', lstar), nrow = 2) 
# also looks fine...

# hmm... so what's the issue...?
test.sim = sim(pars[1,], 0, init.rows = 5 * 1e6)

test.sim %>% group_by(gen) %>% summarise(n = n()) %>% plot(type = 'l')
# hmm...

test.sim %>%
  group_by(gen) %>%
  summarise(zbar = mean(z_i)) # looks fine

test.sim %>%
  group_by(gen) %>%
  summarise(zbar = var(z_i)) %>%
  ggplot(aes(x = gen, y = zbar)) +
  geom_line()
# hmm... randomness or no?
