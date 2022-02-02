# Script for looking at population growth rates (lambda) to get analytical
# expectation vs.simulation results.

library(ggplot2)

rm(list = ls())

##### Source code
source('model_source/sim_model1_functions.R')

##### Look at just one simulation run

pars = data.frame(
  n.pop0 = 100, s.max = 0.9, r = (1.2 / (0.9)) - 1, wfitn = 2,
  sig.a = sqrt(0.5), sig.e = sqrt(0.5), alpha = 0.0001,
  timesteps = 30
)

set.seed(2090100)

test.sim = sim(params = pars, theta.t = 2.5, init.rows = 1e5)

table(test.sim$gen)

test.sim %>%
  ggplot(aes(x = gen, y = (b_i - theta_t))) +
  geom_point(position = position_jitter(width = 0.25)) +
  geom_line(aes(group = i), size = 0.1)

# Population growth rate is lambda column

test.sim.lambda = test.sim %>%
  group_by(gen) %>%
  summarise(n = n()) %>%
  mutate(lambda = exp(c(diff(log(n)), NA))) %>%
  filter(!is.na(lambda))

ggplot(test.sim.lambda, aes(x = gen, y = lambda)) +
  geom_line() +
  geom_point()

# Think that expectation should be s (1 + r)

test.sim.sr = test.sim %>%
  group_by(gen) %>%
  summarise(s2r = mean(s_i * (1 + r_i)))

merge(test.sim.lambda, test.sim.sr, all.y = FALSE) %>%
  ggplot(aes(x = s2r, y = lambda)) +
  geom_point() +
  geom_segment(aes(x = 0.5, xend = 1.3, y = 0.5, yend = 1.3),
               size = 0.1)
# maybe a little biased high...

# (n.b. is nonlinear averaging occurring here?)

test.sim.sr = test.sim %>%
  mutate(s1r = s_i * (1 + r_i)) %>%
  group_by(gen) %>%
  summarise(s2r = mean(s_i * (1 + r_i)),
            s1r = mean(s1r))

ggplot(test.sim.sr) +
  geom_segment(aes(x = 0.5, xend = 1.5, y = 0.5, yend = 1.5),
               size = 0.1, linetype = 2) +
  geom_point(aes(x = s1r, y = s2r))

# great - no non-linear averaging

##### Run this on a batch of simulations

li = vector('list', 20)

set.seed(210182197)

for (k in 1:length(li)) li[[k]] = sim(pars, 2.5, 2 * 1e4) %>% mutate(trial = k)

mi = do.call(rbind, li)

test.growth.rates = mi %>%
  group_by(trial, gen) %>%
  summarise(n = n(),
            elam   = mean(s_i * (1 + r_i))) %>%
  mutate(lambda = exp(c(diff(log(n)), NA))) %>%
  filter(!is.na(lambda))

test.growth.rates %>%
  ggplot(aes(x = elam, y = lambda)) +
  geom_point(aes(colour = log(n)), size = 2) +
  geom_segment(aes(x = 0.5, xend = 1.5, y = 0.5, yend = 1.5),
               size = 0.1) +
  scale_color_viridis_c()

test.growth.rates %>%
  ggplot(aes(x = elam, y = elam - lambda)) +
  geom_point(aes(colour = log(n)))

# Hmm... there is a weird-looking pattern in residuals...
# Not sure if this is an artifact or not.
# I think it's mostly due to sampling variation at small population size.

test.growth.rates %>%
  ggplot(aes(x = n, y = elam - lambda)) +
  geom_point(aes(colour = log(n))) +
  scale_x_log10()
# Yes. Possibly biased low at small size but this is probably also just sampling
# variance

### Next up would probably be to look at variance in growth rate?

# Run some simulations with _zero_ genotypic variance

li = vector('list', 20)

set.seed(53528535)

pars = pars %>% mutate(sig.a = 0, sig.e = 0, alpha = 0)

for (k in 1:length(li)) li[[k]] = sim(pars, 0, 2 * 1e4) %>% mutate(trial = k)

mi = do.call(rbind, li)

test.growth.rates = mi %>%
  group_by(trial, gen) %>%
  summarise(n = n(),
            f = sum(fem),
            elam   = mean(s_i * (1 + r_i))) %>%
  mutate(lambda = exp(c(diff(log(n)), NA))) %>%
  filter(!is.na(lambda))

test.growth.rates %>% ggplot(aes(x = gen, y = n, group = trial)) + geom_line() + scale_y_log10()

# hmm... that's not right... is it?

li[[1]]
# good

head(test.growth.rates)

test.growth.rates %>%
  ggplot(aes(x = elam, y = lambda)) +
  geom_point(aes(colour = log(n)), size = 2) +
  geom_segment(aes(x = 1.0, xend = 1.5, y = 1.0, yend = 1.5),
               size = 0.1) +
  scale_color_viridis_c()

### Look at variance in growth rates

set.seed(116127)

init.known.f = init.sim(pars, theta0 = 0)

n    = nrow(init.known.f)
sbar = mean(init.known.f$s_i)
rbar = mean(init.known.f$r_i)
nf   = sum(init.known.f$fem)

# # one test to look at
# n.tp1.known.f = propagate.sim(init.known.f, params = pars, theta = 0)

# n.known.f = vector(mode = 'double', length = 900)
trials.known.f = vector(mode = 'list', length = 900)

for (k in 1:900) trials.known.f[[k]] = propagate.sim(init.known.f, params = pars, theta = 0)

n.known.f = sapply(trials.known.f, nrow)

mean(n.known.f)
# compare with
n * sbar + nf * sbar * 2 * rbar
# fine once accounting for the 2

var(n.known.f)
n * sbar * (1 - sbar) + nf * (2 * rbar) * sbar * (1 - (2 * rbar) * (1 - sbar))
# 2 in the first term would help maybe...

# What is variance in number of offspring?
kids.known.f = sapply(trials.known.f, function(x) nrow(x %>% filter(gen > 0, age < 2)))

var(kids.known.f)
nf * (2 * rbar) * sbar * (1 - (2 * rbar) * (1 - sbar))
# this is actually pretty good

# So the issue may be with the survival?
ents.known.f = sapply(trials.known.f, function(x) nrow(x %>% filter(gen > 0, age > 1)))
ents.known.f
var(ents.known.f)
# no...

# erm, is there covariation? more surviving parents means more kids...
cov(kids.known.f, ents.known.f)
# uh oh
plot(kids.known.f, ents.known.f)

var(kids.known.f) + var(ents.known.f) + 2 * cov(kids.known.f, ents.known.f)
# shoot
# okay kiss our nice analytical solutions goodbye

# covariance estimate:
2 * nf * rbar * sbar * (1 - sbar)
# simulated covariance
cov(kids.known.f, ents.known.f)
# not that far off?

# analytical estimate incl. covariance
n * sbar * (1 - sbar) +
  nf * (2 * rbar) * sbar * (1 - (2 * rbar) * (1 - sbar)) +
  2 * nf * (2 * rbar) * sbar * (1 - sbar)
# that's actually pretty good.

### Try with varying $N_f$ just to see if conditioning on N_f works

n.trials = 4900

liszt = vector('list', n.trials)

pars = pars %>% mutate(timesteps = 1)

set.seed(20093)

for (k in 1:n.trials) liszt[[k]] = sim(params = pars, theta = 0, init.rows = 600) %>% mutate(trial = k)

all.trials = do.call(rbind, liszt)

head(all.trials)
tail(all.trials)

all.sum = all.trials %>%
  group_by(trial) %>%
  summarise(f0 = sum(fem[gen < 1]),
            n0 = 100,
            n1 = sum(gen > 0)) %>%
  ungroup()

table(all.sum$f0)

# let's try variances...

r = 2 * pars$r
s = pars$s.max

var.f = all.sum %>%
  group_by(f0) %>%
  summarise(expvar = n0*s*(1-s) + f0*s*r*(1+r*(1-s)) + 2 * f0*s*r*(1-s),
            obsvar = var(n1),
            n = n()) %>%
  distinct(f0, .keep_all = TRUE)

# am I a fucking idiot? why is it duping rows...
# anyway

var.f %>%
  ggplot(aes(x = f0)) +
  geom_line(aes(y = expvar), colour = 'blue') +
  geom_point(aes(y = obsvar, alpha = log(n)), colour = 'black')
# erm... kinda looks right...

# hmm... kinda looks right?

var.f %>%
  ggplot(aes(x = f0)) +
  geom_line(aes(y = 0), colour = 'blue') +
  geom_point(aes(y = obsvar - expvar, alpha = log(n)), colour = 'black')

# looks good to me!

# very cool... of course there will be variation around rbar, etc.

