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
# maybe a little biased low?

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
  geom_point() +
  geom_segment(aes(x = 0.5, xend = 1.5, y = 0.5, yend = 1.5),
               size = 0.1)

test.growth.rates %>%
  ggplot(aes(x = elam, y = elam - lambda)) +
  geom_point()

# Looks fine to me. Think this is a reliable unbiased analytical mean

### Next up would probably be to look at variance in growth rate?

