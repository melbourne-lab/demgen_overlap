##### Script for testing simulation code (model 1 - single stage
#####

##### Source code
source('model_source/sim_model1_functions.R')

##### First - does it run at all

pars = data.frame(
  n.pop0 = 100, s.max = 0.9, r = (1.2 / (0.9)) - 1, wfitn = 2,
  sig.a = sqrt(0.5), sig.e = sqrt(0.5), alpha = 0.0001,
  timesteps = 10
)

set.seed(2090100)

# Try one simulation
test.sim = sim(params = pars, theta.t = 2, init.rows = 1e5)
print(paste('Plain simulation run for', pars$timesteps, 'generations and',
            length(unique(test.sim$i)), 'individuals'))

### Look at outputs of this simulation

## Look at population size over time
table(test.sim$gen)
## Look at one longer-lived individual
test.sim %>%
  filter(i %in% 9)
# yes - fecundity changes in each time step! nice.
# s_i changes because of changing population size
## Look at variance in r_i to see if it's actually updating
test.sim %>%
  group_by(i) %>%
  filter(n() > 1, fem) %>%
  summarise(var.r = var(r_i))
# Okay - some variance. But I suppose a lot of reproductive failure.
## Check sex ratios
test.sim %>%
  group_by(gen) %>%
  summarise(p.fem = mean(fem))
# interestingly biased high... possibly due to high survival that makes uneven
# sex ratios persistent? this would be interesting.
## Check reproductive output in each generation
test.sim %>%
  group_by(gen) %>%
  summarise(offspr = sum(r_i),
            mean.r = mean(r_i))
# not sure why it's higher in the first generation - sampling error? sex ratio?
## Look at mean genotype over time
test.sim %>%
  group_by(gen) %>%
  summarise(b.bar  = mean(b_i),
            b.bar0 = mean(b_i[age == 1]),
            b.bar1 = mean(b_i[age >  1]))
# Adaptation is present!

##### Run a simulation with changing theta

pars = data.frame(
  n.pop0 = 100, s.max = 0.9, r = (1.2 / (0.9)) - 1, wfitn = 2,
  sig.a = sqrt(0.5), sig.e = sqrt(0.5), alpha = 0.0001,
  timesteps = 10
)

# Here: theta fluctuates randomly, no autocorrelation
theta.t = 2 + rnorm(11, 0, 0.4)

set.seed(970200)

test.sim = sim(params = pars, theta.t = theta.t, init.rows = 1e5)

print(paste('Rnorm theta simulation run for', pars$timesteps, 'generations and',
            length(unique(test.sim$i)), 'individuals'))

### Some tests

## Population size over time
table(test.sim$gen)
## Is theta updating over time?
test.sim %>%
  distinct(gen, theta_t)
# it is good
## What's happening to the rate of adaptation
test.sim %>%
  group_by(gen) %>%
  summarise(b.bar = mean(b_i),
            d.bar = mean(b_i - theta_t))
# yes - adaptation still occurring although in fits and starts
## What is happening to survival?
test.sim %>% filter(i %in% 51)
# does appear to be changing - cool!

##### Try very strong density dependence on survival

pars = data.frame(
  n.pop0 = 100, s.max = 0.9, r = (1.2 / (0.9)) - 1, wfitn = 2,
  sig.a = sqrt(0.5), sig.e = sqrt(0.5), alpha = 0.005,
  timesteps = 10
)

# Here: theta fluctuates randomly, no autocorrelation
theta.t = 2

set.seed(534337)

test.sim = sim(params = pars, theta.t = theta.t, init.rows = 1e5)

print(paste('Strong NDD simulation run for', pars$timesteps, 'generations and',
            length(unique(test.sim$i)), 'individuals'))

###  Tests
# Population size over time
table(test.sim$gen)
# Incredibly strong! Very unrealistic influences of density dependence no?
test.sim %>%
  group_by(gen) %>%
  summarise(n = n(),
            sbar = mean(s_i))
# yes - survival is incredibly low for the first several timesteps!
test.sim %>%
  group_by(gen) %>%
  summarise(bbar = mean(b_i))
# although adaptation does appear to happen quickly! surprising...
test.sim %>%
  group_by(gen) %>%
  summarise(r = mean(s_i * r_i * 2))

