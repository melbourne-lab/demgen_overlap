##### Script for testing simulation code (model 1 - single stage

### Source code
source('model_source/sim_model1_functions.R')

### First - does it run at all

pars = data.frame(
  n.pop0 = 100, s.max = 0.9, r = (1.2 / (0.9)) - 1, wfitn = 3,
  sig.a = sqrt(0.5), sig.e = sqrt(0.5), alpha = 0.0001,
  timesteps = 10
)

test.sim = sim(pars, theta_t = 3, init.rows = 1e5)

# still buggy lol curse you namespace!