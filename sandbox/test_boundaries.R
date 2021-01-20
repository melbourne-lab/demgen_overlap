library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyselect)

rm(list = ls())

source('model_source/sim_functions_annual.R')
source('model_source/sim_functions_biennial.R')

### Define global params

pars = data.frame(
  end.time = 10,
  init.row = 1e4,
  n.loci = 25, 
  n.pop0 = 40,
  w.max = 1.1,
  wfitn = sqrt(1 / 0.14),
  sig.e = 0.5
)

### Define sim params

theta.stat = 2
# Has mean initial fitness is 0.628

sigma.fluc = 2
# Has mean fitness of 0.843

# (could also us 1 and 1.7)

### Define sims

n.trials = 100

an.stat.list = vector('list', n.trials)
bi.stat.list = vector('list', n.trials)

an.fluc.list = vector('list', n.trials)
bi.fluc.list = vector('list', n.trials)

set.seed(505)

for (trial in 1:n.trials) {
  
  an.stat.list[[trial]] = sim1(params = pars, theta_t = theta.stat)
  bi.stat.list[[trial]] = sim2(params = pars, theta_t = theta.stat)
  
  theta.fluc = rnorm(pars$end.time, 0, sd = sigma.fluc)
  
  an.fluc.list[[trial]] = sim1(params = pars, theta_t = theta.fluc)
  bi.fluc.list[[trial]] = sim2(params = pars, theta_t = theta.fluc)
  
  
  print(trial)
  
}

# oh noes! got an error at trial 33
#  Error: Internal error: `arg_match()` expects a symbol
#  traceback seems to fail at 




