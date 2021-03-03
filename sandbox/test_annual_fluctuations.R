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

sigma.fluc = 0.5
# Has mean fitness of 0.843

# (could also us 1 and 1.7)

### Define sims

n.trials = 400

an.stat.list = vector('list', n.trials)

an.fluc.list = vector('list', n.trials)

set.seed(505)

for (trial in 1:n.trials) {
  
  an.stat.list[[trial]] = sim1(params = pars, theta_t = theta.stat)
  
  theta.fluc = rnorm(pars$end.time, theta.stat, sd = sigma.fluc)
  
  an.fluc.list[[trial]] = sim1(params = pars, theta_t = theta.fluc)
  
  
  print(trial)
  
}


anstat = unroller(an.stat.list)
anfluc = unroller(an.fluc.list)

anstat %>%
  group_by(trial, gen) %>%
  summarise(n = n()) %>%
  ggplot() +
  geom_line(aes(x = gen, y = n,  group = trial),
            size = 0.25) +
  scale_y_log10()

ansizes = rbind(
  anstat %>% mutate(fluc = FALSE),
  anfluc %>% mutate(fluc = TRUE)
) %>%
  group_by(trial, gen, fluc) %>%
  summarise(n = n(),
            th = theta[1]) %>%
  ungroup()

ansizes %>%
  filter(trial %in% 1:40) %>%
  ggplot(aes(x = gen, y = n)) +
  geom_line(aes(group = interaction(trial, fluc), colour = fluc), size = 0.5) +
  facet_wrap(~ fluc)

anmeans = ansizes %>%
  arrange(trial, desc(gen), fluc) %>%
  distinct(trial, fluc, .keep_all = TRUE) %>%
  filter(gen < pars$end.time) %>%
  uncount(weights = pars$end.time - gen) %>%
  group_by(trial, fluc) %>%
  mutate(gen = gen + (1:(pars$end.time-gen[1])),
        n = 0) %>%
  ungroup() %>%
  rbind(ansizes) %>%
  group_by(gen, fluc) %>%
  summarise(nbar = mean(n),
            nvar = var(n))

anmeans %>%
  ggplot(aes(x = gen, y = nbar)) +
  geom_line(aes(group = fluc, colour = fluc)) +
  geom_ribbon(aes(ymin = nbar - 2 * sqrt(nvar/n.trials), 
                  ymax = nbar + 2*sqrt(nvar/n.trials),
                  fill = fluc),
              alpha = 0.2)
