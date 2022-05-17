# After talking with Dan Doak (5/13), I wondered about the effects of selection
# on the genetic variance over time.
# Lynch and Walsh's demonstration that additive genetic variance and breeding
# value variance are the same may have been predicated on the absence of
# selection. Selection should narrow phenotypic variance (and phenotypic
# variance is restored by mutations).
# Here, I am simulating breeding value and phenotypic variance over time in the
# absence of selection (w = infinity) to see if the variance changes over time.
# SN - 17 May 2022

# Load packages
library(ggplot2)
library(dplyr)
library(tidyr)

# Clear namespace
rm(list = ls())

# Get model source code
source('model_source/sim_model1_functions.R')

# Get parameters for simulation
pars = data.frame(trial = 1:1000) %>%
  mutate(n.pop0 = 1000, 
         # No selection occurring
         wfitn = Inf,
         s.max = 0.7,
         sig.a = sqrt(0.25), 
         sig.e = sqrt(0.25),
         # No density dependence
         alpha = 0,
         timesteps = 20,
         # Because there is no variance in survival, set this to have expected
         # growth rate 1 without accounting for heterogeneity
         r = (1 / s.max) - 1)

# Run simulations
# For each trial, store breeding value and phenotypic variance for
# offspring (age 0) and adults (age > 0)

liszt = vector('list', nrow(pars))

set.seed(5521015)

for (trial in pars$trial) {
  
  liszt[[trial]] = sim(params = pars[trial,], theta.t = 0, init.rows = 1e5) %>%
    mutate(age = ifelse(!gen, 0, age)) %>%
    group_by(gen, adult = age > 1) %>%
    summarise(varb = var(b_i),
              varz = var(z_i)) %>%
    mutate(trial.no = trial)
  
  print(trial)
  
}

sapply(liszt, nrow)

all.trials = do.call(rbind, liszt)

# Here: get all variances (turn to long form)
all.vars = all.trials %>%
  gather(key = bz, value = var.val, -c(trial.no, gen, adult))

# Visualize a couple of trials
all.vars %>%
  filter(trial.no < 26) %>%
  ggplot(aes(x = gen, y = var.val)) +
  geom_line(aes(colour = bz, linetype = adult)) +
  facet_wrap(~ trial.no)
# Looks stable (with expected degree of variance),
# adults are more stable than offspring, it seems
# (which makes sense - should be autocorrelation here)
# also looks like no consistent gap betwen parents and adults
  
# Get a summary of variances across trials
all.vars.summ = all.vars %>%
  group_by(gen, adult, bz) %>%
  summarise(var.bar = mean(var.val),
            var.var = var(var.val),
            n = n())

head(all.vars.summ)

# Plot variances
all.vars.summ %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = var.bar,
      colour = adult
    )
  ) +
  geom_ribbon(
    aes(
      ymin = var.bar - 2 * sqrt(var.var / n),
      ymax = var.bar + 2 * sqrt(var.var / n),
      fill = adult
    ),
    alpha = 0.25
  ) +
  facet_wrap(~ bz)
# Looks about stable over time, no gap

# Here - explicitly plot a gap
# (paired adults and offspring gaps within a trial)
all.gaps = all.vars %>%
  filter(gen > 0) %>%
  mutate(adult = ifelse(adult, 'adult', 'offsp')) %>%
  spread(adult, var.val) %>%
  mutate(gap = adult - offsp)

# Summarize over all trials
all.gaps.summ = all.gaps %>%
  group_by(gen, bz) %>%
  summarise(gap.bar = mean(gap),
            gap.var = var(gap),
            n = n())

# Plot
all.gaps.summ %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = gap.bar,
      colour = bz
    )
  ) +
  geom_ribbon(
    aes(
      ymin = gap.bar - 2 * sqrt(gap.var / n),
      ymax = gap.bar + 2 * sqrt(gap.var / n),
      fill = bz
    ),
    alpha = 0.25
  ) +
  facet_wrap(~ bz)

# Hmm... *might* be above zero, but this could be sampling variance
# also seems to return to below zero eventually

# Okie dokie - this makes sense Seems like the result of breeding value variance
# equaling segregational variance (or being proportional - one half?) holds if
# there is no selection and breeding value variance will change with selection
