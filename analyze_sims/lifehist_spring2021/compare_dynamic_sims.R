

library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyselect)

# Read in wrapper aux functions
source('model_source/sim_aux_handling_functions.R')

# Function for adding zeros (to get mean population size)
add.zeros = function(n.summary, mxgen, group) {
  
  # note grouping variables involved
  
  n.summary %>%
    select(gen, n, trial, age, varn) %>%
    arrange(trial, desc(gen), age, varn) %>% 
    distinct(trial, age, varn, .keep_all = TRUE) %>%
    filter(gen < mxgen) %>%
    uncount(weights = mxgen - gen) %>%
    mutate(n = 0) %>%
    group_by(trial, age, varn) %>%
    mutate(gen = gen + (1:(mxgen - gen[1]))) %>%
    ungroup() %>%
    rbind(n.summary %>% select(gen, n, trial, age, varn)) %>%
    arrange(trial, age, gen, varn)
  
}

# Read in data
dynamic.fast = rbind(
  read.csv('run_sims/dynamic/dynamic_novar_out.csv') %>% mutate(varn = 0),
  read.csv('run_sims/dynamic/dynamic_lovar_out.csv') %>% mutate(varn = 1),
  read.csv('run_sims/dynamic/dynamic_hivar_out.csv') %>% mutate(varn = 2)
)

# Check... is the degree of variation actually different?
dynamic.fast %>%
  distinct(trial, gen, varn, .keep_all = TRUE) %>%
  group_by(trial, varn) %>%
  summarise(env.var = var(theta)) %>%
  group_by(varn) %>%
  summarise(sig2.bar = mean(env.var),
            sig2.var = var(env.var)) # k that's weird but w/e there appears to be a noticeable effect

# Summarise population size

dyn.f.n = add.zeros(dynamic.fast, 30)

dyn.f.n.summ = dyn.f.n %>%
  group_by(gen, age, varn) %>%
  summarise(nbar = mean(n),
            nvar = var(n),
            n = n()) %>%
  mutate(varn = factor(varn))

dyn.f.n.summ %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = nbar,
      group = interaction(age, varn),
      colour = factor(age),
      linetype = varn
    )
  )

dyn.f.n.summ %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = nbar,
      group = interaction(age, varn),
      colour = factor(age),
      linetype = varn
    )
  ) +
  geom_ribbon(
    aes(
      ymin = nbar - 2 * sqrt(nvar / n),
      ymax = nbar + 2 * sqrt(nvar / n),
      fill = factor(age),
      group = interaction(age, varn)
    ),
    alpha = 0.2
  ) +
  facet_wrap(~ varn, ncol = 1)
   
# Rate of adaptation

# compare only the populations that we have all trials for

dyn.f.evo = dynamic.fast %>%
  mutate(d = gbar - theta) %>%
  group_by(gen, age, varn) %>%
  filter(n() %in% 1000) %>%
  summarise(dbar = mean(d),
            dvar = var(d),
            n = n())

dyn.f.evo %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = dbar,
      colour = factor(age),
      linetype = factor(varn)
    )
  )

# Wow that's neat. The rate of adaptation doesn't seem to matter...

dyn.f.fit = dynamic.fast %>%
  mutate(r = wbar * exp(-0.005 * n),
         w = wbar) %>%
  group_by(gen, age, varn) %>%
  filter(n() %in% 1000) %>%
  summarise(rbar = mean(r),
            rvar = var(r),
            wbar = mean(w),
            wvar = var(w),
            n = n())

dyn.f.fit %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = rbar,
      colour = factor(age),
      linetype = factor(varn)
    )
  )

dyn.f.fit %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = wbar,
      colour = factor(age)
    )
  ) +
  geom_line(
    aes(
      y = rbar,
      colour = factor(age)
    ),
    linetype = 2
  ) +
  facet_wrap(~ varn, ncol = 1)
# I'm not sure this plot is helpful...

# Extinction probabilities, times to extinction.

dyn.f.extinct = dynamic.fast %>%
  group_by(trial, age, varn) %>%
  filter(max(gen) < 30) %>%
  summarise(ext.gen = max(gen))

dyn.f.extinct %>%
  ggplot() +
  geom_histogram(
    aes(
      x = ext.gen,
      fill = factor(age)
    ),
    position = 'identity',
    alpha = 0.4,
    binwidth = 1
  ) +
  facet_wrap(~ varn, ncol = 1)

# Histogram too...