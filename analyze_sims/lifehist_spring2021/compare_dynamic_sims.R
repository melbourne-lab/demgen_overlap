

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
    select(gen, n, trial, age, varn, speed) %>%
    arrange(trial, desc(gen), age, varn, speed) %>% 
    distinct(trial, age, varn, speed, .keep_all = TRUE) %>%
    filter(gen < mxgen) %>%
    uncount(weights = mxgen - gen) %>%
    mutate(n = 0) %>%
    group_by(trial, age, varn, speed) %>%
    mutate(gen = gen + (1:(mxgen - gen[1]))) %>%
    ungroup() %>%
    rbind(n.summary %>% select(gen, n, trial, age, varn, speed)) %>%
    arrange(trial, age, gen, varn, speed)
  
}

# Read in data
dynamic = rbind(
  read.csv('run_sims/dynamic/dynamic_novar_out.csv') %>% mutate(varn = 0, speed = 2),
  read.csv('run_sims/dynamic/dynamic_lovar_out.csv') %>% mutate(varn = 1, speed = 2),
  read.csv('run_sims/dynamic/dynamic_hivar_out.csv') %>% mutate(varn = 2, speed = 2),
  read.csv('run_sims/dynamic/dynamic_slow_novar_out.csv') %>% mutate(varn = 0, speed = 1),
  read.csv('run_sims/dynamic/dynamic_slow_lovar_out.csv') %>% mutate(varn = 1, speed = 1),
  read.csv('run_sims/dynamic/dynamic_slow_hivar_out.csv') %>% mutate(varn = 2, speed = 1)
  
)

# Check... is the degree of variation actually different?
dynamic %>%
  distinct(trial, gen, varn, speed, .keep_all = TRUE) %>%
  group_by(trial, varn, speed) %>%
  summarise(env.var = var(theta)) %>%
  group_by(varn, speed) %>%
  summarise(sig2.bar = mean(env.var),
            sig2.var = var(env.var)) # k that's weird but w/e there appears to be a noticeable effect

# Summarise population size

dyn.n = add.zeros(dynamic, 30)

dyn.n.summ = dyn.n %>%
  group_by(gen, age, varn, speed) %>%
  summarise(nbar = mean(n),
            nvar = var(n),
            n = n()) %>%
  ungroup() %>%
  ungroup() %>%
  mutate(age.fac = factor(age, labels = c('Annual', 'Short-lived perennial', 'Long-lived perennial')),
         env.var = factor(varn, labels = c('No environmental variance',
                                           'Low environmental variance', 
                                           'High environmental variance')),
         speed.f = factor(speed, labels = c('Slow environmental change', 
                                            'Fast environmental change')))

dyn.n.summ %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = nbar,
      group = interaction(age, varn, speed),
      colour = factor(age),
      linetype = varn
    )
  )

dyn.n.summ %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = nbar,
      group = interaction(age, varn, speed),
      colour = factor(age.fac)
    ),
    size = 1.2
  ) +
  geom_ribbon(
    aes(
      ymin = nbar - 2 * sqrt(nvar / n),
      ymax = nbar + 2 * sqrt(nvar / n),
      fill = factor(age.fac),
      group = interaction(age, varn, speed)
    ),
    alpha = 0.2
  ) +
  scale_y_log10() +
  facet_wrap(speed.f ~ env.var) +
  guides(colour = guide_legend(''), fill = guide_legend('')) +
  labs(x = 'Time step', y = 'Population size') +
  theme(
    panel.background = element_rect(fill = 'white'),
    panel.grid = element_line(colour = 'gray88'),
    legend.position = 'none'
  ) +
  ggsave('analyze_sims/lifehist_spring2021/t2_popsize.pdf',
         width = 9, height = 7)
   

# Rate of adaptation

# compare only the populations that we have all trials for

dyn.f.evo = dynamic %>%
  mutate(g = gbar) %>%
  group_by(gen, age, varn, speed) %>%
  filter(n() %in% 1000) %>%
  summarise(gbar = mean(g),
            gvar = var(g),
            n = n(),
            theta = mean(theta)) %>%
  ungroup() %>%
  mutate(age.fac = factor(age, labels = c('Annual', 'Short-lived perennial', 'Long-lived perennial')),
         env.var = factor(varn, labels = c('No environmental variance',
                                           'Low environmental variance', 
                                           'High environmental variance')),
         speed.f = factor(speed, labels = c('Slow environmental change', 
                                            'Fast environmental change')))

dyn.f.evo %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = theta
    ),
    size = 2
  ) +
  geom_line(
    aes(
      y = gbar,
      group = age,
      colour = age.fac
    ),
    size = 1.5
  ) +
  labs(x = 'Generation', y = '') +
  facet_wrap(speed.f ~ env.var) +
  theme(
    legend.position = 'none',
    panel.background = element_rect(fill = 'white'),
    panel.grid = element_line(colour = 'gray88')
  ) +
  ggsave('analyze_sims/lifehist_spring2021/t2_evolvin.pdf',
         width = 9, height = 7)  

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

dyn.f.extinct = dyn.n %>% 
  group_by(varn, speed, age, trial) %>%
  summarise(extinct = any(!n)) %>%
  group_by(varn, speed, age) %>%
  summarise(p = mean(extinct)) %>%
  ungroup() %>%
  mutate(age.fac = factor(age, labels = c('Annual', 'Short-lived perennial', 'Long-lived perennial')),
         env.var = factor(varn, labels = c('No environmental variance',
                                           'Low environmental variance', 
                                           'High environmental variance')),
         speed.f = factor(speed, labels = c('Slow environmental change', 
                                            'Fast environmental change')))

dyn.f.extinct %>%
  ggplot() +
  geom_point(
    aes(
      x = age.fac,
      y = p,
      colour = factor(age)
    ),
    size = 3
  ) +
  geom_segment(
    aes(
      x    = age.fac,
      xend = age.fac,
      y    = p - 2 * sqrt(p * (1-p) / 1000),
      yend = p + 2 * sqrt(p * (1-p) / 1000),
      colour = factor(age)
    ),
    size = 1.2
  ) +
  facet_wrap(speed.f ~ env.var) +
  labs(x = '', y = 'Probability of extinction') +
  facet_wrap(speed.f ~ env.var) +
  theme(
    legend.position = 'none',
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.background = element_rect(fill = 'white'),
    panel.grid = element_line(colour = 'gray88')
  ) +
  ggsave('analyze_sims/lifehist_spring2021/t2_extinct.pdf',
         width = 9, height = 7)
