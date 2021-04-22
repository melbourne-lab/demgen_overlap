

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
    select(gen, n, trial, age, harsh) %>%
    arrange(trial, desc(gen), age, harsh) %>% 
    distinct(trial, age, harsh, .keep_all = TRUE) %>%
    filter(gen < mxgen) %>%
    uncount(weights = mxgen - gen) %>%
    mutate(n = 0) %>%
    group_by(trial, age, harsh) %>%
    mutate(gen = gen + (1:(mxgen - gen[1]))) %>%
    ungroup() %>%
    rbind(n.summary %>% select(gen, n, trial, age, harsh)) %>%
    arrange(trial, age, gen, harsh)
  
}

# Read in data
fluct = rbind(
  read.csv('run_sims/stationary/stationary_out.csv') %>% 
    mutate(theta = 2.5) %>%
    select(gen, n, gbar, zbar, wbar, theta, p.fix.pos, p.fix.neg, v, trial, age) %>%
    mutate(harsh = 0),
  read.csv('run_sims/fluctuating/fluctuating_less_harsh_out.csv') %>% 
    mutate(harsh = 1),
  read.csv('run_sims/fluctuating/fluctuating_harsh_out.csv') %>% 
    mutate(harsh = 2)
)

### Summarise population size

fluct.n = add.zeros(fluct, 30)

fluct.n.summ = fluct.n %>%
  group_by(gen, age, harsh) %>%
  summarise(nbar = mean(n),
            nvar = var(n),
            n = n()) %>%
  ungroup() %>%
  mutate(age.fac = factor(age, labels = c('Annual', 'Short-lived perennial', 'Long-lived perennial')),
         env.var = factor(harsh, labels = c('No variance', 'Low variance', 'High variance')))
  
fluct.n.summ %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = nbar,
      group = interaction(age, harsh),
      colour = factor(age),
      linetype = factor(harsh)
    )
  ) + 
  scale_y_log10()

# More serious plot:
fluct.n.summ %>%
  ggplot(aes(x = gen)) +
  geom_ribbon(
    aes(
      ymin = nbar - 2 * sqrt(nvar / n),
      ymax = nbar + 2 * sqrt(nvar / n),
      group = interaction(age, harsh),
      fill = age.fac
    ),
    alpha = 0.2
  ) +
  geom_line(
    aes(
      y = nbar,
      group = interaction(age, harsh),
      colour = age.fac
    )
  ) +
  scale_y_log10() +
  facet_wrap(~ env.var) +
  guides(colour = guide_legend(''), fill = guide_legend('')) +
  labs(x = 'Time step', y = 'Population size') +
  theme(
    panel.background = element_rect(fill = 'white'),
    panel.grid = element_line(colour = 'gray88'),
    legend.position = 'bottom'
  )

# ope... looks like an old and outdated version of the novr script is hres

### Rate of adaptation

fluct.evo = fluct %>%
  mutate(g = gbar) %>%
  group_by(gen, age, harsh) %>%
  filter(n() %in% 1000) %>%
  summarise(gbar = mean(g),
            gvar = var(g),
            n = n(),
            theta = mean(theta))

fluct.evo %>%
  ggplot(aes(x = gen)) +
  # geom_line(
  #   aes(
  #     y = theta
  #   )
  # ) +
  geom_line(
    aes(
      y = gbar,
      group = age,
      colour = factor(age)
    )
  ) +
  facet_wrap(~ harsh, ncol = 1)

fluct.ext = fluct.n %>%
  group_by(harsh, age, trial) %>%
  summarise(extinct = any(!n)) %>%
  group_by(harsh, age) %>%
  summarise(p = mean(extinct))

fluct.ext %>%
  ggplot() +
  geom_point(
    aes(
      x = factor(age),
      y = p,
      colour = factor(age)
    ),
    size = 3
  ) +
  geom_segment(
    aes(
      x    = factor(age),
      xend = factor(age),
      y    = p - 2 * sqrt(p * (1-p) / 1000),
      yend = p + 2 * sqrt(p * (1-p) / 1000),
      colour = factor(age)
    )
  ) +
  facet_wrap(~ harsh) +
  theme(legend.position = 'none')
