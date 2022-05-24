library(ggplot2)
library(dplyr)
library(tidyr)

rm(list = ls())

source('model_source/sim_model1_functions.R')

pars = expand.grid(
  trial.no = 1:20,
  sig.a    = c(sqrt(0.25), sqrt(0.5)),
  sig.e    = c(0, sqrt(0.5)),
  s.max    = c(0.7, 0.9)
) %>%
  mutate(
    n.pop0 = 2000,
    wfitn  = 2,
    alpha  = 0,
    timesteps = 25,
    r = (1 / s.max) * sqrt((wfitn^2 + sig.a^2 + sig.e^2) / wfitn^2) - 1
  )

liszt = vector('list', nrow(pars))

set.seed(50522105)

for (trial in 1:nrow(pars)) {
  liszt[[trial]] = sim(params = pars[trial,], init.rows = 1e5, theta.t = 0) %>%
    group_by(gen) %>%
    summarise(
      n = n(),
      var.b = var(b_i),
      var.z = var(z_i)
    ) %>%
    mutate(trial.no = trial)
  print(trial)
}

all.trials = merge(
  x = do.call(rbind, liszt),
  y = pars %>% mutate(trial.no = 1:nrow(.)) %>% select(trial.no, sig.a, sig.e, s.max)
)

head(all.trials)

all.trials %>%
  ggplot(aes(x = gen, y = n)) +
  geom_line(
    aes(
      group = trial.no,
      colour = factor(sig.a),
      )
    ) +
  facet_wrap(sig.e ~ s.max)
# I'm confused by this...
# ah... the variance change influences the growth rates (dur)

all.trials %>%
  ggplot(aes(x = gen, y = var.b)) +
  geom_line(
    aes(
      group = trial.no,
      colour = factor(sig.a),
    )
  ) +
  facet_wrap(sig.e ~ s.max)
# seems like equilibration is happening by ~gen 5

all.trials %>%
  ggplot(aes(x = gen, y = var.z)) +
  geom_line(
    aes(
      group = trial.no,
      colour = factor(sig.a),
    )
  ) +
  facet_wrap(sig.e ~ s.max)

var.sum = all.trials %>%
  select(sig.a, sig.e, s.max, gen, var.b, var.z) %>%
  gather(key = vartype, value = var, -c(sig.a, sig.e, s.max, gen)) %>%
  group_by(s.max, sig.a, sig.e, gen, vartype) %>%
  summarise(var.bar = mean(var),
            var.var = var(var),
            n = n())

head(var.sum)

var.sum %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = var.bar,
      colour = s.max,
      group = interaction(sig.a, sig.e, s.max, vartype)
    )
  ) +
  geom_ribbon(
    aes(
      ymin = var.bar - 2 * sqrt(var.var / n),
      ymax = var.bar + 2 * sqrt(var.var / n),
      group = interaction(sig.a, sig.e, s.max, vartype),
      fill = s.max
    ),
    alpha = 0.1
  ) +
  facet_wrap(sig.e + sig.a ~ vartype, nrow = 4) +
  theme(legend.position = 'bottom')

# that's very small s.e.!!!

# what is our prediction??

# pars %>%
#   distinct(sig.a, sig.e, s.max) %>%
#   mutate(pred.var = )