library(ggplot2)
library(dplyr)
library(tidyr)

rm(list = ls())

source('model_source/sim_model1_functions.R')

pars = expand.grid(
  trial.no = 1:40,
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
      colour = factor(s.max),
      group = interaction(sig.a, sig.e, s.max, vartype)
    )
  ) +
  geom_ribbon(
    aes(
      ymin = var.bar - 2 * sqrt(var.var / n),
      ymax = var.bar + 2 * sqrt(var.var / n),
      group = interaction(sig.a, sig.e, s.max, vartype),
      fill = factor(s.max)
    ),
    alpha = 0.1
  ) +
  facet_wrap(paste0(sig.e^2, ';', sig.a^2) ~ vartype, nrow = 4) +
  theme(legend.position = 'bottom')

# that's very small s.e.!!!

# summation wrapper
expr.sum = function(s, sig2, w2) {
  runsum = 0
  for (i in 1:1000) runsum = runsum + ( ((s^(i-1)) * (1 - s))^2 * (w2 * sig2/ (w2 + sig2*(i+1))))
  return(runsum)
}

expr.sum(0.7, 0.5, 2^2)
expr.sum(0.9, 0.5, 2^2)
expr.sum(0.7, 0.25, 2^2)
expr.sum(0.9, 0.25, 2^2)

# what is our prediction??

pred.var.sum = merge(
  x = var.sum,
  y = pars %>%
    distinct(sig.a, sig.e, s.max) %>%
    mutate(pred.var = expr.sum(s.max, sig.a^2, pars$wfitn[2]^2))
)

pred.var.sum %>%
  filter(gen < 8, vartype %in% 'var.b') %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = var.bar,
      group = interaction(sig.a, sig.e, s.max, vartype)#,
#      colour = vartype
    )
  ) + 
  geom_line(
    aes(
      y = pred.var,
      group = interaction(sig.a, sig.e, s.max)
    ),
    linetype = 2
  ) +
  facet_wrap(paste0(sig.e^2, ';', sig.a^2) ~ s.max, nrow = 4) +
  theme(legend.position = 'bottom')

# even after correcting, these are too low...

# hmm.... I guess there is covariance??

set.seed(214339)

test.sim = sim(params = pars[1,], init.rows = 1e5, theta.t = 0)

nrow(test.sim)

test.sim %>%
  # filter(i > pars$n.pop0[1]) %>%
  group_by(gen, age) %>%
  summarise(varb = var(b_i)) %>%
  spread(age, varb)

test.vars.age = test.sim %>%
  filter(i > pars$n.pop0[1]) %>%
  group_by(gen, age) %>%
  summarise(varb = var(b_i))

test.vars.age %>%
  ggplot(aes(x = age, y = varb)) +
  geom_line(aes(colour = gen, group = gen)) +
  scale_colour_viridis_c()

# oh dear god... is there negative covariance??

test.vars.age %>%
  filter(gen > 1) %>%
  ggplot(aes(x = age, y = varb)) +
  geom_line(aes(group = gen)) +
  facet_wrap(~ gen)

test.vars.age %>%
  ggplot(aes(x = age, y = varb)) +
  geom_line(aes(group = gen)) +
  facet_wrap(~ gen)

test.vars.age %>%
  mutate(cohort = gen - age) %>%
  ggplot(aes(x = cohort, y = varb)) +
  geom_line(aes(group = gen, colour = gen)) +
  scale_color_viridis_c(option = 'B') +
  theme(panel.background = element_rect('gray66'),
        panel.grid = element_blank()) #+
  #facet_wrap(~ gen)

test.vars.age %>%
  mutate(cohort = gen - age) %>%
  ggplot(aes(x = gen, y = varb)) +
  geom_line(aes(group = cohort)) +
  # scale_color_viridis_c(option = 'B') +
  # theme(panel.background = element_rect('gray66'),
  #       panel.grid = element_blank()) +
  facet_wrap(~ cohort)

# maybe look at age breakdown?
test.age.dist = test.sim %>%
  group_by(gen, age) %>%
  summarise(n = n()) %>%
  group_by(gen) %>%
  mutate(p = n / sum(n)) %>%
  ungroup()

test.age.dist %>%
  ggplot(aes(x = age, y = p)) +
  geom_line(aes(group = gen))

test.age.dist.off = test.sim %>%
  filter(i > pars$n.pop0[1]) %>%
  group_by(gen, age) %>%
  summarise(n = n()) %>%
  group_by(gen) %>%
  mutate(p = n / sum(n)) %>%
  ungroup()

test.age.dist.off %>%
  filter(gen > 1) %>%
  ggplot(aes(x = age, y = p)) +
  geom_line(aes(group = gen)) +
  geom_line(aes(y = (pars$s.max[1]^(age-1)) * (1 - pars$s.max[1])),
            colour = 'blue', linetype = 2)

#facet_wrap(~ gen)

### Try a long sim to look at the age distribution...

set.seed(214339)

test.long = sim(params = pars[1,] %>% mutate(timesteps = 100), 
                init.rows = 1e5, theta.t = 0)

test.long %>%
  filter(i > pars$n.pop0[1]) %>%
  group_by(gen, age) %>%
  summarise(n = n()) %>%
  group_by(gen) %>%
  mutate(p = n / sum(n)) %>%
  ungroup() %>%
  filter(gen > 1) %>%
  ggplot(aes(x = age, y = p)) +
  geom_line(aes(group = gen)) +
  geom_line(aes(y = (pars$s.max[1]^(age-1)) * (1 - pars$s.max[1])),
            colour = 'skyblue', linetype = 2, size = 3)
# great - looks good to me
# but!! need to do some division to ensure that we don't include an "age 0"
  
test.vars.age.long = test.long %>%
  filter(i > pars$n.pop0[1]) %>%
  group_by(gen, age) %>%
  summarise(varb = var(b_i))

test.vars.age.long %>%
  filter(gen > 1) %>%
  ggplot(aes(x = gen, y = varb)) +
  geom_line(aes(group = gen - age)) #+

test.long %>%
  #filter(i > pars$n.pop0[1]) %>%
  group_by(age) %>%
  summarise(varb = var(b_i)) %>%
  ggplot(aes(x = age, y = varb)) +
  geom_line() +
  geom_line(aes(y = pars$sig.a[1]^2 * pars$wfitn[1]^2 / (pars$wfitn[1]^2 + ((age+1) * pars$sig.a[1]^2))),
            colour = 'blue', linetype = 2)
# hmm... what's up with all of this?????
# still not tota

# are these fluctuations normal...?

list2 = vector('list', 10)

set.seed(482)

for (trial in 1:9) {
  list2[[trial]] = sim(params = pars[1,], init.rows = 1e5, theta.t = 0) %>%
    group_by(gen, age) %>%
    summarise(varb = var(b_i),
              n = n()) %>%
    mutate(trial.no = trial)
  print(trial)
}

vars2 = do.call(rbind, list2)

vars2 %>%
  ggplot(aes(x = age, y = varb)) +
  geom_line(aes(group = gen, colour = gen)) +
  scale_colour_viridis_c() +
  facet_wrap(~ trial.no)

vars2 %>%
  filter(gen > 0) %>%
  ggplot(aes(x = age, y = varb)) +
  geom_line(aes(group = trial.no, colour = factor(trial.no))) +
  scale_colour_brewer(palette = 'Set1') +
  facet_wrap(~ gen, nrow = 5)

### longer simulations

list.long = vector('list', 50)

set.seed(482)

for (trial in 1:50) {
  list.long[[trial]] = sim(params = pars[1,] %>% mutate(timesteps = 100), 
                           init.rows = 1e5, theta.t = 0) %>%
    group_by(age) %>%
    summarise(varb = var(b_i),
              n = n()) %>%
    mutate(trial.no = trial)
  print(trial)
}

vars.long = do.call(rbind, list.long)

vars.long %>%
  ggplot(aes(x = age, y = varb)) +
  geom_line(aes(group = trial.no))

vars.long %>%
  group_by(age) %>%
  filter(n() > 1) %>%
  summarise(varb.bar = mean(varb),
            varb.var = var(varb),
            n = n()) %>%
  mutate(varb.pred = pars$sig.a[1]^2 * pars$wfitn[1]^2 / (pars$wfitn[1]^2 + ((age+1) * pars$sig.a[1]^2))) %>%
  ggplot(aes(x = age)) +
  geom_line(aes(y = varb.bar)) +
  geom_ribbon(
    aes(
      ymin = varb.bar - 2 * sqrt(varb.var / 9),
      ymax = varb.bar + 2 * sqrt(varb.var / 9)
    ),
    alpha = 0.1
  ) +
  geom_line(aes(y = varb.pred), colour = 'blue', linetype = 2)

# oh... this looks correct to me???

### So I guess the issue is covariance between age groups...

