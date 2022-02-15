library(ggplot2)
library(dplyr)
library(tidyr)

rm(list = ls())

source('model_source/sim_model1_functions.R')

lambda.max = 1.1

pars = data.frame(s.max = rep(c(0.5, 0.7, 0.9), each = 50)) %>%
  mutate(n.pop0 = 500, 
         wfitn = 2,
         sig.a = sqrt(0.5), 
         sig.e = sqrt(0.5),
         timesteps = 50) %>%
  # Determine density dependence (depends on N0)
  mutate(alpha = log(lambda.max) / n.pop0) %>%
  # Determine fecundities
  mutate(r = (lambda.max / s.max) * sqrt((wfitn^2 + sig.a^2 + sig.e^2) / wfitn^2) - 1) %>%
  mutate(trial = 1:nrow(.))

# Try it out with theta = 0

li = vector('list', nrow(pars))

set.seed(1906)

for (k in 1:nrow(pars)) li[[k]] = sim(pars[k,], theta.t = 0, init.rows = 1e5) %>% mutate(trial = k)

mm = merge(x = do.call(rbind, li), y = pars %>% select(trial, s.max))

size.summ = mm %>%
  group_by(gen, trial) %>%
  summarise(s.max = s.max[1],
            n = n()) %>%
  group_by(gen, s.max) %>%
  summarise(nbar = mean(n),
            nvar = var(n),
            nn = n()) %>%
  mutate(s.max = factor(s.max))

size.summ %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = nbar,
      group = s.max,
      colour = s.max
    )
  ) +
  geom_ribbon(
    aes(
      ymin = nbar - 2 * sqrt(nvar / nn),
      ymax = nbar + 2 * sqrt(nvar / nn),
      group = s.max,
      fill = s.max
    ),
    alpha = 0.2
  )

# hmm... shit? fuck?

varb.summ = mm %>%
  group_by(gen, trial) %>%
  summarise(s.max = s.max[1],
            varb = var(b_i)) %>%
  group_by(gen, s.max) %>%
  summarise(varb.bar = mean(varb),
            varb.var = var(varb),
            nn = n()) %>%
  mutate(s.max = factor(s.max))

varb.summ %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = varb.bar,
      group = s.max,
      colour = s.max
    )
  ) +
  geom_ribbon(
    aes(
      ymin = varb.bar - 2 * sqrt(varb.var / nn),
      ymax = varb.bar + 2 * sqrt(varb.var / nn),
      group = s.max,
      fill = s.max
    ),
    alpha = 0.2
  )

# what on Earth...
# fuck
# okay, genetic variation *increases*

varz.summ = mm %>%
  group_by(gen, trial) %>%
  summarise(s.max = s.max[1],
            varz = var(z_i)) %>%
  group_by(gen, s.max) %>%
  summarise(varz.bar = mean(varz),
            varz.var = var(varz),
            nn = n()) %>%
  mutate(s.max = factor(s.max))

varz.summ %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = varz.bar,
      group = s.max,
      colour = s.max
    )
  ) +
  geom_ribbon(
    aes(
      ymin = varz.bar - 2 * sqrt(varz.var / nn),
      ymax = varz.bar + 2 * sqrt(varz.var / nn),
      group = s.max,
      fill = s.max
    ),
    alpha = 0.2
  )

# erm - but phenotypic variance is decreasing???

merge(varb.summ, varz.summ) %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = varb.bar,
      group = s.max,
      colour = s.max
    ),
    linetype = 2,
  ) +
  geom_ribbon(
    aes(
      ymin = varb.bar - 2 * sqrt(varb.var / nn),
      ymax = varb.bar + 2 * sqrt(varb.var / nn),
      group = s.max,
      fill = s.max
    ),
    alpha = 0.2
  ) +
  geom_line(
    aes(
      y = varz.bar,
      group = s.max,
      colour = s.max
    )
  ) +
  geom_ribbon(
    aes(
      ymin = varz.bar - 2 * sqrt(varz.var / nn),
      ymax = varz.bar + 2 * sqrt(varz.var / nn),
      group = s.max,
      fill = s.max
    ),
    alpha = 0.2
  )

# what on earth...

# Look at the 0.9 cases

mm %>%
  filter(trial %in% 101) %>%
  group_by(gen) %>%
  summarise(varb = var(b_i), barz = var(z_i))

mm %>%
  filter(trial %in% 101) %>%
  ggplot(aes(x = gen, y = b_i)) +
  geom_line(aes(group = i), size = 0.05) +
  geom_point(size = 0.1)

mm %>%
  filter(trial %in% 101) %>%
  ggplot(aes(x = gen, y = z_i)) +
  geom_line(aes(group = i), size = 0.05) +
  geom_point(size = 0.1)

mm %>%
  filter(trial %in% 101, gen %in% 0:24) %>%
  #distinct(i, .keep_all = TRUE) %>%
  ggplot(aes(x = b_i, y = z_i)) +
  geom_point(aes(colour = gen))
  #facet_wrap(~ gen)

dvarb.dn = mm %>%
  group_by(gen, trial) %>%
  summarise(n    = n(),
            varb = var(b_i),
            s.max = s.max[1]) %>%
  group_by(trial) %>%
  mutate(dn = c(diff(n), NA),
         dv = c(diff(varb), NA))

# The biggest change is happening in the first generation - let's just look at those

mm %>%
  filter(gen < 2) %>%
  group_by(gen, trial) %>%
  summarise(s.max = s.max[1],
            varb = var(b_i),
            varz = var(z_i)) %>%
  ggplot(aes(x = gen, y = varz)) +
  geom_line(aes(group = trial), size = 0.2) +
  geom_point() +
  facet_wrap(~ s.max)

# look at breeding value changes

mm %>%
  filter(gen < 2, trial %in% 104) %>%
  group_by(gen, trial) %>%
  summarise(s.max = s.max[1],
            varb = var(b_i),
            varz = var(z_i))

tr104 = mm %>% filter(gen < 16, trial %in% 104)
tr104.l = tr104 %>% filter(age < 2, gen > 0) %>% distinct(gen, .keep_all = TRUE)

tr104 %>%
  ggplot(aes(x = i, y = b_i)) +
  geom_point(aes(colour = age > 1)) +
  geom_segment(data = tr104.l %>% select(-gen),
               aes(x = i, xend = i, y = -2.5, yend = 2.5),
               size = 0.2) +
  facet_wrap(~ gen) +
  scale_color_manual(values = c('blue', 'red')) +
  theme(legend.position = 'none')

# is it because of the extreme parents reproducing?
# (seems like the selection strength should influence equilibrium dist'n...)

tr104 %>%
  ggplot(aes(x = i, y = b_i)) +
  geom_point(aes(colour = age > 1)) +
  geom_segment(data = tr104.l %>% select(-gen),
               aes(x = i, xend = i, y = -2.5, yend = 2.5),
               size = 0.2) +
  facet_wrap(~ gen) +
  scale_color_manual(values = c('blue', 'red')) +
  theme(legend.position = 'none')

tr104 %>%
  ggplot(aes(x = z_i)) +
  geom_density(aes(group = gen, colour = gen), size = 0.5) +
  scale_color_viridis_b(option = 'C')
# wider tail...

trall = mm %>% filter(gen < 3, (trial %% 50) < 20)

trall %>%
  ggplot(aes(x = i, y = b_i)) +
  geom_point(aes(colour = age > 1), alpha = 0.01) +
  facet_wrap(s.max ~ gen, ncol = 3) +
  scale_color_manual(values = c('blue', 'red')) +
  theme(legend.position = 'none')
# okay - this makes some sense...
# offspring just push boundaries...
# phenotypes though...

trall %>%
  ggplot(aes(x = i, y = z_i)) +
  geom_point(aes(colour = age > 1), alpha = 0.01) +
  facet_wrap(s.max ~ gen, ncol = 3) +
  scale_color_manual(values = c('blue', 'red')) +
  theme(legend.position = 'none')

# hmm... phenotypic variance in *offspring* is higher
# guess there is phenotypic variance added by offspring and lost as parents drop out
# and it's the relative rates of these that influences overall variance?
# wvariance is *lost* in parents and parents that remain will have low variance
# also note weighted variance will be weighted by relative size of offspring/parents
# fast LH strategy: lots of offspring, few surviving parents, so 
# slow LH strategy: fewer offspring, many parents (still losing variance)
#
# but then, why do we see different patterns between b_i and z_i???
# I have absolutely no idea and it's freaky.

# Okay... lessons from this.
# Different strategies will, as expected, have both different rates of change in
# variances and also equilibrium levels of variance.
# Weirdly, phenotypic and genetic variances do not necessarily show the same
# patterns and it's not clear why.
# The original point of this script was to look at population dynamics for
# different parameter combinations with the same max growth rate (because this
# depends on phentotypic variance, it turns out). Though I didn't get a chance
# to look at this, annoyingly populations here do show different patterns even
# with the same maximum growth rate).
# I do see the expected pattern where faster LH strategies have lower growth
# rates - this may be due to the variance in the growth rate, but it actually
# probably has more to do with the equilibrium level of phenotypic variance
# (fuck!) although admittedly the phenotypic variance dynamics may also have to
# do with process variance as well (ahhh)

