# Looking at model performance with different smax values
# when theta is stable and at zero (i.e., ideal non-fluctuating environment)
# SN - init 20 Jan 2022

library(ggplot2)
library(dplyr)
library(tidyr)

rm(list = ls())

source('model_source/sim_model1_functions.R')

lambda.max = 1.2

pars = data.frame(s.max = rep(c(0.5, 0.7, 0.9), each = 50)) %>%
  mutate(r = (lambda.max / s.max) - 1) %>%
  mutate(n.pop0 = 250, 
         wfitn = 2,
         sig.a = sqrt(0.5), 
         sig.e = sqrt(0.5), 
         alpha = 0.00073,
         timesteps = 30) %>%
  mutate(trial = 1:nrow(.))

### Theta = 0

l.th0 = vector('list', nrow(pars))

set.seed(3856666)

for (k in 1:nrow(pars)) {
  l.th0[[k]] = sim(pars[k,], theta.t = 0, init.rows = 1e6) %>% mutate(trial = k)
  print(k)
}

d.th0 = do.call(rbind, l.th0)

nrow(d.th0)

table(d.th0$trial)
# nice...?

th0.summary = d.th0 %>%
  group_by(trial, gen) %>%
  summarise(n = n(),
            bbar = mean(b_i))

th0.demo = merge(
  x = th0.summary %>%
    select(-bbar) %>%
    rbind(expand.grid(trial = 1:150, gen = 0:30, n = 0)) %>%
    group_by(trial, gen) %>%
    summarise(n = sum(n)), 
  y = pars %>% select(trial, s.max)) %>%
  group_by(s.max, gen) %>%
  summarise(nbar = mean(n),
            nvar = mean(n),
            nn = n())

th0.demo %>%
  ggplot(aes(x = gen, colour = factor(s.max), fill = factor(s.max))) +
  geom_line(aes(y = nbar)) +
  geom_ribbon(aes(ymin = nbar-2*sqrt(nvar/nn), ymax = nbar+2*sqrt(nvar/nn)),
              alpha = 0.2) +
  scale_y_log10()
# oh...
# why are these declining? and why is there difference in rate of decline?

# also there's the dumb artifact of reaching SSD arg...

th0.summary %>%
  group_by(trial) %>%
  summarise(extinct = max(gen) < 30) %>%
  merge(y = pars %>% select(trial, s.max)) %>%
  group_by(s.max) %>%
  summarise(p.ext = mean(extinct))
# extinction is... pretty rare. There's honestly just more decline in the shorter-lived
# seems bad...

merge(x = th0.summary, y = pars %>% select(trial, s.max)) %>%
  ggplot(aes(x = gen, y = bbar, colour = log(n+1))) +
  geom_point(position = position_jitter(width = 0.5),
             alpha = 0.5) +
  scale_color_viridis_c() +
  facet_wrap(~ s.max)

# might have causality backwards here... may be more variance in b
# because populations are small.
# but why are populations small...

merge(d.th0, pars %>% select(trial, s.max)) %>%
  filter(s.max < 0.7) %>%
  group_by(trial, gen) %>%
  summarise(r.sum = sum(r_i),
            r.bar = mean(r_i)) %>%
  ggplot(aes(x = gen, y = r.bar)) +
  geom_line(aes(group = trial), size = 0.1)

# More breeding -> more outlying genos -> load?

merge(d.th0, pars %>% select(trial, s.max)) %>%
  filter(age %in% 1, gen > 0) %>%
  ggplot(aes(x = gen, y = b_i)) +
  geom_point(position = position_jitter(width = 0.5),
             alpha = 0.05) +
  facet_wrap(~ s.max)
# hmm...

s2a.offspr.summary = d.th0 %>%
  filter(age %in% 1, gen > 0) %>%
  group_by(trial, gen) %>%
  summarise(n = n(),
            bbar = mean(b_i),
            bvar = var(b_i)) %>%
  merge(y = pars %>% select(trial, s.max)) %>%
  group_by(s.max, gen) %>%
  summarise(nbar = mean(n),
            nvar = var(n),
            s2a.bar = mean(bvar),
            s2a.var = var(bvar))

s2a.offspr.summary %>%
  ggplot(aes(x = gen)) +
  geom_line(aes(y = s2a.bar,
                colour = factor(s.max)),
            size = 3) +
  geom_ribbon(aes(ymin = s2a.bar - 2 * sqrt(s2a.var / 150),
                  ymax = s2a.bar + 2 * sqrt(s2a.var / 150),
                  fill = factor(s.max)),
              alpha = 0.2)
# there does look to be more variance for smaller s.max...

merge(d.th0, pars %>% select(trial, s.max)) %>%
  filter(age %in% 1, gen %in% 1:5) %>%
  ggplot(aes(x = b_i)) +
  geom_segment(aes(x = -2, xend = -2, y = 0, yend = 400),
               colour = 'red') +
  geom_segment(aes(x = 2, xend = 2, y = 0, yend = 400),
               colour = 'red') +
  geom_histogram() +
  facet_grid(rows = vars(s.max), cols = vars(gen))

# Not really seeing evidence of a shifted distribution but could be wrong?

merge(d.th0, pars %>% select(trial, s.max)) %>%
  filter(age %in% 1, gen %in% 1:16) %>%
  ggplot(aes(x = b_i, group = factor(s.max))) +
  geom_density(aes(colour = factor(s.max))) +
  facet_wrap(~ gen, nrow = 4)

# These look mostly similar although there is a higher peak (consistently) for the 0.9 no?
# ehhh doesn't look consistent...

# Quantiles and MAD 
qmad = d.th0 %>%
  filter(age %in% 1, gen > 0) %>%
  group_by(trial, gen) %>%
  summarise(q05 = quantile(b_i, 0.05),
            q95 = quantile(b_i, 0.95),
            mad = mad(b_i)) %>%
  merge(y = pars %>% select(trial, s.max)) %>%
  group_by(s.max, gen) %>%
  summarise(m05 = mean(q05),
            v05 = var(q05),
            m95 = mean(q95),
            v95 = var(q95),
            mdd = mean(mad),
            vad = var(mad))

qmad %>%
  ggplot(aes(x = gen, colour = factor(s.max), fill = factor(s.max))) +
  geom_line(aes(y = m05)) +
  geom_ribbon(aes(ymin = m05+2*sqrt(v05/50), ymax = m05-2*sqrt(v05/50)),
              alpha = 0.25)

qmad %>%
  ggplot(aes(x = gen, colour = factor(s.max), fill = factor(s.max))) +
  geom_line(aes(y = m95)) +
  geom_ribbon(aes(ymin = m95+2*sqrt(v95/150), ymax = m95-2*sqrt(v95/150)),
              alpha = 0.25)

qmad %>%
  ggplot(aes(x = gen, colour = factor(s.max), fill = factor(s.max))) +
  geom_line(aes(y = mdd)) +
  geom_ribbon(aes(ymin = mdd+2*sqrt(vad/150), ymax = mdd-2*sqrt(vad/150)),
              alpha = 0.25)

# yeah,... idk, there may be a slight effect for the first few time steps but it disappears

# guess it's worth looking at two more things:
# (1) variance in population growth rate? would decrease growth right...
# (2) maybe also look at equilibration for genotypic variance (issue before)

# variance in growth rate (1) can be done here

th0.summary %>%
  group_by(trial,) %>%
  mutate(log.lambda = c(diff(log(n)), NA)) %>%
  filter(!is.na(log.lambda)) %>%
  merge(y = pars %>% select(trial, s.max)) %>%
  group_by(gen, s.max) %>%
  summarise(var.lam = var(log.lambda)) %>%
  ggplot(aes(x = gen, y = var.lam, colour = factor(s.max))) +
  geom_line()

# Well... what to make of this.
# Stuff happening later - probably due to small population size (vortex!)
# But even earlier there's more variance in the growth rate for small pops
# As theory would expect...
# Do I need to adjust lambda max to account for this? Ugh...
# 
