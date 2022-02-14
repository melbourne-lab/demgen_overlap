# Script for incorporating NDD (survival) and phenotypic variance into STABLE
# population growth rates.
# A prior, related model (Brett's for NSF proposal) showed evidence of declining
# genotypic variance in the first several time steps
# Here, I'm looking at that decline over time.
# I also through some trial and error figure out parameter combinations to
# ensure stable growth rates (and, as it turns out, also phenotypic variance) so
# long as the population os on the phenotypic optimum

# SN - 14 Feb 2022

rm(list = ls())

source('model_source/sim_model1_functions.R')

pars = data.frame(
  sig.a = sqrt(0.5),
  n.pop0 = 500, s.max = 0.9, r = (1.1 / (0.9)) - 1, 
  wfitn = 2,
  sig.e = 0, alpha = 0.000191,
  timesteps = 50
)

set.seed(20999)

j = sim(pars, theta.t = 0, init.rows = 1e5)

j %>%
  group_by(gen) %>%
  summarise(k = var(b_i)) %>%
  ggplot(aes(x = gen, y = k)) +
  geom_line() +
  geom_point()


j %>%
  ggplot(aes(x = gen, y = b_i)) +
  geom_line(aes(group = i), size = 0.1) +
  geom_point(size = 0.1)

j %>%
  ggplot(aes(x = b_i, group = gen, colour = gen)) +
  geom_density() +
  scale_color_viridis_c(option = 'B')

j %>%
  ggplot(aes(x = b_i)) +
  geom_density() +
  facet_wrap(~ gen)

j %>%
  ggplot(aes(x = s_i)) +
  geom_density() +
  facet_wrap(~ gen)

### oh... is there a cyclical dynamic here? shrinking due to NDD, purging tails, tails re-emerging?

n.trials = 200

li = vector(mode = 'list', length = n.trials)

set.seed(88439)

for (k in 1:n.trials) li[[k]] = sim(pars, theta.t = 0, init.rows = 1e5) %>% mutate(trial = k)

mm = do.call(li, what = rbind)

mm %>%
  group_by(trial, gen) %>%
  summarise(varb = var(b_i)) %>%
  group_by(gen) %>%
  summarise(varb.bar = mean(varb),
            varb.var = var(varb),
            n = n()) %>%
  ggplot(aes(x = gen)) +
  geom_line(aes(y = varb.bar)) +
  geom_ribbon(aes(ymin = varb.bar-2*sqrt(varb.var/n), ymax = varb.bar+2*sqrt(varb.var/n)),
              alpha = 0.2)

# Yes - observable decrease in genetic variation, as expected.
# But why?

varbs = mm %>%
  group_by(trial, gen) %>%
  summarise(varb = var(b_i))

varbs %>%
  ggplot(aes(x = gen, y = varb)) +
  geom_point(position = position_jitter(width = 0.25), alpha = 0.5, size = 2)

mm %>% group_by(trial) %>% summarise(maxg = max(gen)) %>% group_by(maxg) %>% summarise(n = n())
# no extinctions...

mm %>%
  group_by(trial, gen) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = gen, y = n)) + 
  geom_point() +
  scale_y_log10()

# ah... mean population size is changing as well.
# are these related?

# two ways to tell:
# one: assess population growth rate against change in genetic variance
# two: manually adjust parameters to avoid initial population decline

# One: look at population growth rate vs. change in geno variation

dvarb.dn = mm %>%
  group_by(gen, trial) %>%
  summarise(n    = n(),
            varb = var(b_i)) %>%
  group_by(trial) %>%
  mutate(dn = c(diff(n), NA),
         dv = c(diff(varb), NA))

dvarb.dn %>%
  filter(!is.na(dn)) %>%
  ggplot(aes(x = dn, y = dv)) +
  geom_point(aes(colour = gen)) +
  scale_color_viridis_c()

dvarb.dn %>%
  filter(!is.na(dn)) %>%
  ggplot(aes(x = dn, y = dv)) +
  geom_point(aes(colour = gen)) +
  scale_color_viridis_c() +
  facet_wrap(~ (gen > 10))

with(dvarb.dn %>% filter(!is.na(dn)), cor(dn, dv))

# correlation for sure

dvarb.dn %>%
  filter(!is.na(dn)) %>%
  filter(trial %in% 1:25) %>%
  ggplot(aes(x = dn, y = dv)) +
  geom_point(aes(colour = gen), size = 3) +
  facet_wrap(~ trial) +
  scale_color_viridis_c()

aa = log(4/4.5) / (500)

exp(-aa * 500) * 4 / 4.5

# Try again with a weaker density dependence that should keep pop size constant

pars = data.frame(
  sig.a = sqrt(0.5),
  n.pop0 = 500, s.max = 0.9, r = (1.1 / (0.9)) - 1, 
  wfitn = 2,
  sig.e = 0, alpha = 0.000117783,
  timesteps = 50
)

n.trials = 200

li = vector(mode = 'list', length = n.trials)

set.seed(22985)

for (k in 1:n.trials) li[[k]] = sim(pars, theta.t = 0, init.rows = 1e5) %>% mutate(trial = k)

mm = do.call(li, what = rbind)

mm %>%
  group_by(trial, gen) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = gen, y = n)) + 
  geom_point() +
  scale_y_log10()
# this also looks like it's shrinking, although by not as much...

mm %>%
  group_by(trial, gen) %>%
  summarise(n = n()) %>%
  group_by(gen) %>%
  summarise(nbar = mean(n),
            nvar = var(n),
            n = n()) %>%
  ggplot(aes(x = gen)) +
  geom_line(aes(y = nbar)) +
  geom_ribbon(aes(ymin = nbar-2*sqrt(nvar/n), ymax = nbar+2*sqrt(nvar/n)),
              alpha = 0.2) +
  scale_y_log10()

# still declining...

demo.summ = mm %>%
  group_by(trial, gen) %>%
  summarise(b = mean(b_i),
            s = mean(s_i),
            r = mean(r_i)) %>%
  group_by(gen) %>%
  summarise(bbar = mean(b),
            bvar = var(b),
            sbar = mean(s),
            svar = var(s),
            rbar = mean(r),
            rvar = var(r),
            n = n())

demo.summ %>%
  ggplot(aes(x = gen, y = sbar)) +
  geom_line()
# what on earth happened in that first time step?
# also survival still at 0.8... is this what we want?

demo.summ %>%
  ggplot(aes(x = gen, y = rbar)) +
  geom_line() +
  geom_ribbon(aes(ymin = rbar-2*sqrt(rvar/n), ymax = rbar+2*sqrt(rvar/n)),
              alpha = 0.2)

demo.summ %>%
  ggplot(aes(x = gen, y = sbar * (1 + rbar))) +
  geom_line()

head(demo.summ)

# ah... max growth rate should also account for phenotypic variance
# and we can determine alpha just from lambda max and N0...

pars = data.frame(
  sig.a = sqrt(0.5),
  n.pop0 = 500, s.max = 0.9, r = (1.1  / (0.9)) * sqrt((2^2 + 0.5) / 2^2) - 1, 
  wfitn = 2,
  sig.e = 0, alpha = log(1.1) / 500,
  timesteps = 50
)

n.trials = 200

li = vector(mode = 'list', length = n.trials)

set.seed(22985)

for (k in 1:n.trials) li[[k]] = sim(pars, theta.t = 0, init.rows = 1e5) %>% mutate(trial = k)

mm = do.call(li, what = rbind)

# Cool! This is good... seems like so long as population size is stable, genetic variance will remain stable.
# 
# A problem, though, is that our populations will necessarily be shrinking...
# which means that genotypic variance will not likely be stable...

