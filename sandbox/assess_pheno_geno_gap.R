# Here: looking at phenotypic/genotypic divergence
# I noticed in simulations (a couple of weeks ago) that, over time, mean
# population phenotypes and mean population genotypes (breeding values)
# diverged, despite them being initialized with the same mean.
# Interetingly, and intuitively (but unexpectedly) this occurred in adults but
# not in offspring, meaning that it isn't a problem with initializing these.
# A hypothesis is that this is because selection acts on phenotypes and not
# directly on breeding values. Here I ran some simulations to investigate this.
# Mar 16 2021

# Packages
library(ggplot2)
library(dplyr)
library(tidyr)

# Clear namespace
rm(list = ls())

# Source model code
source('model_source/sim_model1_functions.R')

# Initialize parameter table
# We'll run simulations with both changing adult survival (longevity) and
# different levels of phenotypic variance.
# Note that increasing phenotypic variance decreases heritability.
# The sig.e = 0 case corresponds to perfect heritability.

# Growth rate (used to calculate fecundity)
lambda.max = 1.1

# Parameters
# sig.e = environmental variance,
# s.max = max survival
pars = expand.grid(sig.e = c(0, sqrt(0.25), sqrt(0.5)),
                   s.max = c(0.5, 0.7, 0.9),
                   parm.trial = 1:50) %>%
  mutate(n.pop0 = 1000, 
         wfitn = 2,
         sig.a = sqrt(0.25), 
         timesteps = 75) %>%
  # Determine density dependence (depends on N0)
  mutate(alpha = log(lambda.max) / n.pop0) %>%
  # Determine fecundities
  mutate(r = (lambda.max / s.max) * sqrt((wfitn^2 + sig.a^2 + sig.e^2) / wfitn^2) - 1) %>%
  mutate(trial = 1:nrow(.))

# Initialize object for storage
liszt = vector('list', nrow(pars))

# Environmental change
# linearly changing environment w/ no noise
theta.start = 0
theta.trend = 0.06
#theta.sigma = 0.50
theta.sigma = 0

# Length of time to run sims for
t.end = pars$timesteps[1]

set.seed(8820019)

for (k in 1:nrow(pars)) {
  
  this.theta = theta.start + theta.trend*(0:t.end) + rnorm(t.end+1, 0, theta.sigma)
  
  liszt[[k]] = sim(params = pars[k,], theta.t = this.theta, init.rows = 1e6) %>%
    group_by(gen, adult = age > 1) %>%
    summarise(
      n = n(),
      # Summarise 
      dzbar = mean(theta_t - z_i),
      dbbar = mean(theta_t - b_i),
    ) %>%
    ungroup() %>%
    mutate(trial = k)
  
  print(k)
  
}

# Collect all data and bind it with parameter values
all.data = merge(
  x = do.call(rbind, liszt),
  y = pars %>% select(trial, sig.e, s.max)
)

# Plot all phenotypic and genotypic differences
all.data %>%
  ggplot(aes(x = gen, y = dzbar - dbbar, group = interaction(trial, adult))) +
  geom_line(aes(colour = adult), size = 0.1) +
  facet_wrap(s.max ~ sig.e)

# Negative: (th_t - z_t) - (th_t - b_t) < 0 --> 
# for positive theta, the breeding value lag is larger than the phenotypic lag
# (breeding value is lagging behind phenotype)
# More negative --> larger lag

# Initial impressions:
# No heritability means no gap (duh)
# Increasing phenotypic variance increases the gap such that breeding values lag
# behind phenotypes.
# Gap is larger for longer-lived organisms
# (note that longer lived organisms means adults are larger proportion of
# population is adult, so whole-population gap is larger for adults?)
#
# What is the effect that this has on the rate of adaptation? 

# Check for extinctions
all.data %>%
  group_by(trial) %>%
  summarise(last.gen = max(gen)) %>%
  filter(last.gen < 75)
# exactly one - nice.

all.data %>%
  filter(gen < 75) %>%
  group_by(s.max, sig.e, gen, adult) %>%
  summarise(diff.mean = mean(dzbar - dbbar),
            diff.var  = var(dzbar - dbbar),
            n = n()) %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = diff.mean,
      colour = factor(s.max),
      linetype = !adult
    )
  ) +
  geom_ribbon(
    aes(
      ymin = diff.mean - 2 * sqrt(diff.var / n),
      ymax = diff.mean + 2 * sqrt(diff.var / n),
      fill = factor(s.max),
      group = interaction(s.max, sig.e, adult)
      ),
    alpha = 0.2
  ) +
  facet_wrap(~ sig.e)
# Yep - confirms trend above. 

# Now look at population size across groups - does this actually influence
# any demographic patterns?

all.data %>%
  filter(trial < 75) %>%
  group_by(trial, gen, sig.e, s.max) %>%
  summarise(n = sum(n)) %>%
  group_by(gen, sig.e, s.max) %>%
  summarise(nbar = mean(n),
            nvar = var(n),
            nn = n()) %>%
  ggplot(aes(x = gen)) +
  geom_ribbon(
    aes(
      ymin = nbar - 2 * sqrt(nvar / nn),
      ymax = nbar + 2 * sqrt(nvar / nn),
      fill = factor(sig.e)
    ),
    alpha = 0.25
  ) +
  geom_line(
    aes(
      y = nbar,
      colour = factor(sig.e)
    )
  ) +
  facet_wrap(~ s.max)

# Weak effects if anything... advantage of longer longevity just appears to be
# due to initial boost from variance equilibrating...

# Population growth rates would remove this effect

all.data %>%
  filter(trial < 75) %>%
  group_by(trial, gen, sig.e, s.max) %>%
  summarise(n = sum(n)) %>%
  group_by(trial, sig.e, s.max) %>%
  mutate(r_t = c(diff(log(n)), NA)) %>%
  filter(!is.na(r_t)) %>%
  group_by(gen, sig.e, s.max) %>%
  summarise(rtbar = mean(r_t),
            rtvar = var(r_t),
            nn = n()) %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = rtbar,
      group = interaction(s.max, sig.e),
      colour = factor(sig.e)
    )
  ) +
  facet_wrap(~ s.max)
# okay basically no differences in observed growth rate...

### Let's visualize a single simulation instance

set.seed(20090)

this.theta = theta.start + theta.trend*(0:t.end) + rnorm(t.end+1, 0, theta.sigma)

test.sim = sim(params = pars[8,] %>% mutate(timesteps = 25),
               theta.t = this.theta[1:26],
               init.rows = 1e6)

test.sim %>%
  filter(gen %in% 1:16) %>%
  group_by(i) %>%
  mutate(first.gen = min(gen)) %>%
  ungroup() %>%
  ggplot(aes(x = first.gen, y = b_i - theta_t)) +
  geom_point(aes(colour = b_i - z_i)) +
  scale_colour_gradient2(high = 'red', low = 'skyblue',
                         mid = 'white', midpoint = 0) +
  facet_wrap(~ gen, nrow = 4) +
  theme(panel.background = element_rect(fill = 'black'),
        panel.grid = element_blank())

# needs analysis

pheno.luck = test.sim %>%
  group_by(i) %>%
  mutate(first.gen = min(gen)) %>%
  group_by(first.gen, gen) %>%
  summarise(pheno.luck = mean(z_i - b_i),
            age = gen - first.gen,
            n = n()) %>%
  distinct(first.gen, gen, .keep_all = TRUE)

pheno.luck %>%
  ggplot(aes(x = age, y = pheno.luck)) +
  geom_line(aes(group = first.gen, colour = first.gen)) +
  scale_colour_viridis_c()

pheno.luck %>%
  filter(n > 1) %>%
  group_by(age) %>%
  summarise(luck.bar = mean(pheno.luck),
            luck.var = var(pheno.luck),
            n = n()) %>%
  ggplot(aes(x = age)) +
  geom_line(
    aes(
      y = luck.bar
    )
  ) +
  geom_ribbon(
    aes(
      ymin = luck.bar - 2 * sqrt(luck.var / n),
      ymax = luck.bar + 2 * sqrt(luck.var / n)
    ),
    alpha = 0.2
  )
# yes - so "luck" does increase with age, even under this environmental regime!

# Brett calls this the 'selection debt' akin to the "extinction debt"
# could also be "phenotypic extinction debt" or "maladaptation debt"?

### Maybe simulate this out over a bunch of parameters?

# Store results in a new list (use same parameter combos as before)
luck.liszt = vector('list', nrow(pars))

# Simulations same as before (same theta trend - linear, no noise)
set.seed(8820019)

for (k in 1:nrow(pars)) {
  
  this.theta = theta.start + theta.trend*(0:t.end) + rnorm(t.end+1, 0, theta.sigma)
  
  luck.liszt[[k]] = sim(params = pars[k,], theta.t = this.theta, init.rows = 1e6) %>%
    group_by(i) %>%
    mutate(first.gen = min(gen)) %>%
    group_by(first.gen, gen) %>%
    summarise(pheno.luck = mean(z_i - b_i),
              age = gen - first.gen,
              n = n()) %>%
    distinct(first.gen, gen, .keep_all = TRUE) %>%
    mutate(trial = k)
  
  print(k)
  
}

# Collect all data and bind it with parameter values
luck.data = merge(
  x = do.call(rbind, luck.liszt),
  y = pars %>% select(trial, sig.e, s.max)
) %>%
  filter(sig.e > 0)

weighted.luck = luck.data %>%
  group_by(age, sig.e, s.max) %>%
  mutate(n.weight = n / sum(n)) %>%
  summarise(luck.bar = sum(n.weight * pheno.luck),
            luck.var = var(n.weight * pheno.luck),
            n = n()) %>%
  filter(n > 1)

weighted.luck %>%
  ggplot(aes(x = age, y = luck.bar)) +
  geom_line(aes(group = interaction(sig.e, s.max),
                linetype = factor(sig.e),
                colour = factor(s.max)))
# Hmm... weighting by individuals observed in each age class (biasing?),
# there doesn't appear to be an effect of age on the phenotypic lag,
# just how far back it goes...

mean.luck = luck.data %>%
  group_by(age, sig.e, s.max) %>%
  summarise(luck.bar = mean(pheno.luck),
            luck.var = var(pheno.luck),
            n = n()) %>%
  filter(n > 1)

mean.luck %>%
  ggplot(aes(x = age, y = luck.bar)) +
  geom_line(aes(group = interaction(sig.e, s.max),
                linetype = factor(sig.e),
                colour = factor(s.max))) #+
  # geom_ribbon(aes(ymin = luck.bar - 2 * sqrt(luck.var / n),
  #                 ymax = luck.bar + 2 * sqrt(luck.var / n),
  #                 group = interaction(sig.e, s.max),
  #                 fill = factor(s.max)),
  #             alpha = 0.2)

# Not weighting by population size does show a subtle effect of longevity but
# not that much.
# Which is better to use - weighted or unweighted?

### Does weight actually influence the rate of adaptation?

adapt.list = vector('list', 30)

set.seed(295245)

# Run the same parameter combos except perhaps for varying sig.e
for (k in 1:30) {
  
  this.theta = theta.start + theta.trend*(0:t.end) + rnorm(t.end+1, 0, theta.sigma)
  
  adapt.list[[k]] = sim(params = pars[8 + (k%%2),], theta.t = this.theta, init.rows = 1e6) %>%
    group_by(i) %>%
    mutate(first.gen = min(gen)) %>%
    group_by(first.gen, gen) %>%
    summarise(zbar = mean(z_i),
              bbar = mean(b_i),
              pheno.luck = zbar - bbar,
              age = gen - first.gen,
              n = n()) %>%
    distinct(first.gen, gen, .keep_all = TRUE) %>%
    group_by(gen) %>%
    mutate(zbar = weighted.mean(zbar, n),
           bbar = weighted.mean(bbar, n)) %>%
    ungroup() %>%
    mutate(trial = 8 + (k%%2))
  
  print(k)
  
}

# Okay... now how to actually assess luck?
# Is there a way to compress it down to one metric?
# idk...
