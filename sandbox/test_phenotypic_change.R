### Looking at phenotypic change over time...
# In Jun 2022 I found what I thought was an expression for phenotypic change over time
# Here I am testing those out against simultions to see if they are/look correct.
# 

library(ggplot2)
library(tidyr)
library(dplyr)

rm(list = ls())

source('model_source/sim_model1_functions.R')

pars = data.frame(
  n.pop0 = 10000, r = (1.02 / 0.8) - 1, wfitn = 3, s.max = 0.8,
  sig.a = sqrt(0.25), sig.e = sqrt(0.5),
  alpha = 0, mu    = 0, sig.m = 0, timesteps = 15
)

n.trials = 50

liszt = vector('list', n.trials)

set.seed(201)

for (k in 1:n.trials) liszt[[k]] = sim(pars, theta = 1, init.rows = 1e5) %>% mutate(trial = k)

all.data = do.call(rbind, liszt)

nrow(all.data)

all.data %>%
  group_by(trial, gen) %>%
  summarise(n  = n()) %>%
  ggplot(aes(x = gen, y = n)) +
  geom_line(aes(group = trial)) +
  scale_y_log10()
# No extinctions

# Here: try to look at expected phenotypic shift versus observed
some.stuff = all.data %>%
  mutate(d_i = theta_t - z_i) %>%
  group_by(trial, gen) %>%
  summarise(
    # Proportion of population that is adult
    p.adult = mean(age > 1),
    # Breeding value variance
    sig2a = var(b_i),
    # Phenotypic variance
    sig2z = var(z_i),
    # Mean phenotypic gap (Lande's z_t)
    dt = mean(d_i),
    # Phenotypic change after selection (Lande's z_wt)
    dw = dt * (pars$wfitn^2 / (sig2z + pars$wfitn^2))
  )

# Proportion of population that is adult over time

some.stuff %>%
  filter(gen > 0) %>%
  ggplot(aes(x = gen, y = p.adult)) +
  geom_line(aes(group = trial))
# looking good

some.stuff %>%
  filter(gen > 0) %>%
  group_by(gen) %>%
  summarise(pbar = mean(p.adult), pvar = var(p.adult), n = n()) %>%
  ggplot(aes(x = gen)) +
  geom_line(aes(y = pbar)) +
  geom_ribbon(aes(ymin = pbar- 2 * sqrt(pvar / n), ymax = pbar + 2 * sqrt(pvar / n)),
              alpha = 0.1)
# hmm... well slightly below 0.8 but still could consider this constant-ish over time

# Compare observed with expectations
some.stuff %>%
  group_by(trial) %>%
  mutate(exp.d = p.adult * dw + (1 - p.adult) * ((sig2a/sig2z) * dw + ((sig2z-sig2a)/sig2z) * dt)) %>%
  mutate(exp.d = c(NA, exp.d[-10])) %>%
  ungroup() %>%
  ggplot(aes(x = exp.d, dt)) +
  geom_segment(aes(x = 0.7, xend = 1, y = 0.7, yend = 1),
               linetype = 2) +
  geom_point(aes(colour = gen)) +
  scale_color_viridis_c() +
  coord_cartesian() +
  theme(legend.position = 'bottom')
# Looks very close but not quite right...
# Not sure why first time step is so different
# (I guess the discrepancy is growing larger with time - 
# approach to SSD?)
# why is the relationship linear in some places...?
# suggests bias in my estimate

# Going back... does the proportion of adults influence the rate of adaptation?
some.stuff %>%
  group_by(trial) %>%
  mutate(obs.k = c(exp(diff(log(dt))), NA)) %>%
  filter(!is.na(obs.k), gen > 0) %>%
  ggplot(aes(x = p.adult, y = obs.k)) +
  geom_point() +
  stat_smooth(method = 'lm')
# Hmm... if at all it's very weak, so not likely

# Look at phenotypes and genotypes in adults/offspring
age.vars = all.data %>%
  filter(gen > 0) %>%
  group_by(trial, gen, adult = age > 1) %>%
  summarise(
    b.bar = mean(b_i),
    b.var = var(b_i),
    z.bar = mean(z_i),
    z.var = var(z_i),
    n = n()
  )

# Look at mean geno/phenotypes for different age groups
age.vars %>%
  select(-c(b.var, z.var, n)) %>%
  gather(vartype, mean.val, -c(trial, gen, adult)) %>%
  mutate(adult   = ifelse(adult, 'Adult', 'Offspring'),
         vartype = ifelse(vartype %in% 'z.bar', 'phenotype', 'genotype')) %>%
  ggplot(aes(x = gen)) +
  geom_point(
    aes(
      y = mean.val, 
      colour = paste0(adult, '-', vartype)
    ),
    position = position_dodge(width = 0.5)
  ) +
  scale_color_manual(values = c('pink', 'red', 'lightblue', 'blue'), '')

# So there is a growing phenotypic gap but only in adults
# no (or very small) gap for genotypes

# Visualize cohort dynamics
all.cohorts.gen = all.data %>%
  mutate(cohort = gen - (age - 1)) %>%
  group_by(gen, trial, cohort) %>%
  #filter(n() > 1) %>%
  summarise(
    mean.b = mean(b_i),
    mean.z = mean(z_i),
    n = n()
  )

all.cohorts.age = all.data %>%
  mutate(cohort = gen - (age - 1)) %>%
  group_by(age, trial, cohort) %>%
  #filter(n() > 1) %>%
  summarise(
    mean.b = mean(b_i),
    mean.z = mean(z_i)
  )

# Plot some of these
all.cohorts.age %>%
  filter(trial < 17) %>%
  ggplot(aes(x = age)) +
  geom_line(aes(y = mean.b, group = cohort, colour = cohort)) +
  scale_colour_viridis_c() +
  facet_wrap(~ trial, nrow = 4)

all.cohorts.gen %>%
  filter(trial < 17) %>%
  group_by(trial, gen) %>%
  mutate(n_t = sum(n)) %>%
  group_by(trial, gen, cohort) %>%
  mutate(p_cohort = n / n_t) %>%
  ungroup() %>%
  ggplot(aes(x = gen, y = mean.b)) +
  geom_line(aes(group = cohort, colour = cohort)) +
  geom_point(aes(size = p_cohort, colour = cohort)) +
  scale_colour_viridis_c(option = 'A') +
  facet_wrap(~ trial, nrow = 4)

all.cohorts.gen %>%
  #filter(trial %in% 1) %>%
  group_by(trial, gen) %>%
  mutate(n_t = sum(n)) %>%
  group_by(trial, gen, cohort) %>%
  mutate(p_cohort = n / n_t) %>%
  group_by(gen, cohort) %>%
  summarise(p_cohort = mean(p_cohort),
            mean.b = mean(mean.b)) %>%
  ggplot(aes(x = cohort, y = mean.b)) +
  geom_line(aes(group = gen, colour = gen)) +
  geom_point(aes(size = p_cohort), shape = 21) +
  scale_fill_viridis_c()

p.cohorts = all.cohorts.gen %>%
  #filter(trial %in% 1) %>%
  group_by(trial, gen) %>%
  mutate(n_t = sum(n)) %>%
  group_by(trial, gen, cohort) %>%
  mutate(p_cohort = n / n_t) %>%
  group_by(gen, cohort) %>%
  summarise(p_cohort = mean(p_cohort),
            mean.b = mean(mean.b),
            mean.z = mean(mean.z))

p.cohorts %>%
  ggplot(aes(x = cohort, y = mean.z)) +
  geom_line(aes(group = gen, colour = gen)) +
  geom_point(aes(size = p_cohort), shape = 21)

p.cohorts %>%
  mutate(age = 1 + (gen - cohort)) %>%
  gather(vartype, var.mean, -c(gen, cohort, p_cohort, age)) %>%
  ggplot(aes(x = cohort, y = var.mean, colour = gen)) +
  geom_line(aes(group = interaction(gen, vartype))) +
  geom_point(aes(fill = age, size = p_cohort), shape = 21) +
  # scale_colour_gradient(low = 'black', high = 'skyblue') +
  scale_fill_viridis_c(option = 'D') +
  facet_wrap(~ vartype)
# interesting...

ggsave("~/Desktop/lh_adaptation_across_cohorts.png")

p.cohorts %>%
  mutate(age = 1 + (gen - cohort)) %>%
  gather(vartype, var.mean, -c(gen, cohort, p_cohort, age)) %>%
  ggplot(aes(x = gen, y = var.mean, colour = cohort)) +
  geom_line(aes(group = interaction(cohort, vartype))) +
  geom_point(aes(fill = age, size = p_cohort), shape = 21) +
  # scale_colour_gradient(low = 'black', high = 'skyblue') +
  scale_fill_viridis_c(option = 'D') +
  facet_wrap(~ vartype)

ggsave("~/Desktop/lh_adaptation_across_time.png")

# I actually really like the above graphics
# (could perhaps differentiate colours and cohorts...?)

###

# Run some sims looking at two different LHs (fast, slow)


rm(list = ls())

source('model_source/sim_model1_functions.R')

pars = data.frame(s.max = c(0.5, 0.8)) %>%
  mutate(
  n.pop0 = 1000, r = (1.02 / s.max) - 1, wfitn = 2, 
  sig.a = sqrt(0.25), sig.e = sqrt(0.5),
  alpha = 0, mu    = 0, sig.m = 0, timesteps = 10
)

n.trials = 200

liszt = vector('list', n.trials)

set.seed(201)

for (k in 1:n.trials) liszt[[k]] = sim(pars[1 + k %% 2, ], theta = 1, init.rows = 1e5) %>% mutate(trial = k)

# Turn into data frame
all.data = do.call(rbind, liszt) %>%
  mutate(row = 1 + trial %% 2) %>%
  merge(y = pars %>% mutate(row = 1:2) %>% select(s.max, row))

# Get variance by age class (adult, offspring)
age.vars = all.data %>%
  filter(gen > 0) %>%
  group_by(trial, gen, adult = age > 1) %>%
  summarise(
    b.bar = mean(b_i),
    b.var = var(b_i),
    z.bar = mean(z_i),
    z.var = var(z_i),
    s.max = s.max[1],
    n = n()
  )

# Visualise - are there differences in rate of genoptyic change?
age.vars %>%
  ggplot(aes(x = gen)) +
  geom_point(aes(y = b.bar, colour = adult),
             position = position_dodge(width = 0.5)) +
  facet_wrap(~ s.max)
# yes - slower genotypic change for slower LH

# Phenotypic change visualization
age.vars %>%
  ggplot(aes(x = gen)) +
  geom_point(aes(y = z.bar, colour = adult),
             position = position_dodge(width = 0.5)) +
  facet_wrap(~ s.max)
# faster phenotypic change in adults, larger gap for slow LH

# Summarise these across all trials
vars.summ = age.vars %>%
  group_by(s.max, gen, adult) %>%
  summarise(
    b..bar = mean(b.bar),
    b..var = var(b.bar),
    z..bar = mean(z.bar),
    z..var = var(z.bar),
    n = n()
  )

# Visualize breeding value variances aggregated across trials
vars.summ %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = b..bar, 
      colour = adult, 
      linetype = factor(s.max)
    )
  ) +
  geom_ribbon(
    aes(
      ymin = b..bar - 2 * sqrt(b..var / n),
      ymax = b..bar + 2 * sqrt(b..var / n),
      fill = adult,
      group = interaction(adult, s.max)
    ),
    alpha = 0.2
  )
  
# yes - slower change with slower LH
# but is this due to survival directly or through variance?

# Look at mean phenotypic, genotypic change across age groupings
vars.summ %>%
  select(-c(b..var, z..var)) %>%
  gather(key = vartype, value = val, -c(s.max, gen, adult, n)) %>%
  ggplot(aes(x = gen, y = val)) +
  geom_line(
    aes(
      group    = interaction(adult, vartype),
      linetype = adult,
      colour   = vartype
    )
  ) +
  labs(x = 'Generation', y = 'Mean') +
  scale_linetype_manual(values = 2:1, labels = c('Offspring', 'Adult'), '') +
  scale_colour_manual(values = c('blue', 'red'), labels = c('Genotype', 'Phenotype'), '') +
  facet_wrap(~ s.max)
# Yes - looks to me like in aggregate:
# - faster phenotypic change with higher survival
# - slower genotypic change (in both age groups)
# - in aggregate a much larger gap! (env. variation...)

ggsave('~/Desktop/geno_pheno_change_by_age_group.png')

# Look at all trials
age.vars %>%
  select(-c(b.bar, z.bar)) %>%
  gather(vartype, var, -c(trial, gen, adult, s.max, n)) %>%
  ggplot(aes(x = gen, y = var)) +
  geom_line(
    aes(
      group = interaction(trial, adult, vartype),
      linetype = vartype,
      colour = adult
      ),
    size = 0.1
    ) +
  facet_wrap( ~ s.max)
# Variances:
# - slower loss of genotypic variance for slower LH, 
# perhaps larger loss of adult phenotypic variance (due to selection) in adults with slow LH

# 
age.vars %>%
  group_by(gen, adult, s.max) %>%
  summarise(b.var.bar = mean(b.var), z.var.bar = mean(z.var)) %>%
  gather(vartype, var, -c(gen, adult, s.max)) %>%
  ggplot(aes(x = gen, y = var)) +
  geom_line(
    data = all.data %>% 
      group_by(trial, gen) %>% 
      summarise(b.var = var(b_i),
                z.var = var(z_i)) %>%
      mutate(s.max = 0.5 + 0.3 * (trial %% 2)) %>%
      group_by(s.max, gen) %>%
      summarise(b.var = mean(b.var),
                z.var = mean(z.var)) %>%
      gather(vartype, var, -c(s.max, gen)),
    aes(
      group = vartype
    ),
    size = 2
  ) +
  geom_line(
    aes(
      linetype = adult,
      colour = vartype
    )
  ) +
  labs(x = 'Generation', y = 'Variance') +
  scale_linetype_manual(values = 2:1, labels = c('Offspring', 'Adult'), '') +
  scale_colour_manual(values = c('blue', 'red'), labels = c('Genotype', 'Phenotype'), '') +
  facet_wrap(~ s.max)

# faster erosion of genetic variance with more survival (more selection)
# slightly faster genotypic change with low survival
# (attributable to phenotyes?)

##### Conclusions

# - yes - there is faster phenoytpic change with higher adult survival
# BUT this does not necessarily mean that there is faster *genotypic* change
# in fact there is fslower genotypic change occurring wiht slower LHs
