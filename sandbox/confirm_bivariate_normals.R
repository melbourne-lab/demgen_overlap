# summer 2022 (july-august) script demonstrating that breeding values and
# environmental variance components are in fact distributed bivariate normally
# (has not been touched since august)

library(ggplot2)
library(dplyr)
library(tidyr)

rm(list = ls())

source('model_source/sim_model1_functions.R')

pars = data.frame(
  n.pop0 = 5000, r = (1.02 / 0.8) - 1, wfitn = 1.2, s.max = 0.8,
  sig.a = sqrt(0.25), sig.e = sqrt(0.5),
  alpha = 0, mu    = 0, sig.m = 0, timesteps = 1
)

n.trials = 50
liszt = vector('list', n.trials)

set.seed(9941)

for (k in 1:n.trials) liszt[[k]] = sim(pars, theta.t = 3, init.rows = 2*1e5) %>% mutate(trial = k)

all.sims = do.call(rbind, liszt)

all.sims %>%
  mutate(e_i = b_i - z_i,
         age = factor(age)) %>%
  group_by(trial, gen, age) %>%
  summarise(be.cov = cov(b_i, e_i)) %>%
  ggplot(aes(x = gen, y = be.cov, colour = age)) +
  geom_point(position = position_dodge(width = 0.1)) +
  scale_colour_manual(values = c('blue', 'red'))
# yes... there is positive covariance
# surely this means non-normality of phenotypes... hmm... poop)

all.sims %>%
  filter(trial < 5) %>%
  mutate(e_i = b_i - z_i) %>%
  arrange(age, gen) %>%
  ggplot(aes(x = b_i, y = e_i)) +
  geom_point(aes(colour = paste0(gen, ';', age)), alpha = 0.5) +
  scale_colour_manual(values =c('skyblue', 'blue', 'red')) +
  facet_wrap(~ trial) +
  theme(legend.position = 'none')

all.sims %>%
  filter(trial < 5) %>%
  mutate(e_i = z_i - b_i) %>%
  group_by(trial, gen, age) %>%
  mutate(b.cen = b_i - mean(b_i),
         e.cen = e_i - mean(e_i)) %>%
  arrange(age, gen) %>%
  ggplot(aes(x = b.cen, y = e.cen)) +
  geom_point(aes(colour = paste0(gen, ';', age)), alpha = 0.5) +
  scale_colour_manual(values =c('skyblue', 'blue', 'red')) +
  facet_wrap(~ trial) +
  theme(legend.position = 'none')

all.sims %>%
  mutate(e_i = z_i - b_i) %>%
  group_by(trial, gen, age) %>%
  mutate(b.cen = b_i - mean(b_i),
         e.cen = e_i - mean(e_i)) %>%
  # arrange(age, gen) %>%
  ggplot(aes(x = b.cen, y = e.cen)) +
  geom_point(alpha = 0.1) +
  facet_wrap(~ paste0(gen, ';', age))

all.sims %>%
  group_by(trial, age, gen) %>%
  mutate(e_i = z_i - b_i) %>%
  summarise(b.bar = mean(b_i),
            z.bar = mean(z_i),
            e.bar = mean(e_i)) %>%
  gather(var, varval, -c(trial, age, gen)) %>%
  mutate(age = factor(age)) %>%
  filter(gen > 0) %>%
  ggplot(aes(x = var, y = varval, colour = age)) +
  geom_point(position = position_dodge(width = 0.1)) +
  geom_point(data = pars %>%
               mutate(b.bar = 3 * (1 - (sig.e^2 + wfitn^2) / (sig.e^2 + wfitn^2 + sig.a^2)),
                      z.bar = 3 * (1 - wfitn^2 / ((sig.e^2 + wfitn^2 + sig.a^2))),
                      e.bar = 3 * sig.e^2 / (sig.e^2 + wfitn^2 + sig.a^2) - 0) %>%
               select(b.bar, z.bar, e.bar) %>%
               gather(var, varval),
             position = position_dodge(width = 0.1), colour = 'black')

all.sims %>%
  group_by(trial, age, gen) %>%
  mutate(e_i = z_i - b_i) %>%
  summarise(b.bar = mean(3 - b_i),
            z.bar = mean(3 - z_i),
            e.bar = mean(3 - e_i)) %>%
  gather(var, varval, -c(trial, age, gen)) %>%
  mutate(age = factor(age)) %>%
  filter(gen > 0) %>%
  ggplot(aes(x = var, y = varval, colour = age)) +
  geom_point(position = position_dodge(width = 0.1)) +
  geom_point(data = pars %>%
               mutate(b.bar = 3 * ((sig.e^2 + wfitn^2) / (sig.e^2 + wfitn^2 + sig.a^2)),
                      z.bar = 3 * (wfitn^2 / ((sig.e^2 + wfitn^2 + sig.a^2))),
                      e.bar = 3 * (1 - sig.e^2 / (sig.e^2 + wfitn^2 + sig.a^2))) %>%
               select(b.bar, z.bar, e.bar) %>%
               gather(var, varval),
             position = position_dodge(width = 0.1), colour = 'black')
  
pars %>%
  mutate(b.bar= 3 * (sig.e^2 + wfitn^2) / (sig.e^2 + wfitn^2 + sig.a^2),
         z.bar = 3 * wfitn^2 / ((sig.e^2 + wfitn^2 + sig.a^2)),
         e.bar = 3 * sig.e^2 / (sig.e^2 + wfitn^2 + sig.a^2)) %>%
  select(b.bar, z.bar, e.bar) %>%
  gather(var, varval)

  
all.sims %>%
  filter(gen > 0, age > 1) %>%
  group_by(trial) %>%
  summarise(
    b.var = var(b_i),
    e.var = var(z_i - b_i),
    z.var = var(z_i)
  ) %>%
  summarise(b.var = mean(b.var),
            e.var = mean(e.var),
            z.var = mean(z.var))
  
all.sims %>%
  filter(gen > 0, age > 1) %>%
  group_by(trial) %>%
  summarise(
    b.var = var(b_i),
    e.var = var(z_i - b_i),
    z.var = var(z_i)
  ) %>%
  gather(vartype, varval, -trial) %>%
  ggplot(aes(x = vartype, y = varval)) +
  geom_point(position = position_jitter(width = 0.25), colour = 'blue') +
  geom_point(data = pars %>% 
               mutate(b.var = wfitn^2 * sig.a^2 / (sig.a^2 + wfitn^2),
                      e.var = wfitn^2 * sig.e^2 / (sig.e^2 + wfitn^2),
                      z.var = wfitn^2 * (sig.a^2 + sig.e^2) / (wfitn^2 + sig.a^2 + sig.e^2)) %>%
               select(b.var, e.var, z.var) %>%
               gather(vartype, varval),
             colour = 'red')
             
with(pars, wfitn^2 * sig.a^2 / (wfitn^2 + sig.a^2))
with(pars, wfitn^2 * sig.e^2 / (wfitn^2 + sig.e^2))
# why is this working lol...
with(pars, wfitn^2 * (sig.a^2 + sig.e^2) / (wfitn^2 + sig.a^2 + sig.e^2))

### 

par2 = data.frame(
  n.pop0 = 5000, r = (1.02 / 0.8) - 1, wfitn = 3, s.max = 0.8,
  sig.a = sqrt(0.25), sig.e = sqrt(0.5),
  alpha = 0, mu    = 0, sig.m = 0, timesteps = 2
)

n.trials = 50
lisz2 = vector('list', n.trials)

set.seed(9941)

for (k in 1:n.trials) lisz2[[k]] = sim(par2, theta.t = 2, init.rows = 2*1e5) %>% mutate(trial = k)

all.sim2 = do.call(rbind, lisz2)

all.sim2 %>%
  mutate(e_i = b_i - z_i) %>%
  group_by(trial, gen, age) %>%
  mutate(b.cen = b_i - mean(b_i),
         e.cen = e_i - mean(e_i)) %>%
  # arrange(age, gen) %>%
  ggplot(aes(x = b.cen, y = e.cen)) +
  geom_point(alpha = 0.1) +
  facet_grid(rows = vars(age), cols = vars(gen))

all.sim2 %>%
  mutate(e_i = z_i - b_i) %>%
  group_by(trial, gen, age) %>%
  summarise(b.bar = mean(b_i),
            z.bar = mean(z_i),
            e.bar = mean(e_i))

# yes... greater covariance with time
# (probably also greater with selection strength as well!)

vcv.sum = all.sim2 %>%
  mutate(e_i = z_i - b_i) %>%
  group_by(trial, gen, age) %>%
  summarise(
    b.var = var(b_i),
    e.var = var(e_i),
    z.var = var(z_i),
    becor = cor(b_i, e_i)# ,
    # ezbar = b.var + e.var + 2 * becov
  ) %>% 
  gather(vartype, varvalu, -c(trial, gen, age))

vcv.sum %>%
  ggplot(aes(x = gen, y = varvalu)) +
  geom_segment(aes(x = 0, xend = 2, y = 0, yend = 0), 
               linetype = 2) +
  geom_point(aes(colour = factor(age)),
             position = position_dodge(width = 0.25)) +
  facet_wrap(~ vartype)

vcv.sum %>%
  group_by(gen, age, vartype) %>%
  summarise(
    var.bar = mean(varvalu),
    # var.var = var(varvalu),
    # n = n()
  ) %>%
  ggplot(aes(x = gen, y = var.bar)) +
  geom_segment(aes(x = 0, xend = 2, y = 0, yend = 0), 
               linetype = 2) +
  geom_point(aes(colour = factor(age)),
             position = position_dodge(width = 0.25)) +
  facet_wrap(~ vartype)
# Very little erosion of b_i variance, but erosion of e_i variance?

all.sim2 %>%
  filter(trial %in% 2) %>%
  ggplot(aes(x = b_i)) +
  geom_histogram(binwidth = 0.1) +
  facet_grid(rows = vars(gen), cols = vars(age))

vcv.sum %>%
  filter(age %in% 2, gen %in% 1) %>%
  group_by(vartype) %>%
  summarise(
    var.bar = mean(varvalu),
    var.var = var(varvalu),
    n = n()
  )

with(par2, wfitn^2 * sig.a^2 / (wfitn^2 + sig.a^2)) # yes
with(par2, wfitn^2 * sig.e^2 / (wfitn^2 + sig.e^2)) # yes

with(par2,
     -sqrt(sig.a^2 * sig.e^2 / ((sig.e^2 + wfitn^2) * (sig.a^2 + wfitn^2)))
        )

# erm... hmm... okay this seems plausible? at least the first one.
# wait nope it's wrong lol


# fuck.

library(ggplot2)
library(dplyr)
library(tidyr)

# rm(list = ls())

# source('model_source/sim_model1_functions.R')

pars = data.frame(
  n.pop0 = 5000, r = (1.02 / 0.8) - 1, wfitn = 3, s.max = 0.8,
  sig.a = sqrt(0.25), sig.e = sqrt(0.5),
  alpha = 0, mu    = 0, sig.m = 0, timesteps = 2
)

n.trials = 500

liszt = vector('list', n.trials)

set.seed(59)

for (k in 1:n.trials) liszt[[k]] = sim(params = pars, theta.t = 3, init.rows = 5000 * 4) %>% mutate(trial = k)

# twostep = do.call(rbind, liszt)

twostep.summ = do.call(rbind, liszt) %>% #twostep %>%
  # filter(age %in% 1:2) %>%
  mutate(cohort = gen - (age - 1)) %>%
  group_by(trial, gen, cohort) %>%
  summarise(
    b.bar = mean(theta_t - b_i),
    b.var = var(theta_t - b_i),
    z.bar = mean(theta_t - z_i),
    z.var = var(theta_t - z_i)
  ) %>%
  gather(vartype, varval, -c(trial, gen, cohort)) %>%
  group_by(gen, cohort, vartype) %>%
  summarise(
    varbar = mean(varval),
    varvar = var(varval),
    n = n()
  ) 

twostep.summ %>%
  filter(vartype %in% 'b.bar') %>%
  mutate(cohort = factor(cohort)) %>%
  ggplot(aes(x = gen)) +
  geom_point(aes(y = varbar, colour = cohort)) +
  geom_ribbon(
    aes(
      ymin = varbar - 2 * sqrt(varvar / n),
      ymax = varbar + 2 * sqrt(varvar / n),
      fill = cohort
    ),
    alpha = 0.1
  ) +
  geom_line(aes(y = varbar, colour = cohort)) + 
  geom_point(
    data = pars %>% 
      select(wfitn, sig.e, sig.a) %>%
      mutate(b.1 = 3 * (wfitn^2 + sig.e^2) / (wfitn^2 + sig.a^2 + sig.e^2)),
    aes(x = 1, y = b.1)
  ) +
  geom_point(
    data = pars %>% 
      select(wfitn, sig.e, sig.a) %>%
      mutate(b.2 = 3 * (wfitn^2 + sig.e^2) / (wfitn^2 + sig.a^2 + sig.e^2) * (wfitn^2 + sig.a^2) / (wfitn^2 + 2*sig.a^2)),
      # mutate(b.2 = 3 * (wfitn^2 + sig.e^2) / (wfitn^2 + sig.a^2 + sig.e^2) * 
    aes(x = 2, y = b.2)
  )
# ah... okay looks like there's a mistake
# (although this does look close to the expression for the next cohort...
# but the next cohort has no correlation...?)

twostep.summ %>%
  mutate(age = gen - cohort + 1) %>%
  filter(vartype %in% 'b.var') %>%
  mutate(cohort = factor(cohort)) %>%
  ggplot(aes(x = age)) +
  geom_point(aes(y = varbar, colour = cohort)) +
  geom_ribbon(
    aes(
      ymin = varbar - 2 * sqrt(varvar / n),
      ymax = varbar + 2 * sqrt(varvar / n),
      fill = cohort
    ),
    alpha = 0.1
  )  +
  geom_line(aes(y = varbar, colour = cohort)) # + 
  # geom_point(
  #   data = pars %>% 
  #     select(wfitn, sig.e, sig.a) %>%
  #     mutate(b.1 = sig.a^2 * (wfitn^2 + sig.e^2) / (wfitn^2 + sig.a^2 + sig.e^2)),
  #   aes(x = 1, y = b.1)
  # ) 
  # geom_point(
  #   data = pars %>% 
  #     select(wfitn, sig.e, sig.a) %>%
  #     mutate(b.2 = sig.a^2 * wfitn^2 / (wfitn^2 + 2*sig.a^2)),
  #   aes(x = 2, y = b.2)
  # )

# analytic stuff. messy and shitty
# 

twostep1 = do.call(rbind, liszt) %>% #twostep %>%
  # filter(age %in% 1:2) %>%
  mutate(cohort = gen - (age - 1)) %>%
  group_by(trial, gen, cohort) %>%
  summarise(
    b.bar = mean(theta_t - b_i),
    z.bar = mean(theta_t - z_i)
  ) %>%
  gather(vartype, varval, -c(trial, gen, cohort))

twostep1 %>%
  group_by(trial, vartype, cohort) %>%
  mutate(dd = c(exp(diff(log(varval))), NA)) %>%
  filter(!is.na(dd)) %>%
  ggplot(aes(x = gen, y = dd)) +
  geom_line(
    aes(
      group = trial,
      colour = factor(cohort)
    ),
    size = 0.1
  ) +
  facet_wrap(~ vartype)

twostep2 = twostep1 %>%
  group_by(trial, vartype, cohort) %>%
  mutate(dd = c(exp(diff(log(varval))), NA)) %>%
  group_by(gen, cohort, vartype) %>%
  summarise(var.bar = mean(dd),
            var.var =  var(dd),
            n = n())

twostep2 %>%
  filter(gen < 3) %>%
  ggplot(aes(x = gen)) +
  geom_point(
    aes(
      y = var.bar,
      colour = cohort
    )
  ) +
  geom_line(
    aes(
      y = var.bar,
      colour = cohort,
      group = cohort
    )
  ) +
  facet_wrap(~ vartype)

# but how do we know if this is due to variance loss?
# 
# next? try looking at correlation vs. rate of adaptation?
# (large number of trials, different param combos)

twostep.covs = do.call(rbind, liszt) %>% #twostep %>%
  # filter(age %in% 2) %>%
  filter(gen > 0, age == gen + 1) %>%
  mutate(e_i = z_i - b_i) %>%
  group_by(trial, age) %>%
  summarise(
    b.bar = mean(theta_t - b_i),
    b.var = var(theta_t - b_i),
    e.bar = mean(e_i),
    e.var = var(e_i),
    becov = cov(b_i, e_i),
    becor = cor(b_i, e_i),
    z.bar = mean(theta_t - z_i),
    z.var = var(theta_t - z_i)
  )

head(twostep.covs)

twostep.covs %>%
  select(trial, age, becor) %>%
  ggplot(aes(x = becor)) +
  geom_histogram(binwidth = 0.01) +
  facet_wrap(~ age)

pars %>%
  select(wfitn, sig.a, sig.e) %>%
  mutate(rho   = -sqrt(sig.a^2 * sig.e^2 / ((sig.a^2 + wfitn^2) * (sig.e^2 + wfitn^2))),
         sig.x = sqrt(sig.a^2 * (wfitn^2 + sig.e^2) / (wfitn^2 + sig.e^2 + sig.a^2)),
         sig.y = sqrt(sig.e^2 * (wfitn^2 + sig.a^2) / (wfitn^2 + sig.a^2 + sig.e^2))) %>%
  mutate(gamma = ((wfitn^2*rho) - (1-rho^2)*sig.x*sig.y) / sqrt((wfitn^2 + (1-rho^2)*sig.x^2) * (wfitn^2 + (1-rho^2)*sig.y^2)),
         sig.l = sig.x^2 * ((wfitn^2 + (1-rho^2)*sig.y^2) / (wfitn^2 + sig.x^2 + sig.y^2 + 2*rho*sig.x*sig.y)),
         sig.b = sig.y^2 * ((wfitn^2 + (1-rho^2)*sig.x^2) / (wfitn^2 + sig.x^2 + sig.y^2 + 2*rho*sig.x*sig.y)),
         bbvar = (wfitn^2 / (wfitn^2 + (1-rho^2)*sig.x^2)) * 3 * (wfitn^2 + sig.e^2) / (wfitn^2 + sig.e^2 + sig.a^2))

twostep.covs %>% group_by(age) %>% summarise(rho.obs = mean(becor))
# oh shit
# I got it???
# ohhhhh fuuuuuck yeah

twostep.covs %>% group_by(age) %>% summarise(sig.alph = mean(b.var), sig.beta = mean(e.var))
# yes! these look right

twostep.covs %>% group_by(age) %>% summarise(bobs2 = mean(b.bar), eobs2 = mean(e.bar))
# breeding value at least looks correct!
