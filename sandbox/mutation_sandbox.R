library(parallel)
library(ggplot2)
library(dplyr)
library(tidyr)

source('model_source/sim_model1_functions.R')

n.trials = 2000

pars.list = data.frame(
  n.pop0 = 1000, s.max = 0.9, r = (1.1 / (0.9)) - 1, wfitn = sqrt(10),
  sig.a = 1, sig.e = sqrt(4), alpha = 0.0000,
  mu    = 0, sig.m = sqrt(0),
  timesteps = 1
) %>%
  mutate(trial = n.trials) %>%
  uncount(trial) %>%
  mutate(trial = 1:nrow(.)) %>%
  split(f = 1:nrow(.))

set.seed(409)

outp = mclapply(pars.list, 
                function(par.row) {
                  sim(params = par.row, theta.t = 0, init.rows = 3 * 1e4) %>%
                    mutate(e_i = z_i - b_i) %>%
                    group_by(gen, cohort = gen - age + 1) %>%
                    summarise(
                      n = n(),
                      bbar = mean(b_i),
                      zbar = mean(z_i),
                      ebar = mean(e_i),
                      bvar = var(b_i),
                      zvar = var(z_i),
                      evar = var(e_i),
                      ecov = cov(b_i, e_i)
                    ) %>%
                    mutate(k = par.row$trial)
                },
                mc.cores = 6)

outd = do.call(rbind, outp)

head(outd)

outd %>%
  filter(!cohort) %>%
  ggplot(aes(x = gen, y = bvar)) +
  geom_point() +
  geom_line(aes(group = k))

outd %>%
  filter(!cohort) %>%
  ggplot(aes(x = gen, y = zvar)) +
  geom_point() +
  geom_line(aes(group = k))

outd %>%
  filter(!cohort) %>%
  group_by(gen) %>%
  summarise(bvar.bar = mean(bvar),
            bvar.var = var(bvar),
            n = n ())

with(pars.list[[1]], sig.a^2 * ((wfitn^2 + sig.e^2) / (wfitn^2 + sig.e^2 + sig.a^2)))
# yes

outd %>%
  filter(!cohort) %>%
  group_by(gen) %>%
  summarise(zvar.bar = mean(zvar),
            zvar.var = var(zvar),
            n = n ())

with(pars.list[[1]], (sig.a^2 + sig.e^2) * (wfitn^2 / (wfitn^2 + sig.e^2 + sig.a^2)))
# yes (for this cohort at least!)

outd %>%
  filter(!cohort) %>%
  group_by(gen) %>%
  summarise(ecov.bar = mean(ecov),
            ecov.var = var(ecov),
            n = n ())
# hmm... okay sot his covariance is here no matter what
# actually I guess this makes sense...?
# but... what does it even look like to have this covariance when there is no selection??

outd %>%
  filter(!cohort) %>%
  group_by(gen) %>%
    summarise(ebar.bar = mean(ebar),
              ebar.var = var(ebar),
              n = n())
# I mean... these look about the same...

##

set.seed(4729)

sim.trial = sim(pars.list[[1]] %>% mutate(timesteps = 20), theta.t = 0, init.rows = 8 * 1e4)

sim.trial %>%
  filter(age == gen + 1) %>%
  mutate(e_i = z_i - b_i) %>%
  ggplot(aes(x = b_i, y = e_i)) +
  geom_point() +
  facet_wrap(~ gen)
# yeah... def looks like there's something emerging
# wow... so even without any of the pheno variance this is happening
# wait holy shit... wouldn't this influence mean breeding values?

sim.trial %>%
  filter(age == gen + 1) %>%
  mutate(e_i = z_i - b_i) %>%
  ggplot(aes(x = gen, y = b_i)) +
  geom_point(position = position_jitter(width = 0.1))

sim.trial %>%
  filter(age == gen + 1) %>% 
  group_by(gen) %>%
  summarise(bbar = mean(b_i),
            zbar = mean(z_i),
            n = n())
# wow... what the absolute fuck lmao  
# maybe it's a fluke??

sim.trial %>%
  group_by(gen, cohort = gen - age + 1) %>%
  summarise(bbar = mean(b_i),
            # zbar = mean(z_i),
            n = n()) %>%
  ggplot(aes(x = gen - cohort + 1, y = bbar)) +
  geom_point(aes(colour = cohort)) +
  geom_line(aes(group = cohort, colour = cohort)) +
  scale_colour_viridis_c()

sim.trial %>%
  group_by(gen, cohort = gen - age + 1) %>%
  summarise(bbar = mean(b_i),
            zbar = mean(z_i)) %>%
  gather(key = vartype, value = varval, -c(gen, cohort)) %>%
  ggplot(aes(x = gen, y = varval)) +
  geom_point(aes(colour = cohort)) +
  geom_line(aes(group = cohort, colour = cohort)) +
  scale_colour_viridis_c() +
  facet_wrap(~ vartype)

out2 = mclapply(pars.list, 
                function(par.row) {
                  sim(params = par.row %>% mutate(timesteps = 8), theta.t = 0, init.rows = 3 * 1e4) %>%
                    mutate(e_i = z_i - b_i) %>%
                    group_by(gen, cohort = gen - age + 1) %>%
                    summarise(
                      n = n(),
                      bbar = mean(b_i),
                      zbar = mean(z_i),
                      ebar = mean(e_i),
                      bvar = var(b_i),
                      zvar = var(z_i),
                      evar = var(e_i),
                      ecov = cov(b_i, e_i)
                    ) %>%
                    mutate(k = par.row$trial)
                },
                mc.cores = 6)

out3 = do.call(rbind, out2)

out3 %>%
  group_by(gen, cohort) %>%
  summarise(bbar = mean(bbar),
            zbar = mean(zbar),
            ebar = mean(ebar)) %>%
  gather(key = vartype, value = varval, -c(gen, cohort)) %>%
  ggplot(aes(x = gen, y = varval)) +
  geom_point(aes(colour = cohort)) +
  geom_line(aes(group = cohort, colour = cohort)) +
  scale_colour_viridis_c() +
  facet_wrap(~ vartype)

outm %>%
  group_by(gen, cohort) %>%
  summarise(bvarbar = mean(bvar) / pars.list[[1]]$sig.a^2,
            zvarbar = mean(zvar) / (pars.list[[1]]$sig.a^2 + pars.list[[1]]$sig.e^2),
            evarbar = mean(evar) / pars.list[[1]]$sig.e^2) %>%
  gather(key = vartype, value = varval, -c(gen, cohort)) %>%
  ggplot(aes(x = gen, y = varval)) +
  geom_point(aes(colour = cohort)) +
  geom_line(aes(group = cohort, colour = cohort)) +
  scale_colour_viridis_c() +
  # scale_y_log10() + # (though this would be linear on a log scale) * no wait no it isn't duh
  facet_wrap(~ vartype) 
# happy channuka!!
# tyhis actually is looking pretty cool!
# and it makes sense!
# although note the breeding value loss over time.

### well... anyway nice

# if each generation is losing some fraction of its variance each turn
# can restore it with mutations...
# 




# hmm... not really sure what to make of this. looks like it's all small?
# (maybe cohorts are small and noisy by the end?)

### okay anyway

replicate(1000, 
  rnorm(100, 0, sd = rbinom(100, 1, 0.05) * sqrt(1 / (5 + 10) / ((0.05 * 100)^2))) %>%
    var(),
  
) %>%
  mean()
# actually... maybe this is okay?
# target is: 1 / (5 + 10) = .06
# ah... wait shit nvm I guess lmao

replicate(10000, 
          rnorm(100, 0, sd = sqrt(1 / (5 + 10))) %>%
            var(),
          
) %>%
  mean()
# right on this is what I want right

replicate(10000, 
          c(rnorm(5, 0, sd = 1), rep(0, 95)) %>%
            var(),
          
) %>%
  mean()

replicate(10000, 
          rnorm(100, 0, sd = rbinom(100, 1, 0.05)) %>%
            var(),
          
) %>%
  mean()
# okay... this works same as above
# maybe I just have the wrong sigma...?

replicate(10000, 
          rnorm(100, 0, sd = rbinom(100, 1, 0.05) * (0.05 * 100)) %>%
            var(),
          
) %>%
  mean()
# okay so this does *not* give you 1
# what's up with that lmao what am I missing?
# some jensen's thing, perhaps?

replicate(10000, 
          rnorm(100, 0, sd = sqrt(rbinom(100, 1, 0.05) / 0.05) ) %>%
            var(),
          
) %>%
  mean()
# ayy lmao duh

# or is it a law of total variance? uncertainty not just in the normal
# variable but also in the binomial...
# or is it due to sampling?

replicate(1000, 
          rnorm(100, 0, sd = rbinom(100, 1, 0.05) * sqrt(((0.05 * 100)^2) / (5 + 10))) %>%
            var(),
          
) %>%
  mean()
# eh... who knows.


pars.mute = data.frame(
  n.pop0 = 1000, s.max = 0.9, r = (1.1 / (0.9)) - 1, wfitn = sqrt(10),
  sig.a = 1, sig.e = sqrt(4), alpha = 0.0000,
  timesteps = 1
) %>%
  mutate(
    mu = 0.05, 
    sig.m = sqrt((sig.a^2)^2 / (wfitn^2 + sig.a^2 + sig.e^2) / mu)
) %>%
  mutate(trial = n.trials) %>%
  uncount(trial) %>%
  mutate(trial = 1:nrow(.)) %>%
  split(f = 1:nrow(.))

set.seed(420)

outm = mclapply(pars.mute, 
                function(par.row) {
                  sim(params = par.row, theta.t = 0, init.rows = 3 * 1e4) %>%
                    mutate(e_i = z_i - b_i) %>%
                    group_by(gen, cohort = gen - age + 1) %>%
                    summarise(
                      n = n(),
                      bbar = mean(b_i),
                      zbar = mean(z_i),
                      ebar = mean(e_i),
                      bvar = var(b_i),
                      zvar = var(z_i),
                      evar = var(e_i),
                      ecov = cov(b_i, e_i)
                    ) %>%
                    mutate(k = par.row$trial)
                },
                mc.cores = 6)

outm = do.call(rbind, outm)

nrow(outm)

outm %>%
  group_by(gen, cohort) %>%
  summarise(
    bvar.bar = mean(bvar),
    bvar.var = var(bvar)
  )
# hooo I got it lmao willya looka that

set.seed(420)

outm = mclapply(pars.mute, 
                function(par.row) {
                  sim(params = par.row %>% mutate(timesteps = 4), theta.t = 0, init.rows = 5 * 1e4) %>%
                    mutate(e_i = z_i - b_i) %>%
                    group_by(gen, cohort = gen - age + 1) %>%
                    summarise(
                      n = n(),
                      bbar = mean(b_i),
                      zbar = mean(z_i),
                      ebar = mean(e_i),
                      bvar = var(b_i),
                      zvar = var(z_i),
                      evar = var(e_i),
                      ecov = cov(b_i, e_i)
                    ) %>%
                    mutate(k = par.row$trial)
                },
                mc.cores = 6)

outm = do.call(rbind, outm)

nrow(outm)

outm %>%
  group_by(gen, cohort) %>%
  summarise(
    bvar.bar = mean(bvar),
    bvar.var = var(bvar)
  )

outm %>%
  filter(abs(gen - cohort) < 2) %>%
  group_by(cohort, gen) %>%
  summarise(
    bvar.bar = mean(bvar),
    bvar.var = var(bvar)
  )

# okay... still some losses happening here
# rats.
# but why...?

set.seed(4729)

sim.trial = sim(pars.mute[[1]] %>% mutate(timesteps = 5), theta.t = 0, init.rows = 8 * 1e4)

sim.trial %>%
  ggplot(aes(x = i, y = b_i)) +
  geom_point(aes(colour = age > 1)) +
  facet_wrap(~ gen)

# yeah idk
# (but also think about how fucking baller a mutant could be when adapting?)

pars.mute = data.frame(
  n.pop0 = 1000, s.max = 0.9, r = (1.1 / (0.9)) - 1, wfitn = sqrt(10),
  sig.a = 1, sig.e = sqrt(4), alpha = 0.0000,
  timesteps = 1
) %>%
  mutate(
    mu = 0.1, 
    sig.m = sqrt((sig.a^2)^2 / (wfitn^2 + sig.a^2 + sig.e^2) / mu)
  ) %>%
  mutate(trial = n.trials) %>%
  uncount(trial) %>%
  mutate(trial = 1:nrow(.)) %>%
  split(f = 1:nrow(.))

set.seed(420)

outm = mclapply(pars.mute, 
                function(par.row) {
                  sim(params = par.row %>% mutate(timesteps = 4), theta.t = 0, init.rows = 3 * 1e4) %>%
                    mutate(e_i = z_i - b_i) %>%
                    group_by(gen, cohort = gen - age + 1) %>%
                    summarise(
                      n = n(),
                      bbar = mean(b_i),
                      zbar = mean(z_i),
                      ebar = mean(e_i),
                      bvar = var(b_i),
                      zvar = var(z_i),
                      evar = var(e_i),
                      ecov = cov(b_i, e_i)
                    ) %>%
                    mutate(k = par.row$trial)
                },
                mc.cores = 6)

# problem still appears...
# (presumably because new generation's breeding values is mixed with older, less-diverse gene pool)
# (yeah actually I think this has to be it)

# what are the long term dynamics?

outm = mclapply(pars.mute[1:200], 
                function(par.row) {
                  sim(params = par.row %>% mutate(timesteps = 40), theta.t = 0, init.rows = 3 * 1e4) %>%
                    mutate(e_i = z_i - b_i) %>%
                    group_by(gen, cohort = gen - age + 1) %>%
                    summarise(
                      n = n(),
                      bbar = mean(b_i),
                      zbar = mean(z_i),
                      ebar = mean(e_i),
                      bvar = var(b_i),
                      zvar = var(z_i),
                      evar = var(e_i),
                      ecov = cov(b_i, e_i)
                    ) %>%
                    mutate(k = par.row$trial)
                },
                mc.cores = 6)

#############
### real shit
#############

trials.per = 250

params.mute = expand.grid(
  mu = c(1, 10, 50) / 100, 
  s.max = c(0.5, 0.9),
  trial = 1:trials.per
  ) %>%
  mutate(
    n.pop0 = 1000, r = (1.1 / (s.max)) - 1, wfitn = sqrt(10),
    sig.a = 1, sig.e = sqrt(4), alpha = 0.0000,
    timesteps = 10,
    sig.m = sqrt((sig.a^2)^2 / (wfitn^2 + sig.a^2 + sig.e^2) / mu)
  ) %>%
  mutate(trial = 1:nrow(.)) %>%
  split(f = 1:nrow(.))

set.seed(31676)

outy = mclapply(params.mute, 
                function(par.row) {
                  sim(params = par.row, theta.t = 0, init.rows = 1e5) %>%
                    mutate(e_i = z_i - b_i) %>%
                    group_by(gen, cohort = gen - age + 1) %>%
                    summarise(
                      n = n(),
                      bbar = mean(b_i),
                      zbar = mean(z_i),
                      ebar = mean(e_i),
                      bvar = var(b_i),
                      zvar = var(z_i),
                      evar = var(e_i),
                      ecov = cov(b_i, e_i)
                    ) %>%
                    mutate(k = par.row$trial)
                },
                mc.cores = 6)

outy = merge(
  x = do.call(rbind, outy),
  y = do.call(rbind, params.mute) %>% select(trial, s.max, mu),
  by.x = 'k', by.y = 'trial'
)

outy %>%
  filter(s.max > 0.5) %>%
  group_by(gen, cohort, s.max, mu) %>%
  summarise(bvarbar = mean(bvar) / params.mute[[1]]$sig.a^2,
            zvarbar = mean(zvar) / (params.mute[[1]]$sig.a^2 + params.mute[[1]]$sig.e^2),
            evarbar = mean(evar) / params.mute[[1]]$sig.e^2) %>%
  gather(key = vartype, value = varval, -c(gen, cohort, s.max, mu)) %>%
  ggplot(aes(x = gen, y = varval)) +
  geom_point(aes(colour = cohort)) +
  geom_line(aes(group = cohort, colour = cohort)) +
  scale_colour_viridis_c() +
  facet_wrap(vartype ~ mu) 
# looks the same lmao

outy %>%
  filter(s.max < 0.9) %>%
  group_by(gen, cohort, s.max, mu) %>%
  summarise(bvarbar = mean(bvar) / params.mute[[1]]$sig.a^2,
            zvarbar = mean(zvar) / (params.mute[[1]]$sig.a^2 + params.mute[[1]]$sig.e^2),
            evarbar = mean(evar) / params.mute[[1]]$sig.e^2) %>%
  gather(key = vartype, value = varval, -c(gen, cohort, s.max, mu)) %>%
  ggplot(aes(x = gen, y = varval)) +
  geom_point(aes(colour = cohort)) +
  geom_line(aes(group = cohort, colour = cohort)) +
  scale_colour_viridis_c() +
  facet_wrap(vartype ~ mu) 
# hmm.

# so it would appear that there is no effect of mutation rate on variances at least
# let's just look at pop size briefly

with(outy, table(k, gen)) %>% (function(x) sum(!x))
# looks good to me

outy %>%
  group_by(k, gen, s.max, mu) %>%
  summarise(n = sum(n)) %>%
  group_by(gen, s.max, mu) %>%
  summarise(nn = mean(n),
            se = sqrt(var(n) / n())) %>%
  ggplot(aes(x = gen, y = nn)) +
  geom_ribbon(aes(ymin = nn-2*se, ymax = nn+2*se,
                  group = interaction(s.max, mu)), 
              alpha = 0.1) +
  geom_line(aes(linetype = factor(s.max), colour = factor(mu)))

# oh it also doesn't influence population size...
