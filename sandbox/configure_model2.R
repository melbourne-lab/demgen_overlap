library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyselect)
library(parallel)

rm(list = ls())

source('model_source/sim_model2_functions.R')

pars = expand.grid(
    l.ratio = c(0.85, 0.9, 0.95),
    h2 = c(0.2, 0.5, 1.0),
    s  = c(0.1, 0.5, 0.9)
  ) %>%
  mutate(
    n.pop0 = 1000,
    l.max = 2,
    lstar = l.max * l.ratio,
    r.max = (l.max / s) - 1,
    sig.p = sqrt(1 - ((lstar - s) / (s*r.max))^2),
    # p0    = r.max/sqrt(1+sig.p^2) / (1 + r.max/sqrt(1+sig.p^2)),
    sig.a = sqrt(sig.p^2 * h2),
    sig.e = sqrt(sig.p^2 * (1-h2)),
    sig.m = sqrt(sig.a^4 / (1 + sig.e^2 - sig.a^2)), # ah - needs to be scaled by p0
    mu    = 1,
    wfitn = 1,
    kceil = 3000,
    alpha = log(lstar) / n.pop0,
    gbar0 = 0,
    timesteps = 20
  ) 

set.seed(401)
test.sim = sim2(pars[1,], theta.t = 0, init.rows = 1e4 * 20)
# hmm...

test.sim %>%
  group_by(gen) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = gen, y = n)) +
  geom_line() +
  scale_y_log10()
# wait what...?

test.sim %>%
  ggplot(aes(x = gen, y = b_i)) +
  geom_point(position = position_jitter(width = 0.25), alpha = 0.5)
# ehhhhhnhgn

test.sim %>%
  group_by(gen) %>%
  summarise(v = var(b_i)) %>%
  ggplot(aes(x = gen, y = v)) +
  geom_line()

set.seed(50922)

test.sims = mclapply(
  pars[c(7, 16, 25),] %>% 
    mutate(timesteps = 1) %>% 
    uncount(1000) %>% 
    mutate(try.no = 1:(nrow(.))) %>% 
    split(.$try.no),
  function(pars) {
    sim2(pars, theta.t = 0, init.rows = 3000 * 10) %>%
      group_by(gen, new = (!age | !gen)) %>%
      summarise(n = n(), bvar = var(b_i)) %>%
      mutate(trial = pars$try.no, s = pars$s)
  }
  ) %>%
  do.call(what = rbind)

test.sims %>%
  filter(new) %>%
  group_by(gen, s) %>%
  summarise(bvarbar = mean(bvar), bvarvar = var(bvar)) %>%
  ggplot(aes(x = gen)) +
  geom_line(aes(y = bvarbar, group = s, colour = s)) +
  geom_ribbon(
    aes(
      ymin = bvarbar - 2*sqrt(bvarvar / 1000),
      ymax = bvarbar + 2*sqrt(bvarvar / 1000),
      group = s, fill = s
    ),
    alpha = 0.2
  )

# 
test.sims %>%
  mutate(s = factor(s)) %>%
  filter(new) %>%
  group_by(trial, s) %>%
  mutate(dbvar = diff(bvar)) %>%
  group_by(s) %>%
  summarise(bvarbar = mean(dbvar), bvarvar = var(dbvar)) %>%
  ggplot(aes(x = s)) +
  geom_point(aes(y = bvarbar, group = s)) +
  geom_segment(
    aes(
      xend = s,
      y    = bvarbar - 2*sqrt(bvarvar / 1000),
      yend = bvarbar + 2*sqrt(bvarvar / 1000),
      colour = s
    ),
  )
# honestly this is good enough for me I think

set.seed(333)

# test.sims = mclapply(
#   pars[c(7,16),] %>% mutate(timesteps = 1) %>% uncount(400) %>% mutate(try.no = 1:(nrow(.))) %>% split(.$try.no),
#   function(pars) {
#     sim2(pars, theta.t = 0, init.rows = 3000 * 10) %>%
#       group_by(gen, age) %>%
#       summarise(n = n()) %>%
#       mutate(trial = pars$try.no, s = pars$s)
#   }
# ) %>%
#   do.call(what = rbind)
# 
# test.sims %>% 
#   group_by(trial, gen) %>%
#   mutate(p.age = n / sum(n)) %>%
#   ggplot(aes(x = age, y = p.age)) +
#   geom_line(aes(group = trial), size = 0.1) +
#   facet_wrap(gen ~ s)
#   
# test.sims.age = test.sims %>% 
#   group_by(trial, gen) %>%
#   mutate(p.age = n / sum(n)) %>%
#   group_by(s, gen, age) %>%
#   summarise(p.age.bar = sum(p.age) / 400)
# 
# test.sims.age %>%
#   ggplot(aes(x = age, y = p.age.bar, colour = gen, group = gen)) +
#   geom_line() +
#   geom_point() +
#   facet_wrap(~ s)
# # age distribution looks spot on (good...)
# 
# test.sims.age %>% pivot_wider(names_from = gen, values_from = p.age.bar)

set.seed(229359)

test.sims = mclapply(
  pars[12,] %>% mutate(timesteps = 1) %>% uncount(10) %>% mutate(try.no = 1:(nrow(.))) %>% split(.$try.no),
  function(pars) {
    sim2(pars, theta.t = 0, init.rows = 3000 * 10) %>% mutate(trial = pars$try.no, s = pars$s)
  }
) %>%
  do.call(what = rbind)

test.sims %>%
  mutate(e_i = z_i - b_i, new = (!age & gen > 0)) %>%
  ggplot(aes(x = e_i, group = interaction(trial, gen, new))) + 
  geom_density(aes(colour = gen, linetype = new)) +
  facet_wrap(~ trial)
# uhh... these look the same to me...

test.sims %>% group_by(trial, gen) %>% summarise(evar = var(z_i - b_i)) %>%
  pivot_wider(names_from = gen, values_from = evar)

trys.per = 10
set.seed(4523)

sim.out1 = mclapply(
  pars %>% uncount(trys.per) %>% mutate(try.no = 1:(nrow(.))) %>% split(.$try.no),
  function(pars) {
    sim2(pars, theta.t = 0, init.rows = 3000 * 20) %>%
      mutate(e_i = z_i - b_i) %>%
      group_by(gen) %>%
      summarise(
        n = n(),
        kbar = mean(age),
        bbar = mean(b_i),
        bvar = var(b_i),
        evar = var(e_i),
        zbar = mean(z_i),
        zvar = var(z_i)
      )  %>%
      mutate(
        trial = pars$try.no, 
        s     = pars$s,
        lstar = pars$lstar,
        h2    = pars$h2
      )
  },
  mc.cores = 4
) %>%
  do.call(rbind, .)

sim.out1 %>%
  mutate(h2 = factor(h2)) %>%
  ggplot(aes(x = gen, y = n)) +
  geom_line(aes(colour = h2, group = trial), linewidth = 0.25) +
  facet_wrap(s ~ lstar)

# argh...

sim.out1 %>%
  mutate(h2 = factor(h2)) %>%
  ggplot(aes(x = gen, y = bvar)) +
  geom_line(aes(colour = h2, group = trial), linewidth = 0.25) +
  facet_wrap(s ~ lstar)

sim.out1 %>%
  mutate(h2 = factor(h2)) %>%
  ggplot(aes(x = gen, y = zvar)) +
  geom_line(aes(colour = h2, group = trial), linewidth = 0.25) +
  facet_wrap(s ~ lstar)
# rats

sim.out1 %>%
  mutate(h2 = factor(h2)) %>%
  ggplot(aes(x = gen, y = evar)) +
  geom_line(aes(colour = h2, group = trial), linewidth = 0.25) +
  facet_wrap(s ~ lstar)

# okay not sure why this is not working
# one candidate: maybe I did something wrong with the order of survival...

parz = expand.grid(
  l.ratio = c(0.85, 0.9, 0.95),
  h2 = c(0.2, 0.5, 1.0),
  s  = c(0.1, 0.5, 0.9)
) %>%
  mutate(
    n.pop0 = 10000,
    l.max = 2,
    lstar = l.max * l.ratio,
    wfitn = 1,
    r.max = (l.max / s) - 1,
    sig.p = sqrt(1 - ((lstar - s) / (s*r.max))^2),
    sig.a = sqrt(sig.p^2 * h2),
    sig.e = sqrt(sig.p^2 * (1-h2)),
    sig.m = sqrt(sig.a^4 / (1 + sig.e^2 - sig.a^2)),
    kceil = 30000,
    alpha = log(lstar) / n.pop0,
    gbar0 = sqrt(-2 * log((.55 - s/lstar)/(1 - s/lstar))),
    timesteps = 20
  ) 

set.seed(9033)

sim.out2 = mclapply(
  parz %>% uncount(20) %>% mutate(try.no = 1:(nrow(.))) %>% split(.$try.no),
  function(pars) {
    sim2(pars, theta.t = 0, init.rows = 3000 * 20) %>%
      group_by(gen) %>%
      summarise(
        n = n(),
        kbar = mean(age),
        bbar = mean(b_i),
        bvar = var(b_i),
        zbar = mean(z_i),
        zvar = var(z_i)
      )  %>%
      mutate(
        trial = pars$try.no, 
        s     = pars$s,
        lstar = pars$lstar,
        h2    = pars$h2
      )
  },
  mc.cores = 4
) %>%
  do.call(rbind, .)

sim.out2 %>%
  mutate(s = factor(s)) %>%
  ggplot(aes(x = gen, y = n)) +
  geom_line(aes(colour = s, group = trial), linewidth = 0.25) +
  # scale_y_log10() +
  facet_wrap(paste('heritability', h2) ~ paste0('eqm. growth', lstar))
# yep - looks like... 

merge(
  sim.out2,
  parz %>% select(s, lstar, h2, gbar0),
  by = c('s', 'lstar', 'h2')
) %>%
  mutate(bbar = bbar / gbar0) %>%
  mutate(s = factor(s)) %>%
  ggplot(aes(x = gen, y = bbar)) +
  geom_line(aes(colour = s, group = trial), linewidth = 0.25) +
  facet_wrap(h2 ~ lstar)

sim.out2 %>%
  mutate(s = factor(s)) %>%
  ggplot(aes(x = gen, y = zbar)) +
  geom_line(aes(colour = s, group = trial), linewidth = 0.25) +
  facet_wrap(h2 ~ lstar)

sim.out2 %>%
  mutate(s = factor(s)) %>%
  ggplot(aes(x = gen, y = zvar)) +
  geom_line(aes(colour = s, group = trial), linewidth = 0.25) +
  facet_wrap(h2 ~ lstar)

test.sim = sim2(parz[22,] %>% mutate(gbar0 = 2), theta.t = 0, init.rows = 3000 * 1e4)

test.sim %>%
  group_by(gen) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = gen, y = n)) +
  geom_line() +
  scale_y_log10()

test.sim %>%
  group_by(gen) %>%
  summarise(bbar = mean(b_i)) %>%
  ggplot(aes(x = gen, y = bbar)) +
  geom_line()

test.sum.coh = test.sim %>%
  group_by(gen, cohort = gen - age) %>%
  summarise(bbar = mean(b_i), n = n()) %>%
  group_by(gen) %>%
  mutate(p = n / sum(n),
         b.cont = p * bbar / sum(p * bbar)) 

test.sum.coh %>%
  ggplot(aes(x = gen, y = bbar)) +
  geom_line(aes(group = cohort)) +
  geom_point(aes(fill = p), size = 3, shape = 21) +
  scale_fill_viridis_c()

test.sum.coh %>%
  ggplot(aes(x = gen, y = b.cont)) +
  geom_line(aes(group = cohort)) +
  geom_point(aes(fill = p), size = 3, shape = 21)

test.sum.age = test.sim %>%
  group_by(gen, age) %>%
  summarise(bbar = mean(b_i), n = n()) %>%
  group_by(gen) %>%
  mutate(p = n / sum(n),
         b.cont = p * bbar / sum(p*bbar))

# eh
test.sum.age %>%
  ggplot(aes(x = gen, y = bbar)) +
  geom_line(aes(group = age)) +
  geom_point(aes(fill = p), size = 3, shape = 21) +
  scale_fill_viridis_c()

test.sum.age %>%
  ggplot(aes(x = age, y = p)) +
  geom_line(aes(group = gen, colour = gen)) +
  geom_point(aes(colour = gen)) +
  scale_colour_viridis_c()

test.sum.age %>%
  ggplot(aes(x = gen, y = p)) +
  geom_line(aes(group = age, colour = age > 0)) +
  geom_point(aes(colour = age > 0))

test.sum.age %>%
  ggplot(aes(x = gen, y = b.cont)) +
  geom_line(aes(group = age, colour = age > 0)) +
  geom_point(aes(colour = age > 0))
