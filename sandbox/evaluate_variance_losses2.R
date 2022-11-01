### Looking again at variance losses
# this time: mostly thinking about no mutation, sig.e = 0,
# 0-2 generations, mostly focusing on first cohort
# includes the geometric initial age structure incorporated into variances
# but something still isn't quite working...
# SN - 1 Nov 2022

library(parallel)
library(ggplot2)

rm(list = ls())

source('model_source/sim_model1_functions.R')

pars.list = data.frame(
  n.pop0 = 1000, s.max = 0.9, r = (1.1 / (0.9)) - 1, wfitn = sqrt(10),
  sig.a = 1, sig.e = 0, alpha = 0.0000,
  kceil = 10000,
  timesteps = 20
) %>%
  mutate(
    mu = 1, 
    sig.m = sqrt((sig.a^2)^2 / (wfitn^2 + sig.a^2 + sig.e^2) / mu)
  ) %>%
  uncount(n.trials) %>%
  mutate(trial.no = 1:n.trials) %>%
  split(f = 1:nrow(.))

set.seed(409)

simmy = sim(pars.list[[1]], theta.t = 0, init.rows = 5 * 1e5)

cohort.vars = simmy %>%
  group_by(i) %>%
  mutate(cohort = ifelse(!min(gen), 0, gen - age + 1)) %>%
  group_by(gen, cohort) %>%
  summarise(bvar = var(b_i))

cohort.vars %>%
  ggplot(aes(x = gen, y = bvar)) +
  geom_line(aes(group = cohort, colour = cohort)) +
  geom_point(aes(colour = cohort)) +
  scale_colour_viridis_c()

cohort.vars %>%
  ggplot(aes(x = gen - cohort, y = bvar)) +
  geom_line(aes(group = cohort, colour = cohort)) +
  geom_point(aes(colour = cohort)) +
  scale_colour_viridis_c()
# losing variance even with mutation
# how bad is it without mutation?

# re-run but with mutation rate 0

simme = sim(pars.list[[1]] %>% mutate(mu = 0), theta.t = 0, init.rows = 5 * 1e5)

cohort.vars = simme %>%
  group_by(i) %>%
  mutate(cohort = ifelse(!min(gen), 0, gen - age + 1)) %>%
  group_by(gen, cohort) %>%
  summarise(bvar = var(b_i))

cohort.vars %>%
  ggplot(aes(x = gen, y = bvar)) +
  geom_line(aes(group = cohort, colour = cohort)) +
  geom_point(aes(colour = cohort)) +
  scale_colour_viridis_c()

cohort.vars %>%
  ggplot(aes(x = gen - cohort, y = bvar)) +
  geom_line(aes(group = cohort, colour = cohort)) +
  geom_point(aes(colour = cohort)) +
  scale_colour_viridis_c()
# WOOOF
# okay so mutations... are good for re-establishing variance

#####---------------------------------------------------------------------------

# okay... age structure is geometric w/ prob r/(1+r)
# let's try this...

# I think this is wrong... see different version below
try.loss.expr = function(w2, s2, r) {
  sumo = 0
  for (k in 0:100){ 
    sumo = sumo +
      (r/(1+r))^k * ((w2 + k*s2) * (w2 + (k+1)*s2))^(-1)
  }
  sumo = sumo * (s2^2 * w2) / (1 + r)
}

pars.list = data.frame(
  n.pop0 = 1000, s.max = 0.9, r = (1.1 / (0.9)) - 1, wfitn = sqrt(10),
  sig.a = 1, sig.e = 0, alpha = 0.0000,
  kceil = 10000,
  timesteps = 20
) %>%
  mutate(
    mu = 1, 
    sig.m = sqrt(try.loss.expr(wfitn^2, sig.a^2, r))
  )

simmy = sim(pars.list, theta.t = 0, init.rows = 5 * 1e5)

cohort.vars = simmy %>%
  group_by(i) %>%
  mutate(cohort = ifelse(!min(gen), 0, gen - age + 1)) %>%
  group_by(gen, cohort) %>%
  summarise(bvar = var(b_i))

cohort.vars %>%
  ggplot(aes(x = gen, y = bvar)) +
  geom_line(aes(group = cohort, colour = cohort)) +
  geom_point(aes(colour = cohort)) +
  scale_colour_viridis_c()

cohort.vars %>%
  ggplot(aes(x = gen - cohort, y = bvar)) +
  geom_line(aes(group = cohort, colour = cohort)) +
  geom_point(aes(colour = cohort)) +
  scale_colour_viridis_c()

# what is happening with variances in the first cohort?

cohort0.agevars = simmy %>%
  group_by(i) %>%
  filter(!min(gen)) %>%
  group_by(gen, age.cohort = age - gen) %>%
  summarise(bvar = var(b_i))

cohort0.agevars %>%
  filter(!is.na(bvar)) %>%
  ggplot(aes(x = gen, y = bvar)) +
  geom_line(aes(group = age.cohort, colour = age.cohort)) +
  scale_color_viridis_c()
# not sure this is helpful... ugh

# ---------------------------------------------------------------

n.trials = 100

pars.list = data.frame(
  n.pop0 = 1000, s.max = 0.9, r = (1.1 / (0.9)) - 1, wfitn = sqrt(10),
  sig.a = 1, sig.e = 0, alpha = 0.0000,
  kceil = 10000,
  timesteps = 20
) %>%
  mutate(
    mu = 1, 
    sig.m = sqrt(try.loss.expr(wfitn^2, sig.a^2, r))
  ) %>%
  uncount(n.trials) %>%
  mutate(trial.no = 1:n.trials) %>%
  split(f = 1:nrow(.))

set.seed(409)

out1 = mclapply(pars.list,
                function(pars) {
                  sim.output = sim(params = pars, theta.t = 0, init.rows = 5 * 1e5) %>%
                    group_by(gen) %>%
                    summarise(bvar = var(b_i) * (1 - 1/n())) %>%
                    mutate(trialno = pars$trial.no)
                }
) %>%
  do.call(rbind, .)

out1 %>%
  ggplot(aes(x = gen, y = bvar, group = trialno)) +
  geom_line(size = 0.1)

out1 %>%
  group_by(gen) %>%
  summarise(
    bvarbar = mean(bvar),
    bvarvar = var(bvar),
    n = n()
  ) %>%
  ggplot(aes(x = gen)) + 
  geom_line(aes(y = bvarbar)) +
  geom_ribbon(
    aes(
      ymin = bvarbar - 2 * sqrt(bvarvar / n),
      ymax = bvarbar + 2 * sqrt(bvarvar / n)
    ),
    alpha = 0.1
  )
# ugh shit fuck

###---------------------------

# okay what's going on?
# could be two things:
# (1) the expression for variance lost is wrong
# (2) the expression for variance lost is correct,
# but not being put into the model correctly 
# (e.g., needing to be counterbalanced somehow by age structure?)

# here - just see if this expression gets THE VERY FIRST generation's variance lost correct

pars.list = data.frame(
  n.pop0 = 1000, s.max = 0.9, r = (1.1 / (0.9)) - 1, wfitn = sqrt(10),
  sig.a = 1, sig.e = 0, alpha = 0.0000,
  timesteps = 2
) %>%
  mutate(
    mu = 0, 
    sig.m = 0
  ) %>%
  uncount(n.trials) %>%
  mutate(trial.no = 1:n.trials) %>%
  split(f = 1:nrow(.))

set.seed(409)

out1 = mclapply(pars.list,
                function(pars) {
                  sim.output = sim(
                    params = pars, 
                    theta.t = 0, 
                    init.rows = 5 * 1e5
                  ) %>%
                    group_by(i) %>%
                    mutate(cohort = ifelse(!min(gen), 0, gen - age)) %>%
                    group_by(gen, cohort) %>%
                    summarise(bvar = var(b_i) * (1 - 1/n())) %>%
                    mutate(trialno = pars$trial.no)
                }
) %>%
  do.call(rbind, .)

out1 %>%
  ggplot(aes(x = gen, y = bvar)) +
  geom_line(aes(group = interaction(trialno, cohort),
                colour = cohort))

out1 %>%
  filter(cohort < 2) %>%
  ggplot(aes(x = gen, y = bvar, group = trialno)) +
  geom_line(size = 0.1) +
  facet_wrap(~ cohort)

# okay look at variance lost JUST in the 0th cohort
aaa = out1 %>%
  group_by(trialno, cohort) %>%
  mutate(vardiff = c(diff(bvar), NA)) %>%
  ungroup() %>%
  filter(!is.na(vardiff)) %>%
  group_by(cohort, gen) %>%
  summarise(vdmean = mean(vardiff))
# hmmmmmmmm

# okay - think I shouldn't worry about the vardiff changing between gens
# because there is no mutation
# if I nail the mutation rate, then maybe they will end up being the same?

(with(pars.list[[1]], try.loss.expr(wfitn^2, sig.a^2, r))) # okay - this is not right
# (empirical value is -0.059, this is -0.087)

# still wrong!!!
try.loss.expr = function(w2, s2, r) {
  sumo = 0
  for (k in 0:100) { 
    sumo = sumo +
      (r/(1+r))^k * ((w2 + k*s2)^(-1) - (w2 + (k+1)*s2)^(-1))
  }
  sumo = sumo * s2^2 * w2 / (1 + r)
}

# ugh... doesn't look like an algebra problem.

testy = data.frame(age = 0:30) %>%
  mutate(sig.a = with(pars.list[[1]], sqrt(sig.a^2 * wfitn^2 / (wfitn^2 + age*sig.a^2)))) %>%
  mutate(propn = with(pars.list[[1]], (r / 1+r)^age / (1+r))) %>%
  mutate(propsiga = propn^2 * sig.a^2)

sum(testy$propsiga)
# why is this not what I observed.... as bvar_0??
# feel like I'm going crazy!!!!!

out1 %>% filter(cohort < 1, gen < 1) %>% summarise(bvar = mean(bvar))
# what the fuck??

#####----------------------

# is the issue that I'm not getting cohort variances right?
# two timesteps only here - allow me to look at first couple of gens
# for cohort 0

pars.list = data.frame(
  n.pop0 = 10000, s.max = 0.9, r = (1.1 / (0.9)) - 1, wfitn = sqrt(10),
  sig.a = 1, sig.e = 0, alpha = 0.0000,
  timesteps = 2
) %>%
  mutate(
    mu = 0, 
    sig.m = 0
  ) %>%
  uncount(n.trials) %>%
  mutate(trial.no = 1:n.trials) %>%
  split(f = 1:nrow(.))

set.seed(409)

out3 = mclapply(pars.list,
                function(pars) {
                  sim.output = sim(
                    params = pars, 
                    theta.t = 0, 
                    init.rows = 5 * 1e5
                  ) %>%
                    group_by(i) %>%
                    mutate(cohort = ifelse(!min(gen), 0, gen - age + 1)) %>%
                    group_by(gen, cohort, age) %>%
                    summarise(bvar = var(b_i)) %>%
                    mutate(trialno = pars$trial.no)
                }
) %>%
  do.call(rbind, .)

out3 %>%
  filter(!cohort, gen < 1) %>%
  group_by(age) %>%
  summarise(bvarbar = mean(bvar), n = n())
# yeahhhhh lmao something is wrong here
# these variances don't match what I thought I was constucting???
# (look below... yes they are)

testy %>% mutate(bvarexp = round(sig.a^2, 3)) %>% slice(1:10)
# okay actually this looks fine lmao I was just forgetting to square it

merge(
  out3 %>%
    filter(!cohort, gen < 1) %>%
    group_by(age) %>%
    summarise(bvarbar = mean(bvar), n = n()),
  testy %>%
    mutate(bvarexp = sig.a^2) %>%
    select(age, bvarexp)
) %>%
  ggplot(aes(x = bvarexp, y = bvarbar)) +
  geom_point() +
  geom_segment(aes(x = 0.4, y = 0.4, xend = 1, yend = 1))
# yeah that looks bout right!
# okay... so the age-by-age stuff looks right???
# in the initial cohort, the variance for each age groups is about what is expected

# but then why the heck do they not add up to what I expect as the overall variance...?

out3 %>%
  filter(!cohort, gen < 2) %>%
  group_by(age, gen) %>%
  summarise(bvarbar = mean(bvar), n = n()) %>%
  pivot_wider(names_from = gen, values_from = bvarbar) %>%
  ungroup() %>%
  mutate(expo  = `0` * with(pars.list[[1]], wfitn^2 / (wfitn^2 + sig.a^2)),
         expot = c(NA, expo[-nrow(.)]))
# compare expot (expectation based on `0`) to `1` (i.e. t = 0 to t = 1)
# they look... pretty close?

out3 %>%
  filter(!cohort, gen < 2) %>%
  group_by(age, gen) %>%
  summarise(bvarbar = mean(bvar), n = n()) %>%
  pivot_wider(names_from = gen, values_from = bvarbar) %>%
  ungroup() %>%
  mutate(expo  = `0` * with(pars.list[[1]], wfitn^2 / (wfitn^2 + sig.a^2)),
         expot = c(NA, expo[-nrow(.)])) %>%
  ggplot(aes(x = expot, y = `1`, colour = age)) +
  geom_point() +
  geom_segment(aes(x = 0.2, xend = 0.9, y = 0.2, yend = 0.9),
               inherit.aes = FALSE) +
  scale_colour_viridis_c()

# okay so I think I'm getting this right? mostly??
# hmm... very slightly off for middle range but good at the poles??? why? does that matter??

# okay so what the fuck is wrong???

#######------------------

simto = sim(pars.list[[1]] %>% mutate(timesteps = 1), theta.t = 0, init.rows = 1e4)

nrow(simto)

simto %>% filter(!(!age & gen %in% 1))

simto %>%
  filter()

#######--------------

# okay... I think I'm initializing correctly, but getting overall variances
# that are different than what I'm expecting????
# why is this?????
# in cohort 0, the age-variance breakdown w/in the cohort *looks* right
# but when I try to convert it to an overall variance for that cohort, 
# it's pretty wrong and I have no idea why
