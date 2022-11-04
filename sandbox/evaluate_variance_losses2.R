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
  mutate(propn = with(pars.list[[1]], (r / (1+r))^age / (1+r)),
         propo = with(pars.list[[1]], (1 / (1+r))^age / (1 + (1/r)))) %>%
  mutate(propsiga = propn^2 * sig.a^2,
         proosiga = propo^2 * sig.a^2)

sum(testy$propsiga)
sum(testy$proosiga) # AHHH WHAT IS GOING ON... the proportions is right but weighted sum is dogshit

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
                    mutate(cohort = ifelse(!min(gen), 0, gen - age)) %>%
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
    filter(!cohort, !gen) %>%
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

# out3 %>%
#   filter(!cohort, gen < 2) %>%
#   group_by(age, gen) %>%
#   summarise(bvarbar = mean(bvar), n = n()) %>%
#   pivot_wider(names_from = gen, values_from = bvarbar) %>%
#   ungroup() %>%
#   mutate(expo  = `0` * with(pars.list[[1]], wfitn^2 / (wfitn^2 + sig.a^2)),
#          expot = c(NA, expo[-nrow(.)])) %>%
#   ggplot(aes(x = expot, y = `1`, colour = age)) +
#   geom_point() +
#   geom_segment(aes(x = 0.2, xend = 0.9, y = 0.2, yend = 0.9),
#                inherit.aes = FALSE) +
#   scale_colour_viridis_c()

set.seed(409)

out4 = mclapply(pars.list,
                function(pars) {
                  sim.output = sim(
                    params = pars, 
                    theta.t = 0, 
                    init.rows = 5 * 1e5
                  ) %>%
                    group_by(i) %>%
                    mutate(cohort = ifelse(!min(gen), 0, gen - age)) %>%
                    group_by(gen, cohort) %>%
                    summarise(bvar = var(b_i)) %>%
                    mutate(trialno = pars$trial.no)
                }
) %>%
  do.call(rbind, .)

out4 %>%
  filter(!cohort, gen < 2) %>%
  group_by(gen) %>%
  summarise(bvarbar = mean(bvar), n = n())

# 0.752 - with(pars.list[[1]], try.loss.expr(wfitn^2, sig.a^2, r))
# IT WORKS
# VAR LOSS EXPRESSION WORKS




# okay so I think I'm getting this right? mostly??
# hmm... very slightly off for middle range but good at the poles??? why? does that matter??

# okay so what the fuck is wrong???

#######------------------

simto = sim(pars.list[[1]] %>% mutate(timesteps = 1), theta.t = 0, init.rows = 1e4)

nrow(simto)

simto %>% filter(!(!age & gen %in% 1))

# just looking at cohort 0...
cohort0.agevars = simto %>%
  filter(!(!age & gen %in% 1)) %>%
  group_by(gen, age.cohort = age - gen) %>%
  summarise(bvar = var(b_i)) 

cohort0.agevars %>%
  filter(!is.na(bvar)) %>%
  ggplot(aes(x = age.cohort, y = bvar)) +
  geom_point(aes(colour = gen))

cohort0.agevars %>%
  filter(!is.na(bvar)) %>%
  pivot_wider(names_from = gen, values_from = bvar) %>%
  ggplot(aes(x = `0`, y = `1`)) +
  geom_segment(aes(x = 0, xend = 0.9, y = 0, yend = 0.9),
               inherit.aes = FALSE) +
  geom_point(aes(fill = age.cohort), shape = 21, size = 4) +
  scale_fill_viridis_c()

# is the sum(p^2 sig.a^2) formula just the wrong way to do variance??

cohort0.agevars = simto %>%
  # Isolate cohort 0
  filter(!(!age & gen %in% 1)) %>%
  # Get prob(age) for each age in cohort in each time step
  group_by(gen) %>%
  mutate(n.coh.gen = n()) %>%
  group_by(gen, age.cohort = age - gen) %>%
  summarise(
    n.age = n(),
    p.age = n.age / n.coh.gen,
    bvar = var(b_i)
  ) %>%
  # idk why this is needed lmao
  distinct(gen, age.cohort, .keep_all = TRUE)

cohort0.agevars %>%
  ungroup() %>%
  filter(!gen) %>%
  summarise(should.be.1 = sum(p.age))
# thank god

cohort0.agevars %>%
  filter(!gen) %>%
  mutate(exp.p = with(pars.list[[1]], ((1 / (1 + r))^age.cohort) / (1 + (1/r)))) %>%
  ggplot(aes(x = age.cohort)) +
  geom_line(aes(y = p.age)) +
  geom_line(aes(y = exp.p), colour = 'blue')
# ohhhh lmao I've been calculating the var wrong LMAO LMAO LMAO
# TPROB OF BEING IN AGE IS PROPORTIONAL TO 1/(1+r), NOT r/(1+r)

cohort0.agevars %>%
  filter(!gen) %>%
  ungroup() %>%
  mutate(xxx = p.age^2 * bvar) %>%
  summarise(xxy = sum(xxx, na.rm = TRUE))
# okay yeah... yeah... fuck this is fucked
# not even on the same order of magnitude
# what the fuck...

# wait... 
cohort0.agevars %>%
  filter(!gen) %>%
  ungroup() %>%
  mutate(xxx = p.age * bvar) %>%
  summarise(xxy = sum(xxx, na.rm = TRUE))
# this looks right...? just using p as weight instead of p^2?
# HAVE I BEEN MISUNDERSTANDING VARIANCE THIS WHOLE TIME???

simto %>%
  slice(1:pars.list[[1]]$n.pop0) %>%
  summarise(varby = var(b_i))
# mother of christ

try.loss.expr = function(w2, s2, r) {
  sumo = 0
  for (k in 0:100) { 
    sumo = sumo +
      (1/(1+r))^k * ((w2 + k*s2)^(-1) - (w2 + (k+1)*s2)^(-1))
  }
  sumo = sumo * s2 * w2 * r / (1 + r)
}

(with(pars.list[[1]], try.loss.expr(wfitn^2, sig.a^2, r)))
# hmm...
# WAIT THIS IS VERY CLOSE TO WHAT WE SAW ABOVE???

######--------

# TRY ONE MORE TIME BUT WITH THE LOSS EXPRESSION CORRECTED
# (for the actual age stuff)

try.loss.expr = function(w2, s2, r) {
  sumo = 0
  for (k in 0:100) { 
    sumo = sumo +
      (1/(1+r))^k * ((w2 + k*s2)^(-1) - (w2 + (k+1)*s2)^(-1))
  }
  sumo = sumo * s2 * w2 * r / (1 + r)
}

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

# shit... thought that would really work



###### ----------
# NEXT CHECKS
# - is my fixed try.var.loss() or whatever fun actually correct? looks closer at least
# - if it is, is there an implementation error? (error 2 instead of error 1, above)


# var.loss expression works!
# trying one mutation rate...

muto.fun = function(w2, s2, r) {
  sumo = 0
  for (k in 0:1000) {
    sumo = sumo +
      (1 / (1+r))^k * ( (w2 + k*s2)^(-1) - (w2 + (k+1)*s2)^(-1) )
  }
  sumo = sumo * s2 * w2
  return(sumo)
}

pars.list = data.frame(
  n.pop0 = 1000, s.max = 0.9, r = (1.1 / (0.9)) - 1, wfitn = sqrt(10),
  sig.a = 1, sig.e = 0, alpha = 0.0000,
  kceil = 10000,
  timesteps = 20
) %>%
  mutate(
    mu = 1, 
    sig.m = sqrt(muto.fun(wfitn^2, sig.a^2, r))
  ) %>%
  uncount(n.trials) %>%
  mutate(trial.no = 1:n.trials) %>%
  split(f = 1:nrow(.))

set.seed(409)

simmy = sim(pars.list[[1]], theta.t = 0, init.rows = 5 * 1e5)

cohort.vars = simmy %>%
  group_by(i) %>%
  mutate(cohort = ifelse(!min(gen), 0, gen - age)) %>%
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

simmy %>%
  group_by(gen) %>%
  summarise(bvar = var(b_i)) %>%
  ggplot(aes(x = gen, y = bvar)) + 
  geom_line()

simmy %>%
  group_by(gen) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = gen, y = n)) +
  geom_line() +
  scale_y_log10()

# hmmm... still losing some ahhhh

out5 = mclapply(pars.list,
                function(pars) {
                  sim.output = sim(params = pars, theta.t = 0, init.rows = 5 * 1e5) %>%
                    group_by(gen) %>%
                    summarise(bvar = var(b_i) * (1 - 1/n())) %>%
                    mutate(trialno = pars$trial.no)
                }
) %>%
  do.call(rbind, .)

out5 %>%
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
# noooooooooooooooooo
# really thought I had it...

cohort.vars %>%
  # mutate(r = pars.list[[1]]$r) %>%
  # mutate(cohr.weighted.bv = (r+1)/r * (bvar * (1/(1+r)))^(gen-cohort)) %>%
  # group_by(gen) %>%
  # mutate(mean.weighted.bv = sum(cohr.weighted.bv)) %>%
  ggplot(aes(x = gen, y = bvar)) +
  geom_line(aes(group = cohort, colour = cohort)) +
  # geom_line(aes(y = cohr.weighted.bv), size = 2,
  #           data = . %>% distinct(gen, .keep_all = TRUE)) +
  geom_point(aes(colour = cohort)) +
  geom_line(
    data = simmy %>% group_by(gen) %>% summarise(bvar = var(b_i)),
    size = 2
  ) +
  scale_colour_viridis_c()
# hmm......... the decline is really not that big... but it's still there...

### Newton's method to solve for initial?

# f we are solving
f_fun = function(x, w, s, r) {
  sumo = 0
  for (k in 0:1000) {
    sumo = sumo + 
      (1/(1+r))^k * (w^2 * x^2) / (w^2 + k * x^2)
  }
  sumo = sumo * (r / (1 + r))
  return(sumo - s^2)
  # NOTE sumo aka x is sigma, NOT sigma^2
}

# derivative of f wrt x
f_prm = function(x, w, r) {
  sumo = 0
  for (k in 0:1000) {
    nasty.quot = (2*x*w^2 * (w^2 + k*x^2) - (x^2 * w^2) * (2*k*x)) / (w^2 + k*x^2)^2
    sumo = sumo + ((1/(1+r))^k * nasty.quot)
  }
  sumo = sumo * (r / (1 + r))
  return(sumo)
}

newt.tol = function(x0, tol, w, s, r) {
  # xold = x0
  # xnew = xold - (f_fun(xold, w2, s2, r) / f_prm(xold, w2, r))
  # while (abs(xnew - xold) > tol) {
  #   xold = xnew
  #   xnew = xold - (f_fun(xold, w2, s2, r) / f_prm(xold, w2, r))
  # }
  # return(xnew)
  
  xold = x0
  
  while(abs(f_fun(xold, w, s, r)) > tol) xold = xold - (f_fun(xold, w, s, r) / f_prm(xold, w, r))
  
  return(xold)
  
}

(lala = newt.tol(1, 0.000001, pars.list[[1]]$wfitn, 1, pars.list[[1]]$r))

f_fun(lala, pars.list[[1]]$wfitn, 1, pars.list[[1]]$r) + 1

# Mutation rate again
muto.fun = function(w, s, r) {
  sumo = 0
  for (k in 0:1000) {
    sumo = sumo +
      (1 / (1+r))^k * ( (w^2 + k*s^2)^(-1) - (w^2 + (k+1)*s^2)^(-1) )
  }
  sumo = sumo * s^2 * w^2
  return(sumo)
  # sumo is sigma_m^2, not sigma_m
}

###### --------------

n.trials = 400

pars.list = data.frame(
  n.pop0 = 1000, s.max = 0.9, r = (1.1 / (0.9)) - 1, wfitn = sqrt(10),
  sig.a = 1,
  sig.e = 0, alpha = 0.0000,
  kceil = 10000,
  timesteps = 20
) %>%
  mutate(
    sig.a = newt.tol(1, 1e-6, wfitn, sig.a, r)
  ) %>%
  mutate(
    mu = 1, 
    sig.m = sqrt(muto.fun(wfitn, sig.a, r))
  ) %>%
  uncount(n.trials) %>%
  mutate(trial.no = 1:n.trials) %>%
  split(f = 1:nrow(.))

set.seed(909202)

out5 = mclapply(pars.list,
                function(pars) {
                  sim.output = sim(params = pars, theta.t = 0, init.rows = 5 * 1e5) %>%
                    group_by(gen) %>%
                    summarise(bvar = var(b_i) * (1 - 1/n())) %>%
                    mutate(trialno = pars$trial.no)
                }
) %>%
  do.call(rbind, .)

set.seed(409)

simmy = sim(pars.list[[1]], theta.t = 0, init.rows = 5 * 1e5)

cohort.vars = simmy %>%
  group_by(i) %>%
  mutate(cohort = ifelse(!min(gen), 0, gen - age)) %>%
  group_by(gen, cohort) %>%
  summarise(bvar = var(b_i))

cohort.vars %>%
  ggplot(aes(x = gen, y = bvar)) +
  geom_line(aes(group = cohort, colour = cohort)) +
  geom_line(data = simmy %>% group_by(gen) %>% summarise(bvar = var(b_i)),
            aes(x = gen, y = bvar),
            linetype = 2) +
  geom_point(aes(colour = cohort)) +
  scale_colour_viridis_c()

cohort.vars %>%
  ggplot(aes(x = gen - cohort, y = bvar)) +
  geom_line(aes(group = cohort, colour = cohort)) +
  geom_point(aes(colour = cohort)) +
  scale_colour_viridis_c()

simmy %>%
  filter(gen < 16) %>%
  group_by(i) %>%
  mutate(cohort = ifelse(!min(gen), 0, gen - age)) %>%
  ggplot(aes(x = age)) +
  geom_density(aes(colour = gen, group = gen),
               size = 0.5) +
  scale_colour_viridis_c() +
  facet_wrap(~ cohort)
# wait... what??
# age distribution... is geometric but only past age 0?
# I guess age structure/fitness structure influences # offspring...
# wait... no? number of offspring should be the same throughout?

out5 %>%
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
# this is very, very close... like so unbelievably close!


#######--------------

# well... thought I figured something out but maybe I didn't
# at least it looks like I learned corrected these previously-wrong misunderstandings:
#   - variance really is sum over p*sig^2, not p^2*sig^2
#   - probability of being in class k is proportional to 1/(1+r)^k, not (r/(1+r))^k

# update 11/3/22 - I think it's working now???
# or at least much much closer than ever before (in this simple case)
# 

