##########
# Simulated experiment
# Here: trying to capture age structure for populations
# SN - init 6 Mar 2023
##########

# Packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(parallel)

# Clear namespace
rm(list = ls())

# Load source code
source('model_source/sim_model1_functions.R')

### Load in parameters

# Trials per parameter combo
trys.per = 200

# Parameters
pars = expand.grid(
  # (Equilibrium) size of initial cohort
  s.max = c(0.1, 0.5, 0.9),
  # Heritability of fitness
  h2    = c(.25, .5, 1)
) %>%
  # Demographic rates
  mutate(
    # Maximum expected lifetime fitness
    w.max = 3,
    # Gamma squared (pheno variance / sel pressure)
    sig.z = sqrt(.4),
    # Equilibrium lifetime fitness
    wstar = w.max * (1 - s.max) / (sqrt(1 + sig.z^2) - s.max),
    # Mean fecundity
    r     = w.max * (1 - s.max) / s.max,
    # Equilibrium population growth rate
    lstar = (s.max + w.max * (1 - s.max)) / (s.max + (w.max/wstar) * (1 - s.max)),
    # Initial population size
    n.pop0 = 20000,
    # Strength of density dependence
    alpha = log(lstar) / n.pop0,
    # Ceciling-type carrying capacity just in case
    kceil = 30000,
    p0    = (w.max * (1 - s.max)) / (w.max * (1 - s.max) + s.max)
  ) %>%
  # Genetic info
  group_by(lstar, s.max, h2, p0) %>%
  mutate(
    # Gamma-parameterization
    # wfitn = 1 in gamma parameterization
    wfitn = 1,
    # Phenotypic standard deviation in new cohorts
    sig.0 = sqrt(newt.method.g1(.1, 1e-8, s.max / lstar, r)),
    # Breeding value standard deviation in new cohorts
    sig.a = sqrt(h2 * sig.0^2),
    # Non-inherited standard dxeviation in new cohorts
    sig.e = sqrt((1-h2) * sig.0^2),
    # Population-wide breeding value standard deviation
    sig.p = sqrt(gamma.a.calc(sig.a^2, s.max / lstar, r, sig.e^2)),
    mu    = 1,
    sig.m = sqrt(wfitn^2 * (sig.a^2 - (sig.p^2 - p0*sig.a^2)/(1-p0))),
    gbar0 = 2
  ) %>%
  ungroup() %>%
  # Other junk
  mutate(timesteps = 50)

# Run simulations
set.seed(4523)

sim.out2 = mclapply(
  pars %>% uncount(trys.per) %>% mutate(try.no = 1:(nrow(.))) %>% split(.$try.no),
  function(pars) {
    sim(pars, theta.t = 0, init.rows = 50 * 30000) %>%
      mutate(
        e_i = z_i - b_i
      ) %>%
      group_by(gen, age) %>%
      summarise(
        n = n(),
        r = sum(r_i),
        bbar = mean(b_i),
        bvar = var(b_i),
        ebar = mean(e_i),
        evar = var(e_i),
        zbar = mean(z_i),
        zvar = var(z_i)
      )  %>%
      mutate(
        trial = pars$try.no, 
        p0    = pars$p0,
        lstar = pars$lstar,
        h2    = pars$h2
      )
  },
  mc.cores = 12
) %>%
  do.call(rbind, .)

# nrow(sim.out2)

write.csv(
  sim.out2,
  file = 'run_sims/out/sim_results_m1_ages.csv',
  row.names = FALSE
)

# # Age description/summary
# 
# testo = sim.out2 %>%
#   group_by(trial) %>%
#   mutate(max.gen = max(gen)) %>%
#   group_by(h2, p0, lstar) %>%
#   filter(gen <= min(max.gen)) %>%
#   select(-max.gen)
# 
# testa = testo %>%
#   group_by(h2, p0, lstar, trial, gen) %>%
#   mutate(n.total = sum(n),
#          p       = n / n.total) %>%
#   group_by(h2, p0, lstar, age, gen) %>%
#   summarise(p = sum(p) / trys.per) %>%
#   merge(
#     testo %>%
#       distinct(h2, p0, lstar, age) %>%
#       mutate(p.add = 0),
#     all.x = TRUE, all.y = TRUE
#   ) %>%
#   mutate(p = p + p.add) %>%
#   select(-p.add)
# 
# testa %>%
#   filter(lstar < 1.8) %>%
#   ggplot(aes(x = age, y = p, group = gen)) +
#   geom_line(aes(colour = gen)) +
#   scale_colour_viridis_c() +
#   facet_wrap(p0 ~ h2)
# 
# testb = testa %>%
#   # filter(lstar < 1.8) %>%
#   arrange(h2, p0, lstar, gen, age) %>%
#   group_by(h2, p0, lstar, gen) %>%
#   mutate(p.age = cumsum(p)) %>%
#   mutate(p.sta = c(0, p.age[-n()])) %>%
#   ungroup() %>%
#   mutate(age2 = ifelse(age > 4, 5, age)) 
# 
# testb %>%
#   filter(lstar < 1.8) %>%
#   ggplot(aes(x = gen)) +
#   geom_ribbon(
#     aes(
#       ymin = p.sta,
#       ymax = p.age,
#       group = age,
#       fill = age2
#     )
#   ) +
#   scale_fill_viridis_c() +
#   facet_wrap(h2 ~ p0)
# 
# testb %>%
#   filter(p0 < 0.6, lstar < 1.8) %>%
#   ggplot(aes(x = gen)) +
#   geom_ribbon(
#     aes(
#       ymin = p.sta,
#       ymax = p.age,
#       group = age,
#       fill = age2
#     ),
#     colour = 'black', linewidth = 0.05
#   ) +
#   scale_fill_viridis_c(option = 'A') +
#   facet_wrap(h2 ~ ., ncol = 3) +
#   theme(
#     panel.background = element_blank()
#   )
# 
# besta = testo %>%
#   group_by(h2, p0, lstar, age, gen) %>%
#   summarise(
#     across(contains('bar'), mean),
#     n = n()
#   ) %>%
#   pivot_longer(
#     cols = contains('bar'),
#     names_to = 'var',
#     values_to = 'val'
#   )
# 
# head(besta)
# 
# besta %>%
#   filter(lstar < 1.8, h2 < .5) %>%
#   group_by(p0, age) %>%
#   filter(all(n == trys.per)) %>%
#   ungroup() %>%
#   ggplot(aes(x = gen, y = val)) +
#   geom_point(aes(colour = age)) +
#   geom_line(
#     aes(
#       group = age,
#       colour = age
#     )
#   ) +
#   scale_colour_viridis_c() +
#   facet_wrap(var ~ p0)
# 
# cesta = testo %>%
#   mutate(cohort = gen - age) %>%
#   group_by(h2, p0, lstar, cohort, gen) %>%
#   summarise(
#     across(contains('bar'), mean),
#     n = n()
#   ) %>%
#   pivot_longer(
#     cols = contains('bar'),
#     names_to = 'var',
#     values_to = 'val'
#   )
# 
# head(cesta)
#   
# cesta %>%
#   filter(lstar < 1.8, h2 < .5) %>%
#   group_by(p0, cohort, var) %>%
#   filter(n() > 1) %>% 
#   #filter(all(n == 5)) %>%
#   ungroup() %>%
#   ggplot(aes(x = gen, y = val)) +
#   geom_line(
#     aes(
#       group  = cohort,
#       colour = cohort 
#     )
#   ) +
#   scale_colour_gradient2(low = 'skyblue', 
#                          high = 'green', mid = 'black', 
#                          midpoint = 0) +
#   # scale_colour_viridis_c() +
#   facet_wrap(var ~ p0)
# 
# cesta %>%
#   filter(lstar < 1.8, h2 < .5) %>%
#   group_by(p0, cohort, var) %>%
#   filter(n() > 2) %>% 
#   ungroup() %>%
#   ggplot(aes(x = gen, y = val)) +
#   geom_point(
#     aes(
#       colour = cohort,
#       shape = var,
#       alpha = n
#     )
#   ) +
#   scale_shape_manual(values = c(1, 19)) +
#   scale_colour_gradient2(low = 'purple',
#                          high = 'orange', mid = 'black',
#                          midpoint = 0) +
#   facet_wrap( ~ p0)
