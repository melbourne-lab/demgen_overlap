# Script a 

library(parallel)
library(ggplot2)
library(dplyr)
library(tidyr)

rm(list = ls())

source('model_source/sim_model1_functions.R')

lambda.max = 1.2

pars = expand.grid(sig.a = sqrt(c(2, 4)), s.max = c(0.5, 0.9), local.k = 1:50)  %>%
  mutate(
    n.pop0 = 5000, 
    wfitn = 10,
    mu    = 0, sig.m = sqrt(0),
    timesteps = 75,
    k = 1:nrow(.)
  ) %>%
  mutate(
    # params dependent on other params...
    sig.e = sqrt(10 - sig.a^2),
    r     = (lambda.max / s.max) - 1,
    alpha = (log(2)/2 - log(lambda.max)) / n.pop0
  ) %>%
  split(f = 1:nrow(.))

set.seed(409)

sissy.time = system.time(
  outp <- mclapply(pars, 
                function(par.row) {
                  sim(params = par.row, theta.t = 8, init.rows = 1e6) %>%
                    group_by(gen) %>%
                    summarise(
                      n = n(),
                      b = mean(theta_t - b_i),
                      z = mean(theta_t - z_i),
                      e = mean(z_i - b_i),
                    ) %>%
                    mutate(k = par.row$k)
                },
                mc.cores = 6)
)

sissy.time

outd = do.call(rbind, outp)

head(outd)
nrow(outd)
table(outd$k)
tail(outd)

outp = merge(
  x = outd,
  y = do.call(rbind, pars) %>% select(k, s.max, sig.a)
)

outp %>%
  ggplot(aes(x = gen, y = n)) +
  geom_line(aes(group = k, colour = factor(s.max), linetype = factor(sig.a^2)), size = 0.1) +
  scale_y_log10()
# results seem similar to before...


outn = outp %>%
  select(k, gen, n, s.max, sig.a) %>%
  rbind(
    expand.grid(
      k = 1:length(pars), gen = 0:75, n = 0 
    ) %>%
    merge(y = do.call(rbind, pars) %>% select(k, s.max, sig.a))
  ) %>%
  group_by(k, gen, s.max, sig.a) %>%
  summarise(n = sum(n)) %>%
  group_by(gen, s.max, sig.a) %>%
  summarise(
    nbar = mean(n),
    n.se = sqrt(var(n) / n())
  )

outn %>%
  ggplot(
    aes(
      x = gen, 
      colour = factor(s.max), 
      fill   = factor(s.max),
      linetype = factor(sig.a^2))) +
  geom_line(
    aes(
      y = nbar
    )
  ) +
  geom_ribbon(
    aes(
      ymin = nbar - 2 * n.se,
      ymax = nbar + 2 * n.se
    ),
    alpha = 0.2
  )

# hell mf yeah
# look at that... longevity does not affect time to reach fitness = 1 (fuck)
# but does affect depth reached (variance?)
# heritability appears to influence both depth reached (possibly through pop size?) 
# and also time until ...


##
##
## old (when sims were not summarized after being run)
##
##
outd %>%
  group_by(k, gen) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = gen, y = n)) +
  geom_line(aes(group = k, colour = factor(s.max), linetype = factor(sig.a^2))) +
  scale_y_log10()

outd %>%
  group_by(k, gen) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = gen, y = n)) +
  geom_line(aes(group = k, colour = factor(s.max), linetype = factor(sig.a^2))) +
  scale_y_log10()

