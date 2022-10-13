# First pass script for brute-forcing our way to back-engineered mutation rates
# to get mutation-selection balance
# For now just do mutation, assuming that h^2 = 1, or that sig.e = 0
# SN - 13 October 2022

library(ggplot2)
library(dplyr)
library(tidyr)
library(parallel)

rm(list = ls())

source('model_source/sim_model1_functions.R')

# Initial parameters
sm.init = c(0, 2.5)

ntrials = 500 # trials per parameter
target  = 1
hi.mean = 0
hi.se   = 0
lo.mean = 0
lo.se   = 0

sm.try = sm.init

while (
  (abs(hi.mean - target) > hi.se) &
  (abs(lo.mean - target) > lo.se)
) {

  cat("\nRates: ")
  cat(sm.try)
  
  pars.list = expand.grid(sig.m = sm.try, n.trials = 1:ntrials) %>%
    mutate(
      n.pop0 = 1000,
      l.max = 1.2,
      sig.a = 1, sig.e = 0, alpha = 0,
      s.max = 0.9, r = (l.max / (0.9)) - 1, 
      wfitn = sqrt( (sig.a^2 + sig.e^2) / ((l.max^2) - 1) ),
      mu    = 1,
      kceil = 5000,
      timesteps = 20
    ) %>%
    mutate(trial = 1:nrow(.)) %>%
    split(f = .$trial)

  out1 = mclapply(pars.list,
                  function(x) {
                    sim(params = x, theta.t = 0, init.rows = 5000 * 31) %>%
                      group_by(gen) %>%
                      summarise(
                        n = n(),
                        bbar = mean(b_i),
                        bvar = var(b_i) * (1 - 1/n)
                      ) %>%
                      mutate(
                        trial.no = x$trial
                      )
                  },
                  mc.cores = 6
                ) %>%
    do.call(what = rbind) %>%
    merge(
      y = do.call(rbind, pars.list) %>% select(trial, sig.m),
      by.x = 'trial.no', by.y = 'trial'
    ) %>%
    group_by(gen, sig.m) %>%
    summarise(
      n = n(),
      bvarbar = mean(bvar),
      bvarvar = var(bvar)
    ) %>%
    filter(gen > pars.list[[1]]$timesteps / 2)

  time.avg = out1 %>%
    mutate(
      bvarste = sqrt(bvarvar / n)
    ) %>%
    group_by(sig.m) %>%
    summarise_at(
      vars(bvarbar, bvarste),
      mean
    )
  
  hi.mean = time.avg %>%
    filter(sig.m %in% sm.try[2]) %>%
    pull(bvarbar)
  hi.se   = time.avg %>%
    filter(sig.m %in% sm.try[2]) %>%
    pull(bvarste)
  cat(paste("\nhigh mean and se", round(hi.mean, 2), round(hi.se, 2)))
  lo.mean = time.avg %>%
    filter(sig.m %in% sm.try[1]) %>%
    pull(bvarbar)
  lo.se   = time.avg %>%
    filter(sig.m %in% sm.try[1]) %>%
    pull(bvarste)
  cat(paste("\nlow mean and se", round(lo.mean, 2), round(lo.se, 2)))
  
  cat("")
  
  # 
  if (
    all(out1$bvarbar[out1$sig.m %in% sm.try[1]] < target) &
    all(out1$bvarbar[out1$sig.m %in% sm.try[2]] > target)
    ) {
      if (
        mean(abs(out1$bvarbar[out1$sig.m %in% sm.try[1]] - target)) > 
        mean(abs(out1$bvarbar[out1$sig.m %in% sm.try[2]] - target))
      ) {
        sm.try[1] = sm.try[1] + (1/2) * abs(diff(sm.try))
      } else {
        sm.try[2] = sm.try[2] - (1/2) * abs(diff(sm.try))
      }
    } else if (
      all(out1$bvarbar[out1$sig.m %in% sm.try[1]] > target)
      ) {
        sm.try[1] = sm.try[1] - (1/2) * abs(diff(sm.try))
    } else if (
      all(out1$bvarbar[out1$sig.m %in% sm.try[2]] < target)
    ) {
        sm.try[2] = sm.try[2] + (1/2) * abs(diff(sm.try))
    } else {
      cat("\noh god oh fuck")
      break
    }
  
}

# out1 %>%
#   ggplot(aes(x = gen, y = n)) +
#   geom_line(
#     aes(
#       group = trial.no,
#       colour = factor(sig.m^2)
#     ) 
#   ) +
#   scale_y_log10()
# # cool... makes sense?
# # ah probably the breeding value variance is gigantic lmoa
# 
# out1 %>%
#   ggplot(aes(x = gen, y = bvar)) +
#   geom_line(
#     aes(
#       group = trial.no,
#       colour = factor(sig.m^2)
#     ) 
#   )
# 
# out1 %>%
#   group_by(gen, sig.m) %>%
#   summarise(
#     n = n(),
#     bvarbar = mean(bvar),
#     bvarvar = var(bvar)
#   ) %>%
#   ggplot(aes(x = gen)) +
#   geom_line(
#     aes(
#       y = bvarbar,
#       colour = factor(sig.m^2)
#     )
#   ) +
#   geom_ribbon(
#     aes(
#       ymin = bvarbar - 2 * sqrt(bvarvar / n),
#       ymax = bvarbar + 2 * sqrt(bvarvar / n),
#       fill = factor(sig.m^2)
#     ),
#     alpha = 0.1
#   )

