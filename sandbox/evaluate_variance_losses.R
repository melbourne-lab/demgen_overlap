### Sandbox
### Script looking at losses in genetic, non-genetic comp., and phenotypic variance over time
### trying to find a way to offset loss of variance through mutation
### SN

library(parallel)
library(ggplot2)

source('model_source/sim_model1_functions.R')

n.trials = 100

pars.list = data.frame(
  n.pop0 = 1000, s.max = 0.9, r = (1.1 / (0.9)) - 1, wfitn = sqrt(10),
  sig.a = 1, sig.e = sqrt(4), alpha = 0.0000,
  timesteps = 20
) %>%
  mutate(
    mu = 1, 
    sig.m = sqrt((sig.a^2)^2 / (wfitn^2 + sig.a^2 + sig.e^2) / mu)
  ) %>%
  mutate(trial.no = 1:nrow(.)) %>%
  uncount(n.trials) %>%
  split(f = 1:nrow(.))

set.seed(409)

out1 = mclapply(pars.list,
                
                function(pars) { 
                  
                  sim.output = sim(params = pars, theta.t = 0, init.rows = 5 * 1e5) %>%
                    mutate(e_i = z_i - b_i)
                  
                  s1 = sim.output %>%
                    group_by(gen) %>%
                    summarise(
                      glob.n = n(),
                      glob.bbar = mean(b_i),
                      glob.zbar = mean(z_i),
                      glob.ebar = mean(e_i),
                      glob.bvar = var(b_i) * (1 - 1/glob.n),
                      glob.zvar = var(z_i) * (1 - 1/glob.n),
                      glob.evar = var(e_i) * (1 - 1/glob.n),
                      glob.ecov = cov(b_i, e_i)
                    )
                  
                  s2 = sim.output %>%
                    group_by(gen, cohort = gen - age + 1) %>%
                    summarise(
                      n = n(),
                      bbar = mean(b_i),
                      zbar = mean(z_i),
                      ebar = mean(e_i),
                      bvar = var(b_i) * (1 - 1/n),
                      zvar = var(z_i) * (1 - 1/n),
                      evar = var(e_i) * (1 - 1/n),
                      ecov = cov(b_i, e_i)
                    )
                    
                  return(
                    merge(s1, s2, by = 'gen') %>% mutate(k = pars$trial.no)
                  )
              
                }
) %>%
  do.call(rbind, .)

out1 %>%
  group_by(gen, k) %>%
  summarise(zvar.bar = mean(zvar)) %>%
  ggplot(aes(x = gen, y = zvar.bar)) +
  geom_line(aes(colour = factor(k), group = k))
# hmm... what the fuck do I do about this...

out1 %>%
  group_by(gen, k) %>%
  summarise(bvar.bar = mean(bvar)) %>%
  ggplot(aes(x = gen, y = bvar.bar)) +
  geom_line(aes(colour = factor(k), group = k))

out1 %>%
  group_by(gen, k) %>%
  summarise(evar.bar = mean(evar)) %>%
  ggplot(aes(x = gen, y = evar.bar)) +
  geom_line(aes(colour = factor(k), group = k))
# well... at least this isn't impacted
# but still. wtf.

out1 %>%
  filter(k %in% 1, cohort %in% 1) %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = bvar#,
      #group = cohort
    )
  )

# lmao damn wtf fuck

sim.output = sim(params = pars.list[[1]], theta.t = 0, init.rows = 5 * 1e5) %>%
  mutate(e_i = z_i - b_i)

s1 = sim.output %>%
  group_by(gen) %>%
  summarise(
    glob.n = n(),
    glob.bvar = var(b_i) * (1 - 1/glob.n),
    glob.zvar = var(z_i) * (1 - 1/glob.n),
    glob.evar = var(e_i) * (1 - 1/glob.n),
    glob.ecov = cov(b_i, e_i)
  )

s2 = sim.output %>%
  group_by(gen, cohort = gen - age + 1) %>%
  summarise(
    n = n(),
    bvar = var(b_i) * (1 - 1/n),
    zvar = var(z_i) * (1 - 1/n),
    evar = var(e_i) * (1 - 1/n),
    ecov = cov(b_i, e_i)
  )

s3 = merge(s2, s1, by = 'gen')

head(s3)
tail(s3)
nrow(s3)
# this... looks right...?

s3 %>%
  ggplot(aes(x = gen)) +
  geom_line(aes(y = bvar, group = cohort)) +
  geom_line(aes(y = glob.bvar), linetype = 2, size = 1.2)
# why does this work,..?

s3 %>%
  ggplot(aes(x = gen)) +
  geom_line(aes(y = evar, group = cohort)) +
  geom_line(aes(y = glob.evar), linetype = 2, size = 1.2) +
  geom_point(aes(y = evar), size = 0.8)

s3 %>%
  mutate(age = gen - cohort + 1) %>%
  ggplot(aes(x = gen)) +
  geom_line(aes(y = evar, group = age)) +
  geom_line(aes(y = glob.evar), linetype = 2, size = 1.2) +
  geom_point(aes(y = evar), size = 0.8)


# hmm...
# this kinda looks (sorta) like the val in gen 2 is it...
# (wait... but then survival wouldn't matter??)

# what if I just... run a bunch of sims, see the equilibrium rate, and just init a population at some inverse of that proportion?

rm(list = ls())

source('model_source/sim_model1_functions.R')

trials.per = 50

pars.list = expand.grid(
  sig.e = sqrt(2^(1:4) - 1),
  trial = 1:trials.per
) %>%
  mutate(
    n.pop0 = 1000, s.max = 0.9, r = (1.2 / (0.9)) - 1, wfitn = sqrt(10),
    sig.a = 1, alpha = 0.0000,
    timesteps = 30
) %>%
  mutate(
    mu = 1, 
    sig.m = sqrt((sig.a^2)^2 / (wfitn^2 + sig.a^2 + sig.e^2) / mu)
  ) %>%
  mutate(trial.glob = 1:nrow(.)) %>%
  split(f = 1:nrow(.))

set.seed(333)

out1 = mclapply(pars.list,
                
                function(pars) { 
                  
                  sim(params = pars, theta.t = 0, init.rows = 1e5) %>%
                    mutate(e_i = z_i - b_i) %>%
                    group_by(gen) %>%
                    summarise(
                      n = n(),
                      bbar = mean(b_i),
                      zbar = mean(z_i),
                      ebar = mean(e_i),
                      bvar = var(b_i) * (1 - 1/n),
                      zvar = var(z_i) * (1 - 1/n),
                      evar = var(e_i) * (1 - 1/n),
                      ecov = cov(b_i, e_i)
                    ) %>%
                    mutate(trial.no = pars$trial.glob[1])

                },
                mc.cores = 6
) %>%
  do.call(rbind, .)

head(out1)
nrow(out1)

out2 = out1 %>%
  pivot_longer(n:ecov, names_to = 'vartype', values_to = 'varval') %>%
  merge(y = do.call(rbind, pars.list) %>% 
          select(trial.glob, sig.e) %>%
          mutate(sig.e = round(sig.e, 2)) %>%
          rename(trial.no = trial.glob),
        by = 'trial.no') %>%
  group_by(sig.e, gen, vartype) %>%
  summarise(
    var.mean = mean(varval),
    var.se   = sd(varval),
    n.obs    = n()
  )

out2 %>%
  filter(vartype %in% 'evar') %>%
  mutate(sig.e = factor(sig.e)) %>%
  ggplot(aes(x = gen, group = sig.e)) +
  geom_line(aes(y = var.mean, colour = sig.e)) +
  scale_y_log10()

out2 %>%
  filter(vartype %in% 'evar') %>%
  mutate(sig.e = factor(sig.e)) %>%
  group_by(sig.e) %>%
  mutate(scaled.mean = var.mean / var.mean[1]) %>%
  ggplot(aes(x = gen, group = sig.e)) +
  geom_line(aes(y = scaled.mean, colour = sig.e))

# CONCLUSION
# (with mutations at least), 
# the equilibrium portion of initial ENVIRONMENTAL variance depends on sig.e
# (heritability maybe? or phenotypic var...)

# try now with phenotypic var constant, heritability varying

pars.list = expand.grid(
  sig.a = sqrt(1:4),
  trial = 1:trials.per
) %>%
  mutate(
    sig.e = sqrt(10 - sig.a),
    n.pop0 = 1000, s.max = 0.9, r = (1.2 / (0.9)) - 1, wfitn = sqrt(10),
    alpha = 0.0000,
    timesteps = 30
  ) %>%
  mutate(
    mu = 1, 
    sig.m = sqrt((sig.a^2)^2 / (wfitn^2 + sig.a^2 + sig.e^2) / mu)
  ) %>%
  mutate(trial.glob = 1:nrow(.)) %>%
  split(f = 1:nrow(.))

set.seed(333)

out1 = mclapply(pars.list,
                
                function(pars) { 
                  
                  sim(params = pars, theta.t = 0, init.rows = 1e5) %>%
                    mutate(e_i = z_i - b_i) %>%
                    group_by(gen) %>%
                    summarise(
                      n = n(),
                      bbar = mean(b_i),
                      zbar = mean(z_i),
                      ebar = mean(e_i),
                      bvar = var(b_i) * (1 - 1/n),
                      zvar = var(z_i) * (1 - 1/n),
                      evar = var(e_i) * (1 - 1/n),
                      ecov = cov(b_i, e_i)
                    ) %>%
                    mutate(trial.no = pars$trial.glob[1])
                  
                },
                mc.cores = 6
) %>%
  do.call(rbind, .)

head(out1)
nrow(out1)

out2 = out1 %>%
  pivot_longer(n:ecov, names_to = 'vartype', values_to = 'varval') %>%
  merge(y = do.call(rbind, pars.list) %>% 
          select(trial.glob, sig.e) %>%
          mutate(sig.e = round(sig.e, 2)) %>%
          rename(trial.no = trial.glob),
        by = 'trial.no') %>%
  group_by(sig.e, gen, vartype) %>%
  summarise(
    var.mean = mean(varval),
    var.se   = sd(varval),
    n.obs    = n()
  )

out2 %>%
  filter(vartype %in% 'evar') %>%
  mutate(sig.e = factor(sig.e)) %>%
  ggplot(aes(x = gen, group = sig.e)) +
  geom_line(aes(y = var.mean, colour = sig.e)) +
  scale_y_log10()

out2 %>%
  filter(vartype %in% 'evar') %>%
  mutate(sig.e = factor(sig.e)) %>%
  group_by(sig.e) %>%
  mutate(scaled.mean = var.mean / var.mean[1]) %>%
  ggplot(aes(x = gen, group = sig.e)) +
  geom_line(aes(y = scaled.mean, colour = sig.e))

# okay - here, the pattern is still present

### Try phenotypic variance changing, heritability same?

pars.list = expand.grid(
  sig.z = sqrt(((1:4)/2) * 10),
  trial = 1:trials.per
) %>%
  mutate(
    sig.a = sqrt((.25) * sig.z),
    sig.e = sqrt((.75) * sig.z),
    n.pop0 = 1000, s.max = 0.9, r = (1.2 / (0.9)) - 1, wfitn = sqrt(10),
    alpha = 0.0000,
    timesteps = 25
  ) %>%
  mutate(
    mu = 1, 
    sig.m = sqrt((sig.a^2)^2 / (wfitn^2 + sig.a^2 + sig.e^2) / mu)
  ) %>%
  mutate(trial.glob = 1:nrow(.)) %>%
  split(f = 1:nrow(.))

set.seed(333)

out1 = mclapply(pars.list,
                
                function(pars) { 
                  
                  sim(params = pars, theta.t = 0, init.rows = 1e5) %>%
                    mutate(e_i = z_i - b_i) %>%
                    group_by(gen) %>%
                    summarise(
                      n = n(),
                      bbar = mean(b_i),
                      zbar = mean(z_i),
                      ebar = mean(e_i),
                      bvar = var(b_i) * (1 - 1/n),
                      zvar = var(z_i) * (1 - 1/n),
                      evar = var(e_i) * (1 - 1/n),
                      ecov = cov(b_i, e_i)
                    ) %>%
                    mutate(trial.no = pars$trial.glob[1])
                  
                },
                mc.cores = 6
) %>%
  do.call(rbind, .)

head(out1)
nrow(out1)

out2 = out1 %>%
  pivot_longer(n:ecov, names_to = 'vartype', values_to = 'varval') %>%
  merge(y = do.call(rbind, pars.list) %>% 
          select(trial.glob, sig.e) %>%
          mutate(sig.e = round(sig.e, 2)) %>%
          rename(trial.no = trial.glob),
        by = 'trial.no') %>%
  group_by(sig.e, gen, vartype) %>%
  summarise(
    var.mean = mean(varval),
    var.se   = sd(varval),
    n.obs    = n()
  )

out2 %>%
  filter(vartype %in% 'evar') %>%
  mutate(sig.e = factor(sig.e)) %>%
  ggplot(aes(x = gen, group = sig.e)) +
  geom_line(aes(y = var.mean, colour = sig.e)) +
  scale_y_log10()

out2 %>%
  filter(vartype %in% 'evar') %>%
  mutate(sig.e = factor(sig.e)) %>%
  group_by(sig.e) %>%
  mutate(scaled.mean = var.mean / var.mean[1]) %>%
  ggplot(aes(x = gen, group = sig.e)) +
  geom_line(aes(y = scaled.mean, colour = sig.e))

# okay - same pattern as before... maybe even larger

# so there is no one variable that controls this...

#########
##### What about if we simplify... take out mutations

pars.list = expand.grid(
  sig.z = sqrt(((1:4)/2) * 10),
  trial = 1:trials.per
) %>%
  mutate(
    sig.a = sqrt((.25) * sig.z),
    sig.e = sqrt((.75) * sig.z),
    n.pop0 = 1000, s.max = 0.9, r = (1.2 / (0.9)) - 1, wfitn = sqrt(10),
    alpha = 0.0000,
    timesteps = 20,
    mu = 0, sig.m = 0
  ) %>%
  mutate(trial.glob = 1:nrow(.)) %>%
  split(f = 1:nrow(.))

set.seed(333)

out1 = mclapply(pars.list,
                
                function(pars) { 
                  
                  sim(params = pars, theta.t = 0, init.rows = 1e5) %>%
                    mutate(e_i = z_i - b_i) %>%
                    group_by(gen) %>%
                    summarise(
                      n = n(),
                      bbar = mean(b_i),
                      zbar = mean(z_i),
                      ebar = mean(e_i),
                      bvar = var(b_i) * (1 - 1/n),
                      zvar = var(z_i) * (1 - 1/n),
                      evar = var(e_i) * (1 - 1/n),
                      ecov = cov(b_i, e_i)
                    ) %>%
                    mutate(trial.no = pars$trial.glob[1])
                  
                },
                mc.cores = 6
) %>%
  do.call(rbind, .)

# okay - here although there is possibly modest erosion still happening,
# trend persists
# what to do...

trials.per = 25

pars.list = expand.grid(
  sig.a = sqrt(1:4),
  trial = 1:trials.per
) %>%
  mutate(
    sig.e = sqrt(10 - sig.a),
    n.pop0 = 1000, s.max = 0.9, r = (1.2 / (0.9)) - 1, wfitn = sqrt(10),
    alpha = 0.0000,
    timesteps = 20,
    mu = 0, sig.m = 0
  ) %>%
  mutate(trial.glob = 1:nrow(.)) %>%
  split(f = 1:nrow(.))

set.seed(333)

out1 = mclapply(pars.list,
                
                function(pars) { 
                  
                  sim(params = pars, theta.t = 0, init.rows = 1e5) %>%
                    mutate(e_i = z_i - b_i) %>%
                    group_by(gen) %>%
                    summarise(
                      n = n(),
                      bbar = mean(b_i),
                      zbar = mean(z_i),
                      ebar = mean(e_i),
                      bvar = var(b_i) * (1 - 1/n),
                      zvar = var(z_i) * (1 - 1/n),
                      evar = var(e_i) * (1 - 1/n),
                      ecov = cov(b_i, e_i)
                    ) %>%
                    mutate(trial.no = pars$trial.glob[1])
                  
                },
                mc.cores = 6
) %>%
  do.call(rbind, .)

# ugh... fucking christ
# not sure how the fuck I am supposed to do this if I can't get motherfucking equilibrium values
# fuck

##### Well... let's maybe for now just make sig.e = 0 and see what happens

trials.per = 30

pars.list = expand.grid(
  sig.a = sqrt(1:4),
  trial = 1:trials.per
) %>%
  mutate(
    sig.e = 0,
    n.pop0 = 1000, s.max = 0.9, r = (1.2 / (0.9)) - 1, wfitn = sqrt(10),
    alpha = 0.0000,
    timesteps = 20,
    mu = 0, sig.m = 0
  ) %>%
  mutate(trial.glob = 1:nrow(.)) %>%
  split(f = 1:nrow(.))

set.seed(333)

out1 = mclapply(pars.list,
                
                function(pars) { 
                  
                  sim(params = pars, theta.t = 0, init.rows = 1e5) %>%
                    mutate(e_i = z_i - b_i) %>%
                    group_by(gen) %>%
                    summarise(
                      n = n(),
                      bbar = mean(b_i),
                      zbar = mean(z_i),
                      ebar = mean(e_i),
                      bvar = var(b_i) * (1 - 1/n),
                      zvar = var(z_i) * (1 - 1/n),
                      evar = var(e_i) * (1 - 1/n),
                      ecov = cov(b_i, e_i)
                    ) %>%
                    mutate(trial.no = pars$trial.glob[1])
                  
                },
                mc.cores = 6
) %>%
  do.call(rbind, .)

out2 = out1 %>%
  pivot_longer(n:ecov, names_to = 'vartype', values_to = 'varval') %>%
  merge(y = do.call(rbind, pars.list) %>% 
          select(trial.glob, sig.a) %>%
          mutate(sig.a = round(sig.a, 2)) %>%
          rename(trial.no = trial.glob),
        by = 'trial.no') %>%
  group_by(sig.a, gen, vartype) %>%
  summarise(
    var.mean = mean(varval),
    var.se   = sd(varval),
    n.obs    = n()
  )

out2 %>%
  filter(vartype %in% 'bvar') %>%
  mutate(sig.a = factor(sig.a)) %>%
  ggplot(aes(x = gen, group = sig.a)) +
  geom_line(aes(y = var.mean, colour = sig.a)) #+
  # scale_y_log10()

out2 %>%
  filter(vartype %in% 'bvar') %>%
  mutate(sig.a = factor(sig.a)) %>%
  group_by(sig.a) %>%
  mutate(scaled.mean = var.mean / var.mean[1]) %>%
  ggplot(aes(x = gen, group = sig.a)) +
  geom_line(aes(y = scaled.mean, colour = sig.a)) +
  scale_y_log10()

# okay... so w/ no mutations, perfect heritability... differences in loss of diversity
# (ugh fuck even this is not perfectly geometric... goddamnit)
# maybe start with simplest case of mutation where there is no env. variance?

# mutation rate: if the first cohort loses proportion sig.a^2 / (sig.a^2 + w^2).

trials.per = 25

pars.list = expand.grid(
  sig.a = sqrt(1:4),
  trial = 1:trials.per
) %>%
  mutate(
    l.max = 1.2,
    sig.e = 0,
    n.pop0 = 1000, s.max = 0.9, r = (l.max / (0.9)) - 1, wfitn = sqrt(10),
    alpha = 0.0000,
    timesteps = 20
  ) %>%
  mutate(
    mu = 1, 
    p  = s.max / l.max,
    sig.m = sqrt((sig.a^2 / (1-p)^2) * (1 - p^2 * (sig.a^2 / (sig.a^2 + wfitn^2))))
  ) %>%
  mutate(trial.glob = 1:nrow(.)) %>%
  split(f = 1:nrow(.))

set.seed(401)

out1 = mclapply(pars.list,
                
                function(pars) { 
                  
                  sim(params = pars, theta.t = 0, init.rows = 1e5) %>%
                    mutate(e_i = z_i - b_i) %>%
                    group_by(gen) %>%
                    summarise(
                      n = n(),
                      bbar = mean(b_i),
                      zbar = mean(z_i),
                      bvar = var(b_i) * (1 - 1/n),
                      zvar = var(z_i) * (1 - 1/n),
                    ) %>%
                    mutate(trial.no = pars$trial.glob[1])
                  
                },
                mc.cores = 6
) %>%
  do.call(rbind, .)

out2 = out1 %>%
  pivot_longer(n:zvar, names_to = 'vartype', values_to = 'varval') %>%
  merge(y = do.call(rbind, pars.list) %>% 
          select(trial.glob, sig.a) %>%
          mutate(sig.a = round(sig.a, 2)) %>%
          rename(trial.no = trial.glob),
        by = 'trial.no') %>%
  group_by(sig.a, gen, vartype) %>%
  summarise(
    var.mean = mean(varval),
    var.se   = sd(varval),
    n.obs    = n()
  )

out2 %>%
  filter(vartype %in% 'bvar') %>%
  mutate(sig.a = factor(sig.a)) %>%
  ggplot(aes(x = gen, group = sig.a)) +
  geom_line(aes(y = var.mean, colour = sig.a)) #+
# scale_y_log10()

out2 %>%
  filter(vartype %in% 'bvar') %>%
  mutate(sig.a = factor(sig.a)) %>%
  group_by(sig.a) %>%
  mutate(scaled.mean = var.mean / var.mean[1]) %>%
  ggplot(aes(x = gen, group = sig.a)) +
  geom_line(aes(y = scaled.mean, colour = sig.a)) +
  scale_y_log10()

# well... fuck that's not good

# fuck man what the fuck!!!
# is there covariance or something? jesus fucking christ.

parsi = pars.list[[1]] %>% mutate(timesteps = 3)

set.seed(429)

sim.out = sim(params = parsi, theta = 0, init.rows = 5000)

head(sim.out)
table(sim.out$gen)

sim.out %>%
  mutate(cohort = gen - age + 1) %>%
  group_by(gen, cohort) %>%
  summarise(
    n = n(),
    bvar = var(b_i) * (1 - 1 / n())
  ) %>%
  group_by(gen) %>%
  mutate(
    p = n / sum(n)
  ) %>%
  mutate(
    exp.v = sum(p^2 * bvar)
  )%>%
  ggplot(
    aes(
      x = gen, y = bvar, group = cohort
    )
  ) +
  geom_line()
  
# yeah... there *has* to be covariance among generations
# fuuuuck...
# (joint normal again of b_0, b_1?)
# uhhh

### SN 10/10/22 conclusion:
# it seems like the amount of variance lost is quite complicated, not standard
# across params even after some amount of normalization
# what to do about this??
# well... seems like a good start is just to look at case where sig.e = 0,
# as the loss of this variation will complicate analysis of loss of genetic var
# there's covariance among generations so can't just do a simple summation
# guess I'll just have to try brute-forcing... ugh

