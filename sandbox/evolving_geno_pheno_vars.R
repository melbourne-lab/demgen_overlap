library(ggplot2)

rm(list = ls())

source('model_source/sim_model1_functions.R')

lambda.max = 1.1

pars = data.frame(s.max = rep(c(0.5, 0.7, 0.9), each = 500)) %>%
  mutate(n.pop0 = 500, 
         wfitn = 2,
         sig.a = sqrt(0.5), 
         sig.e = sqrt(0.5),
         timesteps = 10) %>%
  # Determine density dependence (depends on N0)
  mutate(alpha = log(lambda.max) / n.pop0) %>%
  # Determine fecundities
  mutate(r = (lambda.max / s.max) * sqrt((wfitn^2 + sig.a^2 + sig.e^2) / wfitn^2) - 1) %>%
  mutate(trial = 1:nrow(.))

# Sims

li = vector('list', nrow(pars))

set.seed(759200)

for (k in 1:nrow(pars)) li[[k]] = sim(pars[k,], theta.t = 2, init.rows = 1e5) %>% mutate(trial = k)

mm = merge(x = do.call(rbind, li), y = pars %>% select(trial, s.max))

mm %>%
  group_by(trial, gen) %>%
  summarise(s.max = s.max[1],
            n = n()) %>%
  ggplot(aes(x = gen, y = n)) +
  geom_line(aes(group = trial)) +
  scale_y_log10() +
  facet_wrap(~ s.max)

# woke
# keep in mind the pheno variance though...

mmsum = mm %>%
  group_by(trial, gen) %>%
  summarise(s.max = s.max[1],
            n = n(),
            b = mean(b_i),
            z = mean(z_i)) %>%
  group_by(s.max, gen) %>%
  summarise(nbar = mean(n),
            nvar = var(n),
            bbar = mean(b),
            bvar = var(b),
            zbar = mean(z),
            zvar = var(z),
            nn = n())

mmsum %>%
  ggplot(aes(x = gen)) +
  geom_line(aes(y = bbar, group = s.max, color = factor(s.max))) +
  geom_ribbon(aes(ymin = bbar-2*sqrt(bvar/nn), ymax = bbar+2*sqrt(bvar/nn),
                  fill = factor(s.max)),
                  alpha = 0.2)
# hmm... genotypic change is actually the same for the first several generations??
              
mmsum %>%
  filter(gen < 5) %>%
  ggplot(aes(x = gen)) +
  geom_line(aes(y = bbar, group = s.max, color = factor(s.max))) +
  geom_ribbon(aes(ymin = bbar-2*sqrt(bvar/nn), ymax = bbar+2*sqrt(bvar/nn),
                  fill = factor(s.max)),
              alpha = 0.2)
# hmm... 
# here, looks like it's the same for generation 1 an diverges after
# (maybe just need more samples to see the difference... - or all have same pheno variance?)

## Within a group, does the pheno variance influence the 

fast.summ = mm %>%
  mutate(b_i = theta_t - b_i,
         z_i = theta_t - z_i) %>%
  filter(s.max < 0.6) %>%
  group_by(trial, gen) %>%
  summarise(s.max = s.max[1],
            n = n(),
            bbar = mean(b_i),
            bvar = var(b_i),
            zbar = mean(z_i),
            zvar = var(z_i))

fast.summ = fast.summ %>%
  filter(zbar > 0, bbar > 0) %>%
  mutate(d.bbar = c(exp(diff(log(bbar))), NA),
         d.zbar = c(exp(diff(log(zbar))), NA)) %>%
  filter(!is.na(d.bbar))

fast.summ %>%
  #filter(gen < 2) %>%
  ggplot(aes(x = gen)) +
  geom_line(aes(y = bvar, group = trial), size = 0.1)

fast.summ %>%
  filter(gen < 9) %>%
  filter(abs(d.bbar) < 2) %>%
  ggplot(aes(x = bvar, y = 1 - d.bbar)) +
  geom_point(size = 3) +
  facet_wrap(~ gen)
  # scale_color_viridis_c() #+
  # scale_y_log10()

fast.summ %>%
  filter(gen < 8) %>%
  ggplot(aes(x = gen, y = d.bbar)) +
  geom_line(aes(group = trial), size = 0.1)

fast.summ %>%
  group_by(gen) %>%
  summarise(db.bar = mean(d.bbar),
            db.var = var(d.bbar),
            n = n()) %>%
  ggplot(aes(x = gen)) +
  geom_line(aes(y = db.bar)) +
  geom_ribbon(aes(ymin = db.bar-2*sqrt(db.var/n), ymax = db.bar+2*sqrt(db.var/n)),
              alpha = 0.1)
# welp... rate of adaptation is changing over time. not good!
# the increase I think is due to incredibly large outliers

# well, okay, different question
# even with adaptation/genotypic change, does variance equilibrate
# (if it does, then we maybe can at least brute-force a way to get equivalent phenotypic variances?
# could we do this? seems like phenotypic variance depends on other parameters
# fuck...)

fast.summ %>%
  group_by(gen) %>%
  summarise(sig2a     = mean(bvar),
            sig2a_var = var(bvar),
            n = n()) %>%
  ggplot(aes(x = gen)) +
  geom_line(aes(y = sig2a)) +
  geom_ribbon(aes(ymin = sig2a-2*sqrt(sig2a_var/n),
                  ymax = sig2a+2*sqrt(sig2a_var/n)),
              alpha = 0.1)
# that looks *basically* equilibrated...  

# two more questions
# - first, do we see the same patterns in all LH groups
#   (are all LH groups also equilibrating, and does this pattern look similar to
#   that from the theta = 0 cases)
# - second, does the equilibrium vary depending on the new phenotypic optimum??
#   (what else causes it to vary? surely, the initial starting variance...)
  
varsum = mm %>%
  group_by(trial, gen) %>%
  summarise(s.max = s.max[1],
            sigb = var(b_i),
            sigz = var(z_i)) %>%
  group_by(s.max, gen) %>%
  summarise(sigb.m = mean(sigb),
            sigb.v = var(sigb),
            sigz.m = mean(sigz),
            sigz.v = var(sigz),
            n = n())


varsum %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = sigb.m,
      colour = factor(s.max)
    )
  ) +
  geom_ribbon(
    aes(
      ymin = sigb.m - 2 * sqrt(sigb.v/n),
      ymax = sigb.m + 2 * sqrt(sigb.v/n),
      fill = factor(s.max)
    ),
    alpha = 0.2
  )

varsum %>%
  ggplot(aes(x = gen)) +
  geom_line(
    aes(
      y = sigz.m,
      colour = factor(s.max)
    )
  ) +
  geom_ribbon(
    aes(
      ymin = sigz.m - 2 * sqrt(sigz.v/n),
      ymax = sigz.m + 2 * sqrt(sigz.v/n),
      fill = factor(s.max)
    ),
    alpha = 0.2
  )
# okay, well this looks about the same qualitatively, which is good

### Looking at equlibrium 

n.trials = 1800

ell = vector('list', n.trials)

set.seed(486605)

for (k in 1:n.trials) ell[[k]] = sim(pars[1,], theta.t = 1 + (k%%3), init.rows = 1e5) %>% mutate(trial = k)

ll = do.call(rbind, ell)

table(ll$theta_t)

ll %>%
  group_by(trial, gen) %>%
  summarise(n = n(),
            theta_t = theta_t[1]) %>%
  group_by(theta_t, gen) %>% 
  summarise(nbar = mean(n),
            nvar = var(n),
            nn = n()) %>%
  ggplot(aes(x = gen)) +
  geom_line(aes(y = nbar, group = theta_t, colour = theta_t)) +
  geom_ribbon(aes(ymin = nbar - 2 * sqrt(nvar/nn),
                  ymax = nbar + 2 * sqrt(nvar/nn),
                  group = theta_t,
                  fill = theta_t),
              alpha = 0.1) +
  scale_y_log10()

# summaries of variances
varsum.th = ll %>%
  group_by(trial, gen) %>%
  summarise(vb = var(b_i),
            vz = var(z_i),
            theta_t = theta_t[1]) %>%
  group_by(theta_t, gen) %>% 
  summarise(vb.bar = mean(vb),
            vb.var = var(vb),
            vz.bar = mean(vz),
            vz.var = var(vz),
            nn = n())

varsum.th %>%
  mutate(theta_t = factor(theta_t)) %>%
  ggplot(aes(x = gen)) +
  geom_line(aes(y = vb.bar, group = theta_t, colour = theta_t)) +
  geom_ribbon(aes(ymin = vb.bar - 2 * sqrt(vb.var/nn),
                  ymax = vb.bar + 2 * sqrt(vb.var/nn),
                  group = theta_t,
                  fill = theta_t),
              alpha = 0.1) +
  scale_y_log10()

# shit, fuck, etc.
# great! fucking christ
# anyway the equilibrium levels are I guess dependent on the selection strength

varsum.th %>%
  mutate(theta_t = factor(theta_t)) %>%
  ggplot(aes(x = gen)) +
  geom_line(aes(y = vz.bar, group = theta_t, colour = theta_t)) +
  geom_ribbon(aes(ymin = vz.bar - 2 * sqrt(vz.var/nn),
                  ymax = vz.bar + 2 * sqrt(vz.var/nn),
                  group = theta_t,
                  fill = theta_t),
              alpha = 0.1) +
  scale_y_log10()

# what the absolute fuck??

# Lessons learned
# - yep, it does look so far like phenotypic variance 