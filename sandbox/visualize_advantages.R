### Visualizing the "advantage" of overlapping generations
# Here: comapring relative (rate sof adaptation, fitness) of semelparous and iteroparous
# populations.
# Using ratios (iteroparous/semelparous) because it seems more reasonable
# (multiplicative process over time).

# First, try out visualizing rate of adaptation as a function of proportion of
# population surviving
w = 2
e = 2
a = 2

p = (1:9)/10

(1 - (a^2)/(w^2 + e^2 + a^2) - (p)*(e^2)/(w^2 + e^2 + a^2)) %>% plot(type = 'l')
# oh duh it's just linear

# Here: a useful wrapper function for determining rate of adaptation
k1 = function(w, e, a, p = 0) (1 - (a^2)/(w^2 + e^2 + a^2) - p * (e^2) / (w^2 + e^2 + a^2))
# w - inv. selection pressure (low number is higher selection pressure)
# e - s.d. of environmentally-determined phenotypic variance
# a - additive genetic variance
# p - proportion of population that is surviving adults

# Here: visualize selection strength over a wide range of parameters
# (note: closer to 1 means SLOWEr adaptation)

expand.grid(w = (1:25)/5, a = 1, h = (1:99)/100, p = (0:99)/100) %>%
  mutate(e = sqrt(a^2 * (1-h)/h)) %>%
  mutate(k = k1(w,e,a,p)) %>%
  ggplot(aes(x = h, y = p, fill = k)) +
  geom_raster() +
  scale_fill_viridis_c() +
  facet_wrap(~ w, nrow = 5) +
  theme(legend.position = 'bottom')

# in all cases: increasing p increases the rate of adaptation
# but effects seem greatest for very low-h traits
# and in cases where h is high, selection is weak, barely any difference

# Ratios of the rate of adaptation
# (smaller number here means relatively faster adaptation)

expand.grid(w = (1:16)/4, a = 1, h = (1:99)/100, p = (1:99)/100) %>%
  mutate(e = sqrt(a^2 * (1-h)/h)) %>%
  mutate(k.ratio = k1(w,e,a,p) / k1(w,e,a,0)) %>%
  ggplot(aes(x = h, y = p, fill = k.ratio)) +
  geom_raster() +
  scale_fill_viridis_c() +
  facet_wrap(~ w, nrow = 4) +
  theme(legend.position = 'bottom')

# relative effects of perenniality: 
# very large for high selection (makes sense)
# larger (in general) for less heritable traits

# Here: relative fitness
# (plugging the respective rates into a Gaussian function)
# higher number in this case means relatively HIGHER fitness

expand.grid(w = (1:12)/4, a = 1, h = (1:99)/100, p = (1:99)/100) %>%
  mutate(e = sqrt(a^2 * (1-h)/h)) %>%
  mutate(s.ratio = exp(-k1(w,e,a,p)^2 / (2*w^2)) / exp(-k1(w,e,a,0)^2 / (2*w^2))) %>%
  ggplot(aes(x = h, y = p, fill = s.ratio)) +
  geom_raster() +
  scale_fill_viridis_c() +
  facet_wrap(~ w, nrow = 3, ncol = 4) +
  theme(legend.position = 'bottom')
# ah... color ramp here is not great, getting blown out by high advantage for
# *very extreme* selection pressure, high survival, low heritability
# (beyond that - can't visualize much)

# Look at a smaller range of selection pressures
# (starting at 1 - weaker)
expand.grid(w = (2:9)/2, a = 1, h = (1:99)/100, p = (1:99)/100) %>%
  mutate(e = sqrt(a^2 * (1-h)/h)) %>%
  mutate(s.ratio = exp(-k1(w,e,a,p)^2 / (2*w^2)) / exp(-k1(w,e,a,0)^2 / (2*w^2))) %>%
  ggplot(aes(x = h, y = p, fill = s.ratio)) +
  geom_raster() +
  scale_fill_viridis_c() +
  facet_wrap(~ w, ncol = 4) +
  theme(legend.position = 'bottom')
# Yes same issue... seems like the advantage really takes most effect for VERY strong selection
# (and at that only ideal situations with very large proportion of adults)

expand.grid(w = 1:4, a = 1, h = (1:99)/100, p = (1:99)/100) %>%
  mutate(e = sqrt(a^2 * (1-h)/h)) %>%
  mutate(s.ratio = exp(-k1(w,e,a,p)^2 / (2*w^2)) / exp(-k1(w,e,a,0)^2 / (2*w^2))) %>%
  ggplot(aes(x = h, y = p, fill = s.ratio)) +
  geom_raster() +
  scale_fill_viridis_c() +
  facet_wrap(~ w, nrow = 2, ncol = 2) +
  theme(legend.position = 'bottom')


### Lesons:
# Intuitino is right. Largest benefits for very low heritability, very (very) stron gselection
