rm(list = ls())

source('model_source/sim_model1_functions.R')

lcs = data.frame(lstar = (91:110)/100) %>%
  group_by(lstar) %>%
  mutate(g20 = newt.method.g1(.1, 1e-8, .9 / lstar, (1/3)))

expand.grid(
  l.ccc = (91:111)/100,
  age   = 0:15
) %>%
  merge(lcs) %>%
  mutate(
    p.k = (1/3) / (1 + (1/3)) * (.9 / l.ccc)^age / sqrt(1 + age*g20)
  ) %>%
  ggplot(aes(x = age, y = p.k)) +
  geom_line(
    data = data.frame(age = 0:15) %>% mutate(p.k = (1/4) * (3/4)^age),
    linetype = 2, colour = 'gold', linewidth = 1.5
  ) +
  geom_line(aes(colour = l.ccc, group = l.ccc)) +
  scale_color_gradient2(low = 'red', high = 'blue', mid = 'white', midpoint = 1, 'lambda*') +
  theme(panel.background = element_rect(fill = 'black'), panel.grid = element_blank())

ggsave('~/Documents/Research/boulder/other_image_materials/quasi-geometric.pdf')

  #scale_colour_viridis_c()

taco = expand.grid(p0 = (1 + 2*(1:4))/10, l.ccc = (17:23)/20) %>%
  mutate(
    l.max = 1.2,
    s.max = l.max * (1 - p0),
    r     = p0 / (1 - p0)
  ) %>%
  group_by(p0, l.ccc) %>%
  mutate(g20 = newt.method.g1(.1, 1e-8, s.max / l.ccc, r))

taco = taco %>%
  group_by(p0, l.ccc) %>%
  mutate(g2  = gamma.calc(g20, s.max / l.ccc, r))

taco %>%
  mutate(p0 = factor(p0)) %>%
  ggplot(aes(x = l.max / l.ccc, y = g20, colour = p0, group = p0)) +
  geom_line() +
  geom_point(size = 3) +
  labs(x = 'lambda.max / lambda*', y = 'gamma^2_0 (newest cohort variance)') +
  lims(y = c(0, 5)) +
  scale_colour_viridis_d()

ggsave('~/Documents/Research/boulder/other_image_materials/lrat_g0_gt1.pdf')

taco %>%
  mutate(p0 = factor(p0)) %>%
  ggplot(aes(x = l.ccc / l.max, y = g20, colour = p0, group = p0)) +
  geom_line() +
  geom_point(size = 3) +
  labs(x = 'lambda* / lambda.max', y = 'gamma^2_0 (newest cohort variance)') +
  lims(y = c(0, 5)) +
  scale_colour_viridis_d()

ggsave('~/Documents/Research/boulder/other_image_materials/lrat_g0_lt1.pdf')


taco %>%
  mutate(p0 = factor(p0)) %>%
  ggplot(aes(x = l.max / l.ccc, y = g2, colour = p0, group = p0)) +
  geom_line() +
  geom_point(size = 3) +
  labs(x = 'lambda.max / lambda*', y = 'gamma^2 (population-level variance)') +
  lims(y = c(0, 2.5)) +
  scale_colour_viridis_d()

ggsave('~/Documents/Research/boulder/other_image_materials/lrat_g_gt1.pdf')


taco %>%
  mutate(p0 = factor(p0)) %>%
  ggplot(aes(x = l.ccc / l.max, y = g2, colour = p0, group = p0)) +
  geom_line() +
  geom_point(size = 3) +
  labs(x = 'lambda* / lambda.max', y = 'gamma^2 (population-level variance)') +
  lims(y = c(0, 2.5)) +
  scale_colour_viridis_d()

ggsave('~/Documents/Research/boulder/other_image_materials/lrat_g_lt1.pdf')

taco %>%
  mutate(p0 = factor(p0)) %>%
  ggplot(aes(x = l.ccc, y = 1 / (1+g2), colour = p0, group = p0)) +
  geom_line() +
  geom_point(size = 3) +
  scale_colour_viridis_d()

taco %>%
  mutate(p0 = factor(p0)) %>%
  ggplot(aes(x = l.ccc, y = g20, colour = p0, group = p0)) +
  geom_line() +
  geom_point(size = 3) +
  scale_colour_viridis_d()
# cooler! interestingly different from the population-level variance!

taco %>%
  ggplot(aes(x = p0, y = g2, colour = l.ccc, group = l.ccc)) +
  geom_line() +
  scale_colour_viridis_c()
# boring!

taco %>%
  ggplot(aes(x = p0, y = 1/(1+g2), colour = l.ccc, group = l.ccc)) +
  geom_line() +
  geom_point(size = 3) +
  scale_colour_viridis_c()

taco %>%
  mutate(p0 = factor(p0)) %>%
  ggplot(aes(x = l.ccc / l.max, y = g2, colour = p0, group = p0)) +
  geom_line(
    data = data.frame(x = (70:100)/100) %>% mutate(y = x^(-2) - 1),
    aes(x = x, y = y),
    linetype = 2, inherit.aes = FALSE
  ) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  labs(
    x = expression(paste(lambda["*"], '/', hat(lambda))), 
    y = expression(paste("Phenotypic variance, ", gamma^2))
  ) +
  lims(y = c(0, 2)) +
  scale_color_brewer(palette = 'Spectral', expression(p[0]))

ggsave("sandbox/sandbox_figs/eg_gamma.png", width = 8, height = 5)
