# Looking at relationships between phenotypic distributions and distributions of
# frailties (frailties defined as proportional to the squared distance of the
# phenotype from the optimum.)
# SN - 10 May 2023

library(ggplot2)
library(dplyr)
library(tidyr)

data.frame(
  z = (-300:300)/100
) %>%
  mutate(w = exp(-z^2 / 2)) %>%
  mutate(s = 1 - w) %>%
  mutate(f = -log(s)) %>%
  ggplot(aes(x = z, y = 1 / f)) + geom_line()

data.frame(
  z = rnorm(100)
) %>%
  mutate(w = exp(-z^2 / 2)) %>%
  mutate(s = 1 - w) %>%
  mutate(f = -log(s)) %>%
  ggplot(aes(x = f)) + geom_density()

# gamma emerging here... but frailty should be 1/f not f?


data.frame(
  z = rnorm(10000)
) %>%
  mutate(f = z^2) %>%
  ggplot(aes(x = f, after_stat(density))) + geom_histogram(binwidth = 0.1) +
  geom_line(data = data.frame(z = (1:60)/10) %>% mutate(pz = dgamma(z, shape = 1/2, scale = 2)),
            aes(x = z, y = pz), colour = 'red')
# looks good...

data.frame(
  z = rnorm(10000)
) %>%
  mutate(f = z^2 / 2) %>%
  ggplot(aes(x = f, after_stat(density))) + geom_histogram(binwidth = 0.1) +
  geom_line(data = data.frame(z = (1:60)/10) %>% mutate(pz = dgamma(z, shape = 1/2, scale = 1)),
            aes(x = z, y = pz), colour = 'red')
# neat!

data.frame(
  z = rnorm(10000, sd = sqrt(1/3))
) %>%
  mutate(f = (1/3) * z^2 / 2) %>%
  ggplot(aes(x = f, after_stat(density))) + geom_histogram(binwidth = 0.1) +
  geom_line(data = data.frame(z = (1:30)/10) %>% mutate(pz = dgamma(z, shape = 1/2, scale = (1/3) * 2 / 2)),
            aes(x = z, y = pz), colour = 'red')

data.frame(
  z = rnorm(10000, sd = sqrt(1/3))
) %>%
  mutate(f = (1/3) * (z / sqrt(1/3))^2 / 2) %>%
  ggplot(aes(x = f, after_stat(density))) + geom_histogram(binwidth = 0.1) +
  # geom_line(data = data.frame(z = (1:30)/10) %>% mutate(pz = dchisq(z, df = 1)),
  #           aes(x = z, y = pz), colour = 'red') +
  geom_line(data = data.frame(z = (1:30)/10) %>% mutate(pz = dgamma(z, shape = 1/2, scale = (1/3) * 2 / 2)),
            aes(x = z, y = pz), colour = 'blue')
# oh I think I see what's wrong

data.frame(
  z = rnorm(10000, sd = sqrt(1/3))
) %>%
  mutate(f = z^2 ) %>%
  ggplot(aes(x = f, after_stat(density))) + geom_histogram(binwidth = 0.1) +
  geom_line(data = data.frame(z = (1:30)/10) %>% mutate(pz = dgamma(z, shape = 1/2, scale = (1/3) * 2)),
            aes(x = z, y = pz), colour = 'red')
# lmao hell yearh

data.frame(
  z = rnorm(1000)
) %>%
  mutate(w = exp(-z^2 / 2)) %>%
  mutate(f = -log(w)) %>%
  ggplot(aes(x = z, y = f, colour = w)) +
  geom_point() +
  geom_line(data = data.frame(x = (-300:300)/100, y = ((-300:300)/100)^2 / 2), 
            aes(x = x, y = y), colour = 'black') +
  scale_colour_viridis_c()
