# Here - looking to see what happens with iterative gaussian selection on
# normally distributed phenotypes
# SN - 9 May 2022
# (quick script - no comments)

n_0 = 1000
w_2 = 5
s20 = 1

n.tries = 1000

outsput = data.frame(try = 1:n.tries, output.var = 0)

for (trial in 1:n.tries) {
  outsput$output.var[trial] = data.frame(i = 1:n_0, z_i = rnorm(n_0, 0, s20)) %>%
                                mutate(w_i = exp(-z_i^2/(2*w_2))) %>%
                                filter(as.logical(rbinom(n_0, prob = w_i, size = 1))) %>%
                                summarise(var(z_i)) %>%
                                unlist()
}

hist(outsput$output.var)

mean(outsput$output.var) + (c(-1, 0, 1) * 2 * sqrt(var(outsput$output.var) / n.tries) )
# this does include 5/6!

# Well let's look at iterations

time.len = 50

outsput2 = data.frame(iter = 0:time.len, var0 = 0, var1 = 0)
outsput2$var0[1] = s20
outsput2$var1[1] = s20

for (timestep in 1:time.len) {
  outsput2$var0[timestep + 1] = outsput2$var0[timestep] * (w_2 / (w_2 + s20))
  outsput2$var1[timestep + 1] = outsput2$var1[timestep] * (w_2 / (w_2 + outsput2$var1[timestep]))
}             

outsput2 %>%
  ggplot(aes(x = iter)) +
  geom_line(aes(y = var0), colour = 'green') +
  geom_line(aes(y = var1), colour = 'blue') #+
  #scale_y_log10()

# oh god but some of these numbers are so fucking round...
# there must be a way to do this analytically!!

outsput2[1:10,]
  
# but will it ever converge?? hmm...

