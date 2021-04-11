### Auxiliary functions for handling sim results


unroll.sims = function(sim.list) {
  # Wrapper for converting lists into a single data frame
  # Note: requires individual id numbers begin with 1 in all trials
  sim.list %>%
    do.call(what = rbind) %>%
    mutate(trial = cumsum(i %in% 1))
}

unroll.sums = function(summ.list) {
  # Wrapper for converting lists of summaries into a single data frame
  # Note: will fail in cases where populations go extinct in one gen
  summ.list %>%
    do.call(what = rbind) %>%
    mutate(trial = cumsum(c(1, as.numeric(diff(gen) < 0))))

}

summarise.demo = function(sim.out) {
  # Wrapper to get pop size, genotype, fitness from one sim
  sim.out %>%
    group_by(gen) %>%
    summarise(n = n(),
              gbar = mean(g_i),
              zbar = mean(z_i),
              wbar = mean(w_i),
              theta = theta[1])
}

summarise.geno = function(sim.out, pars) {
  sim.out %>%
    select(-c(i, g_i, z_i, w_i, r_i, fem, theta)) %>%
    gather(key = loc.copy, value = val, -gen) %>%
    mutate(locus = gsub('^[ab]', '', loc.copy)) %>%
    group_by(gen, locus) %>%
    summarise(p = mean(val > 0)) %>%
    group_by(gen) %>%
    summarise(p.fix.pos = mean(p == 1),
              p.fix.neg = mean(p == 0),
              v = sum(2 * p * (1 - p)) / pars$n.loci[1])
}