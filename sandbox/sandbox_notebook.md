# Sandbox notes

The sandbox directory is for scripts that play around with simulation code. It's informal analysis. Scripts are messy and not very well commented.

This markdown is for describing scripts (briefly) and noting lessons learned from them.

### Model 1

##### `assess_onestage_growth_rates.R` (January 2022)

Script that looks both at expectation and variance in growth rates of populations with (1) no phenotypic variance and (2) no environmental mismatch (perfectly adapted in expectation) or stochasticity.

Expected growth rate is quite simple to figure out - I was also able to figure out an expression for variance in growth rates in these circumstances. However, the result on variance in growth rates probably falls apart when adding in phenotypic variance (sad). 

##### `geno_var_growth_rates.R` (February 14 2022)

Here, I was trying to figure out what growth rates (and their variance) look like with phenotypic variance. I was also looking at NDD to prevent populations from exploding in size. I was only looking at one parameter set (life history). This was in an effort to find a stable initial growth rate (including NDD parameter that allows for a reasonable carrying capacity to initialize the population at).

Key insight from this is that (1) phenotypic variance needs to be part of the estimate of max growth rate, i.e., max growth rate depends on phenotypic variance (potentially a problem if variance is not stable) and (2) NDD parameter can just be set simply from the initial size and max growth rate. Code for doing this is at the end of the script.

##### `geno_pheno_variance_smax.R` (February 15 2022)

Here I was looking at changing genotypic and phenotypic variation across three LH groups. All of these simulations were run with phenotypic variance, density dependence, and perfectly adapted (in expectation) populations.

Annoyingly, the rate of genotypic change increases for all parameter groups (with different equilibrium) for each group, but phenotypic change shows different patterns across groups (slow loses genotypic variance, fast gains it). I have no idea why this difference is here.

I have some thoughts about this at the end of the script but I don't have a definitive answer. Weighted sums of variance of (surviving) parents and offspring may help? 

This confirms the genotypic variance is important... population size grows/shrinks with phenotypic variance and stabilizes after ~20 generations. Fast populations have a lower equilibrium size (due to larger pheno variance) and faster populaitons have a greater one.

(update, February 16 2022: I think I may have been modelling breeding value-assignment incorrectly
I will re-run some of these scripts soon)

### `evolving_geno_pheno_vars.R` (February 16 2022)

Script for visualizing and analyzing changes in phenotypic and genotypic variances over time in a couple of different contexts.

Originally run with an old, possibly wrong way of assigning breeding value variances.
