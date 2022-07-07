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

##### `evolving_geno_pheno_vars.R` (February 16 2022)

Script for visualizing and analyzing changes in phenotypic and genotypic variances over time in a couple of different contexts.

Originally run with an old, possibly wrong way of assigning breeding value variances.

##### `assess_pheno_geno_gap.R` (March 16 2022)

At some point in the last several weeks I noticed that with increasing survival, in adults mean genotypes diverged from mean phenotypes.

Here I explore that a little bit. I ran simulations across longevity and phenotypic noise combinations, tracking the gap between phenotypes and genotypes (breeding values).

There is no gap with no phenotypic noise (perfect heritability), which makes sense.
Increasing both the phenotypic noise and/or the longevity increases this gap.
It does not appear to alter population growth rates though (although it is hard to disentangle this from the still-unresolved issue with phenotypic variance).

A hypothesis for what is happening here is that selection is acting on *phenotypes* and not genotypes.
And, with age, individuals may be selected for with middling/maladapted genotypes, but very good phenotypes.
I ran a trial to visualize this, looking at ``residual'' phenotypic variance to see if it is biased with age.
And yes, it looks like it is! (although the type of environmental change may influence this).
Older individuals remaining in the population tend to have a mean phenotype closer to the phenotypic optimum than their breeding value is.
(Aaron W. noted a simple binomial test for the sign of this ``residual'' variation could also demonstrate the point).

This seems related to but not necessarily the same as the idea of inertia that motivated this project.
The maladapted individuals do contribute lagging breeding values.
But is there an *extra* lag behind the phenotype? Or is the phenotypic variance masking some of this effect?
Brett noted that this may be similar to an extinction debt... e.g., a selection/phenotypic debt.

I did run multiple trials to look at the shape of a ``luck'' curve as a function of age,
i.e., the mean ``residual'' phenotypic variance as a function of age.
It does exist; there is still a strong age effect in the faster LH groups though.
The shape of this luck curve may not actually vary by longevity class
(it depends on the way you calculate means... needs more thought).

I would like to test to see how this luck affects the rate of adaptation... but I'm not sure how to do this.
Return to the problem later.

##### `assess_var_under_selection.R` (May 9 2022)

A quick and dirty script for looking at what happens under iterative selection on normally distributed phenotypes.
(Assuming that the phenotypic distribution is centered approximately at zero over time)

The idea here is that (at pre-environmental change steady state), a cohort of individuals will be subject to iterative rounds of selection that (should) reduce the phenotypic variance. 
Here I am looking around for an easy way to characterize this variance after an arbitrary number of rounds of selection. 
If there is a way to do this, then perhaps there is a way to get the variance among adults analytically.

Thanks to some luck in selecting fixed values (multiples of 5!) it turns out there is a way to get this.
There's a recursion relationship that I was able to describe using some pen and paper (described on overleaf).

##### `no_selection_variance_gaps.R` (May 17 2022)

Talking with Dan D. (5/13) made me wonder to what degree the changing phenotypic variance is due to selection (instead of some error in code).
Does the assumption of breeding value variance and additive genetic variance being equal rely on no selection?
Here, I run some simulations to see if parent/adult phenotypic/breeding value variance (after selection) and offspring variance are the same.

It looks like they are! At least no consistent gap in breeding values.
What does this mean? Well...
Mating and segregation of gametes is occurring indpendently of selection (selection acts on survival).
So, we could do segregational variance just using the parental breeding value variance in that timestep.
(Note - some variance should be restored using mutations)

A potential problem is what to do about small popuplation effects... I will think about this.

This also doesn't really address the issue of what to do for breeding value initialization, but Dan did provide some ideas here.

##### `test_equilibrium_variance.R` (May 24 2022)

I (thought that I) figured out an analytical expression for the equilibrium variance of an age-structured population. Here, I tested if it was correct (and looked at a few other things).
I did this by looking at simulations with initially large populations (2000 individuals) with a mean growth rate of approximately 1 for ~25 generations.

It turns out my expressions were not correct (even after correcting for some silly things).

It looks like the age distribution is geometric, with rate approximately `s.max` (which is interesting, because with selection the mean survival is not `s.max` but instead something below it - I suppose in a large population this is not an issue.) There was some correction needed for the fact that ages begin at one instead of zero, though - this needs to be included in subsequent steps.

Likewise, it does look like the expression for variance after successive rounds of selection holds, although for the very old cohorts there is a ton of noise as would be expected from the small sample sizes.

So what is different? I think it's the covariances among all of the groups. This I would need to play around with by hand.

An intersting pattern I saw surprisingly often was a sort of periodicity to the phenotypes going from cohort to cohort. That is, looking in one generation at the variance from cohort to cohort, there would be an up-down-up-down pattern around the mean. What is this?

Oh yeah, and what about non-genetic phenotypic variance? Most of what I was looking at was a small number of simulations with perfect heritability. But of course the phenotypic variance is what determines the growth rate, and the rate of adaptation... ugh. Perhaps it's worth only looking at perfect heritability in these simulations? But then there is the geno-pheno gap... 

There are some cool visualizations in this script by the way! Looking at variance across cohorts, generations, ages...

##### `test_phenotypic_change.R` (June 2022)

Here I was looking to see if my expressions for phenotypic changes over time were accurate.
They were not quite but they were in the right ballpark.

This also includes looking at phenotypic and genotypic change in adults/newborns.
It turns out there is a pheno-geno gap in adults but not in offspring.
There are some cool visualizations of this in ehre.

Also, here I ran a couple of scripts looking at rate of (pheno/geno) change with high/lower survival.
The gap widens with greater survival (because selection acts on adult phenotypes).
Very cool visualizations in here (I think at least) to demonstrate this.
Seems like the phenotypic advantage to longevity is larger than the genotypic downside.

But yes it does seem like there is a disadvantage in *genotypic* change to survival - does it matter though?
Can we quantify it?

##### `visualize_advantages.R` (July 6 2022)

Silly little script looking at ratios of rate of adaptation and fitness with iteroparity compared to survival.

Results are about as expected, with larger effects of survival for strong selection and low heritability.
Maybe the biggest surprise was how the effects seemed so much incredibly larger for only such low h, low w.

