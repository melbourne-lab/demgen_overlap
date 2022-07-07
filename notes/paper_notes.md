# Overlapping generations project notes

### Schmid, M., et al. 2020. A tradeoff between robustness to environmental fluctuations and speed of evolution. Pre-print. 

- 

Initial thoughts:

Ugh, okay, so yes very similar to what I had originally intended for this project. Done more thoroughly (although the MS was harder for me to follow at some point - lots of stuff in here).

- The life cycle: one of the key findings here is that there is only one stage that is subject to selection (out of three). The "precociousness" results here seem pretty reliant on this model feature; if selection happened similarly on this early stage then the results may disappear? I guess a lot of organisms do have a pre-reproductive stage. I feel like there's an alternative way to set up the life cycle more like the Tulja 1988 paper... Or perhaps it could be thought of as the time to reproductive maturity?
- The analytical framework here... seems very nice. Is it possible for our model to be done like this? The Barfield et al. 2011 equation seems like it could be easily repurposed for a Leslie model-type life cycle akin to ours
	- One of the nice things about it is that it allows for numeric adjustment that accounts for "max" pop growth rate - I was concerned before about stochasticity but actually I guess the "max" growth rate is one way to account for this (is it effective though??)
- I'm not really sure how much rescue plays a part in this study. Seems like there's a result here about population growth that is independent of environmental change. But a large part of it does seem to be tracking a changing environment - as expected, the "shorter" life history strategies seem to be better at doing this.
	- But, what about a sudden environmental shift? Is there any reason to think this would be different? I guess this depends on underlying differences between a B&L-type model and a G&H model
	- Autocorrelation I guess... usually this seems like a cop out for a next step to explore, but in this case, life history theory (e.g. the Tulja paper) has things to contribute here wrt population growth rates. 
- Any reason to suspect that things would be different with selection on fecundity? That's what it is in our model at least.
- Dan's recommendation about spatial mixing (or temporal mixing too?) and some other feedback I guess
- What is this genetic storage effect?

### Ellner, S., and Hairston, N.G. 1994. Role of overlapping generations in maintaining genetic variation in a fluctuating environment. The American Naturalist.

- Temporally fluctuating selection (relative fitnesses of phenotypes varies over time) proposed and rejected (!?) as a process for maintenance of genetic variation
	- Frank and Slatkin 1990 suggest that spatial variation in selection could maintain polymorphisms although temporal variation likely would not
	- quant gen: temporal variation in phenotypic fitness hardly if at all increases genetic variance
- These models though assume discrete and non-overlapping generations
	- some other studies (e.g. Chesson) have shown that multiple opportunities for reproduction are important for species coexistence (storage effect)
	- asexual species in these models should translate directly to genetic terms and also possibly for sexual diploids
	- see Chesson 1984, Seger and Brockmann 1987, Haldane and Jayakar 1963
	- but does this transfer to continuous phenotypes?
		- n.b. theory also predicts that in fluctuating environments, species should hedge their bets and choose a single optimum phenotype
- Under what conditions will genetically variable populations be evolutionarily stable?
##### Model and analysis
- X_i(t) is abundance of individuals of type $i$ (genotype or phenotype) in year $t$
	- 1 > H > 0 individuals reproduce, per capita fecundity Y_i(t) and survival s_r
	- 1 - H individuals survive with rate s_d
	- gamma = Hs_r + (1-H)s_d
		- this is generational overlap? seems to be just a weighted average of survival...? actually I guess this makes sense?
	- X_i(t+1) = X_i(t) * [HY_i(t) + gamma] 
		- i.e. number of offspring produced plus survival
- Selection acts through Y_i(t) (n.b. assumes all individuals have same H, s_d, s_r)
	- say that relative fitness is a function of individual's distance fom optimum under gaussian selection, i.e.,
		- R(D_i - M_t) = exp(-(D_i-M_t)^2/(2w^2))
- Saturating yield defines absolute fitness:
	- Yi(t) = K * (R(D_i-M_t)) / sum_j [HX_j(t) R(D_j-M_t)]
	- i.e., there are always exactly K offspring produced and the number of offspring produced per individual with phenotype i is just proportional to its relative fitness in the population
	- then, sum_j X_j(t+1) = K + gamma \sum_j X_j(t)
	- population size converges on \bar{X} = K / (1-gamma)
	- selection is frequency dependent but not density dependent
	- supposedly the saturating yield is very similar to dynamics in more realistic models of density dependent competition (e.g. Levin et al. 1984)
- Invasibility criterion: for two different types d and D, d can invade D if rho(d,D) > 0 i.e. if log of geometric mean grwoth rate for a small type d subpopulation invating an established population of D individuals
	- rho(d, D) = E[ ln( (1-gamma) R(d-M_t)/R(D-M_t) + gamma ) ]
		- i.e., it's the log of the fitness ratio of non-surviving fraction of the population (?) plus the surviving one
		- for ESS D*, rho(d, D*) <= 0 for all d neq D, rho(D*, D*) optimizes rho? because the only strategy that can invade D* is D* itself? here partial derivative of rho wrt d is zero
	- using a Taylor series approximation of rho, including second derivatives wrt d and D, the signs of which are important...
	- rho_dd > 0 for all candidate ESSs means any monomorphic population is invasible
- For a non-negative function phi(x) where phi(0) = 0, phi(1) = 1:
	- rho_dd = C (gamma b V - E)
	- rho_DD = C (gamma b V + E)
	- V is variance in phi'(D* - M_t)
	- E is mean in phi''(D* - M_t
	- b > 0 is selection strength, C > 0 as well
	- n.b. V >= 0, and fluctuations in M_t means that V > 0
	- increasing gamma does not change D, E, or V
	- increasing gamma or V moves rho_dd from negative to positive in which case any  D* can be invaded
		- i.e. generational overlap "destabilizes" population equilibria
- Under Gaussian stabilizing selection, phi(x) = x^2, b = 1/(2*w^2), and the only candidate ESS is E(M_t)
	- here, condition for no ESS is: gamma var(M_t) / w^2 > 1 so increasing gamma likewise reduces likelihood of this
	- e.g., say s_r = 0, H = 1, so gamma = 0
		- ideal strategy then is to simply adopt the pheno with highest fitness in the avg. environment
		- but if s_r > 0 and the environment does vary, then a riskier phenotype can invade (or be invaded) by relying on the survivors (e.g., egg bank) to cushion survivors in a bad year
- Assuming random mating and infinite population size, these results extend to a diploid model with an arbitrary, finite number of loci
	- did not follow the argument presented...
- Simulations to confirm that ESS outcompetes a "coalition" of diverse invaders (because analytic work assumes that competition is pairwise and local)
	- confirm results
- Something about mixed strategies as well... oy
##### Discussion/Conclusions
- Generational overlap and environmental variation (or?) produce genetic and phenotypic variation
	- similar to Gillespie and Turelli 1989, so long as one genotype is best under all conditions, it's impossible to have genetic variance maintained by selection
- Frank and Slatkin 1990 demonstrate that spatial heterogeneity creating fitness variation can allow polymorphism; this model shows similar result for time
	- this model is similar to Gillespie's "c-haploid" model where c = 1 - gamma
- Other ways environmental variability can contribute to genetic variability:
	- lower fitness variance in heterozygotes compared to homozygoes (e.g., the heterozygote fitness is at least the arithmetic mean of respective homozygotes)
	- fluctuating selection leading to large fluctuations in the optimal phenotype produces mutation-selection model (Kondrashov)
	- Hairston and Munns 1984, Hairston 1988, Hairston and Dillon 1990 shows variation in a fitness trait of a marine copepod related to interannual variation in direction and intensity of selection
- The mixed strategies result: seems to be that genetic variation can not exist with a mixed-strategy bet hedger... or that mixed strategies are always invasible? something like that
- Note here that overlapping generations is suggested to slow down evolutionary change (e.g., Templeton and Levin 1979, Hairston and de Stasio 1988, Venable 1989) in addition to maintaining genetic variation
- Also note that with sustained directional selection, then individuals may also be removed from the genotype pool by this selection which can also reduce genetic variation

The modeling did not necessarily stick with me (the approach did though - invasibility, looking at relative fitnesses) but the analogy to the storage effect did. Overlapping generations (which really is just higher survival) and increasing environmental variation means the equilibria for monomorphisms become unstable. Overlapping generations means that iondividuals that were unfit in one environment can persist until later dates when they will be favored.

Interesting that this happens while also true that the persistence of this genetic variation will also slow the rate of adaptation. It's worth reading the other papers cited in here.

### Hairston, N.G., and De Stasio, B.T. 1988. Rate of evolution slowed by a dormant propagule pool. Nature.

- Theory suggests that dormant propagules will slow down rates of genetic change by shielding part of the gene pool from selection pressures
- Here: looking at the effect of egg diapause eggs in lake sediments for two populations of copepods
	- Study system: two small lakes in RI
	- *Diaptomus sanguineus*: two generations from November to May producing eggs that hatch immediately, and females produce diapausing eggs in spring that begin hatching in November
	- Fish predators, increasing in spring as water temperature rises
		- spring production of dormant eggs produces more predation:
		- e.g., fish removal due to drought in one pond means the switch to diapause happened later for the copepods
- Diapausing eggs do not all hatch in the season they are laid in
	- can remain viable for up to three years, perhaps (two lines of evidence)
	- but note they are susceptible to predatory flies which became abundant after prior absence in the drought-stricken pond
- Drought occurred in 1983; the 1983 population in this pond had delayed onset of diapausal eggs
	- no such change in the unaffected pond
	- 1986 and 1987 had intermediate dates of diapaused egg laying
		- many of these eggs laid in 1983 may not have hatched until 1985-1986
- After fish kill, selection favors individuals with delayed production of diapausing eggs
	- eggs hatching from the dormant pool slows down the rate of change... 
		- (what is the rate of change? returning to pre-fish kill levels?)
		- (eggs hatching during this period will also have delayed phenotype which is why the return to normal is slow?)
		- two years with no production of eggs (1983 - 1985) means that the only eggs in the pool are ones with delayed... wait or is this egg production at all? ahhh

Ahhhh. Okay this was such a bizarrely written paper. Just not clear, *too* concise. There was a drought, which led to delayed production of eggs because of the disappearance of predators. These eggs persisted in the egg bank for a few years. I think after this drought, should the populations have returned to normal? What was the deal with this mass reproductive failure? This belies the ecological need for the egg bank but not necessarily the evolutionary consequences of it. Why are the evolutionary consequences of this mass reproductive die off not examined? The emergence of the eggs from this year, which had the phenotype for delayed reproduction, should have delayed the return of the population to typical pre-drought levels? But what is this in comparison to?

A re-read may be useful but it probably would be better to just ignore this and focus on other papers!

### Pelletier, F., et al. 2007. The evolutionary demography of ecological change: linking trait variation and population growth. Science.

- Here: how does quantitative trait variation influence population growth in Soay sheep?
- Sheep: studied since 1985, structure and size of population is known for every year
	- birth weight collected each spring, adult body weight and hind leg length (size proxy) collected for ~50% of population each summer during summer catch
	- assigning paternity using genetic markers: ~60% of lambs, assignment done with >80% confidence
	- known additive genetic variance for several of these quantitative traits (age-specific or environment-specific??)
- Need to know how variation in a quantitative trait influences survival and recruitment (and sensitivity of population growth to survival and recruitment)
	- see Lande 1982, Ecology
	- approach: calculate the proportion of variation in contributions to pop growth p_t(i) that is accounted for by traits
		- p_t(i) is difference between observed population growth and population growth if the focal individual is removed (so it's jackknifing?)
		- (s_t(i) - \bar{s_t} + f_t(i) - \bar{f_t}) / (N_t-1)
		- proportion in variance in p_t(i) explained by variation in the trait - is this done with a regression? (no - with GAMs)
	- this was done across age/stage class before aggregating across whole population
- Body weight variation: 4.7% of population growth; Hind leg length variation: 3.2% of growth; birth weight: 1.7% of growth
	- most of the effects were in lambs and yearlings
	- small contributions in trait values occurs when trait variation does not account for much variation in individual contributions
		- in prime-age and senescent females, little variation in contribution explained (is this true?)
		- in adult males there are simply very few individuals so the effects are uncertain but also hold small weight
- Additive genetic variation: weight contributes 0.9% of population growth, hind leg 1.4% and birth weight 0.2%
	- citing Fisher: population growth is mean fitness, these values also give measures of heritability
- In low survival years, variation in body weight acccounted for 9% of variation in population growth rate
	- in higher-survival years, it accounts for only 4%, coming mostly from contributions of lambs
	- similar results for hind leg length
	- similar relationship between contributions and NAO - smaller contributions in low oscillation years, which are good for sheep
	- greatest opportunity for selection during harsh environments?
- Evidence of non-linear selection and age- and sex-specific responses to selection that also vary over time

Interesting approach but results look to be kinda messy and scatter-shot. Can't help but feel like there's a data scarcity thing happening here, but maybe not (3533 individuals over nearly 30 years?). I think this is worth re-reading. I like the question though, and the (sorta-?) straightforward approach to partitioning variation in a trait to an individual's contribution to population growth, and how this is also interpretable as selection strength.

I'm not quite sure what it means to say contribution of a trait's variation to population growth though. I'll need to think about this more. Contribution to variation in population (regression/variance partitioning approach) makes more sense to me. Are these saying the same thing? I don't think so. Good to think about sensitivity of population growth to a vital rate, then some combo of sensitivity of the vital rate to the environment (as Dan suggests) but also how to incorporate trait variation (or, in a sim approach, genotypic variation?) into this.

Basically, selection pressure varies across time/environments, and across age as well as sex. Different demographic stages have different selection pressures and these vary in how much they contribute to population growth.

For our simplified model, it seems like the thing here that is most likely to show up is that there's variation in how much variation contributes to population growth over time. But also thinking about relative effects of even just two vital rates (adult ? Age and sex I think are less important. Think more.

### Ozgul, A.,, et al. 2010. Coupled dynamics of body mass and population growth in response to environmental change. Nature.

- Yellow-bellied marmot near BV Colorado 1976 - 2008, using body type as phenotype because it determines survival during hibernation and reproduction upon emergence
	- earlier emergence from hibernation (phenological shift with cc) means earlier birth and therefore longer growing season
	- noted shifts in mean body size evident over study period
	- population size fluctuated around mean until ~2001, then increased in final seven years of study
	- body mass has a significant effect on multiple vital rates
- Fitted an IPM with this data
	- two IPMs: one for stable pre-2000 period, one for growth period after 2000
	- increase in growth (comparing the two growth periods?) was primarily adult survival and juvenile growth
		- (Fig. 2... change looks p stark for p survival of older adults, esp. for larger ones, while the growth functions are weirder... for growth post-2000 even smaller individuals still grow??
- Are increases in mean adult survival due to a change mean August mass in each age class, or a change in the shape of the relationship itself?
	- very clever approach here by comparing trait distributions before/after and their curves before/after
	- both processes contribute
- Decomposed change in body mass into contributions from selection versus other processes (?)
	- Age-structured Price Equation??? see Coulson and Tuljapurkar 2008 !!
 	- seems like mostly changes in mean growth rates and only small contributions of selection (?)
- Changing phenology means changing trait (phenotype) which means changing demography and population growth
	- however it seems like most of this change was ecological response, removal of a constraint, rather than actual evolution
	 - also the change in the functional dependence is interesting... what's going on here?

Interesting paper. Useful for (1) overall look at how traits and trait changes translate into demographic changes but also (2) how climate change can actually increase population growth through changes in growing season length. Of course (as authors acknowledge) the change in the growing season length may coincide with drought though, which may reverse some of these demographic gains.

Would be cool to look at growth (maturation) as a response! This does not necessarily need to be genetic, or at least I don't think so. It is not genetic (nor does it even appear to change) in the Schmid et al. paper.

### Rodriguez-Caro, R., et al. 2021. The limits of demograhpic buffering in coping with environmental variation. Oikos.

- Continuum of strategies for dealing with environmental variation: demographic buffering vs. demographic lability
	- labile species persist in stochastic environments by easily tracking the environment (Koons et al. 2009, Jongejans et al. 2010)
		- they do this by letting the most sensitive vital rates vary, responding to selection or environmental cues
	- demographically buffered species limit temporal variation in vital rates that population growth is most sensitive to (see Boyce et al. 2006, Hilde et al. 2020)
- Chelonians (incl. long-lived tetrapods like Testudinidae) evidence extreme demographic buffering
	- terrestrial tortoises have high variance in reproductive output but this barely influences lambda
	- tortoises typically have near constant through time adult survival rates
- Here: population viability of a tortoise in Spain
	- hypothesis: increases in adverse conditions will lower persistence odds, but the reproduction-survival trade-off will slow declines
	- Rodriguez-Caro et al. 2016 suggests population growth rates are ~1
	- ~decade of data
- Environmental effects: monthly temperature, precip, vegetation productivity
	- environmental drivers of population growth through a "moving window" approach - fitting environmental variables in different time windows to probability of reproduction and number of eggs then comparing models with multiple windows per climate variable
	- precipitation influences the probability of reproduction (positive), but no evidence of clutch size being affected by any environmental variables
- Integral projection model approach with two subkernels for traits
	- P subkernel includes a trait value, z, influencing growth conditioned on survival (size?)
	- F subkernal gives contribution of reproductive individuals given trait value (size?)
		- this one has a bunch parameters
	- Simulated population dynamics out for one site over a 100 year period
		- increased with 50% increase in drought conditions (from 1/10 years to 1.5/10 years)
		- also searched for drought rates to see threshold for population growth rate over drought gradient
	- model possible survival/reproduction trade-offs: constant survival, trade-off where survival increases in drought years, and survival decreasing in drought years (in all cases prob. of reproduction is cut nearly in half in drought years)
##### Results
- Under current climate with constant adult survival or survival-reproduction trade-off, population is demographically stable with log(lambda_s) = 0.0000
	- under positive correlation population is unstable log(lambda_s) = -0.0011
	- with increasing drought positively correlated survival has strong ~linear decline while other survival-reproduction regimes have only slight decline
- Under normal conditions, growth has the largest positive effect on population growth (esp. those of adults)
	- similar conditions under drought; under positive correlation the effects of growth are huge
##### Discussion
- Results here suggest that so long as survival decreases in drought years, populations will not be viable under any drought recurrence scenario
	- (viable yes and/or weakly affected if survival is uncorrelated or if there is a trade-off under drought)
- As usual, growth and survival are important for long-lived species
	- growth more important than survival in this case
- Autocorrelations - important to include, missing here

Interesting study although a little bit limited. Here, for a long-lived species, reproduction is what is most obviously sensitive to the environment. Not really mentioned here is that these effects only under incredibly extreme conditions are enough to make populations inviable when survival is held constant. However, when survival decreases as well, then the droughts threaten viability.

What does this mean for this project? Just that survival of longer-lived individuals is in fact very important for their long-term persistence. Growth as well. Other than that not incredibly much.

### Compagnoni, A., et al. 2021. Herbaceous perennial plants with short generation time have stronger responses to climate anomalies than those with longer generation time. Nature communications.

- Climate change: influencing mean and variance in temperatures and precip
	- who will be most vulnerable? well, plants (sessile) that can't demographically buffer...
- Water availability influences NPP, a proxy for population growth
	- important for seed germination, tissue growth, floral induction, seed set
	- temperature modulates water availability
	- more pronounced effects in arid and cold biomes than wet and temperate ones?
		- aridity means more water limitation
		- cold means more frequent extreme cold spells (too cold for tissue growth)
		- confounding here looking across biomes - could it just be different plant functional types in different environments?
- Generation time: how quickly individuals in a population are "substituted"
	- longer-lived species should have less response (in growth rate) to fluctuations
- Here:
	- is population growth rate more strongly associated with precipitation than temperature?
	- is population growth in water-limited biomes more responsive to precip anomalies?
	- is population growth in cold biomes more responsive to temperature anomalies (than pops in warmer biomes)
	- are longer-lived species less responsive to anomalies?
##### Results
- precipitation anomalies had larger effects on population growth rates than temperature anomalies did (and interactive effects were weak/ns)
	- precipitation with one s.d. over mean increased lambda by ~3%
- no effects (from meta-regression) suggesting that colder or water-limited biomes had stronger respective responses to relevant stressors
- response to both types of climate anomalies was negatively correlated with longevity
- results held when removing graminoids (which, according to Tukey HSD, may have different patterns than herbaceous species)
##### Discussion
- Importance of water on plant performance: forecasts involving precipitation are more uncertain than those of temperature, meaning population projections may also be uncertain
- Results here complement Morris et al. (2008)'s results that longevity can buffer populations from variation in survival and reproduction
- The fact that responses do not change based on biome suggests that populations are (demographically) adapted to cope with climate variation in their respecrtive biomes
	- (what does it mean to be demographically adapted... or rather, not demographically adapted?)
	- maybe, e.g., biomass accumulation is decoupled from demographic processes?
- (of course, possible taxonomic and geographic/biome biases in dataset)
- Mechanistic models featuring microclimate information... (Lembrechts and Lenoir 2020?)
- Pumping the breaks on the generation times and climate sensitivity result
	- this study does not look at density dependence, trophic interactions, anthropogenic drivers, etc. so hard to tell if effects are due to climate or something else
	- also, these results are more relevant to changes in climatic variability than means, and with changing means extrapolation may fail
	- cons bio lit suggests that shorter-lived species may be less (rather than more) subject to climate variability (lower extinction risk?)
##### Methods
- Compadre and Padrino databases, looking for density-independent models with at least six transition matrices
	- 48 species, 144 populations
	- plus some other sp with relationships to climatic drivers? total of 62 sp, 162 populations
	- get lambda for each transition (>3700 matrices, 52 IPMs)
- Climate: 1km2 gridded estimates of monthly min/max temp, total precip from some online dataset (CHELSA?), 1901 - 2016
	- standardized z-scores in 12 months preceding census give strength of anomalies (x - E(x)) / sig(x); abs(z) > 1 is an anomaly
	- skew present in only a small number of populations...
- Autocorrelation modeled in meta-regressions by having an autoregressive error term (ooh)

This is different than I thought... for some reason I assumed this was about evolutionary response?

Validates (on a large, although possibly biased, set of species/populations) the theoretical finding that longer-lived populations will have less response in their population growth rates relative to shorter-lived populations. This is still useful for thinking about population fate in the face of increasing population variability! But not for rates of evolution, at least I don't think so.

Also some interesting stuff here (H1-3!) about specific drivers in specific environments. Namely, that preciptation anomalies matter more than temperature anomalies? Useful for thinking about fates under climate change. Relevant to Thermopsis projects, maybe?

The Cotto paper may be a good read next... I think there was some stuff in there about the rate of adaptation.

### Lanfear, R., Kokko, H., and Eyre-Walker, A. 2014. Population size and the rate of evolution. TREE.

- "Substitution" is new mutation that spreads to fixation, so depends on both mutation rate and fixation rate
	- DNA sequence data allows easy observation of substitution
##### Neutral and nearly neutral mutations
- Some cases where fitness effects are close to zero (s approx 0, Ne|s|<<1) in which case fate of mutation is due to drift
	- weaker drift as Ne increases
	- more neutral/neutral-ish mutations in larger populations is balanced by decreased probability of fixation through drift
	- therefore the substitution rate is just the mutation rate independent of Ne
- N.b. Balloux and Lehmann (2012) showed that size fluctuations and overlapping generations will modify the equality relationship between neutral substitution rate and mutation rates (likely due to changes in mutation rate (?))
	- note also that Ne can be associated with factors like generation time or size fluctuations (larger Ne associated with shorter generation time... wait but wouldn't this also mean that there would be more fluctuations in size?)
##### Non-neutral mutations
- With increasing Ne, natural selection becomres more efficient for producing advantageous mutations
	- theory suggests that the "power" of selection with Ne increases faster than the production of new mutations (i.e., selection gets better at operating at higher Ne, removing mutations, more than it produces more mutations) - Akashi et al. 2012 review on this
		- fewer deleterious substitutions with Ne, more advangateous rates with Ne so long as simplifying assumptions are met
- If mutation rates are high or Ne is large, then whether or not Ne increases or reduces rates of evolution depends on how adaptation influences 	
- [some stuff on selective sweeps etc. that I did not read]
##### Across mutations
- Mutations of different fitness effects have different relationships with Ne
	- "we expect [NeRR - relationship between Ne and rate of adaptaiton] to have a U shape" - large and decreasing for small Ne and increasing when Ne is large?

[did not finish]

I'm not sure this is so relevant because we don't really have Ne (other than N I guess?)... I just want a distribution of fitness effects for mutations.
