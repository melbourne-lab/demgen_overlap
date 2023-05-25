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

### Burger, R., and Lynch, M. 1995. Evolution and extinction in a changing environment: a quantitative-genetic analysis. Evolution.

- Some previous models looking at critical rate of long-term environmental change beyond which extinction is certain: Lynch et al. 1991, Lynch and Lande 1993
	- demonstrate that extinction is certain if environmental change is faster than the maximum sustainable rate of evolution
	- Lynch and Lande 1993: small, sexually reproducing populations, but deterministic analysis suggests that rates of environmental change below the "critical rate" means no extinction (ofc. can't be correct)
	- note that finite populations experience drift from which many generations may be needed to recover
		- these instances of drift may cause large temporary lags behind the environmental change
##### Model
- Randomly mating finite population, discrete generations, density dependent population growth
	- gaussian selection with optimum phenotype theta_t = kt + eps_th (eps_th ~ N(0, sig^_th)
	- intrinsic growth rate is R_t = B \bar{W}_t (for mean fitness \bar{W}_t)
	- ceiling-type density dependence with carrying capacity K
- Lynch and Lande (1993): critical rate of change is k_C is value of k beyond which extinction is certain (because \bar{W}_t eventually falls below 1/B)
##### Analytical work:
- assume that k and sig^2_th are sufficiently small, \bar{W}_t sufficiently large to have approx. constant N_e
- def. s = sig^2_g / (sig^2_g + sig^2_e + w^2) (a measure of strength of selection)
- def. a recursion for the phenotypic distribution Phi(g_{t+1}) in terms of Phi(g_t) and \bar{g}_t
- then:
	- E(g_{t+1}) = E(g_t) + s(kt - E(g_t)) 
		- (i.e., mean genotype is the previous generation's genotype, plus heritable portion of the gap created by environmental shift) (but why is it kt and not k?)
	- V(g_{t+1}) = (\sig^2_g / N_e) + (1 - s)^2V(g_t) + s^2 \sigma^2_theta
		- (i.e., weighted variance of environmental trend-variance (weighted by sq. of selection strength) and prev. generation's variance (sq. of reverse of sel. strength...) plus some amount of drift?)
		- (think about this!)
- Recursion relationships for E(g_t) and V(g_t)..
	- influence of initial genotype wanes as t grows large s/t
		- E(g(t)) -> kt - k/s
		- V(g(t)) -> (w^2 + sig^2_e) / (2N_e) + (sig^2_g \sigma^2_th) / (2(w^2 + sig^2_e))
		- so, mean population genotype does track the change (kt term) but lags behind by amount k/s
			- lag is larger with faster change, less selection ((w^2 + sig^2_e + sig^2_g) / sig^2_g is large)
			- under weak selection rel. to genetic variance, t reaches ~95% of its asymptotic value at around 3/s (why?)
- lambda_t  = (B w) / (sqrt(V_{lambda,t}) * exp(- (E(g_t) - kt)^2 / (2 V_{lambda,t})
	- for V_lambda,t = sig^2_g + sig^2_e + w^2 + V(g_t) + sig^2_th
	- (note that B w / sqrt(V_{lambda,t}) is growth rate if population perfectlyt tracks environment)
	- from this, the critical rate is
		- k_C = s sqrt(2 V_lambda ln(B w / sqrt(V_lambda))
		- that is, increases with selection strength and the log of the growth rate under ideal circumstances?, decreases with sqrt of variances...
##### Extinction times
- two phases of "guaranteed" extinction process: one is time until lambda_t = 1 (population stabilizes) then time when population is decreasing
- assume that V_lambda is approx constant (occurs if genetic variance and variance in mean breeding value is small relative to sig^2_theta, sig^2_e, w^2)
	- let k = kappa * k_C (kappa > 1)
	- lambda_t = (B w / sqrt(V_lambda)) ^ ((1 - kappa^2) * (1 - (1 - s)^t)^2
	- from this it can be demonstrated that the time until lambda_t = 1 is independent of population sie and genetic variance (!)
	- phase 2 is shorter (because of increasingly negative growth rates?)
- time to extinction will increase when including fluctuations in R_t, demographic stochasticity, autocorrelation)
##### Simulations
- Random sampling by pairing of parental breeding values
- No differences or stochasticity in fertility/fecundity (each breeding pair produces 2B offspring)
	- mutations occur at rate n \mu
		- (could be useful... see Lande 1976, Lynch 1988)
- Viability selection before reproduction where survivors of viability are next gen's parents (stochastic process)
###### Sim results
- With large w (weak selection), large lag, and additive variance operates nearly independently of selection strength
	- but with strong selection additive variance declines quickly
	- so, when the rate of env. change is slow, intermediate selection sgrengths are favored? (because with strong selection, variance is quickly lost, but with weak selection, there is a strong lag?)
		- but with a more rapidly changing environment, it's monotonic (strong selection means rapid extinction)
- Larger population size (K) has *weak* effect on risk of extinction
	- with sufficiently rapid environmental change, time to extinction is not so different between large and small populations
	- seems like this is because additive genetic variance does not grow indefinitely with increasing K, k_C asympotically approaches a constant with increasing K
	- with slow environmental change, for sufficiently large populations extinction is rare? little stochasticity?
- Extinction times: distribution is tight for very fast change (sd much smaller than mean) but broad with slower change particularly if below critical density (more variation in extinction times - because extinction relies on stochasticity?)
	- so we can only have confident conclusions in time to extinction for very rapid env. change
	- extinction occurs "deterministically" with rapid change, as populations very quickly have growth rates fall below replacement
	- but for slower k (near or below critical value), mutations allow populations to maintain genetic variation and populations have a small lag
- Additive genetic variance increases for small k and then decreases with increasing k (with larger effects for larger populations)
	- high correlation between genetic variance and extinction time (within param combo) so more genetic variance means longer time to extinction
- With fluctuating environmental optima (k = 0, sig^2_e > 0)
	- monomorphic populations (mu = 0, sigma^2_g = 0, E(g) = V(g) = 0) perform better in these situations than genetically variable populations
		- this is due to the lag in responding to the direction of the environment (no autocorrelation means 50% of moving the wrong way) and variance load
	- with sig^2_e approaching w^2, any population can go extinct suddenly regardless of genetic makeup
##### Discussion
- "if the rate of environmental change is sufficiently slow and the amount of genetic variance for the trait is sufficiently high, the population mean phenotype settles into a quasi-steady state lag behind the environmental optimum"
	- magnitude of the lag influences fitness, influencing population growth
		- lag is influenced by genetic variance, and drift/stochasticity can produce temporary periods of even longer lag
	- with sufficiently high lag, extinction is guaranteed due to fitness below replacement
	- there is a critical value of environmental change that populations can track - typically according to these results as ~10% of one phenotypic standard deviation per generation
		- but factoring in other things, e.g., a more kurtotic distribution of mutational effects, could mean the critical rates could be ~1% of a phenotypic sd
- Results here are similar to those of Huey and Kingsolver (1993) who found that there is an intermediate width of fitness function for which persistence odds are higher
	- here, similar result: with too wide a fitness function, too large of a lag, but with too narrow a fitness function, large load due to slow adaptation
	- "This suggests that in a slowly, but steadily, changing environment, braod generalists and narrow specialists will be most vulnerable to extinction"
- Barrier: theory does not (yet) have good descriptions of extinction times for rates of change near the critical point (mean let alone entire distribution)
	
Super interesting paper - glad I read it. Will need re-reads.

Gist: there's a critical rate of change (k_C) determined by heritability, selection strength, and variance (in environment, phenotype, etc.) above which deterministic population extinction is guaranteed. Below (or at) this critical value, populations will still go extinct due to stochasticity but on longer timescales (more variance).

Populations will (with large t) have a lag proportional to k/s (why?) - higher additive genetic variance rel. to non-genetic phenotypic variance and width of fitness function means less lag

Some cool analytical approximations in here although there are some cases (ignoring various types of stochasticity) where they fail.

As far as previous simulations go - this does explain previous simulation results where populations initially increase (carrying capacity or something) then decrease - an inverse of ghe G&H U. Extinction dynamics here have two parts - first decline to lambda = 1, then lambda below 1. This creates that shape. Also the guidance at the end about speed of change rel. to phenotypic variance is useful. Needs to be very small! E.g., for sigma^2_z = 3, a critical rate of change between 0.03 anmd 0.3?


### Helantera, H., and Uller, T. 2020. Different prospectives on non-genetic inheritance illustrate the versatile utility of the Price equation in evolutionary biology. Phil Trans R Soc B.

- Price equation gives average of a trait between generations
	- (change in mean trait) x (mean fitness) = cov(fitness, trait) + (fitness-weighted change in phenotype between parents and offspring)
	- (change in mean trait) x (mean fitness) = cov(parent fitness, offspring phenotype) + (expected change in trait between generations)
		- first term is selection, second term is transmission
	- (change in mean trait) = (slope of parent-offspring regression) x (selection differntial) + cov(fitness-weighted offspring phenotype, parent phenotype) + (expected change in trait between generations)
		- selection pressure is Cov(rel. fitness of parent, phenotype)
		- breeder's equation is special case of this where cov = 0 and expected change in trait between generations is zero
- Phenothypic change between generations - not trivial! (with multi-cellularity at least)
	- e.g., natal environment, early development environment
	- non-genetic (or "extra-genetic") aspects of inheritance - often beyond parental control, very complex
		- this can vary over time and/or take many forms
		- this complexity means typically it is modeled or conceptualized in only ideal terms
- Here: fitting these non-genetic aspects of inheritance into an evo framework using the Price equation
##### Transmission model of inheritance
- Conceptually similar to how we think about genetic forms of inheritance - some generalizations and transferability of intuition and concepts from, e.g., pop gen
	- but these may get strained or become overly-complex when modeling multiple loci
- Day and Bonduriansky's generalization:
	- 2+ separate channels of inheritance, at least one of which is genetic
	- Variants (that which varies? idk) could be alleles, epi-alleles, quantitative phenotypes, breeding values, maternal resources, etc.
	- modified the Price equation to get change in population mean values of genetic (g) and non-genetic (h - lol!) components
	- (change in mean g (h)) x (mean fitness) = cov(fitness, g (h)) + fecundity-weighted mean transmitted change in g (h) + survival-weighted mean within-lifetime change in g (h)
		- this allows for overlapping generations and within-generation change
	- this can be re-written to say change in g (h) also depends on variance in g (h) and also covariances between the two, multiplied resp. by the selection gradient at present values
		- selection may act on the non-genetic value h both directly and indirectly
	- framework here can explain why covariance arises and how it influences evolutionary trajectories
		- also (of less interest to my model) conditional fitness effects, e.g., fitness effects of g depend on value of h
- Transmission perspective on inheritance (e.g., Day and Bonduriansky) - consequences of any form of inheritance rely on:
	- how the variant affects fitness
	- how the transmission occurs
	- how the phenotype (not the variant?) changes over the lifetime
- Another view - transmission of information between generations
	- e.g., DNA of parents contains some info about the state of the environment that offspring may encounter
##### Phenotypic covariance models
- More in line with the quantitative genetic approach
	- don't track individual variants
	- instead consider how much covariance between phenotypes of parents and offspring is determined by shared genes (with additive effects)
- Assumptions: many loci of small additive effect, joint multivariate normal
	- commonly assumed that the mean phenotype doesn't change between generations unless there is selection or drift
	- these assumptions give the following form of the Price equation:
		- change in mean trait = (regression slope of offspring phenotype on parent phenotype) x (selection differential)
			- selection differential is the covariance between phenotype and fitness
	- heritability here can also be defined as ratio of additive genetic variance to total phenotypic variance
- Note that genes are not the only thing that influence the slope of offspring phenotypes onto genotypes
	- this is one reason why additive genetic variance may be preferred for quantifying heritability
	- in facto total phenotypic variance can be partitioned into the sum of additive genetic variance, the dominance variance, epistatic variance, and environmental variance
		- typically this approach labels non-genetic effects of the parental phenotype as maternal effects (although sometimes maternal effects genotype has an influence as well causing interesting dynamics)
	- the variance partitioning could instead be rewritten as all genetic components fall under Vg (additive genetic, dominance, epistatic), then transmitted non-genetic variation (e.g., epigenetic), and non-transmitted
		- this approach would be useful for, e.g., seeing how non-genetic inheritance affects ability of populations to track environmental change
		- a para in here with possible ways to quantify these (empirically challenging)
##### Developmental models of heredity
- Inheritance does not necessarily need to be conceptualized as transmission and does not *a priori* necessitate more relevance of certain parent-offspring relations (all are equally important? not quite sure what to make of this)
- See models by Rice (multiple) where fitness and offspring phenotypes are random variables (not known at time of reproduction)
	- considers the full distribution of offspring-parent distribution rather than the mean (?) (is this distribution bivariate?)
	- the regression slope or additive genetic variance approaches rely on very specific circumstances, often not met
		- Rice approach takes more information about the adult-offspring relationship (distribution) into account (moments)	
		- (although despite the many published studies and datasets no systematic review of how common non-linearities are)
- Need more work looking at which developmental processes do/don't change the parent-offspring phenotype distribution esp. with a complex genotype-phenotype map

Hmm... okay interesting. Non-genetic inheritance I suppose is more about maternal effects, epigenetics, etc. I guess what I was thinking about is non-inherited traits! I think this is worth revisiting at some point - the quant gen perspective was maybe interesting. Is my quant gen model really modeling additive genetic variance, or really only the sum of inherited components? Is this the difference between broad- and narrow-sense heritability?

The information part also seemed tantalizing. 

### Vander Wal, E., et al. 2013. Evolutionary rescue in vertebrates: evidence, applicaitons, and uncertainty.

- Conservation biology typically preserves neutral genetic variation, not adaptive potential linked to demography (which makes evolutionary rescue possible)
	- is evolutionary rescue prevalent in nature? lab and theory yes
- Documenting rescue: requires evidence of decline, evidence of phenotypic change *through standing variation or mutation* and subsequent population recovery
	- "ecological rescues" may occur through plasticity or other factors
	- obtaining both demographic *and* genetic/evolutionary data is difficult
- Rats and anti-vitamin K resistance
	- anti-vitamin K pesticides used for rodent control (mortality and lethal breeding)
	- individuals from populations exposed to the pesticide fared better than individuals from unexposed populations
	- some combination of standing variation and up to six new mutations
- Australian rabbits exposed to the myxoma virus
	- killed 99% of animals infected and caused population decline across the continent
	- selection:
		- (1) for a less virulent myxoma strain which slowed mortality
		- (2) evolution of resistance to myxoma
	- lab experiments show individuals from exposed populations recover while those in unexposed populations did not
	- the co-evolutionary angle complicates things
- Adaptive tracking: recorded in long-term studies of Galapagos finches, great tits, red squirrels
	- these evidence don't include evidence, though, that failure to track the environment means extinction
- Genetic rescue studies as good examples of fitness-linked traits and their importance for population rescue
- Bell & Gonzalez estimate ~25 generations needed for ER to take place
	- similar or longer times observed in insects in lab settings
	- but 25 generations can be a long time for long-lived organisms
	- more gradual environmental change (relative to the generation time of the species) may alleviate extinction odds

Some stuff here I only skimmed (the phenotype-demography map stuff seems interesting!). I was mainly interested in the claims about generation time, which I have seen cited a few times from this paper.

This actually seems like a hand-waivy argument! It comes from Bell and Gonzalez and their estimation that rescue would take ~25 generations. But Bell and Gonzalez were probably working with a semelparous population, where generation and timestep are conflated. This doesn't seem to be considered here, and instead the relationship between longevity and time to rescue is really just multiplying generation time by generations to rescue. Room for improvement here!

### Lande, R. 2009. Adaptation to an extraordinary environment by evolution of phenotypic plasticity and genetic assimilation. Journal of Evolutionary Biology.

- The "Baldwin effect" (1890s!): plasticity of development can produce a partially-adaptive phenotype following environmental change; natural selection furthers the adaptiveness of the phenotype
	- [from skimming some Wiki articles though perhaps we are more interested in "genetic assimilation" - Waddington - than the Baldwin effect?]
- In evolutionary time, rare anomalies in the environment might necessitate the capacity to accelerate phenotypic adaptation through "transient" evolution of plasticity
	- might be further needed for persistence under widespread global environmental change that is currently occurring
##### Conceptual background
- Reaction norm of a benotype: breeding value as a function of the environment in which offspring develop
	- (with purely additive genetic variance, breeding value and total additive effect on offspring are the same)
	- "elevation" (intercept?) of reaction norm is breeding value
	- environment on x-axis, phenotype on y-axis
	- consider WLOG the initial average environment 0
		- reaction norms are likely to converge around intercepts (elevations) of a certain value, call this phenotype A, for environment 0 (canalization)
		- this canaliation means reduced phenotypic and genotypic variation
		- however, even with a canalization occurring in the normal environment, there may still be considerable variation (among individuals) in the *slopes* of their norms
	- if environment shifts suddenly to state delta, with optimum genotype A+B*delta (such that B is the optimal slope for the norm), adaptation will proceed in two stages: mean population slope approaches B, then intercept/elevation move correspondingly
		- this is because with canalization on the phenotype, there will be more (genetic) variation in the reaction norm slope B compared to the elevation A
- Waddington's "genetic assimilation" as another mechanism for evolution of plasticity
	- Waddington's experiments feature a new phenotype expressed as a plastic response, maintained by selection to the point that it occurs frquently even in the original environment
	- "Waddington apparently misinterpreted his selection experiments on threshold characters, which do not necessarily involve the evolution of plasticity" - not really sure what this means... or at least the second half of the sentence - threshold characters are not subject to plasticity?
	- "Here I define genetic assimilation in an altered environment as the reduction in plasticity and its replacement by genetic evolution, while maintaining the phenotype initially produced by plasticity in the altered environment"
	- reduction in plasticity attributable to the cost of maintaining it
- Models of evolution of plasticity require environmental predictability over time spans relevant to individual development
	- with heterogeneity (in time/space) the environment of development is often different from the environment individuals are selected for
	- this means the slope of the population norm is a fraction of the optimum phenotype's slope
		- (the fraction?) is reduced by correlations between environments of development and selection... [so more correlation means lower slope... I think this makes sense...]
##### Model
- z_t = a + (b eps_{t-tau}) + e
	- phenotype is elevation (a, the breeding value in the ref environment where eps = 0)
	- b is the slope
	- juveniles in generation t are exposed to environment eps_{t - tau}
	- tau is a fraction of a generation before the adult phenotype is expressed and subject to selection
	- e is N(o, sigma_e) noise
	- a, b are bivariate normal with additive variances G_aa, G_bb, covariance G_ab remaining constant
	- as such sigma^2_z = G_aa + 2G_ab eps_{t-tau} + Gbb (eps_{t-tau})^2 + sig^2_e
		- variances will be minimized for eps = -G_ab / G_bb
		- if a population evolves (evolutionary history) at eps = 0, G_ab = 0 must be true to satisfy the above, i.e., slopes and elevations must be uncorrelated
	- change in \bar{a}, \bar{b} is equal to G matrix times selection and selection gradient beta
- evolution of the mean reaction norm (\bar{a}, \bar{b}) in a fluctuating environment with mean of eps = 0, variance sigma^2_eps, autocorrelation rho
	- gaussian selection with width much greater than the phenotypic standard deviation (i.e. gamma >> E[sigma^2_z])
	- upon/before the environmental shift the expectation b_0 is rho_tau B due to adaptation to background stochasticity
		- because for \bar{eps} = 0, E[eps_t | eps_{t - tau}] = rho eps_{t-tau}, and relatedly E[theta_t | epsilon_{t-Tau}] = A + rho B eps_{t-tau} (for optimum theta)
- Environmental change: e_t = U_t delta + xi_t, where U_t = 1 for t > 0 and xi_t is environmental stochasticity
	- expectation after adaptation is completed is E[\bar{a}] = A + (1-rho_tau)B delta, E[\bar{b}] = rho_tau B
	- in this case, expectation E[\bar{b}] = A + B delta as expected
	- for sigma^2_xi << delta^2 (i.e., a very large environmental shift relative to strength of environmental fluctuations), sigma^2_z << omega^2, and G_bb delta^2 approx. equal to G_aa + G_bb delta^2 i.e., G_bb delta^2 >> G_aa (much more additive genetic variance in plasticity than breeding values and/or very large environmental change), 
		- analysis demonstratesthat "phase 2" of genetic assimilation occurs much more slowly than "phase 1" of shift in mean plasticity
##### Discussion
- the phi parameter being close to 1 (i.e., G_bb delta^2 >> G_aa) appears to be important for all of this
	- initial mean plasticity (upon environmental change) is partially adaptive, below the optimum by some fraction rho
	- the adaptive plasticity here (and increased additive genetic variance... hmm...) mean that in the first generation under the new environment, mean fitness decreases but there is transient evolution of increased plasticity
	- "The time scale for the mean phenotype to closely approach the optimum in the new average environment is 1 - phi times shorter than by Darwinian evolution alone (with no plasticity) or by the Baldwin effect (with constant plasticity)"
		- so phi tells us the rate... must be connected somehow to this increase in additive genetic variance (for the trait?)
	- during phase 1 only slight evolution in breeding values
- Phase 2 (genetic assimilation) means reduced plasticity nearly compensated by a changing breeding value as the "small fraction" of adaptive change in mean phenotype is completed
- "Accelerated phenotypic adaptation during phase 1 can alternatively be interpreted as resulting from increased genetic variance in the new environment"
	- also implies an increased genetic correlation between slope and elevation which slows the assimilation (breeding value change)
	- stabilizing selection in the new environment(s) should reduce genetic variance and diminish the correlation eventually
- This might produce the effect of allowing population persistence without rapid phenotypic adaptation allowed by the transient evolution of plasticity

DANG this is cool and a gold mine. I think I can borrow from some of these approaches.

This will require more reading to solidify but so far I think what is going on is: canalization (stabilizing selection) reduces additive genetic variance in breeding values but not necessarily in "slope" of the response (which requires environmental autocorrelations?) (analogous to non-genetic phenotypic variation in our model). I think upon environmental change, this additive variance in plasticity increases and allows rapid change of the phenotype relative to the breeding value (where genetic variation is still low). So there is rapid adaptation of the phenotype followed by a slow catching-up of the breeding value. 

Math approaches here might be useful. There is a variance/covariance matrix involved and I think this can be applied to my model (pheno-geno correlation). The phi parameter I think might relate more to the survival parameter in my model than the heritability but I am not quite sure. I do feel like these models have to be analogous in some way though (although in my model rho of the environment is zero... what maintains the phenotypic variance?) The standing variance in plasticity in the pre-shift environment is a fraction of its optimum, and immediately following the shift it appears that this plasticity can quickly move from the fraction towards its optimum... B in this model clearly seems analogous to e in mine.

Fun fun! Re-read at some point. Not even sure how much additional knowledge there is to be gained by reading Chevin and Lande 2010.
 
### Chapin, F.S., et al. 1993. Evolution of suites of traits in response to environmental stress. Am Nat

- "Stress resistance syndrome" - a suite of traits forming an adaptive strategy for tolerating stressful environments
	- could be linkages of physiological and developmental linkages among traits
	- Contemporary (and older) debate about "large evolutionary changes" (?) and whether they are due to small mutational events or few mutations with a lot of pleiotropy
	- are these stress-related traits (esp. in plants) the result of a small number of genetic changes with large pleiotropy, or the accumulation of many independent evolutionary effects operating in parallel?
- All plants show traits in a low relative growth rate (RGR), low resource acquisition, slow tissue turnover, etc. in low-resource environments
	- some but not all plants are specifically adapted to low resource environments
	- plants adapted to high-resource environments will show plasticity in these traits ("broader" reaction norm) than plants in low-resource environments
- Genetic change in an underlying trait (e.g., a hormone) can turn on the SRS allowing it to be expressed in a broader range of environmental circumstances
	- (turn a gene on to enable a bunch of stress response stuff...)
	- ecological sorting or natural selection for the change could then match genes with habitats
	- in plants, a lot of regulatory hormones are caused by a small number of genes with high heritability
	- for quantitative traits, a few major genes of large effect at low frequency have a larger selection response under *strong* selection than a large number of genes of small effect at intermediate frequency
		- but under a small-to-moderate environmental change, changes in many genes of small effect will have a larger response
	- so, this makes it plausible that there would be a small number of genes of large effect responding to selection
- Connection growth to stress resistance:
	- Slow growth typically means slow leaf turnover (conserving carbon and nutrients lost during senescence)
	- Slow growth also minimizes respiration associated with producing new tissues
	- All-in-all slow growth minimizes dependence on the environment for uptake of new resources
		- this also allows relatively more allocation to other processes which might be helpful under stress resistance (e.g., storage)
- Plant stress response systems are centralized and rely on hormones
	- environmental stresses can trigger changes in hormonal balances
	- this is an preventative, early-acting system that reduces growth and changes allocation before there are severe imbalances (?)
- Seed size influences seedling size and produces a variety of size-related traits in young plants which have cascading effects throughout life
	- large seeds produce seedlings with large leaves, producing self-shading (?) and an overall low photosynthetic rate and high tissue maintenance might produce
		- large seeds are typically associated with low RGR and perhaps this is why
	- experimental reduction of seed size is associated with reduced size and increased RGR, but also decreased reproductive effort
	- likewise removing leaf tissue (experimentally) increases whole-plant RGR and photosynthetic rate in remaining leaves
	- large seeds are adaptive for germination in dry or shaded habitats
- Likewise, tissue nitrogen, photosynthesis, stomatal conductance, and transpiration are related to each other
	- simple genetic changes altering tissue N concentrations could alter a bunch of related traits as well
- Trade-offs often account for (negative?) correlations among traits
	- RGR and defense are typically inversely correlated (low RGR associated with higher defense - makes sense...)
	- (...)
- "Implications for global change" section hypothesizes that shorter-lived species will show faster change due to more turnover
	- cite Antonoics papers and a paper by Corinne's dad (lmao) about intensive selection for high abscissic acid

lmao... another citation for a hypothesis instead of a test/result. Well let's test it! The rest of this stuff was not supre relevant although there was some good stuff about pleiotropy vs. alleles of large effect...


### Cotto, O., et al. 2017. A dynamic eco-evolutionary model predicts slow response of alpine plants to climate warming. Nature Comms.

- SDMs are used for assessing biodiversity loss on long-ish timescales
	-  But, they do not account for eco-evolutionary processes, even though populations may adapt to novel conditions
- Mountain species are especially threatened by climate change
	- Nobody is modeling or making predictions from models including evolutionary dynamics of these species
	- How does local additive genetic variance, life history, landscape structure, dispersal influence rapid evolutionary adaptation?
	- Furthermore most relevant studies to adaptation focus on short-lived zpecies!
		- Theory predicts that long-lived species in stressful mountain environments should have slower evolutionary response (Savolainen et al., 2004, Kuparinen et al., 2010)
- Here: DEEMs (dynaimc eco-evolutionary models) that combine niche-based projections from SDMs with empirical data
	- first, predict current distribution with "static ecological niche models" (SENMs)
	- then use DEEMs to simulate changes in distribution given adaptation driven by different scenarios of environmental change
	- Use four endemic alp species in 15 landscapes of Austrian Alps
	- environmental variables modeled: bedrock carbonates, mean annual temprature, mean annual precipitation
		- model environmental change with IPCC forecasts
	- compare SENMs and DEEMs
	- perform sensitivity analysis of DEEM projectsions to additive genetic variance, strength of selection on survival, adult survival rate
##### Methods
- Hill-Smith analysis (?) to select species to model
- DEEMs fit with an approach like  Nemo
	- Three polygenic quantitative traits, one corresponding to each environmental variable
	- Ten unlinked aditive diploid loci with pleiotropic mutations (at rate mu)
		- mutations drawn from 3D gaussian with fixed variance
		- "continuum-of-allele model" (see ref 49)
		- traits is sum of genotypic value plus random non-genetic component, h2 ~ 0.3
	- Hermaphroditic species with four life stages: seeds, seedlings, pre-reprod adults and reprod adults
		- life cycle: mating, seed production, seed dispersal, aging, seed germination, seed survival in seed bank, clonal reproduction, seedling competition (density-dependent Beverton-Holt), seedling viability selection
		- two-year maturation period before reaching reproductive adult stage
		- adult survival determined by parameter s_a (not necessarily known, tried values 0.7, 0.8, 0.9)
		- seedlings compete against other seedlings and adults previously established in favorable microsites (!!)
		- viability selection is Gaussian, uncorrelated selection between traits
		- migration also occurs (somehow)
	- Simulations in two phases: burn-in allowing species to reach genetic-demographic equilibria and colonization, then 150 years under RCP scenarios
		- 10 reps per combo of species, spatial grid, RCP scenario, mutation rate, selection regime
		- data saved (before selection for each ten year window): individuals per site, individuals per age class, average genotype and phenotype, variance of genotype and phenotype, trait values, average seedling fitness
##### Results
- SENMs predict contraction while DEMs predict early expansion (colonization outside initial range?) and then decrease until final year (extinction debt?)
	- DEEM range size contraction occurs via maladaptation to local environmental change
	- Fig 1: DEEMs show greater change in final relative occupancy as well as a longer phase of decline
	- DEEMs show a decrease in habitat suitability over time (Fig. 2)
	- DEEMs show consistent decline in seedling survival but eventual steady-state-ish (lower than original; ~0.75 -> ~0.2)
- Decrease in range sizes occurs more quickly when population size is reduced (isn't this obvious? is order backwards?)
	- variation in population size correlated (?) with climate specialization and range fragmentation
	- range loss higher where there is larger temporal shift relative to traits (faster population decline)
	- declines associated with reduced strength of selection through drift and more demographic stochasticity
		- extinction vortex!
		- "quick" ability to rebound and adapt to local conditions for those populations that escape the vortex
- Effects of climate change is modulated by strength of selection on seedling survival and baseline adult survival
	- Stronger selection on seedlings means faster decline
	- Higher adult survival means longer persistence of adults, slowing near-term rates of adaptation in species range but slower turnover
	- period of decline of relative population size, particularly adults, then plateau (at ~-0.8) after ~90 years
- Age structure: pre-change "stable" age structure is mostly reproductive adults, then seedlings, then pre-reproductive adults
	- maladaptation decreased frequency of pre-reproductive adults (this is because initially there are fewer seeds surviving selection and thus there is less recruitment)
	- Long-lived adults (not subject to selection!?) restrict recruitment through competition and increasing maladaptation to the new climate
		- because new individuals were not replacing senescing adults, adult frequency decreased as population size decreased
	- Environmental stabilization allows population size to equilibrate and frequency of pre-reproductive adults (adapted to novel environment) to increase	
- Local adaptation vs. global adaptation: simulations where individuals had one genotype (average environment in occupied patches) that did not evolve but did disperse (similar to non-evolving control?)
	- This lack of local adaptation meant faster decline
	- But, if local conditions ameliorate towards the center of the niche-space, phenotypic diversity that local adaptation maintains increases maladaptation
##### Discussion
- Long lifespan limits adaptive capacity but does allow long-term persistence in unsuitable sites
	- in this study long lifespan and limited dispersal means that less specialized and more abundant species have slower declines but very rapid loss of adaptive potential
- Results suggest that species range shifts and local dynamics might be somewhat temporally decoupled, with local demographic and genetic processes occurring more quickly than elsewhere
	- local extinction vortex, meaning that small local population sizes might cause rapid extirpation that deterministic expectations might not capture
- A model with only SENM + demography found slower decrease in range size compared to SENM-only projections
- On model assumptions: unlike in model, bioclimatic variables are often correlated and selection on phenotypic traits may also be correlated
	- (see ref 37)
	- model also has independence of the three traits but genetic correlations are likely
	- other variables might be important - the simulation approach makes it possible to incorporate these
	- Also, key assumption: selection from climate is exerted mostly on seedling survival (refs herein)
		- apparently, simulations with selection influencing adult fecundity did not have the same influence
		- of course in a demographic model of longer-lived species population growth is more sensitive to survival than recruitment rates
		- it's possible that the survival function is not gaussian (e.g., could be truncation instead)
		- community composition might be a better predictor of adult fecundity than environmental variables
		- COMPADRE may be useful for figuring out which vital rates in certain environments are most vulnerable to selection via climatic variables

Super cool stuff. Downside is no analytic tractability, but there is some good stuff in here.

Their result is caused by the fact that, as with the Schmid et al. paper, selection is only acting on one vital rate. The influence of longevity is because adult survival is not subject to selection (only seedling viability, occurring once).

Kinda cool to see the age-distribution-related analyses here, particularly the sort of ripple effects (low seedling survival meant fewer pre-reproductive and then reproductive adults, with a large portion of the population held in the seed/seedling stage).

This NEMO stuff sounds interesting!


### Lindstrom, J., and Kokko, H. 2002. Cohort effects and population dynamics. Ecology Letters.

- Cohort effects (mean differences among cohorts) often arise due to conditions during early development (e.g., Lindstrom 1999)
	- population density may also have effects
	- Cohort effects are different from maternal effects: cohort effects are typically not passed on to offspring in the same way that maternal effects are
	- theory typically suggests that individual differences provide a stabilizing effect on population dynamics
		- e.g., effects of resource monopolization: if some high quality individuals can secure enough resources in bad conditions then this may buffer populations (Lomnicki and Sedziwy 1989)
##### Model
- Discrete time population growth with overlapping generations
- Approach: compare population dynamics + cohort effects with no environmental stochasticity vs. population dynamics + cohort effects + environmental stochasticity vs. population dynamics + no cohort effects + "non-remembered" environmental fluctuations
	-  here final model (no cohort effects but fluctuations still present) are a control/alternative model such that fluctuations only are important in current time step
- Maynard Smith-Slatkin density dependence:
	- p_{i,t} is survival of an individual from year t to t+1 as a function of birth year i
	- p_{i,t} = p_0 / (1 + a_{i,t} sum_j N_{j,t})^b
	- f_{i,t} is fecundity of individual born in year i in time step t
	- f_{i,t} = f_0 / (1 + c_{i,t} sum_j p_{j,t}N_{j,t})^d
	- thus population size in t+1 is: 
		- N_{i,t+1} = p_{i,t} N_{i,t}
		- N_{t+1,t+1} = sum_i f_{i,t+1} p_{i,t} N_{i,t}
	- a_{i,t} and c_{i,t} are resp. strength of density on survival for cohort i
		- a_{i,t} = a_0 exp(-alpha q_i + (1-alpha)q_t)
		- similar expression for c
		- q_t is quality of year t; q_t ~ N(0,sigma^2)
		- so alpha scales relative importance of birth year and current year (alpha = 1 is extreme cohort effect, alpha = 0 is no cohort effect)
##### Results
- Fig. 1 shows that cohort effects can be stabilizing or destabilizing
- Fig. 2: identical deterministic effects, comparison of stochastic dynamics with/without cohort effects
	- seems like results on CV (of population size?) with/without cohort effects are similar to each other? (at least, effects are similar compared to deterministic model... wait what is CV in a deterministic model...?)
- Fig. 3 shows effects of inter-cohort variation producing different dynamics than near-term environmental stochasticity (b > 1)
	- populations of certain size can be composed of cohorts of different qualities, producing homoegeneity in vital rates (compared to the population as a whole)
		- through non-linear averaging, the average rate in the heterogeneous population will be below the mean in the homogeneous population (Jensen's ineq.)
		- variation producing "shallower" density dependence curves...?
	- also, small population size tends to result from poor past conditions nad large population size tends to result from good years
		- with "memory" of past conditions (i.e., cohort effects) then population size will correlate with individual quality (e.g., good years -> large populations of good individuals)
		- this will strengthen NDD after "good" years and weaken it after "bad" years, producing a (stabilizing?) effect
		- Fig. 2 also shows thish: high p_0 (high baseline survival) means adding stochasticity has a greater dampening effect than when p_0 is low
		- high environmental variability will produce the largest difference between environmental variations-only and cohort effects case
- If original density dependence is steep enough to produce a deterministic fluctuation, then shallowing-effects due to NDD will produce more stable dynamics
	- but if the original slope of the density dependence is shallow and resulting dynamics are stable, then cohort effects that produce a shallower relationship will lessen the stability of the system (less likely to return to equilibrium if perturbed)
		- e.g., Fig. 4: shallower density dependence means slower recovery from fluctuations
##### Discussion
- Cohort effects and individual variation will increase variation in fluctuations when underlying dynamics are stable but can stabilize population size when dynamics are unstable
	- increasing fluctuations when: fluctuations producing cohort effects introduce variation to populations similarly to how environmental fluctuations influence individuals (regardless of birth year)
	- stabilizing effects when: cohort effects maintain individual-level variation that produces buffering
	- in nature: vertebrates tend to low growth rates and stable dynamics so cohort effects may be de-stabilizing
		- inverts though have shorter life spans, less overlap in generations, thus unlikely to have cohort effects at all
	- stabilizing effects more likely when there are strong non-linearities in density dependence, high potential growth rates
		- e.g., Soay sheep apparently have cohort effects? (Coltman et al.. 1999, Forchhammer et al. 2001) and have strong density dependence and inherently unstable dynamics
		- prediction from Soay that this model/analysis suggests: sequence of similar years should produce more pronounced flutuations than under stable environments (not sure why... because the quality and population size will be correlated?)
- This model is conservative: all individuals in a cohort are assumed equal, but more structured variation will also have effects on population growth

Interesting stuff. Interesting idea here is that there are possible positive effects of within-population variability in vital rates and cohort effects are one of several ways to produce or maintain them. 

I like the intuition that Jensen's inequality produces a flattening effect (Fig. 3).

I also like the idea of cohort effects producing correlations between population size and individual quality - seems like potential for de-stabilizing dynamics. I suppose with Markovian environmental change, these might be uncorrelated without cohort effects, but with cohort effects there is a kind of temporal autocorrelation (not sure if that is the right term - "memory" in the intro was a good way to put this).

The stabilizing vs. de-stabilizing stuff... seems like the point is that it will somewhat counteract whatever the underlying dynamics are. When are underlying dynamics unstable? Not sure when these would be; feel like I'd mostly be interested in stable dynamics (low growth rates) in which case the cohort effects might be de-stabilizing...


### Purvis, A., Gittleman, J.L., Cowlishaw, G., and Mace, G.M. 2000. Prediction extinction risk in declining species. Proc R Soc B.

- Predictions about attributes that correlate with extinction vulnerability:
	- Small populations (due to stochasticity, local catastrophe, slower adaptation, "mutational meltdown", inbreeding)
	- Island endemics (small geographic range means typically small population and typically no associated LH enemies)
	- Higher trophic level species (magnified effects of disturbance)
	- "Slow" LH such as small litter, slow growth rates, late maturity ("less able to compensate for increased mortality with increased population extinction" - cites MacArthur & Wilson 1967)
	- Complex social structures (persistence is at higher level than individual)
	- Large home ranges (more susceptibility to habitat loss and fragmentation)
	- Large body size (typically correlated with above characteristics)
- Variables selected against current extinction risk assessments
	- mammals and carnivores
	- assessments come from 1996 IUCN list, treating extinction risk as a scaled measure 0-5
		- these come from criteria listed above??
		- analyses excluded species referenced in ref to these characteristics
	- Also performed analysis on subset of recently-declining species
- Phylogenetic correlation incorporated into error structure... somehow?
	- "phylogenetically independent contrasts" - what does this mean?
- Looked at declining species looking at each predictor independently
	- then compared with multiple regression of contrasts
		- (this allowed factoring out (?) geographic range which was highly significant)
	- also looked at multiple regression with model simplification to find minimally adequate models
		- removed predictor with lowest marginal reduction in variance
- Small geographic range and island endemicity have greatest predictive importance
	- less important but still supported: slow life history (gestation length only?), low population density, diurnal activity
	- incorporating geographic size removes significance of island status, but density, gestation length, and body mass still appear positive
- Higher extinction risk assessments correlate with small geographic range, high trophic level, low population in minimally adequate models
- Models of primates and carnivores separately yield similar results to each other, as do separating out declining species only	
	- in these models, minimally dequate model has significant effect of gestation length (positive) for carnivores but not primates
	- likewise, age at sexual maturity is significant in the carnivore model too (negative)
- Small range, low density, high trophic level and low reproductive rates typically associated with higher extinction risk
	- these appear to override or outweigh effects of contemporary population size
- Island endemics less at risk than expected by theory (after accounting for range size), agreeing with two other recent studies
- It is possible that rates of range decline are underrepresented in widespread species compared to narrow ranges
	- but also range size might correlate to habitat specificity, which likely hinders adaptive potential in face of environmental change
- Model "misses" where model seriously underestimates threats occur in areas with unusually high tropical forest loss
	- species might be robust to change but change is simply too much

Interesting. I don't know how the phylogenetic contrasts are modeled (it isn't really explained). Other possible methodological issues: (1) treating response as linear, (2) data source for DV, (3) not reporting collinearity among predictors. Taking results at face value, this does suggest that slower LH populations are classified as higher risk. The specific proxy supporting this is gestation time, which correlates with several other pace-of-life factors. Notably age at first reproduction does not seem to be significantly correlated except in one case (also, holy multiple testing Batman).

Cite or no? I guess for showing how this is a widely supported hypothesis. Perhaps not pointing to as a result in and of itself.


### Millner, J.M., et al. 1999. Repeated selection on morphometric traits in the Soay sheep on St. Kilda. Jouranl of Animal Ecology

- Monitoring sequential selection events: som examples of natural selection of quantitative traits in birds (refs herein)
	- although phenotypic selection may produce little/no change if the variance is environmental rather than genetic
	- in soay sheep (St. Kilda), 1985 - 1996 measurements of winter mortality in relationship to quantitative traits
		- previous studies have shown selection on morphometric characters
##### Methods
- Population fluctuations (600 - 1825 individuals over ~30 years), where crashes may be due to starvation exacerbated by gastro parasites
- Morphomeasurements (annually): body weight, hindleg length, incisor arcade breadth (bite size, food intake rate)
	- correlations among variables
- Estimate selection differentials (change in population  means before/after selection divided by variance before selection)
	- use t-tests to estimate character means of survivors and non-survivors
	- using selection gradients (multiple regressions) to figure out which traits featured largest selection
		- 1990 - 1996, estimating selection gradients each year
	- estimate heritabilities as well
##### Results
- Selection for larger body size w/ stat. sign. differences in body size between living and dead sheep in 4/6 years
	- significant selection on hindleg length in same years
	- only two years of significant selection in females, one in males, for incisor arcade breadth
- Between-year variation: density dependent s/t not discernible at low density but discernible differences in all three traits at higher density
	- (but also note easier to detect these things at high density - larger sample size) 
	- survival in age was correlated with population size in lambs but weaker in adults
	- incisor breadth - negative direction of selection in most years
- Interestingly, no evidence from a regression (hmm...) of body weight over time despite evidence of consistent selection for larger size
	- heritability of body weight was significant @ 0.05 but low
	- higher heritabilities of other traits
##### Discussion
- Ah interesting... body weight is body size + condition, and hindleg length should correlate with body size but not necessarily with condition
	- so the lack of response of hindleg length suggests that trends are driven in part by condition
	- (connection here between heritability and "repeatability" - interesting point!)
	- it's possible that the repeated selection on body weight has depleted genetic variation for the trait, causing low heritability
	- or, there may be low response to selection because of opposing selection occurring at other points in the life cycle e.g., fecundity

Hmm not useful in the way I was hoping (not looking at the consequences of the *repeated* aspect of selection) but it is still a useful result. Might not be high longevity here though so effects might not be detectable. 

There is some other interesting stuff though - e.g, lack of detectable response due to low heritable variation, or the depletion of additive genjetic variance for body size, some notes about body size as a trait and selective pressures upon it.

### Gadgil, M., and Bossert, W.H. 1970. Life historical consequences of natural selection. Am Nat.

(Took notes on paper copy)

Confusing at times. 
Model had functions for consequences of reproductive investment on survival (decreasing w/ investment), growth (decreasing w/ investment), and reprodcution. Assume that in the final stage, reproductive investment is maximized. 

Main predictions are what is important.

- "Profit" and "cost" of reproduction to lifetime fitness at different stages; optimal reproductive investment at a given age when profit minus cost is maximized. 
	- Semelparity is favored if the optimal investment is 0/1 for all stages
	- Iteroparity is favored if the optimal investment is at intermediate values for several life stages.
		- This will happen at if the profit function is concave and the cost function is convex (phrased differently in this paper but this *has* to be what they mean)
- Reproductive effort should increase with age in the case of iteroparity
- If there is an environmental effect with equal investments on all ages after a certain age (i.e., for all ages j > j'), this will lower age of reproduction for semelparous individuals and increase investment prior to that age for iteroparous individuals
- Soething about reproductive potential relationships with size... had a hard time understanding this one
- Population regulation: predation-limiting populations (influencing mortality) and resource-limiting populations (influence mortality and reproduction)
	- Uniform probabilities in survival odds (associated with predation-limitation) will not influence reproductive effort in either semelparous or iteroparous populations
	- If equilibrium density declines (due to uniform influence on survival), relaxation of resource limitation should lower age of reproduction for semelparous individuals and increase reproductive effort in iteroparous individuals

Hmm... useful at the very least for making points that reproductive effort should change with age (not present in our model). We also do not include size (explicitly influencing future fitness) or growth effects. 


### Vaupel, J.W., Manton, K.G., and Stallard, E. 1979. The impact of heterogeneity in individual frailty on the dynamics of mortality. Demography.

- Define mortality measure mu(x,y,z) is mortality of individual of age x at time y with 'frailty' z
	- frailty is constant over lifespan(z > 0)
	- wlog let mu(x,y,1) be mortality of "standard" individual
	- mu(x,y,z) = z * mu(x,y,1)
		- can express this as mu(z) = z * mu
- H(x,y,z) is cumulative mortality hazard of individual born at age y-x, up until age x (from age 0 to x, from time y-x to y)
	- H(x,y,z) = int_0^x mu(t,y-x+t,z) dt
		- can express this as H(z) = z * H
	- survival s(x,y,z), probability of surviving to age x, is exp(-H)
		- if we let s(x,y,1) = s, then s(z) = s^z
- Mortality within a cohort
	- mean is \bar{mu}(x,y) = int_0^infty mu(x,y,z) p(z) dz
	- mean frailty is \bar{z}(x,y) = int_0^infty z p(z) dz
	- from this it follows that \bar{mu}(x,y) = mu(x,y,1) \bar{z}
		- i.e., \bar{mu} = mu * \bar{z}
- Assume frailties within a cohort are gamma distributed upon birth
	- gamma has parameters k (shape) and lambda (scale)
	- k = 1 gives exponential, k growing large approaches normal
	- at age x, frailty of survivors will remain gamma with same k, now lambda(x) = lambda + H(x)
	- mean frailty is \bar{z}(x) = \bar{z}(0) k / (\bar{z}(0) + k + H(x))
		- so mean frailty tends to zero with time
	- frailty of those dying @ age x is also gamma distributed, same labmda, shape is now k+1

There is other stuff in here too applying this to some datasets. I skimmed only.

Useful! I can't help but feel like there's a link between gamma distributions in frailty and normal distributions in traits under gaussian selection. The results just look too clean - just look at that mean frailty at age x expression!


### Partridge, L., and Harvey, P.H. 1988. The ecological context of life history evolution. Science.

- Consensus is that under simplifying assumptions, optimal life history maximizes r (intrinsic rate of population growth for a population at stable stage distribution)
	- this parameter is determined by age-specific survival and fertility
	- (it is missing frequency dependence, tho)
- Iteroparity and Cole's paradox
	- Letting parental survival be P and offspring survival be Y, then rather than one additional offspring (as in the case of Cole's paradox) the cost of iteroparity is actually P/Y
		- (so, higher adult to juvenile survival is expected to lead to perennial life histories)
	- Law (Am Nat, 1979) grew Poa annua from two different environments in a CG
		- one population was from low-density environment with high mortality
		- one population was from pasture with high density and high juvenile mortality
		- in common conditions, pasture population were more likely to survive and breed again
- Reproductive effort: theory here relies on ideas about reproductive cost
	- theory suggests that individuals should maximize their reproductive value (sensu Fisher) at each age through appropriate resource allocation
		- e.g., forgo reproduction at younger age if it allows for increased survival or growth to (st)age with more reproductive capacity
	- reproductive costs are difficult to measure, making advance there difficult
		- observational studies may have confounding (environment, phenotypic effects, etc.)
	- experiments also have issue with genetic vs. purely-phenotypic manipulations of reproductive rate
		- intersting stuff about Drosophila manipulations... genetic correlations between fertility and longevity appears negative but phenotypic correlation is often positive du to effects of body size
	- some stuff on single species studies, comparative studies...
- Temporal variability: an allele or genotype with zero fitness in *some* environment will go extinct in enough time (bet hedging, phenotypic plasticity)
- Constraints on life history evolution
	- Mutational variance: Stearns suggesting that fixed LH traits in a lineage while closely related lineages in similar habitats do not have the trait is a result of shortage of mutations?
	- Body size may limit LH evolution in some ways, although other LH components often correlate more closely with each other than they do with body size (e.g., gestational length and age at weaning)
		- size may not necessarily be a fixed constraint determining LH differences and may instead be evolving in response to LH selection pressures
	- age, selection intensity, mutation accumulation
		- Williams (1957) - expected reduction (?) in ferrtility with advancing age?
		- detected reduction in fertiltiy with age in natural populations of long-lived species?
		- reproductive senesence may be a constraint in maintenance/repair
		- Medawar: evolution of aging is due to selection acting on mutations (or any allele?) with age-dependent effects on fecundity or survival
			- from this we can make the prediction that mutation (positive or negative) with effects early in life will be subject to stronger selection than mutations with affects later in life
- Population density can influence how growth rates can be maximized (r- vs. k-selection)
	- investment for high reproduction allows for high growth at low densities (e.g., when colonizing new environments) - r selection
	- cf investment for competitive ability when frequently at high density - K selection
	- K-selection may maximize numbers in the "critical age group" where DD occurs (Charlesworth book)
	- predictions:
		- can't simultaneously have high r and K
		- growth rate patterns should reflect density
			- see experiments on Drosophila that seem to confirm this
	- drawback: selection at high density favors improved competitive ability, which reduces popualtion size (ref)

Great overview! Won't cite here but good for intro to some interesting questions.

### Coulson, T., and Tuljapurkar, S. 2008. The dynamics of a quantitative trait in an age-structured population living in a variable environment. Am Nat.

- Probability of survival and successful reproduction can be reliatn on a phenotype, incl. quantitative traits
	- mean demographic rates in a population may depend on trait distribution, trait distribution and rates of change may depend on underlying demography, etc.
	- "demographic rates can therefore be considered as phenotype-by-environment interactions (Coulson et al. 2006)"
	- viability selection can alter trait distributions within an aging cohort, and fertility selection upon traits with heritable variation can alter trait distributions over time
		- (these strengths can also vary with environment...)
	- "it should be possible to write an expression to exactly describe change in the mean value of a quantitative trait in an age-structured population over a time step" (similarly to the way MPMs work)
		- benefit of this would be to look for retrospective decomposition of observed fluctuations in trait values to contributions from survival, reproduction, and other factors like density dependence, climate
- What can cause (rapid) phenotypic change within a population?
	- viability selection removing individuals as they age
	- fertility selection allowing genotypes to differentially contribute to successive cohorts
	- genetic, maternal, and environmental affects on offspring size (and presumably also rates of growth?)
	- phenotypic plasticity within the lifespan
	- fluctuations in the demographic composition of the population
- Previous works on quant trait dynamics (assuming non-overlapping generations)
	- Breeder's equation: Delta(\bar{z}) = (Cov(Z,W) / \bar{W}) (V_A / V_P)
		- selection differential times heritability
		- estimating the mean and covariance of fitness, trait across the population
		- covariance term gives the difference in the mean of the distribution of parents (weighted by number of offspring per parent) and mean of population before selection (?)
	- Assume that the mean of the phenotypic distribution among parents (after selection) and the mean of offspring is the same [n.b. we do the same]
	- some other assumptions (see Bulmer 1980) e.g., independence between breeding values and environmental components, traits in parents and offspring are normally distributed
	- Lande (1982) and Charlesworth (1993) derived expressions to handle how traits will change with age (assuming weak selection, stable age structure)
- Intracohort viability (e.g., Vaupel's frailty) can be combined, perhaps, with Price's equation: \Delta(\bar{z}) = \bar{Z}(a+1) = \bar{Z}(a) = Cov(Z,S) / \bar{S}(a)
	- Z is trait (frailty) and a is age, \bar{S}(a) is a mean survival of individuals at age a
	- Vaupel's frailty is assumed to be constant through an individual's life, although we do know that phenotypic plasticity can occur to change an individual's trait over time, or that trait values may change with age, environment, physiological status
		- approaches will need to handle this!
	- likewise, fluctuations in the demographic structure of a population are important!
- Price's (1970) equation: change in mean value in a trait between generations is the sum of (1) differences in mean trait values selected to be parents and those of whole population and (2) differences in mean trait values between parents and offspring
	- traits can be scalar or vector, fitness W is scalar
	- fitness is defined as lifetime offspring production (assuming non-overlapping gneerations)
		- in an annual setting, let W = R, the number of offspring of an individual
	- Say n individuals in year t, individual i with trait z_i produces r_i offspring
		- \bar{Z}(t) = (1/n) sum_i z_i
		- \bar{R}(t) = (1/n) sum_i r_i (offspring numbers)
	- Let trait value of jth individual offspring of i be y_ij = z_i + d_ij
		- (so, mean of parent(s), plus some deviation unique to j)
	- \bar{Z}(t+1) = (1 / (\bar{R}n)) sum_i=1^n sum_j=1^ri (z_i + d_ij)
		- (1 / (\bar{R}n) sum_i=1^n sum_j=1^ri zi term is just parental phenotypes weighted by number of offspring, which authors claim is equal to \bar{Z}(t) + Cov(Z,R) / \bar{R}
			- (n.b., in our model cov(Z,R) = 0, so this is just \bar{Z}(t))
		- second term is "infidelity in transmission of mean phenotype from aprents to offspring", which is usually assumed to be zero in expectation
			- (1/\bar{R}n) sum_i=1^n r \bar{z}_i = \bar{D} + cov(D,R)/\bar{R}
			- let \bar{d}_i be difference in mean trait value of individual i's offspring and parental value z_i, then averaging over all parents (and ignoring or setting \bar{d}_i = 0 for non-parents - this shouldn't affect expectation but it should affect variance yes?)
	- Combining these terms gives \bar{Z} = Cov(Z,R)/\bar{R} + \bar{D} + Cov(D,R)/\bar{R}
	- Def n_p <= n to be number of parents, let \bar{R}_+ be mean number of offspring *among parents* (i.e., excluding zeros)
		- then, \bar{Z}(t+1) = (1 / \bar{R}_+ n_p) sum_i=1^n_p sum_j=1^r_i (z_i+ d_ij)
			- (I think this works because non-overlapping generations means that non-reproducing parents have no contribution to successive generations?)
		- in above expression, first term is \bar{Z}_+ + Cov_+(Z,R) / \bar{R}_+
			- (but because of non-overlapping generations?) this is equal to above? same for second term)
	- We can combine all of this to rewrite Price's equation as follows:
		- Delta(\bar{z}) = (\bar{Z}_+ - \bar{Z}) + Cov_+(Z,R) / \bar{R}_+ + (\bar{D}_+ + Cov_+(D,R)/\bar{R}_+)
			- term 1 is selection of parents among individuals, middle term is effect of ferility selection among parents, final term is difference in mean trait between parents and offspring
	- the breeder's equation will produce the same results as Price's equation under some restrictive circumstances
- Quite impressive decompositions in expressions 11, 12
	- looking at some specific cases, e.g., no selection
	- interesting one is second example: no fertility selection, but viability selection, assuming that demographic structure doesn't change, demonstrating that \bar{Z}(t+1) is sum of
		- (1) age-weighted sum of current trait values, weighted by proportion of future time step?
		- (2) fertility-related contributions of each cohort to the next time step
		- (3) changes due to viability selection, in which case non-zero Cov(Z,S)  for each/any ages will produce a change in the mean population phenotype
- Application to red deer in Scotland: rut in S/O, births in M/J, mortality predominantly in F/M/A
	- trait is weight: residual birth weight is residual from regerssion line of weight and capture date (at time of birth)
		- individuals followed through time, allowing complete reconstruction of LHs
		- birth weight is fixed throughout life, in which case some terms in the age-structured Price equation are zero (plasticity terms, \bar{G}_+(a,t) = 0)
		- maximum litter size is 1, so covariance terms Cov_+(D,R)(a,t) are zero (not entirely sure why...)
	- viability selection differentials appear to be positive, mean of 42g/year (?)
	- fertility selection differentials are positive but lower - mean of ~10g/year
	- recruitment contributes ~+26g/year to change in mean population trait
	- somehow, "approximately 70% of selection occurs through survival"
	- fluctuations in demographic structure only explain ~2g increase in birth weight
	- in total there is a predicted average increase in 70g/year (~1% increase per year in mean calf birth weights)
		- but this is not achieved because \bar{D}(a,t) is often negative? female calf weights are less than birth weights of their mothers?? (OH SHIT - see note below)
	- Fig. 3a - evidence of viability selection. Fig 3b - this seems important... difference in maternal and offspring birth weights as a funciton of age is unimodal (?!)
##### Discussion
- This approach allows us to decompose observed changes in a quantitative trait within a population (which has been observed on rapid timescales!) into contributions of viability selection, fertility selection, phenotypic plasticity, mean difference between offspring and parental values, and fluctuations in demographic structure
	- \bar{D}_+(a,t) term, difference in mean trait value between parents and offspring, can arise through several factors, e.g., additive genetic variance (?), environmental effects, and resource availability
- Viability selection had a larger effect than fertility selection
	- Clutton-Brock (1988) cross a variety of vertebrates found that longevity was important in explaining variation in lifetime reproductive success
	- however, viability selection is somewhat offset by difference of mean parental and offspring trait values
		- "parents tend to produce offspring that are more similar to the pouplation mean than they are to themselves"
		- "the observation that this means parental birth weights are, on average, greater than offspring birth weights has not been preciously described. The reasons for the difference are not clear but could be the result of developmental or energetic constraints (citations), which may or may not have additive genetic components."
			- (some stuff in here about what is meant by AGV?)
- Future work that would be cool:
	- extinding work to higher moments of trait distributions (Rice 2004 may have already started progress on this)
	- would be cool as well to incorporate how these terms are influenced by environmental variables


AHH holy shit this rules. Great framing here and also some unexpected results here that our model provides insight into!

Love the dichotomy they provide between viability and fertility selection, and the fact that the two may have varying influences on adaptation. That is what longevity controls - the relative effects of these two!

The negative \bar{D}(a,t) term (lower offspring weight than parental weight) is the effect of repeated viability selection - within a cohort, the phenotypic value outpaces the underlying genes under directional selection. We don't see the same unimodal pattern they see (where only "prime age" mothers produce offspring close to the population-mean birth weight), but they do see a divergence happening growing beyond a certain age. This is likely the effect that we see, and simplifying assumptions in our model about growth, environment, etc. produce the pattern they see at lower ages.

It would also be cool to (1) try to see if we can reconcile our math with theirs and (2) apply this decomposition to some simulation results!

Should read the Clutton-Brock 1988 next! (aw shit it's a book...)


### Jones, O.R., et al. 2014. Diversity of ageing across the tree of life. Nature.

- Looking at 46 species (11 mammals, 12 non-mammalian vertebrates, 10 inverts, 12 vascular plants, green algae), trying to find variation among species in age-specific reproduction and mortality
	- standardized demographic trajectires, starting with mean age of reproductive maturity, ending at terminal stage where ~5% of adults are alive (beyond this point sample sizes are very small)
	- very large variations in trajectories!
- Fig. 1 is the money fig: variation in shape, etc. of fertility over time, survival/mortality over time
	- Mortality: there is both increase in mortality with age (e.g., humans) and mortality declining with age (e.g., desert tortoise) up to temrinal age
		- seems though like most common pattern is increasing mortality with age, even if this is non-monotonic (e.g., bathtub-shape)
	- Fertility: typically "bell-shaped" although the width of this can vary considerably (i.e., fertility spread out)
		- some species, not just humans, have post-reproductive lifespans - perhaps this is more common than originally thought?
- Senescence: most species can be put on a spectrum of senescence
	- Fast-slow continuum has been proposed for arranging species, but is does not produce any obvious pattern in Fig. 1
	- perhaps it is more interesting to define senescence by its abruptness (low to high mortality in a short time) rather than simply the presence of an increase
	- one interesting possible measure of the "shape" of a senescence curve: ratio of mortality at terminal age to average level fo adult mortality (this is nice - it is time-invariant)
		- using this measure, perhaps a relationship between longevity and degree of senescence may emerge... but when testing, there is no relationship between longevity and length of life and senescence-abruptness
		- "Hence the data support Baudisch's conjecture that pace and shape may be two orthogonal axes of life histories" - Baudisch ref is 2011, Mol. Ecol. Evol.
	- Survivorship curves (prop. of population surviving to a given age) - if mortality increases with age, the log-survival curve is concave while decreasing mortality will have convex log-survival curves
- Some phylogenetic signal although not super strong

Very cool! Uses for us are (1) lots of diversity in vital rates over the lifespan and (2) these do not seem to map very neatly onto longevity. Useful for saying that our simple model is an adequate starting point for looking at the effects of longevity, in particular noting that while there is variation in LH patterns they do not appear to be intrinsic to ageing.


### Miller, DA., et al. Stochastic population dynamics in populations of western terrestrial garter snakes with divergent life histories. Ecology.

- We can perform ecological inference on variation life history strategies using (among other things), comparative studies within and among species to see if there are differences in mean demographic parameters)
	- e.g., how prediation can be a top-down agent of selection while resource variability can be a bottom-up force of selection
	- theory suggests that temporal variation in vital rates/dynamics should impose some selection pressures on LH/demographic traits and should help explain patterns in LH
- Here: testing some predictions about temporal variation and its influences on life history evolution by looking at (distinct) populations of a garter snake
	- Study populations in N. Cali originated from an ancestral source population that diverged into genetically diverged ecotypes that have maintained distinct morphologies and LH strategies
	- Lakeshore (L-fast) ecotype has fast growth, early maturity, large litter sizes
	- Meadow (M-slow) ecotype has slow growth, late maturity, less frequent breeding
	- CG work (Bronikowski 2000) suggests that these differences among populations have a genetic basis
##### Methods
- Six populations in N. California, three L-fast and three M-slow
	- variation in representation, temporal extent of population in each dataset
- Diet: Anurans (frogs, toads), annually-varying portion of the snake diet and a proxy for environmental conditions (form a high proportion of diet in years with little food)
	- anurans are a greater proportion of th diet in meadows (M-slow), where fish are not available
	- for L-fast populations, diet is fish and leech heavy and more consistent among years
	- Prediction: presence of breeding anurans, which varies across years, will have a larger effect in the M-slow populations (due to the presence of alternative food sources in the L-fast populations)
	- years classified as A-present or A-absent
- Matrix models: life cycle graphs made into Lefkovitch matrices, estimated five vital rates for fitting models, made matrix model for each population (in each obs. year?), then looked at differences in LH among populations in relation to environmental mean fitness and mean stochastic population growth rate
	- L-fast population has three stages, M-slow has five
		- age of reproduction in L-fast is three, five in M-slow (but can remain in this class for several years)
	- for each population, got annual matrices, mean annual matrix
	- estimate deterministic growth rates from A and deterministic sensitivities and elasticities from from \tilde{A}
	- also estimated the sensitivities of lambda_s to standard deviation of annual (co-)variation in (among) vital rates using Doak et al.'s (2005) formulation
- Test for effects of buffering: there should be an inverse relationship between E_v (deterministic elasticity of vital rate v) and CV of a vital rate v (from Pfister 1998 and others)
	- scales of analysis: E_v vs CV for all vital rates from all populations at same time (overall), E_v vs CV for all vital rates within one population (within-pop), and E_v vs CV for all populations for one vital rate (within-rate)
		- prior analyses have focused on overall and within-populationships, but here the same vitals are measured across populations
		- within-rate relationships tell us if buffering is occurring in a comparative context (reasons given, not read)
- LTRE gives retrospective analysis
	- Decompose contribution of annual variation in vital rates to variance in lambda and variance in lambda explained by anurans
	- letting C_ij = S_vi S_vj cov(vi,vj), where S_vk is deterministic sensitivity to vital rate k, then variation in lambda is approximated by the sum of C_ij over all ij
	- also used SLTRE for something...
	- something about another analysis where, for a given rate, populations with observed v_it were compared to hypothetical populations where v_it was set to values where breeding anurans were present (or not? one or the other) (seems like the RIE effects from the Louthan et al. paper that I borrowed)
##### Results
- Vital rate differences among populations:
	- Proportion of gravid females and litter size were higher in L-fast than M-slow
	- Higher survival in M-slow than L-fast
	- in M-slow populations, juvenile survival was as high or higher than adult survival (limited cost for delayed age of first reproduction)
		- in L-fast populations, there was more variation in survival with age? (maybe I am misunderstanding?)
		- Fig. 2: generally increasing survival with age in L-fast populations, more stable-seeming in M-slow populations, with main difference being at neonate stage where L-fast populations have lower survival but they seem more similar at adult stages
- Varying population growth rates, with more variation among L-fast populations
	- among populations, pretty similar deterministic and stochastic growth rates
	- no consistent differences in effects of stochastic variation on labmda_s between the two ecotypes
- Deterministic perturbation analysis suggests that fecundity is more "influential" in L-fast populations and survival more "influential" in M-slow
	- in L-fast populations, vital rates were similar in their deterministic elasticities, M-slow populations typically had more variation due to the fact that some of the survival rates had substantially higher elasticity
	- longer juvenile period in meadow populations (three vs. one year) means greater elasticity for juvenile survival in M-slow populations (compared to juv. survival in L-fast)
	- similar adult survival elasticities in both cases, while certain rates (percent gravid, neonate survival... LS, what is this, litter size?) were higher in L-fast populations than M-slow
- Evidence supporting buffering:
	- overall negative correlation between E_v and CV for overall tests, suggesting no cases where v varied a lot and was highly influential in producing variance in lambda
	- did not find support for differences in correlations among populations; pooling together all populations produced a significant negative correlations among CV and E_v
	- similarly, no consistent differences between E_v and CV among vital rates; pooling all vital rates together produced a very strongly-negative correlation between the two
	- four of five vital rates showed a negative relationship when comparing CV/E_v across the two environments, such that increasing CV decreased E_v (lone exception is survival of neonates)
		- P_G has much higher variation than other vital rates, and had nearly double the CV in M-slow populations than L-fast ones
		- percent of gravid females, then, seems to be important for differentiating the life cycles here?
		- variation in percent gravidity explained 75-98% of the annual variation in lambda between groups! effects of variance was higher in M-slow than L-fast
	- small effects of covariance terms
- SLTRE also suggests that PG had large contribution to annual variation in lambda, with greater contributions of PG-variance to growth rates in M-slow populations
	- this is consistent with greater fitness-related costs of prey availability in M-slow populations, as M-slow pouplations may be more susceptible to variability in prey in determining the proportion of the population that may be gravid
##### Discussion
- Gravid proportion seemed to have a disproportionate effect in both ecotypes, with a larger effect in M-slow
	- variability in prey availability "synchronizes reproduction" (hmm...)
	- annual variance in proportion of gravid females had hte largest effect on growth rates, as measured by stochastic sensitivities and elasticities... despite having high CV/low deterministic elasticities??
	- there was more variability in gravidity in M-slow populations than L-fast, and seemed to have larger effects (according to random LTREs) in M-slow than L-fast
- Novel support for Pfister's buffering hypothesis:
	- LH differences among populations of a single species is consistent with buffering and longevity
	- also, this holds at the level of individual vital rates (4/5 vital rates showed decreasing E_v when increasing CV of the rate, moving from one ecotype to the other)
	- inverse relationship between elasticity and CV holds within ecotypes, and moving between ecotypes preserves this pattern in expected ways
		- average CVs of percent gravid females in M-slow populations are doubled compared to L-fast, while M-slow has half of the elastiicty for proportion of gravid females
		- but even with buffering, absolute effect of variability in percent gravid females does have a larger effect on growth rates in M-slow compared to L-fast (although if L-fast had the same level fo annual variability, there would be larger effects in L-fast)
- This study demonstrates how top-down and bottom-up selective effects can influence life history evolution
	- high mortality across all ages in L-fast population suggests high predation produces more extrinsic mortality here, likely due to avian predators
	- buttom-up processes: present study demonstrates prey availability and its influence on LH evolution through variation in the proportion of gravid females in a population


This is cool, although not relevant as I was hoping. Is great for thinking about vital rate buffering and LH evolution. PG can vary substantially from year to year; in a population where PG varies more (seemingly due to the variable amount of food), it has lower determinstic elasticity (but still a larger contribution to variation in growth rates? kinda makes sense but still a potential discrepancy here to think about). But beyond deterministic elasticities, what exactly is the life history evolution? Improved survival and delayed reproduction?


### Dirzo, R., et al. 2014. Defaunation in the Anthropocene. Science.

- Here: patterns and consequences of contemporary anthropogenic impact on terrestrial animals
	- "anthropogenic defaunation" - similar to deforestation although difficult to deterct (unlike advances due to remote sensing)
	- some interesting refs/statements herein:
		- declines in individuals in a population still have large consequences, influencing other species and ecosystem functioning (refs 8, 10)
		- abundance declines to "functionally extinct" levels can happen quickly (2, 12)
- Extinctions: annually losing 10K-50K animal species per year (of 5-9 mil)
	- ofc. does not account for population extirpations, non-extinct rapid declines
	- vertebrate data indicate mean decline of ~28% in number of individuals across species across recent decades
	- less attention in inverts than vertebrates, although it is likely they are at least as severe as vertebrate declines
		- <1% of described invert species have been assessed by IUCN, with ~40% considered threatened
		- e.g., moths in UK - "substantial portion" of moths have experienced range contraction recently!
		- non-leps seem to have more decline than leps
- Patterns: certain lineages appear to be more susceptible to others, e.g., amphibians more threatened, birds less threatened, mammals and reptiles intermediate threat
	- geographic distributions: per 10K km^2, few to 75 (mammals), 125 (birds) listed as in decline due by IUCN
	- more decline in tropics, even considering the greater sp. richness in tropics
	- statistical models try to find predictors or correlates of extinction status/risk
		- refs 31 - 34 (Cardillo and Meijaard, TREE 2012; Davidson et al., PNAS 2009; Ockinger et al. Ecol Let 2010; Lee et al. Phil B 2011)
		- correlates include low reproductive rates
		- these analyses appear to be at the global scale rather than the population scale, though, possibly producing weaker (measured) correlations
			- e.g., body size is a better predictor of extinction risk in island birds than for mainland birds
	- body mass-risk relationsihps associated with extinction risk:
		- Fig. 3: anthropocene extinct animals typically larger than anthro-threatened, which are larger than anthro-nonthreatened
			- however, species that went extinct in the pleistocene were larger than those going extinct now
- Impacts on ecosystem services:
	- pollinators mostly in decline in abundance and diversity
	- pest control: small vertebrate declines can cause multitrophic cascades with influences on herbivores, plant populations
	- nutrient cycling: mobile species can move nutrients long distances, some interesting stuff about pleistocene extinctions changing influx of phosphorus
	- water quality: loss of amphibians increases algae concentration, detritus biomass, reducing nitrogen uptake, while loss of larger animals may lead to formation of more anoxic zones
- Evolution: more rapid evolutionary changes for short-lived organisms (ref 72: Palumbi, Science, 2001)
	- defaunation typically selects for body size
		- however, smaller-sized species may not be able to (quickly) fill niche or functional roles of larger animals, possibly causing cascades


Useful stuff. At the very least, pointing out that there are a number of animals, particularly large ones (which probably will have delayed reproduction and/or longevity), face considerable risk in extinction at species or population level.

Cool (morbid) ideas in here about selection for smaller body size and cascading effecst there. Also a mention (I didn't take a note on it) for pollinator extinctions possibly causing rapid evolution in plant mating systems.


### Cardillo, M. 2003. Biological determinants of extinction risk: why are smaller species less vulnerable? Animal Conservation.

- Growing body of literature suggesting that body size influences extinction risk, even after accounting for phylogeny, although typically these are univariate analyses that don't provide additional insight into why or how this relationship may exist
	- e.g., size may directly determine a species vulnerability, e.g., by being less conspicuous to predators
	- e.g., body size may be a surrogate for other traits (LH, ecological) which affect vulnerability
		- e.g., smaller species may have higher reproductive output and higher densities, although this may come at cost to mobility and energetic efficiency
		- of course these traits may be environment- or context-dependent
- Here: looking at two factors that correlate with body size and see their effects on extinction risk
	- reproductive output: should generate faster population growth, faster recovery from disturbances (Pimm, Jones, and Diamond, 1988)
	- home range: small home range may require less energy, higher population densities, allowing for persistence in low-quality habitats
	- analysis in terrestrial mammals of Australia, where evidence already suggests a body size-extinction risk relationship
		- in Australian mammals, there are strong correlations between these to variables and body size
	- also: Australia can be divided into different environment types (low-productivity arid zone and high-productivity mesic zone)
- Dataset: extinctions and declines in last 200 years in Australia
	- no islands, no invasives
	- extinction risk assigned, 1-5 (using some published criteria elsewhere - Strahan)
	- reproductive output: mean litter size (does not account for number of litters due to data scarcity)
	- home range size: only four extinct species with data on this, so consider caution
- Randomizations to see if extinct or extinct + endangered species were non-random subsets of the continent's fauna
	- tests of residuals for litter size ~ body size, home range size ~ body size
		- select species from the continental species pool and look at residual statistics
		- this generates null distribution of... what exactly/
- Independent contrasts derived to fit multiple regression models
	- phylo assembled through published phylogenies
	- analyses performed separately for each zone (not possible to do in one model because some species occupy both zones as do the LH contrasts)
- Extinct (5) or extinct + endangered (4+5) are non-random with respect to litter size and home range size
	- residuals of litter size on body size are low (lower than random selection of resids?), suggesting that extinct species have smaller litters than expected by chance
	- home range size residuals also larger than expected
	- Differences among zones: there is a difference, where extinction risk is higher in arid zone
		- no difference in residual patterns wrt litter size between zones
		- home range residual patterns vary by zone, with arid-zone species having larger range size than home range
- Testing for phylogenetic confounding: the two taxa with the largest numbers of extinct and endangered species have smaller litter sizes than would be predicted from their body sizes, and slightly larger home ranges than predicted
	- so, there is a possibility of phylogeny contributing to bias in thse factors among extinct and endangered species
- In multiple regression, litter size seems to be consistently significantly negatively related with extinction risk
	- body size alone (no litter size) shows positive association with extinction risk although p ~ 0.09 so marginal
	- seems like adding home range size weakens this relationship (p < 0.01 w/o home range but p ~ 0.01 with home range)
	- but, when all variables are present, litter size is only significant variable
	- w/o home range size, model explains ~12.5% of variance in extinction risk
	- effects of litter size were stronger in the mesic zone than in the arid zone, where effects were n.s.
- Strong evidence that litter size affects extinction rate, contributing to pattern of body size-extinction risk relationship
	- litter size is only one component of reproductive output, as current analysis does not account for multiple litters in a year
	- other hypotheses to explain body size relationships (not tested here): torpor? 
	- Johnson (2002) finds a similar relationship between litter size and extinction during the quaternary (pleistocene?) megafaunal extinctions
- Some possibly-interesting stuff I only skimmed about home-range size (rocky outcrops are important? neat)

Cool! Possibly limited. Limitations (1) just mammals, just australia (2) still need to think about robustness of residual analysis and (3) does not account for multiple litters. Although now that I thinka bout it (3) is just greater iteroparity is it not, so perhaps this is not so bad? 

Interesting that he attributes the low-litter size to slowed recovery to original size. This is what I am arguing. He cites a Pimm et al. paper that I think would be useful for this. That's a good one to read next (Pimm, Jones, Diamond, 1988: On the risk of extinction; Am Nat)


### Salguero-Gomez, R., et al. 2016. Fast-slow continuum and reproductive strategies structure plant life-history variation worldwide. PNAS.

- Plant kingdom has considerable LH variation, e.g., longevity (ref 12), seedbanks (13), dormancy (14), variable magnitude (15)
- Here: combine demographic, phylogenetic, ecological data from natural populations in 418 plant species
	- What are main axes of variation in plant LH strategies?
	- How do phylogenetic ancestry, habitat, growth fom, size constrain variation?
	- How do positions on axes of variation predict population growth and recovery from disturbance?
- Data: COMPADRE database with 105 families, 825 natural populations, each of which has at least four years of data
##### Methods
- Constructed a phylogenetic tree (they just did it on their own? lol)
- Derived nine LH traits:
	- Population turnover (T)
	- Longevity (H, L_a)
	- Growth (gamma, rho)
	- Reproduction (phi, S, R_0, L_omega)
	- (also, not LH score, but habitat data as well...)
- Traits log-transformed, scaled to standard normal for normality, fed into phylogenetically informed PCA
	- Looked for signf differences in PCA scores w/ 3w ANOVA w/ post hoc Tukey
##### Results
- Two PCA axes explain 55% of variation
	- Axis 1 (34% of variation) is a fast/slow axis, where an increasingly positive value implies increasing allocation to longevity-related traits, decrease in turnover rates at the expense of growth and production of new recruits
	- mix of generation time T, mean sexual reproduction number across life cycle weighted by SSD phi, individual plant growth rate gamma
		- Positive loading for T, negative loading for growth and mean sexual reproduction, supports trade-off
	- Axis also includes shape of survivorship curve H, mean age at maturity L_a
		- H is Keyfitz's entropy where sign wrt 1 gives different survival-type curves, so increasing H produces increasingly-I (human-like)-type curve
	- Axis 2 (21% of variation) has increasingly positive values for increasing frequency of reproduction and less shrinkage?
		- reproductive strategy that is independent of mean sexual reproduction: net reproductive rate R_O, degree of iteroparity S (both positively loaded)
		- rate of shrinkage (rho) also appears, negatively loaded
	- Axis 3, explaining much less variation than the rest, contains mature life expectancy and period between age of sexual maturity and mean life expectancy
		- does have eigenvalue >1 for some groupings (herbs) but not others (trees, shrubs) when procedure is repeated on subsetted data
- Predictors of positioning on axes:
	- Major habitat type is a significant predictor on fast-slow axis (tropical and subtropical species live longer lives than others, due to dominance of long-lived trees) (may also be sampling effect though)
	- Habitat type is only a weak predictor of position on reproductive strategy axis
	- Taller plants tend to have greater fast-slow axis scores (slower?) than smaller plant with apical meristems closer to or below the ground
	- Growth form and reprdouctive strategy axis (phylo stuff idk about)
	- Phylogenetic relationships play a weak role in PCA
		- a handful of exceptions suggest there may be infra-class structure
- Some overlap among growth forms and size, e.g., overlap between short-lived trees and shrubs overlap on fast-slow axis with longer-lived herbaceous perennials
	- Runkiaer growth forms used in here
- "Population performance" assessed through damping ratio, which measures the rate of  return to asymptotic dynamics following a disturbance (ref 22 - Caswell...)
	- Rate of recovery is associated strongly with position on both axes s/t rapid recovery occurs more with fast growth, high reproduction, and shorter generation time or low reproduction and frequent shrinkage
	- higher asymptotic population growth rates for faster-growing, iteroparous, highly reproductive species
		- lower growth rates for delayed maturity, low senesence, frequent shrinkage
##### Discussion
- The fast-slow continuum has considerable support in animals (refs, but most important seems to be 11: Stearns's 1999 book)
	- analyses of fast-slow have found secondary axes of altricial-precocial (req'ing parental care or not) and iteroparity (Gaillard et al.)
	- present analysis in plants does not suggest phylogenetic signal or influences of size, allometry, in contrast with animals
	- studies in animals typically find two axes capture ~80% of variation in animal strategies, but in plants only ~50% is explained
		- plants typically have life cycle complexity, e.g., dormant stages or long-term seedbanks
	- plants are indeterminate growers
- There were portions of the two-axis plane that was unoccupied (low on axis 2 and extreme on axis 1, i.e., high-shrinkage and low iteroparity low R_0 and either very fast or very slow)
	- these combos may be inviable, but also because the data was taken from natural populations, they may reflect current state of habitat as well
	- populations that are high on both axes (high iteroparity/R_0 and long-lived?) appear to be cases of successfully expanding populations, typified by successful invasive species
		- can't partition how much of this is due to habitat quality vs. other conditions


Super cool! Not entirely sure what axis 2 is capturing that axis 1 is not (how is iteroparity not loaded onto axis 1?) but otherwise great stuff. I wish they explained what the values of high/low damping ratio corresponded to; this is useful stuff and at the very least is informative in our case for giving us insight into the rate of return to the stable stage distribution. Probably would be good to look into the animal refs too... don't want to read all of Stearns but other refs in here are useful. Oli 2004, Gaillard dt al. 2005 are possible sources for this.


### Pimm, S.L., Jones, H.L., and Diamond, J. 1988. On the risk of extinction. Am Nat.

- Extinction risk is higher for small populations than large, from which it follows that (1) populations with more temporal variation in density will be at higher risk and (2) populations with low intrinsic rates of increase (r) will have increased extinction risk
- Two extreme cases of environmental change: (1) "solely" by demographic accidents in an unvarying environment, or (2) total environmental destruction
	- In models of (1), models assume constant per capita birth rate up to ceiling, K
		- these models predict rapid increase in time to extinction (T) with increasing K
	- Models of (2) should assume identical time to extinction for all populations (? why?)
	- Leigh's model with "modest" environmental variation and demographic effects expects T to scale with log(K)
	- Define C = (1/T)^(1/K), a "corrected" risk of extinction (why? and what is this?)
		- in the case of demographic accidents being the sole contributor to extinction risk, T = a^K and C is constant over K
		- If the time to extinction increases more rapidly than a^K, then C decreases with K; the reverse is true if T increases more slowly than a^K
		- because C is less sensitive to K than T may be, use C (but why)
- r, longevity, and extinction rates: can't provide direct estimates of longevity and r, but we can with body size, which is correlated with r
	- r and longevity are known to be related to body size
	- increasing r will proportionately decrease longevity (observed in animals... but not cited)
	- body size is associated with low r, which should increase extinction rates
	- large body size is associated with high longevity, producing lower extinction rates
	- possibly opposing effects here on extinction: is large body size helpful through increased longevity, or harmful through decreased r?
	- verbal argument: for a population of size one, large body size means higher longevity so slower extinction; for a large population, a large-bodied population will take longer to grow back following a reduction in numbers
		- "Although the large-bodied species may persist longer at low densities than a small-bodied species, remaining at low levels greatly increases the risk of extinction; a large-bodied, slowly growing species may still be at low levels when the next severe reduction in numbers occurs"
		- so, verbal argument suggests that large-bodied is advantageous at low densities, but a disadvantage at high densities
	- measuring in terms of lifetimes (rather than years), small-bodied species are always at an advantage due to high r (why?)
		- but conservation planning often requires planning on the order of years ("The year is not a biologically arbitrary unit of time")
		- also years are also relevant to scales of dispersal
	- Illustration of the above with Leigh's (1981) model, with exponential growth @ rate r up until K and only demographic accidents
		- Table 1 has numerical results from Leigh's model, showing that r has small effects at low K, but increasing r for larger K has increasingly drastic effects on time to extinction (i.e., sensitivity to r is low at low K high at high K)
	- temporal variability in population size also is a result of r (and correlatees with body size and longevity): large-bodied species should generally tend to survive environmental disturbances, but they also will have slower recoveries - what is the net effect on CV of population size?
		- using British bird data (birds subject to abrupt declines in density following harsh winters), small-bodied birds show greater reductions in size following harsh winters than large-bodied birds (Cawthorne and Marchant 1980)
		- but the majority of temporal variation in temporal densities comes from slowness in speed of recovery - Pimm 1984 shows that slowly recovering species have larger CVs than species that recover quickly
		- CV and r negatively correlate in this case, although correlation appears to be weak
- Predictions from above:
	- (1) Small populations will have higher extinction rates than large
		- (confirmed in many populations - see Diamond 1984)
	- (2) at low densities, small-bodied fast-growing short-lived species will have higher extinction risk
	- (3) at high densities, large-bodied, slow-growing, long-lived species will have higher extinction risk
	- (4) populations with high CVs will have higher extinction risk
- Test of these predictions from bird data from bird data from 16 islands off the coast of Britain
	- yearly counts of nesting pairs for 100 bird species, 355 populations (island-species combos)
	- species are taxonomically and trophically diverse
	- at least some migratory species are philopatric
		- migratory and resident species may differ in extinction susceptibility: wandering may cause local extinctions, split time between mainland and islands providing short-lived island populations, migration is risky but also may avoid hrasher winters?
 	- for each population, record mean number of nesting pairs years present; then for each species, calculate a mean number of nesting pairs per island averaged over all island (N)
	- get CV in number of nesting pairs for each populations when there were 5+ consecutive years of breeding, then estimated a mean CV over all islands with estimates
	- also recorded number of years between colonization and extinction w/ continuous breeding
	- classified birds by body size as large or small
- Analysis on reciprocal of # years each population survived (zero if never extinct), "risk of extinction per year"
	- excluded populations with mean number of pairs >18 (18 was largest mean size for which extinction occurred), leavig 316 pop'ns of 67 species
		- 62 of said species had at least one extinct population, five had no extincitons despite mean size that was smaller than 18
	- assume that, because K can not be truly measured (and is a modeling construct), N (avg. population size) is a correlate of real K
		- use this to estimate C = (1/T)^(1/N)
		- use of (1/N) here corrects for differences in mean population size (ahhh nice)
		- use of C here also produces results that are more approximately normal
	- handling non-extinct species: absence of extinctions may be informative at low densities but is less informative at high ddensities, which why populations were excluded if they exceeded mean size 18; for smaller populations C was treated as zero
	- some other technical considerations here - bullet point 4 I skimmed and did not understand fully
	- Analyses performed: ANCOVAs, with predictor N included?
		- analyses separated out by body size; separated out by migrants vs. residents
		- (some technical considerations about N being present in both dependent and independent variables... but they say it's okay. nothing alarming in residuals, any bias is going to be consistent across N)
- Results:
	- as expected, strong relationship where (1/T) decreases with N
	- body size effects: C ~ body size for large- and small-bodied birds signf. in slope and intercept
		- (does intercept matter here?)
		- small-bodied birds: correcting for population size (using C) removes relationship between C and N (n.s. slope)
		- large-bodied birds: C increases with N (p ~ 0.01), so times to extinction increase more slowly with population size
		- below a certain threshold (~7 pairs), large-bodied species have a smaller extinction risk than large; this reverses above the threshold
			- in the regression, the intercept of the regression line is lower than for small-bodied species, but greater slope (Fig. 4)
		- large-bodied residents are less extinction prone at low densities, but increases more rapidly with increasing density
		- this is in line with predictions (2) and (3)
	- migrants do have higher risk of extinction than residents but patterns with population and body size do not appear, although this may be due to to data scarcity
	- five species without extinction: all residents; three larger bodied with a span of densities (1-10) and two are small-bodied with high densities
	- residuals ~ CV does suggest increasing CV does increase residuals (i.e., extinction risk)
- Present study finds relationships between population size and risk found in other studies, although claims to do so in a better way than others
	- This relationship requires accounting for population size in analysis (and presumably to use C as a response and not just 1/T)
	- Indeed, other studies claiming to find effects of body size on extintion risk did not account for body size (and large-bodied organisms tend to have smaller population sizes)
		- Here, there is a residual effect of body size after accounting for population size
		- large body size is correlated with rate of population growth and longevity, meaning we can likely draw inference about longevity and intrinsic growth rates from results of body size
- Large body size, through correlations with population growth rates and longevity, may have mixed effects on extinction risk
	- here, the net result is that large body size lowers extinction risk at low population densities, but increases extinction risk (rel. to small body size) at higher densities
	- for populations at small densities, "demographic accidents" (stochasticity) are likely the largest cause of variation in size that increases extinction risk
- Differences between migrant and resident bird species - maybe not observed before! Migrants had greater corrected extinction risk than resident species of the same pop size and body size
- These results looked at short timescales and small population densities; for medium-sized populations, effects of population size may wane relative to other influences such as CV, r, or longevity
- Conservation concern: some may ask why care about populations of this size, rather than populations in the 20 - 500 size?
	- Fragmentation may produce very small populations
	- Reintroduction into the wild may face such questions, e.g., one population of size 50 or five populations of size 10
	- Captive breeding programs in zoos typically have small populations
	- (Also, rescue means populations may decline to very small sizes!!)


Really useful and informative stuff. First thing to note is that there is a need to correct for population size, which may possibly correlate with other traits (not measured here, but we'd expect large-bodied organisms to be in smaller populations)

But otherwise, I really love the arguments suggesting and the analysis supporting the density-dependent influences of population size. This is an apparent tension that I have struggled with and not really put to words before! The idea is that at small size, large-body size (or, more appropriately, its correlate of higher longevity) will increase extinction risk (or, relatedly, prolong extinction). However, as body size increases, intrinsic population growth rates (r) decrease, such that there is a slower increase. There's some surprising stuff in here about CV being lower for low-r populations (is this due to the decline or due to the r?), but the CV being associated with extinction risk is also useful.

Some of the stats methods in here are ehhh but I think it's still fine to use on the whole.


### Oli, M.K. 2004. The fast-slow continuum and mammalian life-history patterns: an empirical evaluation. Basic and Applied Ecology.

- r-K selection: early attempt by MacArthur and Wilson 1967, Pianka 1970, to explain diversity of LH patterns but eventually deemed to be insufficient on its own
- Proposed fast-slow continuum:
	- "Fast" end - early maturation, large reproductive rates, short generation times
	- "Slow" end is opposite end
	- "some" empirical support this is qualitative
		- Gaillard et al. (1989) did this with a PCA but did not explore what the first principal component measured
		- Read and Harvey (1989) looked at mortality relative to body weight
- F / a ratios as proposed by Oli and Dobson (2003) are ratio of fertility rate (F) to age of first reproduction
	- low F/a is slow, high F/a is fast, cut-offs provided here
- Another hypothesis: LH patterns are due to adult body size
	- also (not sure why in this para) evidence that LH evolution may be phylogenetically constrained
	- relative influences of these are not clear
- Here: characterize pace-of-life of mammals with F/a ratio, then use population dynamics modeling to look at some other stuff
##### Methods
- LH data for 138 mammal populations, estimated fertilities F and survival probabilities P, built matrix models accordingly
	- caveats with Leslie matrices:
		- variation in LHs can cause sizes of projection matrices to vary in size and therefore they can not be directly compared
		- timing of certain life cycle events are not explicitly modeled in these matrices (e.g., ages of first and last reproduction)
	- divide life histories into a two-stage cycle: pre-reproductive J, post-reproductive A
		- survival to adulthood has probability P_j per time unit with a time steps to do so
		- assume that upon reaching reproductive maturity they reproduce with average fertility F and survive with probability P_a until age of last reproduction omega
		- some work from prior Oli and Zinner papers demonstrating that the partial life cycle approach desribed here has deficiencies but overall does a good job of addressing above-mentioned caveats
	- consequence: life cycle can be characterized by five parameters: a, omega, P_j, P_a, F
	- estimated lambda from characteristic equation of these life cycles (neat!)
		- With lambda obtained, can estimate elasticities
- For each population, estimate F/a, bin them as fast, slow, medium
	- perform AoV to assess differences in parameters vary among three groups:
		- variables testded: log body mass, LH params above, lambda, elasticities of each rate
		- REGW multiple range test (is this post-hoc or what?)
- Phylogenetic constraint hypothesis tested with nested ANOVA (order main effect, family nested within order)
	- analyze residuals to see if fast-slow continuum still appears after controlling for phylogeny
- Allometric constraint hypothesis tested by regressiong LH vars, lambda, elasticities on body mass
	- analysis of residuals: do effects among variables persist after removing body size?
##### Results
- Even split among the three pace of life gropus
	- body mass, all LH vars, lambda differed among groups
		- Table 1 has results
		- higher lambda (>1) for fast than medium/slow (~1)(see Fig. 1F)
		- higher F for fast than medium/slow
		- for P_a: slow > medium > fast
- Elasticities: increasing F/a increased elasticity of lambda to reproductive parameters (a, F) and decreased elasticity to survival parameters
- Analysis of residuals after controlling for body size (allometric constraint hypothesis) shows strong evidence of the residuals even after controlling for body size
	- body size relationships being incorporated into models also did not change relationships between F/a and elasticities
- After removing phylogenetic effects, correlations between F/a and life history variables became weaker but still were sign. non-zero

Very cool. Did not critique demo and stats techniques that much, but I did find the approach to be straightforward! I'm not using age of first reproductiona as a parameter here but results here I still think are applicable or useful to us, particularly as there are variables for survival probabilities and fertilities.

Useful results here: (1) yes, evidence of a robust fast-slow continuum that involes advancing age (even if age at first reproduction) and fertility rates (2) differences among these groups in lambda. The elasticities argument is interesting too although I would need to think about it more, but the idea here is that sensitivity to survival increases with increasing F/a (or, slowing pace of life) which makes sense. Does this influence our results at all? I want to say no because the proportional effects of selection/maladaptation load here do not necessarily have to do with any particular rate, yes? Also seems cool and worthwhile that these results remain if controlling for body size, although I'm not sure this matters in our hypothesis.


### Gaillard, J.-M., et al. 2005. Generation time: a reliable metric to measure life-history variation among mammalian populations.

- Response to Oli and Dobson (2003), who found from analysis of 142 mammal populations that the F/a ratio (as described above) predicts elasticities of vital rates to lambda and provides a good proxy for position on the slow-fast continuum
	- they also did not find strong empirical support for age at first reproduction affecting (...? something?) and did not find relationships between phylogeny and body size predicted in literature
- Critique of Oli and Dobson, presented here:
	- Authors ignore prior theoretical work that demonstrates that elasticities of lambda to vital rates shrink with generation time (Charlesworth 1994)
	- F/a ratio has "no theoretical justification"
	- Work presented here suggests that age of first reproductino is likely a good index of position of mammals on the slow-fast continuum
- Prior work on generation time:
	- Leslie (1966) and Charlesworth (1994) demonstrated the importance of generation time, which could be defined as weighted mean age of child birthing-mothers in a population
	- definition provided in here, but with simplifying assumption of fecundity that is independent of age, very nice expression: a + (s / lambda - s)
	- as such, generation time is a function of all vital rates
	- (Lebreton and Clobert 1991 seems useful)
	- seems likely that generation time would be a reliable measure of populations on the slow-fast continuum although this has not been tested yet
	- re-analysis of a subset of O&D's data: PCA on five LH parameters suggested by O&D (same as above)
		- first axis of PCA accounts for 69% of variation in variables, corresponds to timing of start/stop of reproduction, fecundity, survival
		- generation time as defined here is highly correlated with PC1 (rho = 0.9)
	- note, though, that generation time assumes a stable age structure, which may be violated in temporally-varying environments
- The F/a ratio: not theoretically supported, may not be intepretable as originally claimed, and something about elasticities
	- that second point: O&D correctly estimated F, but this means that F/a includes adult survival and is thus can't be straightforwardly used as a measure of reproductive output to age at maturity
	- however, there is a strong correlation between F/a and PC1 (rho = -0.92), but there are several other variables that correlate with PC1 which would also work well as fast-slow indicators
- Is fast-slow placement truly independent of phylogeny or body size, as claimed?
	- ANCOVA where log-generation time is dependent variable, log body mass is covariate, order is a factor with eight levels
		- common slope to all orders of body mass on generation time, although orders have wildly differing intercepts
		- additive effects of body mass and mammalian order explain ~60% of variation in generation time
	- n.b. that O&D's survival estimates come from life tables, which require strong assumptions for valid and accurate inference (stationary age distribution, equal probability of sampling all individuals)
- Conclusions: generation time is a good metric for assessing importance of life history variables to population growth rates and position of a population on the fast-slow continuum
	- generation time may be just as easy to estimate in the field as other variables

Short and sweet. Love it. Okay so this claims that generation time should work just as well if not better as an indicator for position on the fast-slow continuum (for mammals at least) as F/a. Some arguments herein. For the sake of what I'm interested, it is good to know that generation time is a proxy for LH-placement - I even think the expression for T_b appears in some form in my analysis (if a = 0, then T_b would be s/(lambda-s), which I feel appears somewhere...). This does seem to fit in better with my framing than Oli 2004, although Oli 2004 may be useful for providing more evidence of fast-slow placement with lambda (where medium/slow species have lambda ~ 1 and fast species have greater ones). Use these both!


### Rosenblad, K.C., Baer, K.C., and Ackerly, D.D. 2023. Climate change, tree demography, and thermophilization in western US forests. 

- Global climate change is, in some systems, increasing the relative abundance of heat-tolerant ("thermophilic") taxa
	- widespread, observed across taxa, regions, spatial scales (see refs 6-8)
	- but, rates of thermophilization vary with much unexplained variation
	- Possible influences on the rate of thermophilization
		- rates of warming ofc
		- physiology, e.g., plant water and/or temperature regulation mechanisms
		- in forest and understory communities, changes to canopy structure influence solar radiation penetration and temperature changes
		- topography
	- can measure these in part with community weighted mean of some climate-driven trait, e.g., climatic niche
		- relative weights of each species and their traits depends on abundance, and change in abundance relies on demographic processes
		- can use a decomposition to get relative influences of demographic processes (for each species) to community mean change
			- e.g., thermophilization through enhanced recruitment of thermophile species, or through increased mortality of cool-adapted taxa
- Here: 10 year dataset of ~45K forest subplots using USDA FIA dataset
	- model mean temperature indices of tree communities over time using Bayesian hierarchical models
	- includes effects of long-term average climate, "recent" climate change, topography, disturbance, other predictors
	- Q1: are western forests undergoing thermophilization?
	- Q2: what influences the rate of thermophilization?
	- Q3: what are separate contributions of mortality, growth, recruitment to thermophilization
##### Methods
- Baseline censuses 2000 - 2008, resureys 2010 - 2018, in Phase 2 FIA (>10% tree canopy coverage)
	- plots each have four circular subplots, radius 7.32m, wherein DBH of trees with at least DBH of 12.7 cm measured in each subplot
- 110 species in dataset; assessed "niche means", "modeled niche optima", and "simple niche means" of each
	- measure mean or optimal mean annual temperature across sp range
- Regression models use a community temperature index of one subplot at one time point as a response
	- estimate thermophilization with a binary variable for the difference between baseline and resurvey (resp., T1 and T2)
	- fixed effects: baseline MAT, baseline MAP, baseline mean annual CWD (climatic water defecit), temporal differences in variables
		- topographic heat load, binary effects of fire or insect damage
		- intercept is mean of T1 observations, non-interaction coefficients are effects of variables with T1
		- then included interactions of predictors with temporal changes (i.e., difference between driver's influence in T1 and T2)
	- random effects of subplot, subplot-level spatial random effect
- Demographic processes:
	- recruitment is entry into the 12.7+ size class
	- re-ran analyses described above, except with response for T2 now somehow includes effects of one demographic processes
		- here, T1 values are all unchanged, T2 values reflect only observed effects of one processes with remaining processes removed and said trees hel at their T1 observed value?
- Some cool climate data that might be good to look into
##### Results
- Community mean T1 temeprature ("baseline community temp" correlates very well with baseline annual mean temperature; 
	- in resurvey, mean temperature change of ~1/3 deg C
		- mean thermophilization rate across forest was .0039C/y (compared to .032 C/y)
	- thermophilization was greater in plots with more warming (1C thermophilization per 1C warming)
	- also greater thermophilization in drier areas
- Mortality had greatest contribution to thermophilization (~1.1 deaths/subplot)
	- minimal contributions of recruitment at both 12.7cm and 2.5cm thresholds
		- from subplots where saplings 2.5-12.7cm were recorded, 87% of saplings recorded in the 12.7cm class in T2 were saplings in T1 (i.e., 13% made it from zero or <2.5 cm to 12.7cm between surveys)
- Fig. 2: strong pos. effects on baseline temperature, baseline CWD; weak pos. effects of baseline precip, baseline topographic heat load, weak neg. effect of baseline conifer dominance
- Mortality- and recruitment-driven thermophilization greater in plots experiencing more warming, whereas growth-driven thermophilization was stronger in plots with more CWD, lower precip
	- look at data-rich Fig. 4
	- thermophilization greater on north-facing slopes than southern-facing (mortality-driven)
	- thermophilization greater where there was insect damage observed between surveys (strongest per-sd coefficient effect)
	- net-~zero effect of fire, but this is due to differences in mortality- and drought-driven processes
	- greater effects in conifer-dominated subplots, apparently due to growth
##### Discussion
- Widespread pattern of community weighted mean temperature niche optimum increasing, but only at ~1/10 of the rate of actual temperature changes in the same period
	- projected forwards, the ~3C avg. increase in temperatures may only be accompanied by a <.5C increase in community mean temperature niche, in which case the "climate debt" will go up by ~2.5C in addition to other anthropogenic lags
- Some of these processes may have been operating before the study's timeframe, with possibly complex relationships among blah blah blah...
- Also, the study period coincided with a severe drought, in which case patterns observed may be greater than realized future trends
- Lack of contribution from recruitment to thermophilization suggests that new recruits have thermal niches similar to their baseline communities
	- instead, thermophilization is from mortality of individuals (species) unable to withstand the increasing temperatures
	- loss of vulnerable species...
	- this contrasts with ref 22 that "hot spots" of fecundity and recruitment for tree species in the western US are shifting toward cooler, moister regions at rates commensurate with climate change
		- perhaps this is because present study did not look at seedlings directly
		- if this is correct, then perhaps effects of recruitment-driven thermophilization will appear in future resurveys
		- or, difference in spatial grain of analysis explains this difference... didn't try to understand this lmaooo
- Insect damage increased rates of thermophilization, but conflicting effects of fire
	- bark beetles tend to attack drought-stressed trees in the first place
	- insect damage may also have direct effects on understory by modifying the canopy
	- fire produced greater mortality-thermophilization and weaker growth-thermophilization
		- a prior study suggests that cool-associated tree taxa have poor tolerance, in which case increased temperature may just be killing these species
- Stronger thermophilization on pole-facing, cool hillslopes, possibly because equator-facing sites already had more tree-adapted species at baseline
	- produces possible biotic homogenization in topographically-heterogeneous landscapes, incl. reduced turnover among subplots and reduced landscape-level richness


Ooh this is super interesting. Community-weighted trait means similar to a population-weighted trait means - interesting! Wonder also about community-weighted variance and its effects on thermophilization rates... how easy is that to figure out? Or can one even figure out community-weighted variance in the same way? Or is this really just species diversity? Or maybe related but not a perfect correlate.

But regardless it seems cool. The increase in community means did appear to be driven more by survival than anything else. What I think we can use from this is that with increased warming (and related changes), there is a large mortality that influences trait means... not totally sure if there is a straight line here to population means, but assuming there is variation among species traits that influence survival, there should also be intraspecific trait differences that influence this as well right? If there is heritable variation in these traits (big if) then we should expect adaptation.

Note that this fits in much more with the B&L-mode of change! In fact, if we can treat the communities as similar to a population with trait variation, then this in fact does look like their lagging rates of change. Very interesting!

Also I wonder if there is a way to apply longevity to this? Or at least look at species traits from FIA.


### Lande, R. 1982. A quantitative genetic theory of life history evolution.

- Life history evolution typically assumes some quantity is being maximized...
	- Maynard Smith (1978): optimize expected number of individuals over lifetime
	- Fisher (1958), Charlesworth (1980): in a constant environment, LH evolution should maximize population growth rates or population size itself
	- It's also possible that there are multiple stable equilibria
- Dynamic models of LH evolution should include necessary or relevant constraints (can come from quantitative genetics theory)
##### Model
- Assume same patterns of age-specific seletion on both sexes, no sexual dimorphism
- Life history represented as a vector of traits, $z_n$, with elements properly scaled to be gaussian with mean $\bar{z}$ and vcov matrix $P$
	- each of the elements of $P$ can be decomposed into additive genetic and environmental components, s/t $z = x + e$, and $P = G + E$
	- assume breeding values and environmental effects are multivariate normal
	- any system of mating consistent with MV normality is allowed (totally random, assortative by phenotype, assortative by age)
	- elements of $z$ may be, e.g., individual growth curve parameters
- let individuals with phenotype $z$ have fertility and mortality rates at age $a$ defined by $m_a(z)$ and $mu_a(z)$ (respectively)
	- probability of survival is $l_a(z) = \exp{-\int_0^a mu_y(z) dz}$, i.e., exponentiated negative integral of mortality over age
	- fitness at age $a$ of individuals with LH z is $w_a(z) = l_a(z) m_a(z)$
	- let the expected fitness of a newly-conceived cohort at age $a$ be $\bar{w}_a$
	- with assumptions of maximum age of reproduction or non-zero survival, can define $\int_a^\infty exp{-ra} \bar{w}_a da = 1$ 
		-  $r$ is rate of increase of population in the absence of evolutionary changes when phenotypic variation is non-heritable or population is at equilibrium
	- (some other stuff in here: generation time, stable age distribution, reproductive number...)
- Assume selective forces on a population are weak and change only very slowly over time (rel. to typical demographic processes)
	- population then should (disregarding transient dynamics) track a "stable" age distribution with population size $N$ governed by $dN/dt = rN$
	- under assumptions that mean the rate of evolution of mean life history phenotype in unselected (newly-conceived) individuals is $d\bar{z}/t = G \grad r$
		- $\grad r$ is the gradient of $r$ with respect to z, i.e., vector of partial derivatives of $r$ wrt elements of $z$, i.e., sensitivity of growth rate to each LH trait value
	- above can be expressed also as $d\bar{z}/dt = \sum_j=1^n G_ij dr/d\bar{z}_j$, i.e., sensitivity of growth rate to each LH trait value, weighted by variance/covariances
	- it can be shown (appendix) that selection gradient takes form $\grad r = P^{-1} Cov(w(z), z)),$
		- RHS of this equation is partial regression coefficients of trait on relative fitness (holding all else constant)
	- also in appendix: relative fitness of a LH trait  can be approximated by age-specific fecundity and mortality rates
		- weights on these processes depend on present and future reproductive values
	- conclusion from all this: changes in age-specific fecundity and mortality rates, independent of the phenotype $z$, can still alter relative fitnesses of individuals with different life histories (duh?)
	- defining $\bar{W} = exp{r}$ as the mean absolute fitness per unit time allows $\grad r = \grad \log{\bar{W}},$ allowing for some other things...
	- from this it can also be shown that $dr/dt \geq 0,$ i.e., rate of population growth is always increased by life history evolution at rates proportional to the selection gradient and the degree of additive genetic variance/covariance $G$
	- in a constant environment, fitness (and population growth) will typically be optimized at some local optimum of $r,$ akin to the adaptive topography for phenotypic evolution
##### Discussion
- Some studies showing effects of selection and relating them to variances and covariances in fitness components (often negative)
	- morphological, physiological, behavior characteristics are often under joint stabilizing selection
		- as trait distributions move towards this optimum, additive genetic correlations will typically decrease, depleting additive genetic variance and increasing the proportion of nonadditive genetic variance in fitness (is the depletion due to stabilizing selection??)
	- for characters that are not direct components of fitness, approximations where the phenotypic vcov matrices are held constant typically perform well for a few generations

This is more about the evolution of LH rather than evolutionary consequences of LH. Not sure there's a ton in here that is directly relevant to this project, particularly as I am only working on a single trait.


### Willi, Y., and Hoffmann, A.A. 2009. Demographic factors and genetic variation influence population persistence under environmental change. Journal of Evolutionary Biology.

- Lynch and Lande's (1993) model: first one to merge quantitative genetics with populatin persistence
	- Growth rate of sexual population comes from four loads:
		- Phenotypic variation during development (? reads like variance load)
		- Lag load is how much the population mean differs from the optimum
		- Two other types of load due to genetic drift and stochasticity
			- Genetic drift causes deviation of realized phenotypes from expectation and is (inversely) proportional to N_e
			- Stochasticity in selection also causes departure from the optimum
	- Burger and Lynch (1995) noted that stochasticity in genetic variation, demography, and environment can produce extinctions
- Empirical evidence of load including those not mentioned in L&L93 or B&L95, particularly loads at small population size
	- inbreeding depression, accumulation or fixation of deleterious alleles, Allee effects, demographic variability from stochasticity in vital rates
	- refs herein
- Here: empirical look at these types of load and their relationships to population size in a Drosophila system
##### Methods
- Study organism: *Drosophila birchii* native to tropical rainforests in NE Australia
	- collected from six populations, establishment of six "isolines"
	- 28 experimental populations of size 20 (18 reps), 100 (eight reps), 1000 (two reps), created from 15 randomly chosen isofemale lines per site
		- more populations founded by mixing flies from different origins (>20 lines?)
	- populations maintained in vials w/ 20 individuals/vial
- Estimated natural rate of increase ($r = \log{R_0}$ where $R_0$ is mean number of offspring per individual) beginning at generation 4
	- in proceeding 10 generations, 12 populations at size 20 went extinct
	- assess if population mean $r$ and within-population s.d. in $r$ depended on population size (mixed effects? models)
- After ten generations held at or below initial size, populations exposed to directional selection
	- populations flushed (???) to at least 200 flies for 5-6 populations
	- selection for heat-knockdown resistance for 4-5 generations
		- resistance is time until fly was incapacitated following heat exposure
		- most resistant 50% of females transferred to vials for ovipositioning, mated with males with similar resistance levels
		- compared to lines with >200 individuals
		- estimated selection intensity, selection differential, heritability, additive genetic variance
- Simulation stuff: estimate maximal rate of env. change $k_C$ (minutes of heat resistance per generation) for which populations can maintain $r \geq 0$
	- assume trait is subject to stabilizing selection, with optimum moving deterministically and unidirectionaly
	- let $r^*$ be the growth rate in the absence of selection (but including drift and inbreeding depression)
	- formula for $k_C$ as defined by L&L93
	- also use some analytical estimates from B&L95 to estimate growth rates, mean lag, mean variance
	- simulation stuff to get time until extinction
##### Results
- Populations of size 20 had lower $r$ over time than the other two sizes
	- but no sign. difference between 100 and 1000
	- similar pattern for standard deviation in $r$
	- suggests genetic drift (independent of env. change) in small populations
- Heritability of knockdowns was low, but insign. differences among the size groups
- Small population had significantly lower additive genetic variance than the size-100 populations
	- insignificant differences at other contrast levels but mean $V_A$ in size-1000 populations was smaller than for size-100 populations (!)
	- correspondingly higher $k_C$ for size 100 compared to size 20, with size 1000 intermediate (!)
- Tight relationship between mean reproductive rate, stochasticity, and time to extinction in the small populations
	- small populations had no relationship between simulated time to extinction and $k_C$
	- for low $k_C$, differences in time to extinction among size-100 and size-1000, but for intermediate or fast $k_C$ nearly-identical extinction times
- Fig. 3: simulated times-to-extinction as a function of growth rate, stochasticity in growth rate, and additive genetic variacne for different sizes
	- rapid (<50 gen) time to extinction for size-20 populations, unless $r$ was very large and $V_A$ was small
	- for size-100 populations, increasing $V_a$ increased time to extinction only for high $r$, low stochasticity in $r$
	- for size-1000 populations, for intermediate-to-large $r$, increasing $V_A$ increased time to extinctions rapidly, similar effects for lower stochasticity in $r$
##### Discussion
- Differences in reproduction and additive variance between small-medium populations but not so much for largest ones
	- intermediate-size populations did seem to have some negative impact of stochasticity and drift (support...?)
- Relationship between population size and fitness: size-20 populations seem to have ~50% of the fitness of size-1000 popiulations, in line with other studies (e.g., Reed 2005)
	- seems like an effect of inbreeding depression, or at least that's what Reed (2005) concludes and O'Grady et al. (2006) found
- Yes - higher stochasticity in $r$ for smaller populations
- Additive genetic variance: the small-medium comparison result makes sense, the medium-large is harder to explain
	- may be an effect of drift-induced rise observed in some other settings (e.g., Van Buskirk and Willi 2006, wherein traits with a non-additive basis may have an increase in $V_A$ following a bottleneck)


Hmm... might have been useful for the NDD project and its basics but less so here. The interesting results are simulation-based so not something I'm inclined to cite for empiricism. 




# Incomplete stuff

### Isaac, J.L. 2009. Effects of climate change on life history: implications for extinction risk in mammals. Endangered Species Research.

- 

### Beckerman, A.P., Benton, T.G., Lapsley, C.T., and Koesters, N. 2003. Talkin' 'bout my generation: environmental variability and cohort effects. Am Nat.

- Cohort effects: delayed life-history effects synchronized among groups within a population
- In variable environments, there can be time-delayed effects of environment on performance of individuals at future times (through life history traits)
	- (citation here: Beckerman et al., 2002)
	- within-generation effects (e.g., delayed-life history effects) can produce cohort effects
	- maternal or paternal environment effects are examples of among-generation effects
	- Cohort effects will arise when variance in an LH trait within a group is appreciably *smaller* than variance in the trait among the population such that cohorts are statistically distinct in the trait (think like ANOVA within/among group differences)
	- Some empirical examples cited herein
- Lindstrom and Kokko (2002): cohort effects (variation among cohorts) can be produced by cohort-specific density dependence
	- if the NDD is non-linear curve of performance ~ density, then variation among cohorts will produce non-linear averaging in the average effects of NDD across cohorts (rel. to a population with variation that is not structured)
	- L&K02 also show: so long as "underlying population dynamics are stable" (?), then cohort effects can introduce variation and fluctations akin to effects of environmental variation
		- or, cohort effects can decrease temporal variation when "underlying dynamics are variable"


### Vinton, A.C., Gascoigne, S.J.L., Sepil, I., and Salguero-Gomez, R. 2022. Plasticity's role in adaptive evolution depends on environmental change components. TREE.

- Does adaptive plasticity facilitate adaptive responses in changing environments, and if so when and by how much?
- Contradictory evidence on whether plasticity facilitates or hinders adaptive evolution; few general patterns
	- H1: Plasticity slows phenotypic change by masking genotypic variation (e.g., Bogert effect?)
	- H2: Plasticity facilitates evolution by buffering populations while genetic change occurs (e.g., plasticity-first hypothesis, Baldwin effect)
- Under "moving optimum theory" under (continuously?) (uni-?)directionally changing environments
	- There is a critical rate of change that the mean population phenotype must be capable of tracking (probably from B&L95)
		- During this tracking there will be a phenotypic lag (larger lags increase extinction risk)
		- Magnitude of lag determined by evolution (selection strength, genetic variatioin) and ecologicy (life history, plasticity, population dynamics)
	- Plasticity can contribute to persistence and adaptation in relationship to this lag
- Adaptation requires selection on heritable variation in a trait
	- Population-dynamic traits tend to be quantitative traits (refs: Yamamichi 2022, Hill 2010)
	- Breeder's equation: change in a trait is selection differential x narrow-sense heritability
	- Mean change, variability of change, and temporal autocorrelation can influence heritability, genetic variation, selection
- (Mean) rate of environmental change is important
	- Sufficient environmental change is required to produce strong enough selection for the population to track the moving optimum
	- Slow environmental change produces weak selection (ineffective?)
	- Additive genetic variance and heritability can change (possibly increasing genetic potential)
	- But, change beyond the critical rate can be too fast for the optimum to follow, producing increasingly large lag and eventual extinction
- Environmental variation and autocorrelation (novel and harsh environments)
	- "Moderate" environmental variation can optimize selection
	- Positive temporal autocorrelation can also increase additive genetic variance over time (increasing tracking ability)
	- More variability and less autocorrelation means more exposure to novel, usually unfavorable environments
	- "Temporal refugia" (should arise from positive autocorrelation?)
	- But, exposure to unfavorable environments can also increase additive genetic variance (??) (Hoffman and Merilla 1999)
		- occurs because selection will not remove mutations maladaptive only in rare environments
			- (not removing these mutations, because they occur rarely in novel environments, increases variance)
	- Harsh, novel environments can also decrease genetic variance
		- e.g., environmental conditions may prevent expression of underlying genetically-determined trait benefits, reducing heritability (because genetic influence on the trait is lowered)
- Although usually assumed that variation hinders population growth, there variability can have positive or negative effects on population dynamics
	- E.g., through varying effects of density dependence
	- Might also see transient dynamics, or unexpected dynamics due to deviation from age distribution
	- Nonlinearities or correlations between/among vital rates in differing populations may produce buffering strategies
		- Plasticity of vital rates can also influence buffering
	- Likewise, diferent phylo- histories or life history strategies will have different effects on environmental fluctuations
		- Positive autocorrelation increases likelihood of long stretches of adverse conditions
		- Positive autocorrelation also increases likelihood of positive growth that replenishes size and allows tracking of trait values
		- Paniw et al. 2018: density-independent, stage-structured population models show that pace of life and iteroparity positive correlate (??) with sensitivity to autocorrelation
- Within-generational plasticity: because plasticity can evolve or have different forms, its effects on the phenotypic lag are complex
	- Typically, plasticity increases phenotypic lag
	- Assumption under moving optima: plasticity can buffer decreases in population size but come at an energetic cost
	- Mixed results about environmental variability and plasticity combined effects
	- Plasticity's effects rely in part on unreliable cues, in which case autocorrelation is helpful (inasmuch as autocorrelation is a proxy for predictability)
##### Hypotheses
- Benefits of plasticity on adaptation in response to increasing mean rate of change
	- HA1: Plasticity benefits increase with faster environmental change.
		- Slow change means weak selection, high population growth and high heritability of fitness
			- Then, plasticity has little benefit
		- But, mean environmental changes that are too fast may be assisted if plasticity can help the population "catch up" and maintain genetic diversity
	- HA2: Plasticity benefits decrease with faster environmental change
		- Weak selection -> lag load -> plasticity overcomes
		- Increasing rates of change might incur larger costs than benefits of plasticity
		- Small populations under high rates of environmental change are susceptible to drift, which plasticity may increase by moving the phenotypic average and obscuring genetic variation from selection
	- HA3: Plasticity benefits are highest at an intermediate rate of change
		- Some combination of HA1 and HA2
- Benefits of plasticity on adaptation in response to increasing environmental variation
	- HB1: Benefits of plasticity increases with increasing variation
		- With increasing variability, plasticity can dampen detrimental effects of very large fluctuations and buffer populations from extinction
		- This probably has diminishing returns such that some variation is too extreme
	- HB2: Benefits of plasticity decrease with increasing variation
		- Under more stable environments, plasticity speeds up the process of fixing advantageous traits
		- More (environmental) variation decouples phenotypic and genotpyic selection, possibly producing genetic maladaptation
	- HB3: Benefits of plasticity are highest under low and high environmental variability
		- Ability of the trait to reach the peak of the fitness landscape on its own (without plasticity) might be higest under intermediate rates of change
		- In which case, plasticity would have least marginal benefit under intermediate environmental variability
- Benefits of plasticity on adaptation in response to autocorrelation
	- HC1: Benefits of plasticity increase with increasing autocorrelation
		- Higher autocorrelation means higher reliability of temporal cues and predictablity of future states
		- So, plastic responses will be more accurate and more likely to assist under tracking
	- HC2: Benefit of plasticity to adaptive evolution decreases with autocorrelation
		- Autocorrelation can occur at various scales with varios lag-lengths
			- What if the lag is out of sync with the pace of life history? (e.g., generation time)
		- Autocorrelation may cause populations to sit in unfavorable conditions for long periods of time, reducing genetic variation
			- With reduced genetic variation, it's harder for plasticity to help track the moving optimum (but also plasticity might be needed even more here?)
##### Conclusions
- 


### Chevin, L.-M., and Lande, R. 2010. When do adaptive plasticity and genetic evolution prevent extinction of a density-regulated population? Evolution.

(reading again... this time though I'm focusing on the plasticity part)

- Hendry et al. 2008 suggests that much observed phenotypic change following environmental change is phenotypic plasticity, not evolution
	- Lande (2009) showed that with variance in plasticity (slope of norm), after a change (beyond "background" fluctuations in magnitude), phenotypic adaptation involves rapid increase in plasticity and then slow genetic assimilation of the phenotype
- Model: z is *not* a major component of fitness???
	- z = a + b*eps + e
		- a = breeding value, b = plasticity, eps is environmental state, e ~ N(0, sig^2_e)
		- a and b are bivariate normal, sig_a, sig_b, assumed uncorrelated
		- before selection, given eps:
			- zbar = abar + bbar * eps
			- sig2z = sig2a + sig2b eps^2 + sig2e
		- assume sig_a, sig_b are *constant* across generations (*ooo)
	- z subject to Gaussian selection, exp(-(z-B(eps)^2 / (2*omega^2))
		- B*eps is environmental optimum
	- z also subject to cost of plasticity, C(b) = exp(-b^2 / (2*omega_b^2))
	- selection acts first on b (plasticity costs evaluated) then on z (phenotypic load cost)
	- prior state of environment: eps = 0, abar_0 = 0, bbar_0 = alpha * B
		- alpha is "initial relative plasticity" (there is a similar term in Lande 2009)
	- as with Lande (2009): phase I (given sufficient variance in plasticity) means substantial increase in plasticity that accelerates approach to new optimum, phase II plasticity decreases due to cost and is compensated by increase of "elevation" (bvs) of reaction norms
		- equation for maladaptation during any time step t of phase I
		- "If the initial relative plasticity is moderate to low (if alpha^* < 1/2 [...]), consistent with a cost of plasticity, then most of the evolutionary recovery occurs during phase I of phenotypic evolution"

(return... was not ready to think about the math here as seriously as I think I need to)

