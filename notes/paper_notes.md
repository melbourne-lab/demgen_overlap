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

