# Meeting with Dan Doak (November 4 2021)

Meeting with Dan (over zoom) to discuss the project, differentiation from the Schmid et al. effort. 

Schmid's efforts are interesting but limited. In a sense the results are to be expected the way the model is set up. But there are ways to make the model more realistic that may change the results in unexpected ways.

##### Adult survival

In real life cycles, random environmental variation does impact juveniles, but it also is more likely to influence adults. This influence on adults is more important (long-lived organisms more sensitive to adult survival). E.g., pine trees and beetles - adults are more likely targets of attack.

Consider scenarios where multiple vital rates are changing wrt the environment. 
	- E.g., a simple life history with two vital rates - could possibly do something analytical here
	- Adult and juvenile survival both - could do this with a tuning parameter where the balance of effects on adult/juv survival

Schmid's results may only be true because of the low elasticity of the life cycle to juvenile survival. Longevity makes sense here because relatively less of the life cycle is subject to selection...

(Also note here: Schmid has NDD influencing the juvenile stage only, but the density depends on the number of adults? the adults are suppressing the juvenile stage? interesting)

What to measure? Not only extinction risk itself, but also think about that trade-off between longevity and adaptability. Maybe it disappears?

##### Selection over the life cycle

I had been thinking about this question in terms of "intertia" of having older individuals adapted to prior environments sticking around and contributing to future generations. Others (e.g., Schmid, but also Ellner and Hairston) have been thinking about this in terms of parts of the life cycle being sheilded from selection pressures. These are not identical - my original project was trying my best to not include the latter aspect.

Dan points out that a good way of thinking about this is "integrated selection pressuure" over the life cycle. For an annual it's obvious what this is, but there should be a way to do this for longer-lived individuals as well.

One way is to do this analytically, using a sensitivity approach using the chain rule (sensitivity of lambda to vital rates, sensitivity of vital rates to growth). However, this only gives sensitivity to a single environmental change (typically a small one). One could also use Tuljapurkar's approximation to get sensitivity to variance in lambda; a drawback here though is this only applies to simpler modes of environmental change (and may fail with more realistic modes of change - see below).

How to actually operationalize what this is? Another way to do this is to just compare (could be done with simulations) relative fitness approach - compare lifetime fitness of maladapted individual with max possible fitness.

Has this idea of integrated selection pressure been explored or developed thoroughly elsewhere? Maybe in some of eco-demography models but not a ton. Caswell's book would be a good place to look for this.

Also useful to think of this in terms of the Breeder's equation, which is modeling generational change. Here we'd also want to think about generation time and controlling for that. Think of the problem posed last spring - 30 time steps is 30 generations for an annual but 6 generations for a population where individuals live for five years. Would be good to standardize some of these effects by number of generations to make sure that that isn't explaining dynamics. Do effects go away when doing this? [Well, also keep in mind that under unidirectional change the environmental state will be different after 6 years compared to 30].

##### Fitness versus lambda

These will only be unequal with deviations from the stable stage distribution. But, under realistic modes of environmental change, this is actually very likely! There is analytical work on transient elasticity analysis but it makes simplifying assumptions that may not make it worthwhile. Simulations would be better.

It's possible that fitnesses don't actually change that much (rapid adaptation) but that doesn't matter very much for extinction. [Is extinction the most important response? What we are most interested in? Does it matter if something adapts sufficiently if it goes extinct? What does it even mean to adapt sufficiently if the population still goes extinct?]

Is this distinction between which is more important - population growth or mean fitness - something that's been studied, debated, or definitively answered? [I don't remember Dan's exact answer but my sense is no]

##### Model setup

I should be able to adapt code from prior sims to this life cycle. Don't worry about being able to match Schmid's results exactly. 

Realistic environmental change - lots of stuff showing that this is unidirectional change but with increasing variance. Many models do not account for this! Analytical work often does not for sure. This is a good rationale for simulation.

Constructing the life cycle - it's okay (and often very realistic) to have non-reproductive juveniles and a single fecundity for adults.

Standardizing lambda across treatment groups - it was clever (good) of Schmid et al. to numerically find parameters to standardize the population growth rate. I can do this too with a simple life cycle. Or, I could just follow general rules in a reasonable way (e.g. increasing survival means decreasing fecundity...). Could also use a random search to find workable parameters, although Dan notes that this will occasionally produce non-sensical and unrealistic life histories.

##### Varying other vital rates than survival?

Be weary with fecundity. Dan has noted that fecundity increases with environmental change (!) in some of his plants (Silene?) The net effect is still to decrease population growth due to the low sensitivity of the growth rate to reproduction. But, this would complicate some simulations - have fecundity increase with environmental change? Could be worthwhile, accurate, interesting etc. but it complicates the study. There is enough to do with survival.

Survival and growth? Growth may not behave exactly like survival. Not necessarily as tricky as fecundity would be may have different effects than survival. [I'll note here that similarly to fecundity, environmental change could also be increasing growth...]

##### General

Dan said there's also an interesting question in here about what in an iid (or non-iid) environment favors longevity, or when longevity would not be favored. The model is also suited for that related but different question. I didn't quite understand what he meant and we were close to the end of the meeting so I did not ask. But it would be good to follow up more with him on that.

He is happy to look at an early draft of a chapter proposal! If he does not answer or immediately look at it I can badger him.

