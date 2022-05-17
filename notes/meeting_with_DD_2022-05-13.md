# Meeting notes

### Scott, Dan Doak, May 13 2022 (zoom)

Purpose of this meeting was to to talk about what to do with dynamic genetic variance.

Both to talk about finding a steady state given that the variance changes over time (to different degrees among groups) and also to think about the Lynch and Walsh result that breeding value variance and additive genetic variance are identical.

##### Analytical approaches

Without me providing much detail, Dan said that the approach I was using thinking about cohorts made sense to him.
Over email he said that I could use a discretized approach and borrow from matrix model approaches for SSDs (adaptation is movement among classes).

The analytical approach would only work for a stable environment.
A varying environment and/or density dependence would likely require simulations.

He recommended paper(s) by Benton and Grant for NDD and stage distributions.
He also said that for this case, it might be acceptable to model NDD as a ceiling function.

##### Notes on the trade-off

Dan noted (I think correctly) that the mutation rate should probably vary among groups as part of the trade-off.

Dan also made me think that I could model the trade-off as part of a single trait? 
Have fecundity vary by individuals in a way that also includes the trade-off
(survival is high, fecundity is low, e.g.).
This would be a good way to look at longevity as an adaptation rather than just comparing growth rates/extinction rates.
Dan noted that doing this would require the trade-off to be very strong to see an effect!

##### Lynch and Walsh

How does this result work? Dan encouraged me to go back and re-read the notes - maybe even to send the chapter to him.
What are hte assumptions? Is selection occurring? Are there even fitness differences?

My thoughts during/within the meeting: because selection is not acting on fecundity, and selection has already occurred through survival, I could maybe still use the L&W result for segregational variance.

Dan said this did not strike him as intuitive and if he was a reviewer on a MS like this, he would be skeptical. But I do think it's a good idea.

##### Other stuff

There was something about variance in lifetime fitness vs. phenotypic variance. I don't quite remember what.

Dan said it would be good to perhaps control for *variation* in lifetime fitness. 
I have mean fitness the same right now across groups, but by design not variance.
Should I control for variance even if it means the mean variance is the same?
I think I am accounting for the effects of phenotypic variance in population growth, but not the among-group variance...
Hmm...

Dan will be available after June 10.
