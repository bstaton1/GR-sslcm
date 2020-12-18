# GR-SSLCM Meeting 12/17/2020

**Attendees**: Gibson, Justice, Liermann, Staton

_Post-meeting notes in italics_

## Progress Since Last Meeting

* Unfortunately not much to report here, Staton has been pulled to working on other CRITFC projects
* PR [#55](https://github.com/bstaton1/GR-sslcm/pull/55): Introduced a truncated prior on productivity to prevent large draws
  * _No comments here, group seemed to agree this was necessary and well-justified. Some discussion of the truncation on the capacity parameter was had as well, we may revisit the level of truncation used to see if inference changes._
* PR [#57](https://github.com/bstaton1/GR-sslcm/pull/57): Improved model syntax. Rather than hard-coding `[,1]` when not looping over that dimension, we now have a name for what the 1 means (e.g., `o_nat`, `s_female`, etc.). Makes code easier to read and edit. See issue [#56](https://github.com/bstaton1/GR-sslcm/issues/56) for more details and examples
  * _The group was supportive of this change and saw how it improves model readability_
* Majority of Staton's LCM work since last meeting has been focused on ways to estimate hydrosystem survival, ocean survival, and maturity jointly. See below.

## Unresolved Topics

### How to estimate hydrosystem survival, ocean survival, and maturity?

* We believe we can inform hydrosystem survival by fitting to estimates obtained in other studies.
* Using one candidate estimate set, this now allows us to explore ways to estimate ocean survival in addition to maturity
* However, in a full model that estimates all survival and all maturity terms, these processes are confounded. I.e., there are multiple combinations of parameters that could give the same number of fish returning for observation. Other terms: they are non-separable, non-identifiable
* "Full model" -- everything is modeled hierarchically:
  * Year-specific and origin-specific sex apportionment parameters upon ocean entry
  * Year-specific, sex-specific, origin-specific, and age-specific maturation parameters
  * Year-specific, origin-specific, and age-specific survival parameters
* This model has significant identifiability/convergence problems. 
* **GOAL**: find a reasonable simplification that allows estimation of ocean survival and maturity. 
* _The group understood the issue here. Sentiments were raised about how we don't need to model detailed biological processes in the ocean, our model is focused on freshwater dynamics primarily. We just need something that will move fish from Bonneville as juveniles back to Bonneville as adults at the appropriate age, sex, origin and calendar year. However, it was generally agreed that the better we can do with estimating ocean survival terms, the more compelling of a story we will have in terms of the benefits of integrating all available data sets into one state-space model like this._

#### Models Tried

Each model has an associated RMD output file we can look at to see the problems. 

##### Full Model

* Significant convergence issues
* Some parameters (survival & maturity) highly uncertain

##### Alt 1: Time Constant Maturity

* Excessively slow mixing for maturity parameters
* Fits age composition data poorly
* Unrealistic ocean survival time series

##### Alt 2: Small Hydrosystem Survival SE (EXPLORATORY ONLY)

* Didn't really resolve issues of full model

##### Alt 3: Same O1 -> O2 survival as O2 -> O3

* Resolves some convergence/identifability issues
* Best improvements over full model so far
* I think some theoretical justifiation exists
* Not currently accounting for relationships between brood years - should this be done? E.g., Have a parameter that influences survival from O1 -> O2 for brood year y and from O2 -> O3 for brood year y-1
* _Group thought this simplification was appropriate, but wasn't sure if the brood year lag is necessary or how it would be included. Staton will try to come up with something if we end up with this simplification, otherwise keeping it without treating brood years specially should be okay_

##### Alt 4: Combine Parameters

* Liermann's a, b, c approach (show his markdown doc)
* Major convergence problems for _a_ parameter, particularly for hatchery females
* Model suggests essentially zero hatchery smolt are female?
* _The group was a bit puzzled about why this simplification doesn't work better, but it is almost certainly due to the rarity of age 3 females, particularly for hatchery fish. Staton will try another simplification: where rather than having a = s1\*m1, b=s2\*m2(1 - m1)/m1, and Return1 = Smolts\*a, Return2 = Smolts\*a\*b, etc., we just have a1, a2, a3, where Return1 = Smolts\*a1, Return2 = Smolts\*a2, Return3 = Smolts\*a3. The advantage of this kind of composite is that it can still "be" the full model, but we don't have to estimate all of the terms. The disadvantage is that it is harder to know how much interannual variability is due to maturity and how much is due to ocean survival, thus, if we wish to manipulate ocean survival in forward simulations, coming up with multiplication factors will be more difficult than if we had explicit terms for maturity and survival_

##### Alt 5: Estimate Natural Origin Survival parameters, estimate an adjustment factor to obtain hatchery

* Identifiability problems when trying to estimate age-specific scalers
* Issue reduced with fewer scalers
* _The group liked this approach overall as a candidate simplification_

##### Alt 6: Same as Alt 5, but also introduce Alt 3 simplification

* This one seems most promising
* Survival from O1 -> O2 is the same as O2 - > O3 (still not accounting for brood year indexing)
* Estimate one log-odds ratio that changes the natural origin survival into hatchery origin survival in all years
* _The group agreed this simplification produced the most sensical results of those shown, and that if we move forward with one of those shown here, it should be this one_

## Next Steps

* Decide on final structure for ocean survival/maturity, implement it
  * _Staton has a few more alternative structures to evaluate, but if none look better than Alt 6 above, this will be merged into master_
* Decide on best hydrosystem survival data source
  * _If we recognize that our ocean survival terms are more indicies than rigorous estimates (which we will **need** to if we use a composite approach, and to a lesser extent with a large simplification off of the full model, such as Alt 6 above), then the exact values we use here are of lesser importance. The group seemed comfortable with moving ahead with this data set, but we should follow up with the experts on whether they deem our treatment of the data to be largely appropriate or if we are making a gross error somewhere_
* Track down fishery mortality estimates
  * _The group indicated this is of high importance, especially considering hatchery and natural fish are harvested at different rates likely. Staton will follow up looking for harvest information_
* Should we weight spawners of different sexes, ages, or origins differently in reproduction?
  * _Gibson indicated ODFW has a large data set for fecundity (separated for hatchery and natural origin fish, and Grande Ronde specific fish), which is very encouraging. We will likely turn the reproductive function (spawners -> parr) into spawners_at_age_sex_origin\*fecundity_at_age_sex_origin -> parr_

After these issues are resolved, I believe it is time to integrate the four populations into one model

Grande Ronde Model Watershed has requested an update on the life cycle model. Two possible dates:

* Jan 21, 2021
* Feb 18, 2021
* Staton is happy to develop the slides and present them to this group
  * _The group is okay with suggesting the Feb 18th date, and with having Staton as the lead presenter. Other collaborators should be present for the call if at all possible to help answer questions, and everyone will be included as authors on the slides. The group thinks a high-level talk will be best for this audience, i.e., they won't care too much about modeling details._

**Next Meeting Date/Time**: _January 21, 2021 at 10:00am_

## Post-meeting to-do

_All of these are post-meeting notes, so italics are omitted_

* Submit pull request to `GR-sslcm-data` repo containing the hydrosystem survival estimates Staton has been using
* Try a handful of alternatives to the "full model"
* Contact folks about information on fishery rates
* Make an attempt at including finer-scale upstream survival dynamics (i.e., mark-selective fishery)
* Commit and pull request this branch into the master, close related issue
* Rebase a "final" ocean survival branch on top of the new master
* Merge this "final" ocean survival branch into master
* Schedule and send invite for next meeting