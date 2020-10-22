# GR-SSLCM Meeting 10/22/2020

**Attendees**: Gibson, Liermann, Staton

_Post meeting notes in italics_

## Progress Since Last Meeting

### Bug Fix: weir composition data

* PR [#40](https://github.com/bstaton1/GR-sslcm/pull/40)
* Gibson pointed out a problem with calculation of weir age/sex/origin composition
  * Fish captured, sampled, and released downstream of weir then subsequently recaptured were being double-counted.
  * This PR corrected that problem
* _No noteworthy discussion here - agreement this issue needed to be fixed_

### Feature: pre-spawn survival sample size threshold

* PR [#41](https://github.com/bstaton1/GR-sslcm/pull/41)
* Recall issue of chasing noise in UGR pre-spawn carcass data
* Included ability to discard data if less than X fish sampled
* Didn't fix the cases I thought it would
* _Showed the comparative graphic found at the PR link, and team members had the idea to look at the juvenile data resulting from the spawners in these years. We discovered that the years with very low pre-spawn survival were also years of incredibly low trap abundance (at least at the fall trap). So the model is explaining the low  juvenile abundance with low pre-spawn survival, even if it doesn't have pre-spawn survival data. I would have expected the model to just put all of this into recruitment anomalies, rather than put it in both places. The group did not seem to think this is a problem, and even that it is an advantage of having such an integrated model. **Data issue arose**, however, the calculation of trap numbers does not account for fish that passed during times the trap was inoperable. So we should be viewing all trap counts as under-estimates of actual passage, though the extent varies year to year. We did not come up with a solution, but Gibson indicated that ODFW has plans to improve trap estimation in the future to account for this._

### Feature: devise carcass vs. weir comp correction and fit to both data sets

* PR [#42](https://github.com/bstaton1/GR-sslcm/pull/42)
* Age/sex-specific correction factor. No age/sex interaction included, no correction for hatchery vs. natural fish either
* Allows fitting to both carcass and weir composition data while dealing with sampling bias of carcass surveys towards females/larger fish
* Display fit plots
* _Quite a bit of discussion here, centered around the directionality of the estimated correction factors. For Lostine, the correction factors are exactly in the directions we would expect, older fish more likely to be sampled as carcasses relative to weir, and females more likely than males. Generally pattern holds for other populations, but for UGR especially, the pattern is a bit different. Some ideas were raised about how this may be due to differential broodstock removal strategies among populations. We still think the best approach for the full multi-population model is to use a single correction factor for all populations, **we will revisit this topic at that point if the fit is poor**._

### Feature: track hatchery fish from release to adulthood

* PR [#45](https://github.com/bstaton1/GR-sslcm/pull/45) (not yet merged, waiting for conclusion of this discussion)
* New components:
  * Added origin dimension to most objects in model
  * Survival of released smolt from release in tributary to estuary (independent from natural fish)
  * Sex-apportionment of hatchery fish upon arriving to estuary (independent from natural fish)
  * Maturation of hatchery fish at age 3 and age 4 by sex (independent from natural fish)
* Discarded components:
  * p_HOR expansion. No longer needed since there are now process equations to produce adult returns
* Assumes ocean survival and upstream survival are identical among natural and hatchery fish
* Model appears to fit as well as before, no new convergence issues introduced
* Now calculating SAR rates for hatchery and natural fish: gut check, do these seem realistic?
  * _Gibson indicated the SAR rates are on par with the official ODFW numbers_
* Should this PR be merged as is?
  * _Yes, the group agreed this is the right approach to use moving forward - the PR was merged_.

### Issue: hatchery adults in wrong years

* Issue [#44](https://github.com/bstaton1/GR-sslcm/issues/44): hatchery fish sampled as adults in years before process model says they can exist
* Also, hatchery fish returning to MIN, but process model won't allow this currently
* Approach currently used: throw out all hatchery adults in return years that aren't associated with a non-zero smolt release
* Do we need to model a straying process??
* _We talked at length about this, and are most concerned with the MIN issue. For MIN, hatchery fish are showing up in carcass surveys in most years, and there is no process to generate them -- we have to accommodate them somehow. Staton will post a comment on issue [#44](https://github.com/bstaton1/GR-sslcm/issues/44) with more details about what we talked about and where we "settled" (no great solutions as of yet)._

## Next Steps

### More hatchery vs. natural differences

#### Hydrosystem

* Would very much like to bring in external information on survival from LGD to estuary. Would be origin, year specific, common across populations
* Advantage: further partition survival among hatchery and natural fish
  * Would provide estimated survival from hatchery release to top of LGD
* Advantage: _should_ allow estimating ocean mortality
* Disadvantage: will need to fit new CJS models.
* Staton is starting some exploratory analyses via simulation to see what kind of precision we can expect in the estimates of total survival from LGD to estuary.
  * If estimates are too imprecise, it won't help much in estimating ocean survival
  * Using m-array approach rather than state-space expression (i.e., no individual effects) -- MUCH faster computationally, but a bit more involved JAGS modeling
  * One m-array for each tagging occasion: four m-arrays per year (summer, fall, spring natural tagging, plus hatchery releases)
  * Estimate population/year/origin-specific survival from tagging to LGD.
  * Estimate dam/origin/year-specific survival from dam `j` to dam `j+1` for all dams except the last two
  * Estimate dam/year-specific detection at dam `j` for all dams except the last one.
  * For the last survival/detection, this is confounded, so estimate the product.
  * Estimate survival to the last dam as `product/mean(detection[other_dams])`
  * Estimate survival through the last dam as `mean(survival[other_dams])`
  * Assuming 1000 natural fish per tagging event, 5000 tagged hatchery fish, per population. Realistic?
  * Assuming ~0.2 detection prob at each dam, and ~0.7 survival at each dam. Realistic?
* _The team agrees it would be ideal to bring in external hydrosystem survival to enable ocean survival estimation. However, we feel that the CJS route should be a bit of a last resort. We think that some useable information should exist out there, without having to fit GR-specific models. We do not doubt that it is possible to fit such models, but see it as a non-trivial endeavor (though the issue of barging intimidates us quite a bit) and that we could circumvent a lot of work and having to explain a lot of extra steps in the eventual MS dealing with this component if we could find existing estimates that will suit our needs. Staton will contact Bob Lessard (CRITFC) who has been involved with much mainstem survival work to see if he can point us in the right direction._

#### Upstream survival

* Can we dig up any information on:
  * Predation mortality
  * Fishery mortality
* _We are confident we can connect with other scientists who can provide information about how to structure these components of the model_
  * _Staton will contact Doug Hatch about predation estimates_
  * _Staton will contact CRITFC fish management about fishery mortality estimates_
  * _Liermann will contact Tom Cooney about these topics_
  * _Gibson will ask around ODFW to see if anyone has any insights_

# Post Meeting To-Do

_All of this content is post meeting, so no italics necessary_

## Admin

* Record notes, push to branch, merge to master
* Schedule next meeting: 11/12/2020 at 9:00am
* Merge PR [#45](https://github.com/bstaton1/GR-sslcm/pull/45)
* Add more context to issue [#44](https://github.com/bstaton1/GR-sslcm/issues/44)

## Model

* Figure out solution to issue [#44](https://github.com/bstaton1/GR-sslcm/issues/44)
* Once done, try fitting the MIN model for the first time

## Contact

* Bob Lessard about hydrosystem survival rates
* Doug Hatch about predation information
* CRITFC fisheries management about harvest

