# GR-SSLCM Meeting 10/22/2020

**Attendees**: 

## Progress Since Last Meeting

### Bug Fix: weir composition data

* PR [#40](https://github.com/bstaton1/GR-sslcm/pull/40)
* Gibson pointed out a problem with calculation of weir age/sex/origin composition
  * Fish captured, sampled, and released downstream of weir then subsequently recaptured were being double-counted.
  * This PR corrected that problem

### Feature: pre-spawn survival sample size threshold

* PR [#41](https://github.com/bstaton1/GR-sslcm/pull/41)
* Recall issue of chasing noise in UGR pre-spawn carcass data
* Included ability to discard data if less than X fish sampled
* Didn't fix the cases I thought it would

### Feature: devise carcass vs. weir comp correction and fit to both data sets

* PR [#42](https://github.com/bstaton1/GR-sslcm/pull/42)
* Age/sex-specific correction factor. No age/sex interaction included, no correction for hatchery vs. natural fish either
* Allows fitting to both carcass and weir composition data while dealing with sampling bias of carcass surveys towards females/larger fish
* Display fit plots

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
* Should this PR be merged as is?

### Issue: hatchery adults in wrong years

* Issue [#44](https://github.com/bstaton1/GR-sslcm/issues/44): hatchery fish sampled as adults in years before process model says they can exist
* Also, hatchery fish returning to MIN, but process model won't allow this currently
* Approach currently used: throw out all hatchery adults in return years that aren't associated with a non-zero smolt release
* Do we need to model a straying process??

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

#### Upstream survival

* Can we dig up any information on:
  * Predation mortality
  * Fishery mortality

