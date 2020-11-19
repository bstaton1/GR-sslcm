# GR-SSLCM Meeting 11/19/2020

**Attendees**:

## Progress Since Last Meeting

### Progress Update from C. Justice Re. Capacity Index Analysis

* C. Justice will step us through the work on building fish density vs. habitat metric models
* Recall the motivation behind this analysis is to replace the notion of "pool equivalent units" with an index of total suitable habitat that is more transparent to obtain and use in developing restoration/climate change scenarios with respect to freshwater habitat

### Feature: track dynamics of hatchery fish from release to adulthood

* PR [#45](https://github.com/bstaton1/GR-sslcm/pull/45)
* Discussed during last meeting, agreed it should be merged
* Downside: no survival data (yet) for hatchery release to LGD. Survival to estuary for hatchery fish was from release, unlike for natural fish which have it as LGD to estuary

### Feature: fit to survival of hatchery fish from release to LGD

* PR [#52](https://github.com/bstaton1/GR-sslcm/pull/52)
* Manner in which we calculate these survival data is still in question, see issue below
* However, we have some starting estimates we can use, and future development will need this component in place so this PR adds this step
* This PR also includes a bug fix on the calculation of SE(logit(surv)) from SE(surv)

### Issue: method for aggregating survival of hatchery fish to LGD across multiple raceways

* GR-sslcm-data issue [#9](https://github.com/gibsonpp/GR-sslcm-data/issues/9)
* Survival estimates/SEs are available for each raceway
* But we want a total survival term
* How do we aggregate?
  * Two options: weighted average or Monte Carlo simulation
  * A comparison of the two methods is on a new branch in this repo
* **Other topics related to this issue**
  * Are we okay with ignoring the rare occurrences where parr were released in the summer/fall?
  * Are we okay with ignoring smolt release data prior to brood year ~1995? 

### Issue: some MCMC draws of the productivity parameter can be absurdly large

* No formal issue yet created
* Show diagnostic plot for CAT
* The prior is currently: `log_alpha ~ dnorm(0, 0.001); alpha <- exp(log_alpha)`
* Are we okay with using a truncated prior to exclude huge outcomes?
  * This parameter is the max summer parr produced per spawner (regardless of age/sex/origin, currently)
  * If we assume a female can carry 10,000 eggs (which is high), all spawners are female (not true), and all eggs survive to parr stage (not true), the highest this can be is 10,000. Are we okay with using this as the upper bound in the prior? E.g., `log_alpha ~ dnorm(0, 0.001) T(,log(10000))`?
  * I'm currently doing this for the capacity parameter as well - capacity can't be larger than `exp(15) = 3,269,017`. Are we okay with this as well?

## Next Steps

### Feature: fit to hydropower survival estimates

* Intent is to allow estimation of age-specific ocean survival and interannual variability
* From P. Gibson's work in issue [#49](https://github.com/bstaton1/GR-sslcm/issues/49), it seems we have multiple options for informing these hydropower survival estimates
  * They are LGD to BON. So first ocean survival term would be from BON to spring of following year
* Using the estimates from Table A.1 pasted in the comments on that issue, I attempted to fit a model that estimates origin- and age-specific time-varying ocean survival (constant across sexes)
  * There was an issue with hatchery fish where we were missing the survival from release to LGD. A fixed value of 0.6 was inserted to get around this at the time. PR [#52](https://github.com/bstaton1/GR-sslcm/pull/52) will fix that problem.
  * The model displayed MCMC sampling issues for the first time
  * R. Sharma (and B. Lessard) have warned that we will run into problems here because maturity and survival are confounded.
  * Ideas for things to try to remedy this:
    * Force ocean survival from age 3 to age 4 and from age 4 to age 5 to be the same. Would need to be careful with year/age indexing for random deviates: age 4 from brood year y and age 5 from brood year y-1 should be the same
    * Force time constant maturity schedules, and estimate time varying ocean survival terms?
    * (Exploratory only): vastly reduce uncertainty of hydrosystem survival - if we know with much higher certainty how many fish reach BON, does this help fix the issue? Trying to get at the cause of the problem here.
    * Treat age 3 maturation for females differently. Fix it at a very small value? Make it time constant (even if other age/sex/origin combos are time varying)?
    * OTHERS??
* If we want to estimate ocean survival, we need to figure out:
  * The right hydropower survival estimates to use
  * How/if we should address the barging issue (estimates are for non-barged fish only)
  * A solution to the confounding of maturity and ocean survival