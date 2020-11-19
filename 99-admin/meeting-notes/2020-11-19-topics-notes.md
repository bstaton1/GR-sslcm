# GR-SSLCM Meeting 11/19/2020

**Attendees**: Gibson, Justice, Liermann, Staton, White

_Post meeting notes are in italic_

## Progress Since Last Meeting

### Progress Update from C. Justice Re. Capacity Index Analysis

* C. Justice will step us through the work on building fish density vs. habitat metric models
* Recall the motivation behind this analysis is to replace the notion of "pool equivalent units" with an index of total suitable habitat that is more transparent to obtain and use in developing restoration/climate change scenarios with respect to freshwater habitat
* _Following the presentation of this work, the group discussed the merits and disadvantages of this way of producing a habitat capacity index, but did not discuss much on alternative approaches_
  * _Merits: (a) marginal effect plots are clear and align well with a priori expectations, (b) data are from the same region as well will be applying its estimates, (c) translating restoration/climate scenarios into changes in the metric will be relatively straightforward, and at least highly transparent_
  * _Disadvantages: (a) fish data are very noisy, and there is some concern this model will tell us what we want - i.e., data are a snapshot of density at the time of sampling, which could be influenced more by things like redd distribution than by the habitat metrics used as covariates, (b) AIC analyses were developed for use in relatively small model selection problems, less so when there are possibly 100s of models, (c) would it not be better to assume we know which variables are important, then allow the data to tell us the effect size alone, rather than telling us both which are important and their effects sizes (as T. Cooney did), (d) relationships are uncertain, and this does not include out-of-sample uncertainty -- in application, we will be largely (if not entirely) ignoring this uncertainty_
  * _There was some discussion of trying multiple different habitat indices from different analyses to see whether there is much of a difference among approaches. Examples include T. Cooney's PEU metric, M. Bond's metrics (not sure what these are), and K. See's QRF model. These would come at the stage of fitting the multi-population model to historical data, not for projection I believe._

### Feature: track dynamics of hatchery fish from release to adulthood

* PR #45
* Discussed during last meeting, agreed it should be merged
* Downside: no survival data incorporated (yet, as of the version of PR #45) for hatchery release to LGD. Survival to estuary for hatchery fish was from release, unlike for natural fish which have it as LGD to estuary
* _Not much discussion here, this was mostly a reminder that this was a major development that is not a baked-in feature of the model. PR #52 resolves the downside listed_

### Feature: fit to survival of hatchery fish from release to LGD

* PR #52
* Manner in which we calculate these survival data is still in question, see issue below
* However, we have some starting estimates we can use, and future development will need this component in place so this PR adds this step
* This PR also includes a bug fix on the calculation of SE(logit(surv)) from SE(surv)
* _Not much discussion here, this topic had already been discussed at length in [#9](https://github.com/gibsonpp/GR-sslcm-data/issues/9)_
* _Staton forgot to show the bug fix. Interested team members can view the change by navigating to commit 769100f1ab6b385502de0f8287650fd3e9baae3b (commit message: "bug fix: error in get_logit_se" found in PR #52)_

### Issue: method for aggregating survival of hatchery fish to LGD across multiple raceways

* GR-sslcm-data issue [#9](https://github.com/gibsonpp/GR-sslcm-data/issues/9)
* Survival estimates/SEs are available for each raceway
* But we want a total survival term
* How do we aggregate?
  * Two options: weighted average or Monte Carlo simulation
  * A comparison of the two methods is on a new branch in this repo
  * _Decision: we should stick with the weighted average approach. It is easy to explain, and the Monte Carlo approach illustrated that it did not under-represent the uncertainty in aggregate survival, which was the concern there_
* **Other topics related to this issue**
  * Are we okay with ignoring the rare occurrences where parr were released in the summer/fall?
    * _Originially, we were okay with this, but upon further discussion, we decided it would be best to account for it in a simple way. P. Gibson will upload a revised data set where the parr releases will be included with the smolt releases, and their survivals will be included in the weighted average_
  * Are we okay with ignoring smolt release data prior to brood year ~1995? 
    * _Yes, but perhaps we should calculate a rough estimate of how many adults would have survived and returned from these release events_

### Issue: some MCMC draws of the productivity parameter can be absurdly large

* No formal issue yet created
* Show diagnostic plot for CAT
* The prior is currently: `log_alpha ~ dnorm(0, 0.001); alpha <- exp(log_alpha)`
* Are we okay with using a truncated prior to exclude huge outcomes?
  * This parameter is the max summer parr produced per spawner (regardless of age/sex/origin, currently)
  * If we assume a female can carry 10,000 eggs (which is high), all spawners are female (not true), and all eggs survive to parr stage (not true), the highest this can be is 10,000. Are we okay with using this as the upper bound in the prior? E.g., `log_alpha ~ dnorm(0, 0.001) T(,log(10000))`?
  * I'm currently doing this for the capacity parameter as well - capacity can't be larger than `exp(15) = 3,269,017`. Are we okay with this as well?
  * _We are okay with imposing these constraints. In choosing their value, the justification for alpha (productivity) was well received, but it was agreed that for beta (capacity) the choice is less clear as to how to pick the constraint. The team likes the approach of seeing if the posterior is forced up against the bound, and if it is, then the bound needs to be increased. Side note: the current bound for beta does not seem to be influencing the posterior at all, with the exception of performing its intended purpose: excluding outrageously large values from being sampled_

## Next Steps

### Feature: fit to hydropower survival estimates

* Intent is to allow estimation of age-specific ocean survival and interannual variability
* From P. Gibson's work in issue #49, it seems we have multiple options for informing these hydropower survival estimates
  * They are LGD to BON. So first ocean survival term would be from BON to spring of following year
* Using the estimates from Table A.1 pasted in the comments on that issue, I attempted to fit a model that estimates origin- and age-specific time-varying ocean survival (constant across sexes)
  * There was an issue with hatchery fish where we were missing the survival from release to LGD. A fixed value of 0.6 was inserted to get around this at the time. PR #52 will fix that problem.
  * The model displayed MCMC sampling issues for the first time
  * R. Sharma (and B. Lessard) have warned that we will run into problems here because maturity and survival are confounded.
  * Ideas for things to try to remedy this:
    * Force ocean survival from age 3 to age 4 and from age 4 to age 5 to be the same. Would need to be careful with year/age indexing for random deviates: age 4 from brood year y and age 5 from brood year y-1 should be the same
    * Force time constant maturity schedules, and estimate time varying ocean survival terms?
    * (Exploratory only): vastly reduce uncertainty of hydrosystem survival - if we know with much higher certainty how many fish reach BON, does this help fix the issue? Trying to get at the cause of the problem here.
    * Treat age 3 maturation for females differently. Fix it at a very small value? Make it time constant (even if other age/sex/origin combos are time varying)?
    * OTHERS??
    * _The team seemed to agree that these are all plausible fixes and that they would be on board with using any of them (except the exploratory only one) should it reduce the confounding problem_
    * _M. Liermann suggested that others collapse the maturation and survival parameters into single parameters, and that this may help. B. Staton is a bit unclear about how this would work, and will follow up should the currently proposed solutions fail_
    * _However, M. Liermann will work up an "algebraic proof" to help illustrate whether some parameters are identifiable from the available data and model structure, similar in structure to the appendix from the original Liermann et al. manuscript that illustrated that the proportion of summer parr that become fall migrants is identifiable with the assumption that movement survival in spring from the tributary to LGD is the same between life history strategies. The hope is that this may help us structure a model that is identifiable_
* If we want to estimate ocean survival, we need to figure out:
  * The right hydropower survival estimates to use
  * How/if we should address the barging issue (estimates are for non-barged fish only)
  * A solution to the confounding of maturity and ocean survival
* _The team seems supportive of continuing down the path of trying to estimate both maturity and ocean survival_

### Next Meeting

* Thursday, December 17 at 10:00am