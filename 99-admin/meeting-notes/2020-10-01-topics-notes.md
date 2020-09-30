# GR-SSLCM Meeting 10/1/2020

**Attendees**:

## Progress Since Last Meeting

### Added age/sex/origin-specific broodstock removals

* PR [#28](https://github.com/bstaton1/GR-sslcm/pull/28)
  * Apportions estimate of total fish arriving at weir to age/sex/origin based on the observed proportions sampled at the weir (this corrects for partial detection of fish at the weir)
  * Divides the known number of fish removed in each age/sex/origin category by the apportioned total that arrived at the weir to obtain the proportion of each class that was removed.
  * Had we used just the individual fish data, this would over estimate `p_remove`, since it does not include fish that pass the weir but are not sampled. Hence the first step.

### Exploratory Analysis of Pre-spawn Mortality Data Completed

* Staton fitted exploratory models to investigate how complex the pre-spawn mortality component needs to be. Findings were:
  * Inter-annual variability is low for some populations and high for others (especially UGR) - indicating among year variability should be allowed in model
  * Among-population variability is moderate - indicating different populations should have different pre-spawn mortality terms
  * Lack of strong evidence for age or origin composition being important for explaining why some years have higher pre-spawn mortality than others.
  * Some evidence suggesting covariability in pre-spawn mortality among populations
* Results in Rmd file sent to team on 9/22/2020 - can go through document if there are questions

### Incorporated Estimation of Pre-Spawn Mortality in LCM

* PR [#31](https://github.com/bstaton1/GR-sslcm/pull/31): incorporated estimation of pre-spawn mortality terms.
* Uses same hierarchical structure as other survival components: across year mean and logit-normal sd estimated, year-specific parameters drawn as logit-normal random deviates.
* Used to decrement `Sb` (adults passing weir) to obtain `Sa` (adults on spawning grounds that spawn successfully). Apply a sex, origin, and age constant survival term.
* Estimated by fitting to carcass recovery data: `carcs_spawned[y] ~ dbin(phi_Sb_Sa[y], carcs_sampled[y])`
* Incorporated fit plots into the output summary Rmd file
  * Looks good for CAT and LOS, but for UGR, looks like the model is "chasing noise". E.g., 1995: very few fish sampled, but all were pre-spawn morts, and `phi_Sb_Sa[y]` that year is essentially zero
  * **Question**: should we have some criterion for including an observation in a year? E.g., only include pre-spawn mortality data for a year if at least 10 fish were sampled? Or some percentage of all spawners? This would be an attempt to prevent this issue.

### Incorporated Density-Dependent Overwinter Survival

* PR [#33](https://github.com/bstaton1/GR-sslcm/pull/33): added PEU estimates from the Liermann & Sharma model to serve as habitat index of capacity
  * Can be updated to include estimates from Justices recent analysis of other habitat features. That model is coming along -- Justice may wish to provide an update here?
* PR [#34](https://github.com/bstaton1/GR-sslcm/pull/34): added estimation of logistic model coefficients, structured exactly as in Liermann & Sharma model
  * Added plots to show pattern for each LH-type and population
  * Patterns look similar across populations: more fish at end of summer leads to lower over-winter survival, for each LH-type
* **Question**: random deviates around LH-specific deterministic logistic curves are independent. This could be a good place to bring in covariance for the first time: should we allow year-specific values to covary between LH-types? What about when multiple populations are added? Would we rather use among-population covariance here, or within-population but among LH-type covariance?

### Hatchery Release Data Added to [GR-sslcm-data](https://github.com/gibsonpp/GR-sslcm-data) Repo

* Gibson obtained and added these data.
* We can now discuss ways to model the dynamics of hatchery fish

## Next Steps: continue single-population model development

### Changes to any current components before moving ahead?

* Now would be a good time to investigate other model structures before adding more complexity
* Staton is okay with moving ahead at this time, with possible exception of covariance in overwinter survival

### Consistency of age/sex comp signal between weir and carcass recoveries

* If these are fairly consistent, we can use carcass comp in place of weir comp data where these are missing
* Will help inform us as to what to do for MIN age/sex data
* Some more details found in GR-sslcm-data issue [#4](https://github.com/gibsonpp/GR-sslcm-data/issues/4)

### Structure of Hatchery Fish Dynamics

* Ideally, we would model origin-specific:
  * Movement survival from tributary to LGD
  * Movement survival through hydrosystem
  * Proportion female upon arrival to estuary
  * Ocean mortality by age
  * Pr(mature-at-age/sex|not yet mature)
  * Upstream survival (e.g., different fishery exploitation rates on hatchery and wild fish)?
* Treat releases as perfectly known
* Fit to disaggregated age/sex/origin composition (currently aggregated across origin types)
  * Would help estimate sex comp parameters, maturation parameters, and ocean mortality parameters

### What other kinds of external survival data exist?

#### Hydrosystem survival (downstream)

* Can we fit GR-specific CJS models that move fish through hydrosystem? 
* Ideally, these would provide estimates for natural and hatchery fish separately.
* Or are we comfortable using aggregate estimates from another source?
* Fitting to these data *could* allow estimation of ocean mortality terms

#### Upstream survival

* Information on predation below Bonneville?
* Information on fish ladder successful passage rates?
* Information on origin-specific exploitation rates?