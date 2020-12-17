# GR-SSLCM Meeting 12/17/2020

**Attendees**: 

## Progress Since Last Meeting

* Unfortunately not much to report here, Staton has been pulled to working on other CRITFC projects
* PR [#55](https://github.com/bstaton1/GR-sslcm/pull/55): Introduced a truncated prior on productivity to prevent large draws
* PR [#57](https://github.com/bstaton1/GR-sslcm/pull/57): Improved model syntax. Rather than hard-coding `[,1]` when not looping over that dimension, we now have a name for what the 1 means (e.g., `o_nat`, `s_female`, etc.). Makes code easier to read and edit. See issue [#56](https://github.com/bstaton1/GR-sslcm/issues/56) for more details and examples
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

##### Alt 4: Combine Parameters

* Liermann's a, b, c approach (show his markdown doc)
* Major convergence problems for _a_ parameter, particularly for hatchery females
* Model suggests essentially zero hatchery smolt are female?

##### Alt 5: Estimate Natural Origin Survival parameters, estimate an adjustment factor to obtain hatchery

* Identifiability problems when trying to estimate age-specific scalers
* Issue reduced with fewer scalers

##### Alt 6: Same as Alt 5, but also introduce Alt 3 simplification

* This one seems most promising
* Survival from O1 -> O2 is the same as O2 - > O3 (still not accounting for brood year indexing)
* Estimate one log-odds ratio that changes the natural origin survival into hatchery origin survival in all years

## Next Steps

* Decide on final structure for ocean survival/maturity, implement it
* Decide on best hydrosystem survival data source
* Track down fishery mortality estimates
* Should we weight spawners of different sexes, ages, or origins differently in reproduction?

After these issues are resolved, I believe it is time to integrate the four populations into one model

Grande Ronde Model Watershed has requested an update on the life cycle model. Two possible dates:

	* Jan 21, 2021
	* Feb 18, 2021
	* Staton is happy to develop the slides and present them to this group

**Next Meeting Date/Time**:

