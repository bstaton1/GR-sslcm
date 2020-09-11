# GR-SS-LCM Meeting 9/11/2020

**Attendees**: Gibson, Justice, Liermann, Staton, White

*Post meeting notes in italic*

## Progress Since Last Meeting

### Biological data housed in [GR-sslcm-data](<https://github.com/gibsonpp/GR-sslcm-data>)

* PG, please update group on reasons why we are storing data here
* Show how file organization should be handled in local repositories (to keep relative file paths functional)
* Intent is that all data (environmental, hatchery releases, etc.) will be placed here as well

### Infrastructure created for preparing data for analysis

* PR #[13](https://github.com/bstaton1/GR-sslcm/pull/13): take raw biological data files and generate one large data frame, similar to the `together.csv` file used by the original Liermann & Sharma LCM model. Data frame is called `bio_dat` in code, and is created with `source("00-data/prep-bio-data.R")`
  * Variable names **are not** consistent with the naming conventions used in model (by design)
* PR #[15](https://github.com/bstaton1/GR-sslcm/pull/15): take combined `bio_dat` object and prepare it into a list format for JAGS.
  * Function to create list for one population -- variable names **are** consistent with the naming conventions used in model
  * Function to take multiple one-population lists and combine them into one multi-population list
  * Function to instruct JAGS which elements are `NA` values and should be skipped in likelihood calculations -- this is a slick trick, show if we go through the model code
* Some issues remain
  * Only using weir data for composition, so code is not functional for MIN yet
  * Composition data and adult count data are combined across hatchery and wild fish, model accepts a `p_HOR` variable (assumed known w/o error or age/sex structure) to expand wild returning adults to this scale 
  * Will discuss in greater detail in next steps, below

### First revised version of model created

* PR #[19](https://github.com/bstaton1/GR-sslcm/pull/19); comments there include a narrative of how each life stage is handled currently. **CONTRIBUTORS -- PLEASE REVIEW PRIOR TO MEETING**
* **State-initialization approach**: estimate the first 5 brood year adult recruits, then use hyperparameters of ocean survival and maturity by sex to apportion them to the year/age/sex they would return. Open `simple-sim.xlsx` to illustrate how this works.
* Model code is found in `02-model/model-code.R`, and can currently be fitted to CAT, LOS, and UGR
* Source `02-model/fit-model.R` to fit the model (change the value of the `pop` variable to specify which population to fit to)
  * Writes a file (e.g., `02-model/model-output/CAT-output.rds`), which is ignored by git and has list elements for the input data and the output MCMC samples (in `mcmc.list` format)
* Step through model code? Could take a while.
* Naming conventions updated with slight edits after creating this model (PR #[20](https://github.com/bstaton1/GR-sslcm/pull/20))

### Output summary plots created

* PR #[22](https://github.com/bstaton1/GR-sslcm/pull/22): created `03-post-process/output-plots.Rmd`. Show an example output file
  * Creates MCMC diagnostic summaries (look pretty good at this point)
  * Creates plots showing fit to abundance, survival, and composition data - models fit data really well
  * Creates a plot showing Beverton-Holt relationship between spawners and parr, with an OLS fit for comparison
  * Creates boxplots showing posteriors of hyperparameters: across-year mean and SD for various survival parameters, maturity parameters, etc.
* When knitting, select "Knit with Parameters" to change the population to display
* **Any other information to include in these files at this point?**
  * *None at this time*

## Next Steps: continue single-population model development

> Would be really nice to (mostly) finalize dynamics/observation models with single population models before scaling up to multiple population models
>
> *Group agreed with this sentiment. Also, single population model infrastructure should be retained to allow investigation of what influence the hierarchical structure of the multi-population model has on population-specific estimates*

### Decide on treatment of adult data

* Kept it simple for now, because MIN with no weir complicates things a bit
* Currently, fits to aggregate age/sex composition of hatchery + natural origin, and only from weir data
* Would be really nice to fit to hatchery and natural separately, and to fit to weir and carcass recovery data separately
  * Envisioning one multinomial likelihood within a year, each with 12 outcomes (2 sexes, 3, ages each, 2 origins): this would be for weir data
  * Would need to model hatchery releases and survival to adult-hood
  * Not sure about how to treat carcass data - binomial for spawned/prespawned and estimate annual survival probabilities?
    * *Yes, unless there is strong info in the data to suggest prespawn mortality differs with age or origin*
  * Broodstock removals should be age/sex/origin-specific (currently only a single proportion applies)?
    * *Yes, to start out, can use weir data to calculate this proportion by age/sex/origin. However, because some fish pass the weir that are not sampled, this may overestimate the proportion broodstock removals*
  * There are some gaps in in weir age/sex/origin comp data - can I use carcass data to fill these in?
    * *Yes, unless there are strong differences in the composition signal in the raw data. Some concerns raised that age 3 males may be under-represented in carcass data and that females overall may be over-represented*

### Include density dependent overwinter survival

* Currently density independent
* Use logistic model approach, or Beverton-Holt approach?
  * *We favor the logistic approach, primarily because lognormal process error could produce more smolts in the spring than there were parr in the fall*
* Will need habitat capacity index - start with the PEU variable until CJ's analysis is complete

### Decide on relative reproductive output

* Currently, all spawners are worth 1 unit of reproductive output. Should we search for fecundity measures? Hatchery and natural origin fish separate?
  * *Yes, this should be included. Search for data sources to inform these. They should be age/size specific, and ideally for hatchery and natural fish separate. Ideally from GR basin, but don't have to be.*

# Housekeeping

* Do SW and CJ wish to be added to the repo so they can access these materials? 
  * *Yes, they wish to be added as it makes it easier to share the materials for meetings if they have access*
* Are other repo contributors getting emails about every issue, PR, merge commit, etc. that BS makes? If so, is this annoying or do you prefer to keep seeing this?
  * *No, other repo contributors are not getting bombarded by emails for every issue, PR, commit, etc. This is good. I assume that I can `@mention` someone if I want to get their attention*

# Post-Meeting To-Do

*All of this content was written after the meeting, italics omitted for convenience*

* Changes to model/data not dependent on additional information/exploratory analyses

  * Brood stock removals become age/sex/origin specific
  * Include PEU variable in GR-sslcm repo
  * Make overwinter survival density dependent

* Additional exploratory analyses

  * Do carcass age/sex/origin composition data generally agree with weir age/sex/origin composition data in years/populations that have both data sources? - Gibson volunteered
  * Is there evidence that prespawn mortality differs systematically by origin or age?

* Additional information sources

  * Information on fecundity - Gibson volunteered
  * Consistency in survival from tributary to LGD between natural and hatchery origin smolts. Justice had performed this analysis in 2013, and will provide an update to the group on how similar they are
  * Estimates of hatchery releases - forgot to ask someone to volunteer, should follow up with group.

* Administration

  * Add White and Justice to GR-sslcm and GR-sslcm-data repos

  * Commit, push, and merge the changes to this file

  * Create issues for future changes to model code/data as listed above

  * Schedule next meeting

    * 9:00am Thursday October 1st.

    

