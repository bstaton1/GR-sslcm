# GR-SS-LCM Meeting 9/11/2020

**Attendees**:

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

## Next Steps: continue single-population model development

> Would be really nice to (mostly) finalize dynamics/observation models with single population models before scaling up to multiple population models

### Decide on treatment of adult data

* Kept it simple for now, because MIN with no weir complicates things a bit
* Currently, fits to aggregate age/sex composition of hatchery + natural origin, and only from weir data
* Would be really nice to fit to hatchery and natural separately, and to fit to weir and carcass recovery data separately
  * Envisioning one multinomial likelihood within a year, each with 12 outcomes (2 sexes, 3, ages each, 2 origins): this would be for weir data
  * Would need to model hatchery releases and survival to adult-hood
  * Not sure about how to treat carcass data
  * Broodstock removals would need to be age/sex/origin-specific (currently only a single proportion applies)
  * Fitting to both arrival at weir and carcass data would provide (some of) the info for estimating pre-spawn mortality by age/sex/origin?

### Include density dependent overwinter survival

* Currently density independent
* Use logistic model approach, or Beverton-Holt approach?
* Will need habitat capacity index - start with the PEU variable until CJ's analysis is complete

### Decide on relative reproductive output

* Currently, all spawners are worth 1 unit of reproductive output. Should we search for fecundity measures? Hatchery and natural origin fish separate?

# Housekeeping

* Do SW and CJ wish to be added to the repo so they can access these materials? 
* Are other repo contributors getting emails about every issue, PR, merge commit, etc. that BS makes? If so, is this annoying or do you prefer to keep seeing this? There may be a way 

