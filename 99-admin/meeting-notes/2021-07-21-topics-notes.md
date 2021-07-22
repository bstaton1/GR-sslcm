# GR-SSLCM Meeting 7/21/2021

**Attendees**: Gibson, Justice, Liermann, Staton, White

_Post-meeting notes in italics_

## 1. Items to Finalize Before Integrating Populations

### 1.1 Adult Migration

#### 1.1.1 Sea Lion Predation

* Discussion with D. Hatch (CRITFC)
* See Rmd analysis and issue [#69](https://github.com/bstaton1/GR-sslcm/issues/69) for details
* Good to proceed with this strategy?
* _The team agreed that this approach seems reasonable and had no recommendations for changes_

#### 1.1.2 Harvest Mortality

* Probably best to model harvest in at least two areas: below Bonneville (mark-selective) and Zone-6 (non-mark selective)
* There is reason to think age-3 fish harvested at lower rates; could use a fixed (expert opinion) multiplier to represent selectivity
* S. Ellis and K. Self (CRITFC) have information on this, annual rates are ~5-7% in each location, waiting on exact numbers
  * Not yet sure if this will be time varying (would use blocks like for sea lion mortality) or time constant
  * _Team had no recommendations for changing this aspect, other than the general comment below_

#### 1.1.3 "Other" Mortality

* Things like dam passage failures, heat exhaustion, etc.
* S. Ellis and K. Self (CRITFC) have information on this
* Will almost certainly be time constant

**Go over Rmd analysis that lays out plan and shows results of incorporation using rough numbers**

_General comment on adult migration from BON -> LGR: rather than disaggregating these various sources of mortality, could we not just use one term for "converting" fish at BON to fish above LGR, as informed by PIT tag data? We aren't that interested in evaluating changes in these rates for the simulation analysis, so what is the utility in splitting them apart? These are good questions and we should explore going this route. Gibson will ask around to see if any one knows about pre-summarized survival information, and Staton will bring this up to Ellis and Self when we speak next. If these data are not already summarized, we are not interested in performing a new CJS analysis to estimate this survival. **This branch will be merged into master for now (which uses the strategy presented in the Rmd), recognizing that it may change with more information.**_ 

_Further, there was some talk of in-tributary tribal harvest of fish after they are counted at the weirs -- this would not be captured under this approach, and we may want to look more deeply into how many fish are harvested at this stage to see if it should be accounted for._

### 1.2 Brood Stock Removal Calculation

* See issue [#71](https://github.com/bstaton1/GR-sslcm/issues/71) for details
* Currently, we use a proportion removed, which we estimate from weir data as n_pass/n_arrive by age/sex/origin. In the model, we take total returns times this proportion to get adults past weir
* But because total returns are stochastic in the model, the number of fish removed for broodstock will also be stochastic, which doesn't make sense given we know exactly how many fish were removed.
* It makes more sense to subtract fish removed for broodstock. To prevent possibility of negative spawners, use a `max()` constraint.
* I have this change sitting on a branch ready to merge. The `max()` constraint gets triggered very rarely and for only minor offenses (see update comments on [#71](https://github.com/bstaton1/GR-sslcm/issues/71) thread). **Merge it?**
* _**The team liked this change and had no recommended alterations prior to merging, thus it will be merged as-is**. Liermann raised an idea of using the weir removals more explicitly in the observation model as a means to force the process model to return at least as many adults by age/sex/origin as were removed at the weir. However, Staton is not sure how this would be implemented, and given the `max()` constraint is triggered very rarely and for minor discrepancies (several fish), we decided the proposed approach is sufficient._
* _The team agreed that we will need to devise a removal control rule for the simulation model. Gibson indicated that there are some written rules for this that we can lean on, and the EcoLogical model also had a somewhat sophisticated rule built in, so we may be able to do something similar._
* _Gibson raised the point that not all of the fish removed at the weir are used for broodstock (particularly jacks), so it would be best to start referring to this component as "weir removals"._

### 1.3 Autocorrelation

* Staton conducted an analysis looking for temporal autocorrelation in the input data time series and the process model residuals. See Rmd analysis.
  * Current model assumes no autocorrelation
* Strong evidence of autocorrelation in the raw data time series
  * But not all data series are linked explicitly to a process model component
  * Better to look for autocorrelation in the model residuals
* Less autocorrelation found in model residuals, but some processes consistently have reasonably high evidence  of autocorrelation (Pr(AR1 coef > 0) > 0.8) for multiple populations (at least 2/4):
  * Parr recruitment
  * Smolt migration survival from tributary to LGR (HOR & NOR)
  * Proportion of parr that are fall migrants
  * Hydropower survival (HOR & NOR)
  * First year ocean survival
* Should some of this be modeled with environmental covariates?
  * Looked at NOAA stoplight data set for correlations with 1st year ocean survival
  * PDO had moderately high correlation (~-0.5)
  * But might not be the most useful for forecasting future conditions -- PDO has temperature trends removed, only captures cyclical temperature fluctuations
  * Better to model time dependencies using AR(1) process here?
  * _The team agreed that using AR(1) for Yr1 ocean survival is the right approach rather than using PDO_
* Go over Rmd results of including AR(1) on model components that have evidence of autocorrelation
  * Generally, mean parameters became more uncertain with inclusion of AR(1) -- if autocorrelation is really present, then this may be a better characterization of uncertainty
  * Recruitment relationships became much more uncertain, driven mostly by greater uncertainty in capacity
  * How to proceed? Proposal:
    * Include AR(1) dynamics on all of these components except recruitment for now. Revisit recruitment autocorrelation after population integration
  * _The team seemed to agree that adding AR(1) to the recruitment process and making capacity estimation more difficult should be avoided at this time. Liermann raised the point that perhaps capacity becomes much more uncertain with the inclusion of AR(1) because it makes it more difficult for the model to tease apart density dependence from autocorrrelation._
  * _We also discussed the benefits to incorporating AR(1) complexity, and that what really matters is that we capture short term autocorrelation in adult abundance. If we can accomplish this without adding AR(1) to every model component that shows some historical evidence of it, then that is better than introducing this additional complexity that is difficult to estimate. **At this time, we decided to include AR(1) on Yr1 ocean survival and leave it off the other components for now. After integration, we will re-evaluate the evidence for its presence in model residuals and whether near-term projections capture autocorrelation in the returning adult abundance.**_

### 1.4 Fecundity

* Still wish to express parr recruitment as a function of eggs? (rather than total spawners)
  * Convenient way to reflect (possibly changes in) relative reproductive output by year, sex, age, origin
  * Which variability should we account for?
* Gibson put together analysis of heterogenity in fecundity with primary findings:
  * Definitive evidence of age effect
  * Some evidence for population and origin variability
  * Some evidence for interannual variability -- but does not appear to be trending
* Proposal: account for age effect only. Why not others?
  * Most important explanatory covariate based on Gibson's analysis
  * Including year means we need a value for every year -- imputations and projections
  * Including population means we need to impute values for Minam
  * Simpler inference. If the there were really strong patterns in these factors I could see including them, but based on the analysis, doesn't seem warranted.
  * _The team agreed with this reasoning and agreed that we should account for variability due to age only. **This change will be added to the model before integrating populations into one model.**_

### 1.5 Habitat Metric

* Still using PEU variable
* Is the revised weighted usable habitat metric ready to be inserted?
* _There is still some uncertainty about what to do for the Lostine habitat data, since only about half of the historical Chinook extent has AQI information. We will continue using PEU values until we decide on how to handle this issue._
* _There was some discussion about incorporating the Imnaha population into our model. Gibson indicated that it would be possible to compile much of the same data that other populations use, but would be more difficult since those data are collected by Nez Perce and not ODFW. However, no team members were aware of AQI data in that tributary that could inform its capacity index, and no samples from that tributary are available to inform the habitat capacity regression model Justice developed. For these reasons, we plan to leave out Imnaha and proceed with the four populations captured by T. Cooney's work and the Liermann & Sharma model._

## 2. Miscellaneous

### 2.1 GRMW Meeting

- November 16-19, 2021 GR State of the Science meeting
  - Nov. 16 - Devoted to LCM meeting
  - Nov. 17 - Field trips
  - Nov. 18/19 - SOS
- One day in this window will be reserved for presentations about LCM and a discussion of simulation scenarios
- Presentations:
  - Overview of model and places where future scenarios could be reflected -- Staton leads + whole team contributes
  - Overview of revised weighted usable habitat metric and ways different restoration scenarios could be reflected -- Justice leads + whole team contributes
  - Overview of target values analysis -- White leads (maybe mixed in with Justice presentation)

### 2.2 Model Diagram

- Show diagram
- Intent I had in mind: to show main life stages we account for, their structure, and which components are informed by data
- Strengths - accomplishes the intent (?)
- Weaknesses - kind of busy, too detailed for a PPT slide, doesn't show functional form of relationships, doesn't show which processes have shared parameters across populations, upstream migration needs more detail, no distinction between data we fit to with uncertainty (e.g., juvenile survival) and which we assume known (e.g., sea lion mortality, hatchery smolt releases)
- _The team had positive feedback on the diagram but requested more time to look it over. Staton will create an issue with the image for others to post feedback on. Justice raised the idea that it would be good to show where habitat information feeds in to the model._

## 3. Integration of Populations

**MAKE SURE TO TALK ABOUT 3.2 AT LEAST IF THERE IS NOT TIME FOR THE OTHER TOPICS**

### 3.1 Logistics

* Will be a big change to model code: essentially add `j` dimension for population to all quantities
* From a repo perspective, how would we like to retain the single-population version?
  * Option #1: two files tracked at all times; keep single population model and in a new file develop multi-population model
  * Option #2: make commits that alter the current single pop model; implies this will be the only version to run from here on -- single-pop model will live in the commit history
  * Will we want to re-run the single pop model again, or was this just a (long) developmental step?
  * If we wish to keep both versions, I'd like the only major difference between versions to be that one is multipopulation
  * _Although the team thinks it unlikely that we will want to run the single population model after the multi-population model is working, they think it would be unwise to discard that model prematurely. We will use Option #1 for simplicity, but may also create a branch that serves as the most complex working version of the single population model for us to revisit if we wish._

### 3.2 Hierarchical Structure

* See [#80](https://github.com/bstaton1/GR-sslcm/issues/80)
* Liermann & Sharma model had two levels: across-population mean -> population-specific mean -> population- and year-specific value
* Would like to simplify this by dropping top level if possible. Thoughts?
* _The team agreed that the utility of the top-level is limited for this case and that it can be removed. The question came up about how difficult it would be to add in for future use cases (i.e., other basins), and overall we think that relative to the other changes that may need to be made to the model for it to be applicable to another basin, that adding this back in wouldn't be a major hurdle._
* _Gibson raised the question that if we are not using the top-level, what is the utility in building a multi-population model. The answer is that there are other parts of the model where it may be more useful to share information, such as using identical ocean survival parameters, hydrosystem survival parameters, etc. for all populations. An additional benefit is the ability to capture correlated dynamics among populations._

### 3.3 Covariance Structure

* There is evidence of covariance among populations -- see Rmd on topic.
* Can be modeled using multivariate (logit/log) normal distribution.
  * Must estimate covariance matrix -- not the easiest thing to do
  * Complexity depends on if we allow all population pairs have unique correlations, or if it is okay to assume all populations have the same correlation.
* _We did not have time to discuss this topic with any detail. We will revisit it after Staton builds an early version of the multi-population model._

## Next Steps

_All of this is post-meeting content, so the italics are omitted._

* Rebase and merge existing branches (with some edits in parentheses):
  * These meeting notes
  * New weir removal calculations
  * AR(1) on ocean survival (only, remove the AR(1) components that are included for other life stages)
  * Adult migration submodel
* Add fecundity and alter prior on the `alpha` parameter to reflect new definition - parr recruits per egg rather than per spawner
* Create issue for gathering feedback on model diagram
* Handle existing "quality-of-life" improvements in documented issues:
  * [#78](https://github.com/bstaton1/GR-sslcm/issues/78), [#79](https://github.com/bstaton1/GR-sslcm/issues/79)
* Build an early version of the multi-population model
  * Use common hydropower survival terms, common ocean survival terms, but unique sex apportionment and maturity schedules for now.
  * Omit any covariance structure for now
* Schedule meeting to go over integration and next steps in development




