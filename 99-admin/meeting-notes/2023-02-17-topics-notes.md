# GR-sslcm Meeting 2023-02-17

**Attendees:** Brandt, Feldhaus, Gibson, Justice, Staton

*Post-meeting notes in italics*

```diff
-Changes or potential changes to the model in red
+Questions for ODFW to consider internally in green
```

## Meeting Purpose

The purpose of the meeting was to:

1.  Bring in ODFW staff from the Lower Snake River Compensation Plan project (“LS” hereafter) who collect and manage the data for adult returns and for hatchery juveniles (Feldhaus and Brandt)
2.  Review how those data sets are currently being used in the model
3.  Get feedback from the LS team on the model in general

## 1 Adult Returns

### 1.1 Adult abundance

- Basic data value used here is total return to river \[RTR\] for all fish combined \[all ages (3-5), sexes, and origins\], plus measures of uncertainty for each RTR value.
    
    ```diff
    +ODFW team will work separately on the uncertainty estimates for adult abundance.
    ```
    
- *LS team notes that initial abundance estimates for a given return year continue to be tweaked for several years; but, Feldhaus reports that the magnitude of the estimate is unlikely to change much, and he is likely comfortable using the existing estimates through return year 2022.*
    

### 1.2 Weir Removals

- Fish are sampled for age and origin composition at the weir, then apply weir removals
    
- Assume removals are known without error, subtract from total abundance (origin/age-specific)
    
    - Issue: [#71](https://github.com/bstaton1/GR-sslcm/issues/71)
- Data source: weir database records
    
    - All Lostine weir data from 2001-2008 is excluded
    - Replacing weir records for broodstock fish with hatchery spawning records proved not possible
    
    ```diff
    +ODFW: how substantial were the removals likely to have been in these years? Should we try to impute missing values, rather than treating them as zero?
    ```
    
- *Feldhaus notes that for UGR, where only 50% of hatchery-origin fish are ad-clipped, there is a pattern that a consistent proportion (~20-30%) of the fish that are assigned natural origin are, in fact, hatchery-origin (based on genetic PBT data). This discrepancy is not currently accounted for in any of the existing population estimates in LS or VSP sources.*
    
    - *Return year 2014, in particular, has an unusually low proportion hatchery origin fish at the UGR weir -- maybe due to numerous actual hatchery origin fish being mis-assigned as natural origin? Staton notes also: 2014 CAT and 2013 LOS*
    - *Whatever "correction factors" are obtained here, should they also apply to the weir composition counts, not just the weir removals?*
    
    ```diff
    +ODFW team will review available data for estimating a standard correction factor to apply to origin assignments
    -This is a potential minor change to the model -- Staton prefers to perform these adjustments to the data outside of the model if possible. Ideally this would be provided to Staton as a matrix storing the proportion of the fish assigned each age/origin that should have been assigned another age/origin.
    ```
    
- *Here is an example matrix Staton describes above. Each row shows the fraction of fish assigned a given age/origin that are actually of each age/origin. The example below shows the hypothetical case where 25% of all fish assigned as NOR fish are actually HOR fish. In the example, all fish assigned as HOR are actually HOR fish. There are no aging errors.*
    
    |     | NOR-3 | NOR-4 | NOR-5 | HOR-3 | HOR-4 | HOR-5 |
    | --- | --- | --- | --- | --- | --- | --- |
    | **NOR-3** | 0.75 | 0   | 0   | 0.25 | 0   | 0   |
    | **NOR-4** | 0   | 0.75 | 0   | 0   | 0.25 | 0   |
    | **NOR-5** | 0   | 0   | 0.75 | 0   | 0   | 0.25 |
    | **HOR-3** | 0   | 0   | 0   | 1   | 0   | 0   |
    | **HOR-4** | 0   | 0   | 0   | 0   | 1   | 0   |
    | **HOR-5** | 0   | 0   | 0   | 0   | 0   | 1   |
    
- *Weir efficiency is variable across years and populations. UGR weir typically comes out early, thus potentially misses the late part of the adult run. Lostine weir typically goes in late, thus potentially missing the early part of the adult run. Jacks tend to arrive earlier than age-4 and age-5, so there is a potential for some systematic bias in age sampling at the weirs. However, on the whole, LS team judges that this is not enough of a problem to worry about for the LCM.*
    
    ```diff
    +ODFW team will also review disposition codes and how recaptures and outplants in weir removal data are being accounted for.
    +We should revisit LOS weir removal data for the years in which these data are missing. Could we use simple averages of years that do have good data? The model currently treats these years as zero removals which could be problematic if the removals were actually somewhat substantial.
    ```
    

### 1.3 Prespawn Mortality

- Estimated from carcass recoveries
    
    - Issue: [#30](https://github.com/bstaton1/GR-sslcm/issues/30)
    - See also: `GR-sslcm-data/exploratory-analyses/explore-prespawn-data.html` (emailed 2022-09-22)
    - Currently no threshold for minimum number of carcasses required in order to use these data
        - Issue: [#38](https://github.com/bstaton1/GR-sslcm/issues/38)
        - Pull Request: [#41](https://github.com/bstaton1/GR-sslcm/pull/41)
- *Justice mentioned an idea from M. Kaylor that we could potentially calculate a spearate estimate of prespawn mortality by comparing mark-recapture estimates of above-weir spawner abundance against above -weir redd counts. This would not work in the UGR population (due to problems including low sample size for a mark-recapture population estimate and a lack of access to Vey Meadows), but it would be interesting to look into for Imnaha, Lookingglass, maybe Catherine and Lostine populations.*
    
- *Feldhaus mentioned that the number of carcass recoveries varies with level of spawning ground survey effort, such as when wildfire conditions prevent key surveys. Although variability among years and populations in number of carcass recoveries is not itself a problem for the model, systematic bias in probability of recovering prespawn mortalities is a potential problem.*
    
    ```diff
    +ODFW team will review relevant spawning ground survey records to identify  any years/populations (primarily MIN and UGR) where we judge that missed surveys or similar issues make the associated prespawne mortality estimates too unreliable to include (Pre-spawn data will then be set to `NA` for the relevant records in the revised carcass data set).
    -Additionally the team supports re-implementing a threshold (likely n=10 or n=20) of the minimum number of female carcasses recovered in order to estimate prespawn mortality from the carcass data.
    ```
    
    *Although this may have only limited effect on model fit, it better matches our level of confidence in the data for estimating pre-spawn mortality rates.*
    

### 1.4 Carcass sampling for age/origin

- Apply population-specific correction factors to account for systematic differences in carcass recovery probabilities as a function of age. (Should we consider this being origin-specific as well?)
    
    - See exploratory analysis: [explore-weir-vs-carcass-composition.Rmd](https://github.com/gibsonpp/GR-sslcm-data/blob/master/exploratory-analyses/explore-weir-vs-carcass-composition.Rmd)
    - Issue: [#39](https://github.com/bstaton1/GR-sslcm/issues/39)
    - Issue: [#104](https://github.com/bstaton1/GR-sslcm/issues/104)
- This approach does ignore some known issues with: assigning origin at the weir (see section 1.2 above); low carcass recovery numbers; low/variable weir efficiency (which may be non-random with respect to fish age, see section 1.2 also above)
    
- Sex (proportion female) is currently assigned according to a time-constant and age/population-specific proportion from straight carcass recovery data
    
- *Carcass recovery rates are known to vary according to sex, so this approach risks systematically overestimating the proportion female.*
    
- *Brandt proposes just assuming that all age-4 and age-5 fish are 50% female (and that all age-3 fish are 0% female).*
    
- *The full team supports making this change: assuming 50% female is simple, defensible, and avoids known problems with sex assignment in weir and carcass data. Additionally, the parameters governing egg to parr survival are "tuned" to the scale of total egg production based on spawner abundance. If per capita egg production changes in the model, then the $\alpha_j$ and $\beta_j$ parameters will adjust to accommodate it.*
    
    ```diff
    -Revise the model to use 50% as the proportion female for age-4 and age-5 in all years and populations.
    ```
    

### 1.5 Tributary Harvest

- Tributary harvest numbers (all harvest is below weirs) are included in the total return to river adult abundance values. These harvested fish are currently ignored (i.e., the model allows them to have a chance to spawn, it currently ignores any mortality between LGR and the weir as adults).
    
- *Only Lostine has substantial tributary harvest in most years (other populations typically have little or no harvest), but team agrees that we need to account for these harvested fish in some way.*
    
- *Proposed approach: treat harvested fish like weir removals and subtract them out from the population at this point in the model.*
    
    - *Harvest data are problematic, but it's still what we have to work with. We do have rough estimates for composition of harvested fish by origin and age.*
    - *This approach will require deciding how to integrate harvest information (counts by age/origin) with weir data (individual fish records) -- Staton sees simple solution, see below*.
    - *The Lostine situation is complicated by repeated "recycling" of some fish captured at the weir back to the lower river for additional harvest opportunity. It may be appropriate to treat any fish outplanted to the lower Lostine as being removed for harvest.*
    
    ```diff
    +ODFW will reconsider the Lostine weir data relative to harvest there, and revise the data set(s) appropriately.
    +For data organization, Staton suggests to leave RTR estimates and weir records the same as before, and upload additional data file that includes numbers harvested by year, age, origin, and population; these will be included with fish removed at weir.
    -Change to model/data prep: include numbers harvested by year, age, origin, and population in the `B[year,age,origin,population]` data array.
    ```
    

### 1.6 Hatchery-origin returning adults from brood years prior to the start of official hatchery releases

- For simplicity, treat all these fish as strays in the model
    - Issue: [#44](https://github.com/bstaton1/GR-sslcm/issues/44)
    - Also discussed at LCM meeting (2021-01-20, see the [notes](https://github.com/bstaton1/GR-sslcm/blob/master/99-admin/meeting-notes/2021-01-20-topics-notes.md))
- ...even though we know that in reality some of these adults were the product of early hatchery releases (see also section 3.1, below)
    - See also P. Gibson's [exploratory analysis](https://github.com/bstaton1/GR-sslcm/files/5822358/returning-adults-from-early-hat-releases.pdf), also discussed in the [#44](https://github.com/bstaton1/GR-sslcm/issues/44) thread
- *LS team agrees that this is a reasonable way to deal with the issue. Early hatchery releases were haphazard. There was a large degree of straying from early hatchery releases of Rapid River stock fish; but after the official hatchery releases began in ~2000, straying of adults from other populations into our populations is probably not something we need to worry about for the LCM populations (though it would be more of a concern in the Lookingglass population).*

## 2\. Fecundity

- Total egg production is a function of the number of females spawning in nature
    
    - Fecundity (eggs/female) varies by **age** and by **population**, but not by year, origin, or individual size
    - See P. Gibson's [exploratory analysis](https://github.com/gibsonpp/GR-sslcm-data/blob/master/exploratory-analyses/explore-relative-fecundity.Rmd)
    - See section 1.4 above for how number of females is determined
- Data source: hatchery spawning fecundity records, averaged across years
    
- *LS team expressed reservations about using a single fecundity estimate (by age) for all populations in all years. There is some evidence of a trend toward decreasing size at age. Also, average size of returning adults seems to vary by population.*
    
- *Size (length) is the primary predictor of fecundity -- in LS team's judgment, there is no need to worry about additional potential effects of origin, population, or year if length (by population/year/age) is accounted for.*
    
- *Proposed alternative approach for determining relative reproductive output to use in the model:*
    
    - *Use all hatchery spawning records to calculate an overall relationship between length and fecundity*
    - *Use weir and carcass length data to estimate mean length of female spawners by age, for each population and year. This will require some thought for the Minam population.*
    - *Use the regression relationship to calculate fecundity at age for each population and year.*
- *Staton notes that potential benefit of doing this would be to help explain some negative trends in egg-to-parr survival process noise terms. That these negative trends exist without accounting for time trending reproductive output implies that the model could benefit from accounting for it. The trends are not super strong, and vary in strength across populations. Autocorrelation in these process noise terms is of at least equal concern to Staton currently.*
    
- *Staton and Justice also note that a disadvantage of doing this is the need for a fecundity value to be available in every year/age/population (no `NA` values; don't worry about age-3, since we assign 0% as female), including in future years for the simulation analyses. This will almost certainly include some imputation of missing values in the observed period (which may be tricky?). We will also need to assume how the length/fecundity trend will look in the future for the simulation.*
    
    ```diff
    +ODFW team will consider whether we do in fact have adequate data to implement this approach. If so, then calculate the fecundity at age/year/population values as separate data set that can be fed into the model.
    ```
    

## 3\. Hatchery Smolt Releases

### 3.1 Abundance

- Use estimated number of hatchery smolt released per population per year as known without error (data source: LS CWT releases database)
- Hatchery smolt releases in the model begin with BY 1997 (LOS) or BY 1998 (CAT, UGR)
- We decided to ignore (for modeling purposes) some pre-1995 hatchery smolt/parr releases into UGR and CAT (decision was that including them is not worth the additional complexity that would be required in the model)
    - Issue: [#44](https://github.com/bstaton1/GR-sslcm/issues/44)
    - See also P. Gibson's [exploratory analysis](https://github.com/bstaton1/GR-sslcm/files/5822358/returning-adults-from-early-hat-releases.pdf), also discussed in the [#44](https://github.com/bstaton1/GR-sslcm/issues/44) thread

### 3.2 Survival to LGR

- Data source: LS PIT Tag releases database
- Survival estiamtes and variances are aggregated across raceways to get a single estimate per population/year (three instances of parr releases with survival estimates are included in this aggregation)
    - Issue: [#9](https://github.com/bstaton1/GR-sslcm/issues/9)
    - Discussed at 2020-11-19 meeting (see [meeting notes](https://github.com/bstaton1/GR-sslcm/blob/master/99-admin/meeting-notes/2020-11-19-topics-notes.md))
- *LS team agrees that this approach is appropriate (for both survival and abundance of hatchery smolt).*
- *Currently outmigration survival for wild fish in the model is size-based; we could incorporate size for hatchery fish as well, but LS team agrees that this is not of much interest.*
- *One factor that does impact outmigration survival of hatchery smolt is timing (time of year) of release; we could try to model Julian day of release as a predictor of outmigration survival for hatchery smolt. However, the team concluded that this would not be worth the extra effort and complexity. The model as-is could still support hypothetical increases in outmigration survival for hatchery smolt when it comes to simulation.*

## 4\. Juvenile Outmigration Survival Through the Hydrosystem (LGR $\to$ BON)

- Use Fish Passage Center's annual hydrosystem survival estimates
    - Separate estimates for hatchery ("Catherine Creek AP") vs. natural ("Aggregate Wild Snake R. Chinook") origin smolt.
    - Issue: [#49](https://github.com/bstaton1/GR-sslcm/issues/49)
    - Pull request (data repo): [#12](https://github.com/gibsonpp/GR-sslcm-data/pull/12)
- *LS team agrees that this a good data source to use*