# GR-SSLCM Meeting 1/20/2021

**Attendees**:

### Refresher on Model: key differences with Liermann & Sharma model

> This is primarily for R. Sharma, L. Burns, and M. Kaylor, if they attend

* Single population model (for now)
* Model tracks age 3 adults
* Model tracks sex composition as well (starting upon arrival to estuary)
* Hatchery fish have a dynamics model that tracks them from release through adulthood
* Model fits to survival of hatchery fish from release to LGR
* Model fits to survival of hatchery and natural origin fish from LGR to estuary
* Model fits to adult composition data from both weir sampling and carcass surveys (includes correction factor to allow discrepancies here)
* Model fits to pre-spawn mortality data gathered in carcass surveys
* Model does not weight adult spawners differently in reproductive output (currently)
* Habitat components: C. Justice is working on a more involved regression modeling approach to quantify weighted usable habitat to replace the "PEU" variable. 

## Progress Since Last Meeting (PRs merged)

### Fitting to Hydrosystem Survival Data

* PR [#62](https://github.com/bstaton1/GR-sslcm/pull/62)
* Data (estimates) from CSS draft 2020 report
* Separate survival for hatchery and natural origin fish
* Fit to data was incredibly poor using assumed values of ocean survival

### Estimating Ocean Survival

* PR [#63](https://github.com/bstaton1/GR-sslcm/pull/63)
* Estimate survival from OA0 -> OA1: mean, variability, and year-specific deviates (natural origin only)
* Estimate survival from OA1 -> OA2: mean, variability, and year-specific deviates (natural origin only)
* Apply the same estimates to OA2 -> OA3 as for OA1 -> OA2 (natural origin only)
* Survival is same for both sexes
* For hatchery origin fish, estimate a single, time constant scaler that adjusts these estimates to be from natural-origin to hatchery-origin
* Vastly improved fit to hydrosystem survival data
* Maturity schedule still fully estimated
* Convergence is decent, but some indications the model has a hard time estimating these parameters. My hope is that building in all four populations to inform these shared parameters will help

### Use Time Lag in Shared Ocean Survival Parameters

* PR [#65](https://github.com/bstaton1/GR-sslcm/pull/65)
* Because survival OA1 -> OA2 and OA2 -> OA3 occur in two consecutive calendar years for the same brood year,  and survival conditions probably operate on the calendar year basis, it makes sense to use the same parameter for two consecutive brood years, rather than the same brood year

## Unresolved Topics 

**Led by: P. Gibson**

* Should differential reproductive output be accounted for by age/sex/origin/year? More details in issue [#16](https://github.com/bstaton1/GR-sslcm/issues/16)
* Accounting for juveniles spawned downstream of traps (GR-sslcm-data issue [#15](https://github.com/gibsonpp/GR-sslcm-data/issues/15))
* New information about early hatchery releases, issue [#44](https://github.com/bstaton1/GR-sslcm/issues/44)
* Replace weir removal records with hatchery spawning records, (GR-sslcm-data issue [#7](https://github.com/gibsonpp/GR-sslcm-data/issues/7))

## Next Steps

* Seek input from expert(s) on how we treat hydrosystem survival (e.g., the data source, ignoring barging, etc.): is it adequate
* Contact experts for data/information on upstream migration:
  * Fishery mortality: year, origin specific?
  * Predation mortality: below Bonneville; year, origin specific?
  * Passage success rates
* Presentation to partners (Feb 18 at 8:15am)
* Symbology document could use updating
* Output plots document could use some updates
* Combine four populations into one model

**Next Meeting Date/Time:** 