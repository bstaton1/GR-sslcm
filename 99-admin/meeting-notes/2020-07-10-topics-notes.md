# GR-SSLCM Meeting 7/10/2020 Discussion Topics

Notes are _italic_

## Model/Data-Based Topics

### Based on the spreadsheet (simple-sim.xlsx)

* Primary states to track
* How to initialize modeled states
  * _Idea to initialize the states using adults in the early years rather than smolt arriving at LGD._
  * _Questions of whether this will line up with keeping everything as square time series across populations?_
  * _Staton will try building a replicate tab in simple-sim.xlsx to envision what this would look like_

### Which transition probabilities to estimate vs. which to inform directly with external estimates

* Eggs to parr [density-dependent via BH function; ESTIMATE]
  * _Yes, but question remains whether we should use eggs or some other "effective spawners" multiplier_
* Parr overwinter survival [density-dependent via logistic function; ESTIMATE]
  * _Yes, but Sharma indicated this could be done with another BH function. Used logistic originally for consistency with Cooney's model_
* Smolt movement to LGD [random deviates around a mean; ESTIMATE]
  * _Yes_
* Smolt movement from LGD to estuary [EXTERNAL ESTIMATES??]
  * _Team seems to agree that it would be ideal to feed in information on hydrosystem survival, but have concerns about where to get the estimates to inform it. _
  * _Gibson suggested we could re-fit some CJS models with GR data only and treating dams after LGD differently._
* Marine survival (smolt to age 3, age 3 to age 4, age 4 to age 5) [ESTIMATE??]
  * _Team seems to agree that it would be nice to estimate ocean mortality rather than fix it, but has concerns that survival and maturity will be confounded_
* Maturity from ocean juvenile to returning adult [ESTIMATE]
  * _Idea came up to simplify maturity schedule: give all populations the same vector each year, and/or remove the temporal variation here_
* Survival upstream as adults, prior to reaching weirs (harvest, predation, etc.) [EXTERNAL ESTIMATES??]
  * _Should talk with S. Ellis (CRITFC) about GR-specific exploitation rates_
  * _Won't have good data on predation, but there are some papers we can use to create some rough placeholders for now_
* In-tributary fisheries (LSN mostly)
  * _There are data on these removals: numbers of fish harvested by year_
* Brood stock removals
* Pre-spawn mortality [ESTIMATE??]

### How to deal with hatchery production?
* Hatchery adults return and are counted by the weirs.
* Some are passed above the weir and have a chance to contribute to future natural production.
* These features suggest we need biological processes for their generation within the model allowing them to be observed?
* Are there data on hatchery smolt releases (UGR, CCR, LSN)?
* Do we use different survival/maturity terms for hatchery/natural origin fish through hydrosystem, ocean, and upstream phases?
* _We have two options here_
  * _Model hatchery dynamics following release and move them through other parts of the life cycle: this would be ideal, but would involve a lot more complexity_
  * _Perform a pHOS expansion where the model tracks only natural origin fish until counted as adults_

### At what stage do we introduce sex structure?
  * Maturity?
      * _Team thinks this is good_

## Code/Workflow-Based Topics

### Regular meeting structure

* Objective: keeps us on track, regular discussions allow them to be shorter and more focused, easier to assign tasks/check-in than via email alone
* Propose we meet at a semi-regular day/time every 1-2 weeks
* Is there a day/time that will *in general* work for at least the modeling team
  * _Team thinks Thursdays @10am will be best_
  * _Next meeting: Thursday 07-30-2020 @10am_
    * _Tentative meeting: Staton is in field from 07-22 - 07/29_
    * _Will we have anything new to talk about?_
* MS Teams?
  * _Team agreed that this works well_
* Group of required participants (involved first-hand with data, code, analytical decisions): Gibson, Liermann, Sharma, Staton
* Group of optional participants (have an interest in the topics we discuss, may have key input in some discussions): Justice, Jordan, Sedell, White [OTHERS??]

### Naming conventions: pull up HTML file and go through

  * _Overall like Staton's proposed conventions, will start here and can use global find/replace if necessary later. In particular, we may wish to be more explicit about what symbols mean by replacing them with "surv" or "mat"_.
  * _It is not necessarily the most intuitive at first, but it shoots for consistency_

### R package vs. no package

* Staton proposes no package (for now at least)
* _Team agrees, no package for now, reconsider when things are more finalized_

### GitHub collaboration

* Set of “rules” for moving forward – consistency in work flow, avoid merge conflicts
  * _Staton will continue compiling this list and finalize it ASAP_
* Identify a main repository owner; also will serve as main code reviewer, curator, and merger
  * _Team identified Staton as this person_
* Comfort level with moving forward on model development in a GitHub collaborative process?
  * _Team seems fairly comfortable, though some sentiments raised that there will still be some learning to do_

## Post-Meeting To-Do: Staton

_NOTE_: These are all notes so _italics_ not necessary

* Try initializing with adults rather than smolts at dam – build a new tab in simple-sim.xlsx to visualize how this would work.
* Explore possibility of estimating ocean mortality and informing downstream survival with estimates using the simulation code from a while back - mostly to see if this is even possible
* Send invite for next meeting
* Finalize GitHub rules and send out
* Finalize and distribute naming conventions
* Create github repo
  * Add main folder structure w/placeholders where necessary
  * Add collaborators: Gibson, Liermann, Sharma
  * Submit data issue and assign it to Gibson
  * Once Gibson adds data files, Staton will write code to prep data for fitting
  * Staton will write code to fit model to one population
  * Once we agree on the initial structure here, we will build in other populations