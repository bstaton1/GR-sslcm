# GR-sslcm Meeting 2023-01-12

**Attendees**: Gibson & Staton

*Post Meeting Notes in Italics*

*This summary was built based on an email thread between B. Staton and P. Gibson titled "LCM Data Check-in" (2023-01-10 thru 2023-01-20)*

## Purpose

This meeting was called by P. Gibson after receiving an email from B. Staton inquiring about a variety of data related topics, thinking a call would be easier to discuss.

## 1\. Juvenile Mean Length Data

- P. Gibson had indicated that she was in the process of reworking the spring mean length data (the summer mean length data have been considered finalized for sometime) -- Staton wished to check in on the status of this
- *Via commit [861d1a3](https://github.com/gibsonpp/GR-sslcm-data/commit/861d1a3c99aae9db3db7584bff48b36fd003dfb1), P. Gibson updated all mean length data (2023-01-12).*

## 2\. Updating Data with Recent Brood Years

- In preparation for finalizing the model for manuscript submission (target mid-summer), we would like to have the data sets as current as possible. Staton wished to check in about whether adding RY 2020-2022 adult data and BY 2018-2022 juvenile data would be possible.
- *P. Gibson indicated that this shouldn't be too problematic.*

## 3\. Cite-able Documents for Manuscript

- In preparation for manuscript writing, Staton has been compiling references, and has hoped that the majority of data sources we used will have reports that are easily citable. For the data sources P. Gibson compiled and provided, Staton is hoping she can provide a relevant reference where readers can go to learn more about the methods for data collection/data-level estimation (e.g., CJS methods).
- The data sources compiled and provided by P. Gibson include:
    - Summer/fall/winter/spring PIT tag survival monitoring
    - Summer/spring mean length
    - Fall/spring screw trap passage
    - Hatchery smolt releases
    - Survival of hatchery smolt from release to LGR
    - Abundance of adults reaching their natal tributaries
    - Number of adults removed at weir
    - Age/origin composition of adults arriving at weir
    - Age/origin composition of carcasses
    - Spawn status of female carcasses (for pre-spawn mortality estimation)
    - Age-specific sex composition (carcasses)
    - Fecundity-at-age
- *P. Gibson indicated that getting sources for most of the data sources shouldn't be a problem, and that some will be found in the same annual reports. This is relatively low priority.*

## 4\. Data Availability for Sharing

- In preparation of manuscript submission, I've been thinking about code/data availability again. Let's start a conversation about this.
- *P. Gibson indicated that she would initiate discussions to explore getting approval to share the versions of the data we are using along side the model code.*

## 5\. Addition of ODFW Colleagues to Model/Data Team

- From P. Gibson "I think we need to bring in some of my ODFW colleagues (specifically the people who manage the adult data, in the same way that I manage the juvenile data) to get their buy-in on some of the decisions about how adult data are being used in the model."
- *Two of P. Gibson's colleagues at the La Grande ODFW office are more familiar with the adult data sources, and are the appropriate folks to weigh in on the validity of our treatment of the data in the model. These ODFW staff members are J. Feldhaus and E. Brandt*
- *In bringing Feldhaus & Brandt up to speed with the model and data usage, we should have one or more devoted discussions with them. Perhaps one on the model itself (big picture) and one on the specific use of specific data sets.*
- *Prior to these discussions, it would be useful for P. Gibson if B. Staton could compile a list of how each variable in each adult data set on GitHub is currently used, and which are not used at all.*

## To-Dos

*All of this content is post-meeting, so italics are omitted.*

### P. Gibson

**Shorter Term**

- Work on looping in ODFW colleagues as appropriate
    - Will meet with Feldhaus and Brandt about what role they want to play in the LCM going forward, and whether one of them really wants to be involved to the level of potential authorship. After we have that conversation, will send an email to the rest of the LCM team to announce these additions.
    - Schedule a meeting with at least Staton, Gibson, Feldhaus, and Brandt to present how the relevant data sources are currently used and what decisions were made leading up to this point. Gibson will organize a list of discussion topics informed by information on GitHub in issues, pull requests, and exploratory analyses.
    - Personal to-do: work to go over some of the details of data processing, where relevant with Feldhaus and Brandt.
- Update all data sets to go through brood year 2020 (juvenile data)/return year 2022 (adult data). Adult abundance might be the last/most difficult to update.
- Review the composite uncertainty values currently associated with the adult abundance estiamtes; revise as needed; and improve documentation of that processes. Feedback from Staton may be requested for this process.

**Longer Term**

- Identify basic references/methods documentation (primarily agency reports) for the ODFW data sources
- Maybe draft the piece of the main test (and supplemental) methods section that describes the data sources
- Figure out questions about data sharing: confirm official ODFW policy on whether the ODFW data sets will be able to be made publicly available (vs. by request only \[is that even possible?\]), and ensure that we have full collaborator permission (primarily CTUIR/NPT) for whatever level of data sharing we finally settle on.

*Gibson indicates that most of these items will need to wait until late Feb or March.*

### B. Staton

**Shorter Term**

- Compile document summarizing how each variable from each adult data source is used/not used by the data preparation code and model itself. Send to P. Gibson ASAP to facilitate discussions with Feldhaus & Brandt
- Be available to go over model and data topics with Feldhaus & Brandt, potentially at great length and detail, to get feedback on current usages of adult data.
- Continue SSM development: validation, priors, final outputs to display results
- Document meeting notes on GitHub

**Longer Term**

- Finish SSM development
- Write bulk of SSM manuscript content