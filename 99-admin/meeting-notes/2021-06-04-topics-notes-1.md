# GR-SSLCM Meeting 6/4/2021

## Topic: How to model sea lion mortality

**Attendees**: Gibson, Hatch, Justice, Liermann, Staton

### Summary of Discussion

* Some members of the GR-SSLCM modeling team met with Doug Hatch on June 4, 2021.
* We discussed three citeable works that quantify sea lion mortality on spring Chinook: 
  * [Rub et al. (2019)](https://cdnsciencepub.com/doi/pdf/10.1139/cjfas-2018-0290) - Documents a study where adults were captured and tagged near Astoria and any fish that did not successfully pass Bonneville were assumed to have died, with the predominate cause of death attributed to sea lions. They provide estimates from 2010 -- 2015 of survival from all mortality sources other than harvest for the aggregate run of spring Chinook salmon.
  * [Sorel et al. (2020)](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/1365-2664.13772) - Presents a modeling framework for estimating population-specific survival estimates from Astoria to Bonneville. They use the same tagging data as Rub et al., but integrate it with run timing data from 18 populations (including the four included in the GR-SSLCM). They show earlier migrating populations suffer higher mortality than later running populations.
  * [Tidwell and vander Leeuw (2021)](http://pweb.crohms.org/tmt/documents/FPOM/2010/Task%20Groups/Task%20Group%20Pinnipeds/2020%20Pinniped%20Annual%20Report.pdf) - Presents a long-term monitoring effort to track sea lion abundance and predation events immediately downstream of Bonneville. Provides estimates from 2000 - 2020 of aggregate-run mortality that occurs immediately downstream of Bonneville.
* Doug emphasized two key sources of variability that should be accounted for in modeling historical adult returns:
   * __Variability in sea lion mortality across populations__ resulting from differential overlap in run timing and sea lion peak abundance. Sea lion abundance is generally higher earlier in spring because they are drawn in by other prey sources (eulachon) and are at high abundance by the time the earliest populations show up.
   * __Variability in sea lion mortality over time__ resulting from an inter-annual trend of increasing abundance of sea lions showing up downstream of Bonneville.
* Doug thought an approach where we model temporal variability using a set of 3 - 4 time blocks where the survival is constant within a block but differs among blocks would be appropriate and that we could somehow use the output from Sorel et al. to reflect among-population variability in sea lion mortality.
* The problem: total survival (what we need) is only available currently for 2010 - 2015 (in Rub et al. and Sorel et al.). The longer time series (USACE) only covers the area closest to Bonneville.
  * Potential solution: regress Rub et al. against the USACE estiamtes to obtain a predictive model that can be used to produce estimates on the Rub et al. scale but that go back to 2000.
  * Build time blocks off of this predicted time series to pass to the LCM
  * This time series wouldn't be population-specific, would need to some how relate Rub et al. estimates with the Sorel estimates to get population-specific time series

### Follow-up Steps

* Conduct analysis described above comparing Rub et al. and USACE estimates
* Contact author of Sorel et al. to obtain population/year specific estimates
* Build hypothetical time blocks and values for historical LCM
* Run it by team members and get suggestions for improvements/approval for inclusion