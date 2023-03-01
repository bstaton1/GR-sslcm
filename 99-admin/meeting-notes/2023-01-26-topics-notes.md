# GR-sslcm Meeting 2023-01-26

**Attendees:** Liermann & Staton

*Post-meeting notes in italic*

## Purpose

The model has been giving some non-ideal simulated outcomes for adults and juveniles -- primarily over predicting them, as well as a couple other non-ideal behaviors. This meeting was intended to be a technical model-oriented discussion about how important this behavior is to fix, and options to try for doing so.

## Materials Sent in Advance

Sent to M. Liermann by B. Staton on 2023-01-11

- `output-plots-FOR-SOS-medium.html`
- `sim-vs-obs-FOR-SOS-medium.html`
- `GR-sslcm-math-description.pdf`

## General Items to Catch Up On

- Methods for validating model statistically
    - **Quantile residuals** for all stochastic processes: process and observation.
        - Obtain the CDF of the realized value from the posterior predictive distribution. Residuals have a IID U(0,1) distribution if assumptions met
    - Ratio of **summary statistics** of values obtained during **observed and simulated time periods**
        - Summary stats include time series mean (most important in Staton's opinion), CV, and ACF lag-1
        - Intended to measure how well simulated outcomes from model replicate those of the observed period
- Use power function for parr mean length vs. egg density relationship? (Issue [#159](https://github.com/bstaton1/GR-sslcm/issues/159)) *Liermann agrees with Staton's assessment that this should be tried based on the promise shown in the raw data*
- Are correlation parameters too restrictive? Recently found an [approach that may help](https://www.jstor.org/stable/24306780) and may be worth trying for estimating correlation matrices within unique elements for each population pair. *Liermann seemed to be interested in this, but only if it doesn't blow up the model/take a long time to implement/assess*
- Carcass composition data show evidence of over-dispersion *Liermann suggested that year/population random effects could be added to the carcass recovery correction model to account for over-dispersion.*
- Ocean-age-specific delta parameters? *Liermann agreed that it makes sense to assume that the survival advantage of being NOR should diminish with increasing ocean age if model seems to support it. A potential thing to try to improve performance around ocean processes in general is to estimate a different Yr1 ocean survival time-series for NOR and HOR fish ([#162](https://github.com/bstaton1/GR-sslcm/issues/162), rather than the current approach of assuming they are perfectly correlated with an estimated shift in magnitude up or down (i.e., the delta parameters).*

## Major Problem #1: Adult Return Too Abundant and Biased Away from Age-4

**Source of Evidence for Problem and Logic**: `sim-vs-obs-FOR-SOS-medium.html` and `output-plots-FOR-SOS-medium.html` (process noise diagnostic section)

**Description of Problematic Behavior**

- The total return to estuary is too high: LOS, UGR, HOR for most populations; mean 1.2-1.5 times greater in simulated period
- The adult return is biased young and old, i.e., away from age-4. NOR age-3 return is 1.2-1.7 times greater in simulated period; NOR & HOR 1.1-1.5 times greater in simulated period

**Why So Concerning?**

- The adult return will be a key metric for tracking future population status in scenario simulations. Simulated population improvements without intervention aren't ideal -- too optimistic of predictions for future.
- The adult composition drives per capita egg production. However, curiously given the biased age composition, eggs per spawner mean ratios are ~1

**Logic Chain to Describe Behavior**

- Yr1 ocean survival has tendency to pull higher in simulated vs. observed period
    - Causes more fish to survive and mature as age-3: explains positive bias in age-3 return proportion and over abundance of adult returns to estuary
- Pr(Mature age-3) is pulling high for UGR, which has the high abundance, high age-3 return problem the worst
    - Causes fewer fish to make it to age-4 and to decide to mature at age-3 instead
- Pr(Mature age-4) is pulling low
    - Causes too many to make it to age-5 and decide to mature of returning at age-4
- *Liermann understands and agrees this a likely chain of events contributing to the over abundance of adult returns, but isn't sure what is (a) pulling Yr1 ocean survival up, (b) pulling Pr(age-3) up, and (c) pulling Pr(age-4) down*
- *Liermann seemed to see these items as minor to moderate problems -- not ideal, but we shouldn't spend forever tracking down and fixing the root cause. Something similar to "A lot of modelers would have called it done \[a long time ago\]" was said*

**Potential Mechanisms(s) Causing Behavior**

- Use of `dbeta(1,1)` priors for `mu_phi` and `mu_psi` places too much weight on values near and above 0.5 when we know they can't be this high.
    - **Evidence**: The year-specific values in observed period are getting pulled lower than the mean of the posterior predictive distribution (PPD) to explain the data (see process noise term diagnostics for `phi_O0_O1` and `psi_O1` \-\- they are on the low end of PPD on average, `psi_O2` are on the high end of PPD on average), but the prior pulls the `mu_*` parameters up; with no data to pull simulated years back down, get bias in simulated period relative to observed. (*see [#161](https://github.com/bstaton1/GR-sslcm/issues/161))*
    - Potential solution: use more restrictive `dbeta()` priors
- There are unaccounted time trends in maturation probabilities. `mu_psi_*` parameters are informed by full observed time series simulation is based on `mu`

*Liermann agrees better priors should be used and that the flat ones can serve to pull up the posterior by being unnecessarily vague.*

## Major Problem #2: NOR Juvenile Abundance Too High

*This issue was not discussed much during the call.*

**Source of Evidence for Problem**: `sim-vs-obs-FOR-SOS-medium.html`

**Description of Behavior**

- The juvenile abundance is overall too high; parr mean 1.15-1.3 times greater in simulated period for LOS and UGR
- Smolt at LGR 1.2 - 1.8 times greater includes all populations.
- These problems are all for NOR juveniles.

**Why So Concerning**

- Smolt outmigrants will likely be a metric used in summarizing population status, same concerns as for over-abundant adult returns

**Logic Chain to Describe Behavior**

- Parr recruitment is too high in the same populations that have super over abundant adults (LOS and UGR)
- Eggs per spawner ratio between sim and obs period is ~1, so over abundant adults $\to$ over abundant parr
- parr $\to$ LGR smolt survival higher in simulated period (mean ratio 1.08-1.16)
- So in general, too many adults $\to$ too many parr (for some populations) coupled with over-predicted juvenile survival (for all populations, to varying degrees) $\to$ too many smolts at LGR

**Potential Mechansism(s) Causing Behavior**

- The over prediction of parr abundance should be mostly addressed via fixing major issue #1.
- Over-winter survival and migration survival to LGR are both being over predicted -- related to spuriously high variance terms -- these terms are more variable in simulated vs. observed period.
- Potential solution: tighter prior on variance terms to regularize?
- Can't dive too deep on this issue until major problem #1 (adults) sorted out

## Overall Recommendations/Next Steps

*All of this content is post-meeting, so italics are omitted.*

### Regarding Sim vs. Obs

Liermann liked this approach of comparing simulated vs. observed outcomes as a means to validate the model's ability to reproduce past outcomes. However he indicated that we shouldn't require that every diagnostic looks perfect/great before moving forward. We can use these tools to highlight potential areas where the model may be unreliable and to build correction factors for variance terms where appropriate.

**Key recommendation** was made with respect to issues Staton has experienced with imputing the known values (e.g., weir removals, hatchery smolt inputs) for the simulated period. We should use a "tack-on" method ([#160](https://github.com/bstaton1/GR-sslcm/issues/160)) for the simulation validation: fit model to only BY2000+ data (to avoid period with no HOR inputs), then for the simulated period, simply duplicate the entire time series of known values that were used in the observed period. This at least guarantees that the known quantities have the same mean/variance in both periods, even if the covariance with other outcomes (e.g., weir removals and return abundance are correlated) is not accounted for.

### Regarding Quantile Residuals

Liermann liked this approach of standardizing residuals from random processes that have outcomes on different scales. He used QQ plots previously to identify issues with lack of fit previously, we didn't discuss the benefits/detriments of either approach.

### Regarding Prior Distributions

Liermann & Staton agree that we can do better than using `dbeta(1,1)` on all `mu` parameters. See [#161](https://github.com/bstaton1/GR-sslcm/issues/161).

**Key recommendations** were made that we should use priors that favor values &lt;0.5 (perhaps relatively strongly so) for Yr1 ocean and Pr(age-3) and priors that favor &gt;0.5 for Pr(age-4). We should also revisit other priors in the model, specifically for logit-normal SD terms. There is some reason to think that perceived bias in the simulated period is caused by greater variability in simulated period, which when the expected value is low, would induce many more high outcomes. The current prior for all SD parameters is `dunif()` sometimes allowing values as high as 3. This should be reconsidered, since the use of relatively strong priors on variance components may help regularize the inter-annual variability between simulated and observed periods. Diagnostic figures should be added that enable comparing prior vs. posterior distributions to assess the influence of these priors, and we should have some kind of a sensitivity analysis to present.

Liermann recommended to speak with colleagues to build priors. Staton thinks this would take too long and would be too difficult to explain, and the variability of the prior would be difficult to quantify based on talking to a handful of folks. Staton thinks we are better off stating our (ballpark) prior expectations and defining the purpose: to downweight biologically implausible values, generally not to impose much information on the actual shape.