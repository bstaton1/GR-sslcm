
> This repository stores the code code used to fit the model and compile the output presented in the manuscript *A state-space model to quantify density dependence, demographic heterogeneity, and spatial synchrony in Grande Ronde Basin Chinook salmon populations* by B.A. Staton, P.P. Gibson, M. Liermann, C. Justice, M.J. Kaylor, R. Sharma, and S.M. White, which is currently undergoing peer review.

[![ArticleDOI](https://img.shields.io/badge/Article-PLACEHOLDER-blue?logo=doi&logoColor=f5f5f5)]()  
[![GitHub Repo Archive DOI](https://img.shields.io/badge/GitHub%20Repo%20Archive-PLACEHOLDER-blue?logo=github)]()

## Repo Organization

| Subdirectory       | Description                                                                                                                                                                   |
|:-------------------|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `00-data/`         | Scripts for compiling the data files from their raw format in the [GR-sslcm-data](https://github.com/gibsonpp/GR-sslcm-data) repository for model fitting.                    |
| `01-functions/`    | Scripts housing functions used by scripts throughout the repository, organized by those for data preparation, plotting, MCMC initialization, and miscellaneous utility tasks. |
| `02-model/`        | Scripts to fit the model as well as the JAGS model code. MCMC samples are saved here after model fitting.                                                                     |
| `03-post-process.` | Scripts to analyze the model output, including those to produce Rmarkdown files and the figures/tables presented in the manuscript.                                           |
| `99-admin/`        | Notes from various meetings discussing model development and preliminary results.                                                                                             |

The additional script include `00-packages.R` loads most packages used by other scripts in this repository.

## Reproducibility

The output file from running the model is too large to track in a git repo, so this code will need to be executed locally to reproduce the output presented in the manuscript and supplement.

First, download the contents of this repository and the data repository ([GR-sslcm-data](https://github.com/gibsonpp/GR-sslcm-data)), and unzip them such that they are in the same parent directory (called whatever you like), as follows:

``` bash
├── parent
    ├── GR-sslcm
    ├── GR-sslcm-data
```

Then, navigate to the `GR-sslcm` directory and from the command line, execute:

``` bash
Rscript 02-model/fit.model.R
```

To see the help documentation regarding the various arguments and their settings, execute:

``` bash
Rscript 02-model/fit.model.R --help
```

If you do not have access to a command line, the model can still be executed from within [RStudio](https://posit.co/download/rstudio-desktop/). Open the `GR-sslcm.Rproj` and execute:

``` r
source("02-model/fit-model.R")
```

Either option will fit the state-space model, build the output summary file, and the manuscript figures and tables using the same configurations as presented in the manuscript. The supplement will be a `.html` file in the `03-post-process` subdirectory, and the figures/tables will be saved in the `03-post-process/ms-content` subdirectory.

> *The default is to use 4 CPU cores to run 150,000 MCMC iterations per each of the 4 chains.
> This takes approximately 18 hours for MCMC and an additional several hours for post-processing on a modern laptop computer (RAM: 64GB; Processor: 20 threads, 5.40 GHz Max Turbo, 24 MB Cache).
> Similar results to those we present can be obtained with fewer MCMC iterations but lower convergence, for example:*
>
> ``` bash
> Rscript 02-model/fit-model.R --mcmc medium
> ```
>
> *This should take approximately a third of the time, other options would be*:
>
> ``` bash
> Rscript 02-model/fit-model.R --mcmc short
> ```
>
> *or even*:
>
> ``` bash
> Rscript 02-model/fit-model.R --mcmc vshort
> ```
>
> *The last option should take only several minutes for MCMC and render the output in under an hour, however, because the chains will be far from converged, **it is intended for testing the code only**.*

## Dependencies: JAGS and R

Program JAGS was used to run the MCMC for fitting the state-space model (v4.3.1 at the time of manuscript publication). It can be downloaded [here](https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/) and use all defaults for installation.

Program R was used to perform all analyses (v.4.4.1 at the time of manuscript publication). It can be downloaded [here](https://www.r-project.org/) and use all defaults for installation.

The table below shows all R packages (and their versions at time of publication) that are specifically called (typically with `pkg::fn()`) or loaded (using `library(pkg)`).

| Package                                                       | Version | Description                                                                  |
|:--------------------------------------------------------------|--------:|:-----------------------------------------------------------------------------|
| [`abind`](https://CRAN.R-project.org/package=abind)           |   1.4.5 | Combine Multidimensional Arrays                                              |
| [`argparser`](https://CRAN.R-project.org/package=argparser)   |   0.7.2 | Command-Line Argument Parser                                                 |
| [`bookdown`](https://CRAN.R-project.org/package=bookdown)     |    0.40 | Authoring Books and Technical Documents with R Markdown                      |
| [`corrplot`](https://CRAN.R-project.org/package=corrplot)     |    0.94 | Visualization of a Correlation Matrix                                        |
| [`grDevices`](https://CRAN.R-project.org/package=grDevices)   |   4.4.1 | The R Graphics Devices and Support for Colours and Fonts                     |
| [`jagsUI`](https://CRAN.R-project.org/package=jagsUI)         |   1.6.2 | A Wrapper Around ‘rjags’ to Streamline ‘JAGS’ Analyses                       |
| [`kableExtra`](https://CRAN.R-project.org/package=kableExtra) |   1.4.0 | Construct Complex Table with ‘kable’ and Pipe Syntax                         |
| [`knitr`](https://CRAN.R-project.org/package=knitr)           |    1.48 | A General-Purpose Package for Dynamic Report Generation in R                 |
| [`latex2exp`](https://CRAN.R-project.org/package=latex2exp)   |   0.9.6 | Use LaTeX Expressions in Plots                                               |
| [`matrixcalc`](https://CRAN.R-project.org/package=matrixcalc) |   1.0.6 | Collection of Functions for Matrix Calculations                              |
| [`msdown`](https://CRAN.R-project.org/package=msdown)         |   2.2.0 | A Framework for Writing Manuscripts with the rmarkdown and bookdown Packages |
| [`mvtnorm`](https://CRAN.R-project.org/package=mvtnorm)       |   1.2.6 | Multivariate Normal and t Distributions                                      |
| [`posterior`](https://CRAN.R-project.org/package=posterior)   |   1.6.0 | Tools for Working with Posterior Distributions                               |
| [`postpack`](https://CRAN.R-project.org/package=postpack)     |   0.5.4 | Utilities for Processing Posterior Samples Stored in ‘mcmc.lists’            |
| [`reshape2`](https://CRAN.R-project.org/package=reshape2)     |   1.4.4 | Flexibly Reshape Data: A Reboot of the Reshape Package                       |
| [`rmarkdown`](https://CRAN.R-project.org/package=rmarkdown)   |    2.27 | Dynamic Documents for R                                                      |
| [`scales`](https://CRAN.R-project.org/package=scales)         |   1.3.0 | Scale Functions for Visualization                                            |
| [`shiny`](https://CRAN.R-project.org/package=shiny)           |   1.9.1 | Web Application Framework for R                                              |
| [`stringr`](https://CRAN.R-project.org/package=stringr)       |   1.5.1 | Simple, Consistent Wrappers for Common String Operations                     |
| [`this.path`](https://CRAN.R-project.org/package=this.path)   |   2.5.0 | Get Executing Script’s Path                                                  |

Running the code below will display the packages not already installed on your machine:

``` r
pkgs <- c("abind", "argparser", "bookdown", "corrplot", "grDevices", "jagsUI", "kableExtra", "knitr", "latex2exp", "matrixcalc", "msdown", "mvtnorm", "posterior", "postpack", "reshape2", "rmarkdown", "scales", "shiny", "stringr", "this.path")
pkgs[!pkgs %in% rownames(installed.packages())]
```

## Session Info

<details>
<summary>
<b>Click to Expand/Hide Session Info</b>
</summary>

    ## ─ Session info ───────────────────────────────────────────────────────────────
    ##  setting  value
    ##  version  R version 4.4.1 (2024-06-14 ucrt)
    ##  os       Windows 11 x64 (build 22631)
    ##  system   x86_64, mingw32
    ##  ui       RTerm
    ##  language (EN)
    ##  collate  English_United States.utf8
    ##  ctype    English_United States.utf8
    ##  tz       America/Los_Angeles
    ##  date     2025-05-13
    ##  pandoc   3.1.11 @ C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools/ (via rmarkdown)
    ## 
    ## ─ Packages ───────────────────────────────────────────────────────────────────
    ##  ! package        * version  date (UTC) lib source
    ##    abind          * 1.4-5    2016-07-21 [1] CRAN (R 4.4.0)
    ##    argparser      * 0.7.2    2024-04-04 [1] CRAN (R 4.4.1)
    ##    backports        1.5.0    2024-05-23 [1] CRAN (R 4.4.0)
    ##    bookdown       * 0.40     2024-07-02 [1] CRAN (R 4.4.1)
    ##    checkmate        2.3.2    2024-07-29 [1] CRAN (R 4.4.1)
    ##    cli              3.6.3    2024-06-21 [1] CRAN (R 4.4.1)
    ##    colorspace       2.1-1    2024-07-26 [1] CRAN (R 4.4.1)
    ##    corrplot       * 0.94     2024-08-17 [1] CRAN (R 4.4.1)
    ##    digest           0.6.36   2024-06-23 [1] CRAN (R 4.4.1)
    ##    distributional   0.4.0    2024-02-07 [1] CRAN (R 4.4.1)
    ##    evaluate         0.24.0   2024-06-10 [1] CRAN (R 4.4.1)
    ##    fansi            1.0.6    2023-12-08 [1] CRAN (R 4.4.1)
    ##    fastmap          1.2.0    2024-05-15 [1] CRAN (R 4.4.1)
    ##    generics         0.1.3    2022-07-05 [1] CRAN (R 4.4.1)
    ##    glue             1.7.0    2024-01-09 [1] CRAN (R 4.4.1)
    ##    htmltools        0.5.8.1  2024-04-04 [1] CRAN (R 4.4.1)
    ##    httpuv           1.6.15   2024-03-26 [1] CRAN (R 4.4.1)
    ##    jagsUI         * 1.6.2    2024-01-30 [1] CRAN (R 4.4.1)
    ##    kableExtra     * 1.4.0    2024-01-24 [1] CRAN (R 4.4.1)
    ##    knitr          * 1.48     2024-07-07 [1] CRAN (R 4.4.1)
    ##    later            1.3.2    2023-12-06 [1] CRAN (R 4.4.1)
    ##    latex2exp      * 0.9.6    2022-11-28 [1] CRAN (R 4.4.1)
    ##    lifecycle        1.0.4    2023-11-07 [1] CRAN (R 4.4.1)
    ##    magrittr         2.0.3    2022-03-30 [1] CRAN (R 4.4.1)
    ##    matrixcalc     * 1.0-6    2022-09-14 [1] CRAN (R 4.4.0)
    ##    mime             0.12     2021-09-28 [1] CRAN (R 4.4.0)
    ##    msdown         * 2.2.0    2025-04-25 [1] local
    ##    munsell          0.5.1    2024-04-01 [1] CRAN (R 4.4.1)
    ##    mvtnorm        * 1.2-6    2024-08-17 [1] CRAN (R 4.4.1)
    ##    pillar           1.9.0    2023-03-22 [1] CRAN (R 4.4.1)
    ##    pkgconfig        2.0.3    2019-09-22 [1] CRAN (R 4.4.1)
    ##    plyr             1.8.9    2023-10-02 [1] CRAN (R 4.4.1)
    ##    posterior      * 1.6.0    2024-07-03 [1] CRAN (R 4.4.1)
    ##    postpack       * 0.5.4    2022-12-21 [1] CRAN (R 4.4.1)
    ##    promises         1.3.0    2024-04-05 [1] CRAN (R 4.4.1)
    ##    R6               2.5.1    2021-08-19 [1] CRAN (R 4.4.1)
    ##    Rcpp             1.0.13   2024-07-17 [1] CRAN (R 4.4.1)
    ##    renv             1.0.7    2024-04-11 [1] CRAN (R 4.4.1)
    ##    reshape2       * 1.4.4    2020-04-09 [1] CRAN (R 4.4.1)
    ##    rlang            1.1.4    2024-06-04 [1] CRAN (R 4.4.1)
    ##    rmarkdown      * 2.27     2024-05-17 [1] CRAN (R 4.4.1)
    ##    rstudioapi       0.16.0   2024-03-24 [1] CRAN (R 4.4.1)
    ##    scales         * 1.3.0    2023-11-28 [1] CRAN (R 4.4.1)
    ##    sessioninfo      1.2.2    2021-12-06 [1] CRAN (R 4.4.1)
    ##    shiny          * 1.9.1    2024-08-01 [1] CRAN (R 4.4.1)
    ##    stringi          1.8.4    2024-05-06 [1] CRAN (R 4.4.0)
    ##    stringr        * 1.5.1    2023-11-14 [1] CRAN (R 4.4.1)
    ##    svglite          2.1.3    2023-12-08 [1] CRAN (R 4.4.1)
    ##    systemfonts      1.1.0    2024-05-15 [1] CRAN (R 4.4.1)
    ##    tensorA          0.36.2.1 2023-12-13 [1] CRAN (R 4.4.0)
    ##  D this.path      * 2.5.0    2024-06-29 [1] CRAN (R 4.4.1)
    ##    tibble           3.2.1    2023-03-20 [1] CRAN (R 4.4.1)
    ##    utf8             1.2.4    2023-10-22 [1] CRAN (R 4.4.1)
    ##    vctrs            0.6.5    2023-12-01 [1] CRAN (R 4.4.1)
    ##    viridisLite      0.4.2    2023-05-02 [1] CRAN (R 4.4.1)
    ##    xfun             0.46     2024-07-18 [1] CRAN (R 4.4.1)
    ##    xml2             1.3.6    2023-12-04 [1] CRAN (R 4.4.1)
    ##    xtable           1.8-4    2019-04-21 [1] CRAN (R 4.4.1)
    ##    yaml             2.3.10   2024-07-26 [1] CRAN (R 4.4.1)
    ## 
    ##  [1] C:/Users/bstaton/AppData/Local/R/win-library/4.4
    ##  [2] C:/Program Files/R/R-4.4.1/library
    ## 
    ##  D ── DLL MD5 mismatch, broken installation.
    ## 
    ## ──────────────────────────────────────────────────────────────────────────────

</details>
