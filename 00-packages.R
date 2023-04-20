suppressWarnings({
  suppressPackageStartupMessages({
    library(knitr)          # for rendering Rmarkdown documents; install.packages("knitr")
    library(rmarkdown)      # for rendering Rmarkdown documents; install.packages("rmarkdown")
    library(kableExtra)     # for building nice tables in Rmarkdown; install.packages("kableExtra")
    library(reshape2)       # for reformatting data (e.g., long to wide); install.packages("reshape2")
    library(abind)          # for binding together array objects; install.packages("abind")
    library(jagsUI)         # for calling JAGS from R; install.packages("jagsUI")
    library(postpack)       # for processing posterior samples; install.packages("postpack")
    library(scales)         # for transparent colors in plots; install.packages("scales")
    library(mvtnorm)        # for sampling multivariate normal random variables; install.packages("mvtnorm")
    library(latex2exp)      # for creating nice plot labels; install.packages("latex2exp"))
    library(posterior)      # for calculating updated MCMC diagnostics; install.packages("posterior")
    library(argparser)      # for parsing command line arguments for 02-model/fit-model.R script
    library(stringr)        # for miscellaneous string handling
  })
})
