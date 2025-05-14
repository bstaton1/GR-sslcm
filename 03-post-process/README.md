This subdirectory contains scripts to summarize, analyze, and visualize the state-space model output.
The content in this subdirectory is automatically built when executing:

```bash
Rscript 02-model/fit-model.R --sim TRUE
```

To disable building this content and only run MCMC, instead execute:

```bash
Rscript 02-model/fit-model.R --rmd FALSE --ms FALSE
```

* `output-plots.Rmd`: the main Rmarkdown file that builds the manuscript supplemental file.
* `output-plots-children/`: a subdirectory containing many Rmarkdown files, which are independent of one another, and each creates one section of the supplemental file.
* `compare-sim-to-obs.Rmd`: a Rmarkdown file that compares properties of the observed time period to those of a simulated time period as a means to validate the model. To see the contents of this file, run the model using:
    ```bash
    Rscript 02-model/fit-model.R --sim TRUE
    ```
    
    _This content is not presented or discussed in the manuscript._
* `ms-content/`: a subdirectory containing R scripts for creating the figures and tables presented in the manuscript. The script `in-text-quantities.R` in that subdirectory is used to calculate quantities presented in the manuscript text.
