This subdirectory contains scripts for fitting the state-space model.

* `fit-model.R`: This is the main script that is called to fit the model. Can be executed interactively, but also includes command line arguments. Run:
    ```bash
    Rscript 02-model/fit-model.R --help
    ```
* `model-code.R`: This file contains the code that specifies the state-space model in the JAGS language. Many comments are placed throughout defining what each section does.

After MCMC is done, a subdirectory called `model-output` will be created where the output (e.g., input data, output MCMC samples, miscellaneous info about the model run) will be saved.
