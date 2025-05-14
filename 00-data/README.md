All data files for this project are found in the [GR-sslcm-data](<https://github.com/gibsonpp/GR-sslcm-data>) repo.

The scripts in this subdirectory are:

* `prep-bio-data.R`: performs the necessary formatting to translate the raw data files (from [GR-sslcm-data](<https://github.com/gibsonpp/GR-sslcm-data>)) into a large data frame (called `bio`) that can be more easily manipulated for passing to JAGS.
* `prep-env-data.R`: specifies the values of the weighted usable habitat length metric (WUL). The PEU (pool equivalent units) metric was used in previous versions of the model, and not in the final version.
