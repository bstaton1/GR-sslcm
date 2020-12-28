# ::::::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT TO RUN THE JAGS MODEL FOR ONE POPULATION #
# ::::::::::::::::::::::::::::::::::::::::::::::: #

##### STEP 0: SET UP WORKSPACE #####

# load biological data (also loads packages)
source("00-data/prep-bio-data.R")

# load environmental/habitat data
source("00-data/prep-env-data.R")

# load all necessary functions
invisible(sapply(list.files(path = "01-functions", pattern = "\\.R$", full.names = T), source))

# set the output directory: if it doesn't exist, a new directory will be created
out_dir = "02-model/model-output"

# specify a scenario name
scenario = "relate-hat-nat-age-4-5-same"

# handle command line arguments
# run this script via command line: Rscript 02-model/fit-model.R LOS TRUE
# or if in interactive session, uncomment the pop you wish to fit
args = commandArgs(trailingOnly = T)
pop = args[1]
rmd = as.logical(args[2])
mcmc_length = args[3]
rmd_mcmc_plots = T

# if pop not supplied, set a value and warn
if (is.na(pop)) {
  pop = "CAT"
  # pop = "LOS"
  # pop = "MIN" # don't do this one yet, the data prep steps (for age comp mostly) aren't correct for MIN
  # pop = "UGR"
  cat("\n\n'pop' was not supplied as a command line argument.", pop, " will be used.")
}

if (is.na(rmd)) {
  rmd = T
  # rmd = F
  cat("\n\n'rmd' was not supplied as a command line argument.", rmd, " will be used.")
}

if (is.na(mcmc_length)) {
  mcmc_length = "short"
  cat("\n\n'mcmc_length' was not supplied as a command line argument.", mcmc_length, " will be used.")
}

##### STEP 1: PREPARE DATA FOR JAGS #####

# build JAGS data object
jags_data = create_jags_data_one(pop)
jags_data = append_no_na_indices(jags_data)

# some parameters we are assuming known (for now)
# [1] is natural, [2] is hatchery origin
add_jags_data = list(
  mu_phi_M_O1_assume = c(0.2, 0.2),    # survival from arrival to estuary to next spring (become SWA1)
  mu_phi_O1_O2_assume = c(0.8, 0.8),   # survival from SWA1 to SWA2
  mu_phi_O2_O3_assume = c(0.8, 0.8),   # survival from SWA2 to SWA3
  mu_phi_Rb_Ra_assume = c(0.7, 0.7)    # survival upstream as adults in-river. mortality sources: sea lions, fishery, hydrosystem
)

# some dummy variables for performing weir vs. carcass composition correction
add_jags_data2 = list(
  age3 = c(1,0,0),      # is each k age 3? 1 = yes; 0 = no
  age5 = c(0,0,1),      # is each k age 5?
  male = c(0,1)         # is each s male?
)
add_jags_data = append(add_jags_data, add_jags_data2)

# the years in which strays will be needed
if (pop != "MIN") {
  yrs = as.numeric(names(jags_data$Mb_obs[,2,2]))
  first_brood_release = min(yrs[jags_data$Mb_obs[,2,2] > 0], na.rm = T)
  first_adult_return = first_brood_release + jags_data$kmax
  stray_yrs = 6:(which(yrs == (first_adult_return - 1)))
  not_stray_yrs = max(stray_yrs+1):jags_data$ny
} else {
  stray_yrs = 6:jags_data$ny
  not_stray_yrs = numeric(0)
}

add_jags_data3 = list(
  stray_yrs = stray_yrs,
  not_stray_yrs = not_stray_yrs,
  n_stray_yrs = length(stray_yrs),
  n_not_stray_yrs = length(not_stray_yrs)
)
add_jags_data = append(add_jags_data, add_jags_data3)

# create index names, for when it would help improve readibility of the JAGS code
add_jags_data4 = list(
  i_fall = 1,     # fall migrants are i = 1
  i_spring = 2,   # spring migrants are i = 2,
  o_nat = 1,      # natural origin are o = 1,
  o_hat = 2,      # hatchery origin are o = 2,
  s_female = 1,   # females are s = 1,
  s_male = 2     # males are s = 2
)
add_jags_data = append(add_jags_data, add_jags_data4)

# calculate the upper bound on initial adult recruits and add to data
add_jags_data = append(add_jags_data, list(max_init_recruits = max(jags_data$Ra_obs/(add_jags_data$mu_phi_Rb_Ra) * 1.5, na.rm = T)))
jags_data = append(jags_data, add_jags_data)

##### STEP 2: SPECIFY JAGS MODEL #####

# source the jags model code and write to a text file
source("02-model/model-code.R")
jags_file = "02-model/model.txt"
write_model(jags_model_code, jags_file)

##### STEP 3: SELECT NODES TO MONITOR #####

jags_params = c(
  # reproduction
  "alpha", "beta", "sigma_Pb",
  
  # overwinter survival coefficients
  "gamma0", "gamma1",
  
  # hyperparameters: central tendency
  "mu_pi", "mu_phi_Mb_Ma", "mu_phi_Ma_M",
  "mu_omega", "mu_psi_O1_Rb", "mu_psi_O2_Rb", "mu_phi_Sb_Sa",
  "mu_phi_M_O1", "mu_phi_O1_O2", "mu_phi_O2_O3",
  
  # hyperparameters: inter-annual sd
  "sig_Lpi", "sig_Lphi_Pa_Mb", "sig_Lphi_Mb_Ma", "sig_Lphi_Ma_M",
  "sig_Lomega", "sig_Lpsi_O1_Rb", "sig_Lpsi_O2_Rb", "sig_Lphi_Sb_Sa",
  "sig_Lphi_M_O1", "sig_Lphi_O1_O2", "sig_Lphi_O2_O3",
  
  # year-specific parameters
  "pi", "phi_Pa_Mb", "phi_Mb_Ma", "phi_Ma_M", "omega", 
  "psi_O1_Rb", "psi_O2_Rb", "phi_Sb_Sa",
  "phi_M_O1", "phi_O1_O2", "phi_O2_O3",
  
  # derived survival terms
  "phi_Pb_Ma", "phi_Pa_Ma",
  
  # info for initialization
  "p_init", "mu_init_recruits", "sig_init_lrecruits", "init_recruits",
  
  # states
  "Pb", "Pa", "Mb", "Ma", "M", "O", "Rb",
  "Ra", "Sb", "Sa", "q_Ra", "q_Sa_adj", "Sa_tot", "Ra_tot",
  
  # carcass vs. weir correction
  "z", "carc_adj", "n_stray_tot", "stray_comp",
  
  # misc parameters
  "O_phi_scaler_nat_hat"
)

##### STEP 4: SELECT MCMC ATTRIBUTES #####

jags_dims = list(
  n_post = switch(mcmc_length,  "short" = 1000, "medium" = 24000, "long" = 60000),
  n_burn = switch(mcmc_length,  "short" = 500,  "medium" = 20000,  "long" = 30000),
  n_thin = switch(mcmc_length,  "short" = 1,    "medium" = 8,     "long" = 20),
  n_chain = switch(mcmc_length, "short" = 3,    "medium" = 3,     "long" = 3),
  n_adapt = switch(mcmc_length, "short" = 100,  "medium" = 1000,  "long" = 1000),
  parallel = T
)

##### STEP 5: GENERATE INITIAL VALUES #####

jags_inits = lapply(1:jags_dims$n_chain, gen_initials, jags_data)

##### STEP 6: CALL THE JAGS SAMPLER #####

# print a start message
cat("\n\nRunning JAGS on Population:", pop)
starttime = Sys.time()
cat("\nMCMC started:", format(starttime))

# call JAGS
post = jags.basic(
  data = jags_data,
  inits = jags_inits,
  parameters.to.save = jags_params, 
  model.file = jags_file, 
  n.chains = jags_dims$n_chain, 
  n.adapt = jags_dims$n_adapt, 
  n.iter = jags_dims$n_post + jags_dims$n_burn, 
  n.burnin = jags_dims$n_burn, 
  n.thin = jags_dims$n_thin,
  parallel = jags_dims$parallel,
  verbose = F
)

# print a stop message
stoptime = Sys.time()
cat("\nMCMC ended:", format(stoptime))
cat("\nMCMC elapsed:", format(round(stoptime - starttime, 2)))

# delete the text file that contains the JAGS code
unlink(jags_file)

##### STEP 7: SAVE THE OUTPUT #####

if (!dir.exists(out_dir)) dir.create(out_dir)

# create the output object: stores data, posterior samples, and population name
out_obj = list(
  jags_data = jags_data,
  post = post,
  pop = pop,
  scenario = scenario
)

# create the output file name
out_file = paste0(pop, "-output-", scenario, ".rds")

# save the file
cat("\nSaving rds Output")
saveRDS(out_obj, file.path(out_dir, out_file))

# render the output plots if requested
if (rmd) {
  cat("\nRendering Rmd Output")
  setwd("03-post-process")
  rmarkdown::render(input = "output-plots.Rmd",
                    output_file = paste0(pop, "-output-plots-", scenario, ".html"), 
                    params = list(pop = pop, scenario = scenario, mcmc_plots = rmd_mcmc_plots), 
                    quiet = T
  )
  setwd("../")
  cat("\n\nDone.")
}
