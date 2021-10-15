# :::::::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT TO RUN THE JAGS MODEL FOR ALL POPULATIONS #
# :::::::::::::::::::::::::::::::::::::::::::::::: #

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
scenario = "no-sex"

# handle command line arguments
# run this script via command line: Rscript 02-model/fit-model.R LOS TRUE
# or if in interactive session, uncomment the pop you wish to fit
args = commandArgs(trailingOnly = T)
rmd = as.logical(args[1])
mcmc_length = args[2]

if (is.na(rmd)) {
  rmd = FALSE
  cat("\n\n'rmd' was not supplied as a command line argument.", rmd, "will be used.")
}

if (is.na(mcmc_length)) {
  mcmc_length = "short"
  cat("\n\n'mcmc_length' was not supplied as a command line argument.", mcmc_length, "will be used.")
}

##### STEP 1: PREPARE DATA FOR JAGS #####

# build JAGS data object
jags_data = create_jags_data_mult(c("CAT", "LOS", "MIN", "UGR"))
jags_data = append_no_na_indices(jags_data)

# some parameters we are assuming known (for now)
# harvest rates Below and Above BON by year, age, origin
Ub_45_nat = with(jags_data, c(rep(NA, kmax), rep(0.02, ny - kmax)))
Ub_45_hat = with(jags_data, c(rep(NA, kmax), rep(0.08, ny - kmax)))
Ub = abind(list(cbind(Ub_45_nat * 0.25, Ub_45_nat, Ub_45_nat), cbind(Ub_45_hat * 0.75, Ub_45_hat, Ub_45_hat)), along = 3)
dimnames(Ub) = with(jags_data, list(rownames(Ra_obs), kmin:kmax, c("Nat", "Hat")))

# proportion of spawners by age and population that are female
p_female = array(NA, dim = c(jags_data$nk, jags_data$nj))
dimnames(p_female) = list(jags_data$kmin:jags_data$kmax, colnames(jags_data$Ra_obs))
p_female["3",] = 0
p_female["4","CAT"] = 0.55; p_female["4","LOS"] = 0.53; p_female["4","UGR"] = 0.57 
p_female["5","CAT"] = 0.43; p_female["5","LOS"] = 0.41; p_female["5","UGR"] = 0.48 
p_female["4","MIN"] = round(mean(p_female["4",], na.rm = TRUE), 2)
p_female["5","MIN"] = round(mean(p_female["5",], na.rm = TRUE), 2)

add_jags_data = list(
  Ub = Ub,                  # harvest rate after sea lion mortality below BON
  f = c(1904, 3971, 4846),  # fecundity [female age]
  p_female = p_female       # proportion of spawners by age and population that are female
)

# calculate a metric for overall survival from estuary to tributary
# used in obtaining an upper bound for initial adult recruits
overall_phi_Rb_Ra = with(append(jags_data, add_jags_data), mean(phi_SL, na.rm = T) * mean(LGR_adults[,"Nat"]/BON_adults[,"Nat"], na.rm = TRUE))

# some dummy variables for performing weir vs. carcass composition correction
add_jags_data2 = list(
  age3 = c(1,0,0),      # is each k age 3? 1 = yes; 0 = no
  age5 = c(0,0,1),      # is each k age 5?
  j_z = c(1,2,4),       # which j elements will have z estimated (must have both carcass and weir composition data)?
  nj_z = 3              # how many j elements will have z estimated?
)
add_jags_data = append(add_jags_data, add_jags_data2)

# create index names, for when it would help improve readability of the JAGS code
add_jags_data3 = list(
  i_fall = 1,     # fall migrants are i = 1
  i_spring = 2,   # spring migrants are i = 2,
  o_nat = 1,      # natural origin are o = 1,
  o_hat = 2       # hatchery origin are o = 2,
)
add_jags_data = append(add_jags_data, add_jags_data3)

# calculate the upper bound on initial adult recruits and add to data
p_NOR = apply(jags_data$carc_x_obs, 3, function(y) {rowSums(t(apply(y, 1, function(z) z/sum(z)))[,1:6])})
add_jags_data = append(add_jags_data, list(max_init_recruits = apply(((jags_data$Ra_obs * p_NOR)/overall_phi_Rb_Ra * 1.5), 2, max, na.rm = TRUE),
                                           first_LGR_adults = min(which(!is.na(jags_data$LGR_adults[,1])))))

# append all of this additional content to the data object
jags_data = append(jags_data, add_jags_data)

##### STEP 2: SPECIFY JAGS MODEL #####

# source the jags model code and write to a text file
source("02-model/model-code.R")
jags_file = "02-model/model.txt"
write_model(jags_model_code, jags_file)

# store the model code as an object, to be saved in the output object later
# use this version rather than 'model.txt' because this one retains the comments
jags_model_contents = readLines("02-model/model-code.R")
jags_model_contents[1] = "model {"

##### STEP 3: SELECT NODES TO MONITOR #####

jags_params = c(
  # reproduction
  "alpha", "beta", "sigma_Pb",
  
  # overwinter survival coefficients
  "gamma0", "gamma1",
  
  # hyperparameters: central tendency
  "mu_pi", "mu_phi_Mb_Ma", "mu_phi_Ma_O0",
  "mu_psi_O1_Rb", "mu_psi_O2_Rb", "mu_phi_Sb_Sa",
  "mu_phi_O0_O1", "mu_phi_O1_O2", "mu_phi_O2_O3", "mu_phi_Rb_Ra",
  
  # hyperparameters: inter-annual sd
  "sig_Lpi", "sig_Lphi_Pa_Mb", "sig_Lphi_Mb_Ma", "sig_Lphi_Ma_O0",
  "sig_Lpsi_O1_Rb", "sig_Lpsi_O2_Rb", "sig_Lphi_Sb_Sa",
  "Sig_Lphi_O0_O1", "sig_Lphi_Rb_Ra",
  
  # year-specific parameters
  "pi", "phi_Pa_Mb", "phi_Mb_Ma", "phi_Ma_O0", 
  "psi_O1_Rb", "psi_O2_Rb", "phi_Sb_Sa",
  "phi_O0_O1", "phi_O1_O2", "phi_O2_O3",
  "phi_Rb_Ra",
  
  # derived survival terms
  "phi_Pb_Ma", "phi_Pa_Ma",
  
  # info for initialization
  "p_init", "mu_init_recruits", "sig_init_lrecruits", "init_recruits",
  
  # states
  "Pb", "Pa", "Mb", "Ma", "O0", "O", "Rb",
  "Ra", "Sb", "Sa", "q_Ra", "q_Sa_adj", "Sa_tot", "Ra_tot", "f_tot",
  
  # carcass vs. weir correction
  "z", "carc_adj", "n_stray_tot", "stray_comp",
  
  # misc parameters
  "O_phi_scaler_nat_hat",
  
  # misc derived quantities
  "beta_per_peu", "Pb_per_Sa_tot", "Pb_per_f_tot", "Mb_per_Sa_tot", "Sa_tot_per_Sa_tot", "Ra_per_Ma", "phi_O0_Rb_BON",
  
  # residuals
  "Lpi_resid", "Lphi_Pa_Mb_resid", "Lphi_Mb_Ma_resid",
  "Lphi_Ma_O0_resid", "Lpsi_O1_Rb_resid",
  "Lpsi_O2_Rb_resid", "Lphi_O0_O1_resid", "Lphi_O1_O2_resid",
  "Lphi_Sb_Sa_resid", "lPb_resid", "Lphi_Rb_Ra_resid",
  
  # AR(1) coefficients
  "kappa_phi_O0_O1",
  
  # fit statistics for posterior predictive checks
  "weir_x_obs_dev", "weir_x_obs_new_dev", "carc_x_obs_dev", "carc_x_obs_new_dev",
  "Pa_obs_dev", "Pa_obs_new_dev", "Mb_obs_dev", "Mb_obs_new_dev", 
  "Ra_obs_dev", "Ra_obs_new_dev", "Lphi_obs_Pb_Ma_dev", "Lphi_obs_new_Pb_Ma_dev",
  "Lphi_obs_Pa_Ma_dev", "Lphi_obs_new_Pa_Ma_dev", "Lphi_obs_Mb_Ma_dev", "Lphi_obs_new_Mb_Ma_dev",
  "Lphi_obs_Ma_O0_dev", "Lphi_obs_new_Ma_O0_dev", "carcs_spawned_dev", "carcs_spawned_new_dev",
  "LGR_adults_dev", "LGR_adults_new_dev",
  
  # log posterior predictive density
  "weir_x_obs_lppd", "carc_x_obs_lppd", "Pa_obs_lppd", "Mb_obs_lppd", "Ra_obs_lppd",
  "Lphi_obs_Pb_Ma_lppd", "Lphi_obs_Pa_Ma_lppd", "Lphi_obs_Mb_Ma_lppd",
  "Lphi_obs_Ma_O0_lppd", "carcs_spawned_lppd", "LGR_adults_lppd"
  
)

##### STEP 4: SELECT MCMC ATTRIBUTES #####

jags_dims = list(
  n_post = switch(mcmc_length,  "short" = 2000,  "medium" = 24000,  "long" = 60000),
  n_burn = switch(mcmc_length,  "short" = 1000,  "medium" = 20000,  "long" = 60000),
  n_thin = switch(mcmc_length,  "short" = 3,     "medium" = 8,      "long" = 20),
  n_chain = switch(mcmc_length, "short" = 3,     "medium" = 3,      "long" = 3),
  n_adapt = switch(mcmc_length, "short" = 1000,  "medium" = 1000,   "long" = 1000),
  parallel = TRUE
)

##### STEP 5: GENERATE INITIAL VALUES #####

jags_inits = lapply(1:jags_dims$n_chain, gen_initials, jags_data)

##### STEP 6: CALL THE JAGS SAMPLER #####

# print a start message
cat("\n\nRunning JAGS on Multi-Population Model")
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
  verbose = FALSE
)

# print a stop message
stoptime = Sys.time()
cat("\nMCMC ended:", format(stoptime))
cat("\nMCMC elapsed:", format(round(stoptime - starttime, 2)))

# delete the text file that contains the JAGS code
unlink(jags_file)

##### STEP 7: SAVE THE OUTPUT #####

# create the output directory if it doesn't already exist
if (!dir.exists(out_dir)) dir.create(out_dir)

# create the output object: stores data, posterior samples, and population name
out_obj = list(
  jags_model_code = jags_model_contents,
  jags_data = jags_data,
  jags_inits = jags_inits,
  jags_dims = jags_dims,
  jags_time = c(starttime = format(starttime), stoptime = format(stoptime), elapsed = format(round(stoptime - starttime,2))),
  post = post,
  scenario = scenario
)

# create the output file name
out_file = paste0("output-", scenario, ".rds")

# save the file
cat("\nSaving rds Output")
saveRDS(out_obj, file.path(out_dir, out_file))

# render the output plots if requested
if (rmd) {
  # start a timer, this can take a while
  starttime = Sys.time()
  
  # print a progress message
  cat("\nRendering Rmd Output")
  
  # set working dir to post-processing directory
  setwd("03-post-process")
  
  # file name of rendered output
  rmd_out_file = paste0("output-plots-", scenario, ".html")
  
  # render the output
  render(input = "output-plots.Rmd",
         output_file = rmd_out_file,
         params = list(scenario = scenario),
         quiet = TRUE
  )
  
  # open the rendered file when complete
  file.show(rmd_out_file)
  
  # stop the timer
  stoptime = Sys.time()
  cat("\nRmd elapsed:", format(round(stoptime - starttime, 2)))
  
  # set the working dir back
  setwd("../")
  cat("\n\nDone.")
}
