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

# handle arguments supplied via command line
# run like: Rscript 02-model/fit-model.R base medium TRUE TRUE
args = commandArgs(trailingOnly = TRUE)
scenario_name = args[1]
mcmc_length = args[2]
do_sim_vs_obs = as.logical(args[3])
rmd = as.logical(args[4])

# scenario name: change this
if (is.na(scenario_name)) scenario_name = "siw"

# set mcmc length
if (is.na(mcmc_length)) mcmc_length = "vshort"

# include a forward simulation to compare to observed values as a validation?
if (is.na(do_sim_vs_obs)) do_sim_vs_obs = FALSE

# build RMDs?
if (is.na(rmd)) rmd = TRUE

# include posterior predictive checks in JAGS code?
do_pp_check = TRUE

# include posterior predictive density calculation in JAGS code?
do_lppd = FALSE

# build the scenario name
scenario = paste0(scenario_name, "_", mcmc_length, ifelse(do_sim_vs_obs, "_sim", ""))

##### STEP 1: PREPARE DATA FOR JAGS #####

# build JAGS data object
jags_data = create_jags_data_mult(c("CAT", "LOS", "MIN", "UGR"), first_y = ifelse(!do_sim_vs_obs, 1991, 2000), last_y = 2019)

# remove trap abundance data in BY 2007 for UGR
# estimates impossibly low, causing problems with model
jags_data$Pa_obs["2007","fall-mig","UGR"] = NA
jags_data$sig_Pa_obs["2007","fall-mig","UGR"] = NA
jags_data$Mb_obs["2007","spring-mig","NOR","UGR"] = NA
jags_data$sig_Mb_obs["2007","spring-mig","NOR","UGR"] = NA

# add NA indices in the correct locations
jags_data = append_no_na_indices(jags_data)

# set assumed survival past sea lions to 1 in all years
# past time-trending patterns were causing problems with ocean survival estimation/simulation
jags_data$phi_SL[2:nrow(jags_data$phi_SL),] = 1

# proportion of spawners by age and population that are female
# assume 50% female for all fish aged 4 or 5; 0% for age-3
# constant for all years/populations/origin types
Omega = array(0.5, dim = c(jags_data$nk, jags_data$nj))
dimnames(Omega) = list(jags_data$kmin:jags_data$kmax, colnames(jags_data$Ra_obs))
Omega["3",] = 0

add_jags_data = list(
  f = c(1904, 3971, 4846),  # fecundity [female age]
  Omega = Omega             # proportion of spawners by age and population that are female
)

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
  o_nor = 1,      # natural origin are o = 1,
  o_hor = 2       # hatchery origin are o = 2,
)
add_jags_data = append(add_jags_data, add_jags_data3)

# set first year adult PIT tag data exist to inform BON -> LGR survival
add_jags_data = append(add_jags_data, list(first_x_LGR = min(which(!is.na(jags_data$x_LGR[,1])))))

# set upper boundaries for early Rb states
if (!do_sim_vs_obs) {
  add_jags_data = append(add_jags_data, list(max_Rb_init = c(50, 200, 200)))
} else {
  add_jags_data = append(add_jags_data, list(max_Rb_init = c(200, 1000, 1000)))
}

# add variables to scale and center quantities used in regressions
add_jags_data = append(add_jags_data, list(
  E_scale = 10000,
  L_Pb_center = apply(jags_data$L_Pb_obs, 2, mean, na.rm = TRUE),
  L_Pb_scale = apply(jags_data$L_Pb_obs, 2, sd, na.rm = TRUE),
  L_Mb_center = apply(jags_data$L_Mb_obs, 2, mean, na.rm = TRUE),
  L_Mb_scale = apply(jags_data$L_Mb_obs, 2, sd, na.rm = TRUE)
))

# add hyperparameters of priors on ocean survival parameters
add_jags_data = append(add_jags_data, list(
  alpha_prior =               c(2, 8),
  mu_phi_O0_O1_prior =        c(1, 9),
  mu_phi_O1_O2_prior =        c(60, 40),
  mu_phi_O2_O3_prior =        c(70, 30),
  mu_psi_O1_prior =           c(1, 9),
  mu_psi_O2_prior =           c(8.5, 1.5),
  Tau_Lphi_E_Pb_prior =       c(0.30, 2), 
  Tau_lL_Pb_prior =           c(0.10, 2),
  Tau_Lpi_prior =             c(0.30, 2),
  Tau_Lphi_Pa_Mb_prior =      c(0.30, 2),
  Tau_lDelta_L_Pb_Mb_prior =  c(0.15, 500),  # model fails with Cholesky error if this has low DF
  Tau_Lphi_Mb_Ma_prior =      c(0.30, 2),
  Tau_Lphi_Ma_O0_prior =      c(0.30, 2),
  Tau_Lphi_O0_O1_prior =      c(0.30, 2),
  Tau_Lpsi_O1_prior =         c(0.15, 2),
  Tau_Lpsi_O2_prior =         c(0.30, 2),
  Tau_Lphi_Rb_Ra_prior =      c(0.30, 2),
  Tau_Lphi_Sb_Sa_prior =      c(0.15, 2)
))

# append all of this additional content to the data object
jags_data = append(jags_data, add_jags_data)

### ADD ON QUANTITIES FOR FORWARD SIMULATION ###

if (do_sim_vs_obs) {
  jags_data = append_values_for_sim(jags_data)
}

##### STEP 2: SPECIFY JAGS MODEL #####

# write the jags model source code to a text file
jags_source = "02-model/model-code.R"
jags_file = file.path("02-model", basename(tempfile(pattern = "model-", fileext = ".txt")))
write_model_code(jags_source, jags_file)

# toggle on/off the calculation of pp checks and lppd
toggle_data_diagnostics(do_lppd, do_pp_check, jags_file)

# toggle Rb_init for HOR fish
# (this is only for the sim vs. obs)
if (do_sim_vs_obs) {
  toggle_HOR_Rb_init(jags_file)
}

##### STEP 3: SELECT NODES TO MONITOR #####

# nodes to monitor for any model
jags_params = c(
  # reproduction
  "alpha", "beta", "phi_E_Pb", "lambda", "sig_lbeta", "kappa_phi_E_Pb", 

  # length-related quantities
  "omega0", "omega1", 
  "theta0", "theta1", "tau0", "tau1",
  
  # overwinter survival coefficients
  "gamma0", "gamma1",
  
  # hyperparameters: central tendency
  "mu_pi", "mu_phi_Mb_Ma", "mu_phi_Ma_O0",
  "mu_psi_O1", "mu_psi_O2", "mu_phi_Sb_Sa",
  "mu_phi_O0_O1", "mu_phi_O1_O2", "mu_phi_O2_O3", "mu_phi_Rb_Ra",
  
  # "dot" parameters -- expected values after accounting for relationships (dot), AR(1) (dot2)
  "phi_E_Pb_dot", "phi_E_Pb_dot2", "L_Pb_dot", "phi_Pa_Mb_dot",
  "Delta_L_Pb_Mb_dot", "phi_Mb_Ma_dot", "phi_O0_O1_dot2",
  
  # year-specific parameters
  "pi", "phi_Pa_Mb", "phi_Mb_Ma", "phi_Ma_O0", 
  "psi_O1", "psi_O2", "phi_Sb_Sa",
  "phi_O0_O1", "phi_O1_O2", "phi_O2_O3",
  "phi_Rb_Ra", "Delta_L_Pb_Mb",
  
  # derived survival terms
  "phi_Pb_Ma", "phi_Pa_Ma",
  
  # info for initialization
  "p_init", "mu_init_recruits", "sig_init_lrecruits", "init_recruits",
  
  # states
  "Pb", "Pa", "Mb", "Ma", "O0", "O", "Rb",
  "Ra", "Sb", "Sa", "p_Ra", "p_Sa_prime", "Sa_tot", "Ra_tot", "E", "E_sep", "L_Pb", "L_Mb",
  
  # carcass vs. weir correction
  "z", "zeta", "mu_z", "sig_z",
  
  # straying model
  "G", "p_G",
  
  # misc parameters
  "delta_O0_O1", "delta_O1_O2", "delta_O2_O3", "phi_Pb_Pa", "bad_flag",
  
  # misc derived quantities
  "lambda_pop", "Pb_per_Sa_tot", "Pb_per_E", "Mb_per_Sa_tot", "Sa_tot_per_Sa_tot", "Ra_per_Ma", "phi_O0_Rb_BON",
  "E_per_Sa",
  
  # residuals: process model
  "Lpi_resid", "Lphi_Pa_Mb_resid", "Lphi_Mb_Ma_resid",
  "Lphi_Ma_O0_resid", "Lpsi_O1_resid",
  "Lpsi_O2_resid", "Lphi_O0_O1_resid",
  "Lphi_Sb_Sa_resid", "Lphi_E_Pb_resid", "Lphi_Rb_Ra_resid",
  "lL_Pb_resid", "lDelta_L_Pb_Mb_resid",
  
  # quantile residuals: observation model
  "x_Ra_obs_qresid", "x_Sa_prime_obs_qresid",
  "lPa_obs_qresid", "lMb_obs_qresid", "lRa_obs_qresid",
  "Lphi_obs_Pb_Ma_qresid", "Lphi_obs_Pa_Ma_qresid", "Lphi_obs_Mb_Ma_qresid",
  "Lphi_obs_Ma_O0_qresid", "x_LGR_obs_qresid", "x_carcass_spawned_obs_qresid",
  "lL_Pb_obs_qresid", "lL_Mb_obs_qresid",
  
  # AR(1) coefficients
  "kappa_phi_O0_O1"
)

# nodes for posterior predictive checks
pp_check_params = c(
  "x_Ra_dev", "x_Ra_new_dev", "x_Sa_prime_dev", "x_Sa_prime_new_dev",
  "Pa_obs_dev", "Pa_obs_new_dev", "Mb_obs_dev", "Mb_obs_new_dev", 
  "Ra_obs_dev", "Ra_obs_new_dev", "Lphi_obs_Pb_Ma_dev", "Lphi_obs_new_Pb_Ma_dev",
  "Lphi_obs_Pa_Ma_dev", "Lphi_obs_new_Pa_Ma_dev", "Lphi_obs_Mb_Ma_dev", "Lphi_obs_new_Mb_Ma_dev",
  "Lphi_obs_Ma_O0_dev", "Lphi_obs_new_Ma_O0_dev", "x_carcass_spawned_dev", "x_carcass_spawned_new_dev",
  "x_LGR_dev", "x_LGR_new_dev", "L_Pb_obs_dev", "L_Pb_obs_new_dev", "L_Mb_obs_dev", "L_Mb_obs_new_dev"
)

# nodes for log posterior predictive density
lppd_params = c(
  "x_Ra_lppd", "x_Sa_prime_lppd", "Pa_obs_lppd", "Mb_obs_lppd", "Ra_obs_lppd",
  "Lphi_obs_Pb_Ma_lppd", "Lphi_obs_Pa_Ma_lppd", "Lphi_obs_Mb_Ma_lppd",
  "Lphi_obs_Ma_O0_lppd", "x_carcass_spawned_lppd", "x_LGR_lppd"
)

# nodes for monitoring potentially non-vague prior densities
prior_params = c(
  "alpha_pr", "mu_phi_O0_O1_pr", "mu_phi_O1_O2_pr", "mu_phi_O2_O3_pr",
  "mu_psi_O1_pr", "mu_psi_O2_pr",
  "sig_Lphi_E_Pb_pr", "sig_Lphi_O0_O1_pr", 
  "sig_Lpsi_O1_pr", "sig_Lpsi_O2_pr"
)

# nodes for Tau matrices
Tau_params = c(
  "Tau_Lphi_E_Pb", "Tau_lL_Pb", "Tau_Lpi", "Tau_Lphi_Pa_Mb", "Tau_lDelta_L_Pb_Mb", 
  "Tau_Lphi_Mb_Ma", "Tau_Lphi_Ma_O0", "Tau_Lphi_O0_O1", "Tau_Lpsi_O1", "Tau_Lpsi_O2",
  "Tau_Lphi_Rb_Ra", "Tau_Lphi_Sb_Sa"
)

# add these additional nodes if included in JAGS model
jags_params = c(jags_params, prior_params, Tau_params, paste0(Tau_params, "_pr"))
if (do_pp_check) jags_params = c(jags_params, pp_check_params)
if (do_lppd) jags_params = c(jags_params, lppd_params)

##### STEP 4: SELECT MCMC ATTRIBUTES #####

jags_dims = list(
  n_post  = switch(mcmc_length, "vshort" = 100, "short" = 2000, "medium" = 24000, "long" = 50000, "vlong" = 100000),
  n_burn  = switch(mcmc_length, "vshort" = 10,  "short" = 1000, "medium" = 20000, "long" = 30000, "vlong" = 50000),
  n_thin  = switch(mcmc_length, "vshort" = 1,   "short" = 3,    "medium" = 8,     "long" = 10,    "vlong" = 25),
  n_chain = switch(mcmc_length, "vshort" = 2,   "short" = 4,    "medium" = 4,     "long" = 4,     "vlong" = 4),
  n_adapt = switch(mcmc_length, "vshort" = 100, "short" = 100,  "medium" = 1000,  "long" = 3000,  "vlong" = 5000),
  parallel = TRUE
)

##### STEP 5: GENERATE INITIAL VALUES #####

# set random seed so at least initial values are reproducible
set.seed(1234)

# generate initial values
jags_inits = lapply(1:jags_dims$n_chain, gen_initials, jags_data = jags_data)

##### STEP 6: CALL THE JAGS SAMPLER #####

# print a start message
cat("Running JAGS on Multi-Population Model\n")
cat("  MCMC Length Selected: ", mcmc_length, "\n", sep = "")
cat("  Simulating years for validation: ", ifelse(do_sim_vs_obs, "Yes", "No"), "\n", sep = "")
cat("  Simulating data for posterior predictive check: ", ifelse(do_pp_check, "Yes", "No"), "\n", sep = "")
cat("  Monitoring log posterior predictive density for WAIC: ", ifelse(do_lppd, "Yes", "No"), "\n", sep = "")
starttime = Sys.time()
cat("  MCMC started:", format(starttime), "\n")

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
cat("  MCMC ended:", format(stoptime), "\n")
cat("  MCMC elapsed:", format(round(stoptime - starttime, 2)), "\n")

##### STEP 7: PROCESS COVARIANCE MATRICES #####

cat("Decomposing all Precision Matrices\n")

# a function to get the posterior of the vcov components (sig and rhos) from a precision matrix
process_Tau = function(Tau_name) {
  cat("  ", Tau_name, sep = "")
  suppressMessages({
    vcov_decomp(
      post, param = paste0(Tau_name, "["), invert = TRUE,
      sigma_base_name = stringr::str_replace(Tau_name, "Tau", "sig"),
      rho_base_name = stringr::str_replace(Tau_name, "Tau", "rho")
    )
  })
}

# extract the names of all Tau nodes
Tau_names = match_params(post, "^Tau", type = "base_only")

# apply the above function to each one
Sig_components = lapply(Tau_names, function(x) {out = process_Tau(x); cat("\n"); return(out)})

# combine with other posterior samples
Sig_components = lapply(Sig_components, as.matrix)
Sig_components = do.call(cbind, Sig_components)
post = post_bind(post, Sig_components)

# remove all precision matrix nodes: no longer needed
post = suppressMessages(post_remove(post, "^Tau", ask = FALSE))

# remove all diagonal elements of correlation matrices
# these are all 1 for all MCMC draws
diagonal_rhos = sub_index(sub_index("rho_.+[pop,pop]", pop = 1:jags_data$nj), pop = 1:jags_data$nj)
post = suppressMessages(post_remove(post, diagonal_rhos, ask = FALSE))

##### STEP 8: SAVE THE OUTPUT #####

# create the output directory if it doesn't already exist
if (!dir.exists(out_dir)) dir.create(out_dir)

# create the output object: stores data, posterior samples, and population name
out_obj = list(
  jags_model_code = readLines(jags_file),
  jags_data = jags_data,
  jags_inits = jags_inits,
  jags_dims = jags_dims,
  do_lppd = do_lppd,
  do_pp_check = do_pp_check,
  do_sim_vs_obs = do_sim_vs_obs,
  jags_time = c(starttime = format(starttime), stoptime = format(stoptime), elapsed = format(round(stoptime - starttime,2))),
  post = post,
  scenario = scenario
)

# create the output file name
out_file = paste0("output-", scenario, ".rds")

# save the file
cat("Saving rds Output\n")
saveRDS(out_obj, file.path(out_dir, out_file))

# delete the text file that contains the JAGS code
unlink(jags_file)

##### STEP 9: RENDER RMDs #####

# render the output plots if requested
if (rmd) {
  # start a timer, this can take a while
  starttime = Sys.time()
  
  # print a progress message
  cat("Rendering Rmd Output\n")
  
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

# render the simulation vs. observed time series plots if requested and applicable
if (rmd & do_sim_vs_obs) {
  # start a timer, this can take a while
  starttime = Sys.time()
  
  # print a progress message
  cat("\nRendering Sim vs. Obs Rmd Output")
  
  # set working dir to post-processing directory
  setwd("03-post-process")
  
  # file name of rendered output
  rmd_out_file = paste0("sim-vs-obs-", scenario, ".html")
  
  # render the output
  render(input = "compare-sim-to-obs.Rmd",
         output_file = rmd_out_file,
         params = list(scenario = scenario),
         quiet = TRUE
  )
  
  # open the rendered file when complete
  file.show(rmd_out_file)
  
  # stop the timer
  stoptime = Sys.time()
  cat("\nSim vs. Obs Rmd elapsed:", format(round(stoptime - starttime, 2)))
  
  # set the working dir back
  setwd("../")
  cat("\n\nDone.")
}

