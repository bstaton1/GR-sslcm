# :::::::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT TO RUN THE JAGS MODEL FOR ALL POPULATIONS #
# :::::::::::::::::::::::::::::::::::::::::::::::: #

##### STEP 0: SET UP WORKSPACE W/SETTINGS #####

# load all necessary functions
invisible(sapply(list.files(path = "01-functions", pattern = "\\.R$", full.names = T), source))

# print initial starting message
my_cat("Preparing Objects for JAGS:", "...please wait", first = TRUE)

# load biological data (also loads packages)
source("00-data/prep-bio-data.R")

# load environmental/habitat data
source("00-data/prep-env-data.R")

# start a timer for all calculations
starttime_all = Sys.time()

# load all necessary functions (again, prep-bio-dat clears workspace)
invisible(sapply(list.files(path = "01-functions", pattern = "\\.R$", full.names = T), source))

# set the output directory: if it doesn't exist, a new directory will be created
out_dir = "02-model/model-output"

# the list of acceptable mcmc length settings
accepted_mcmc = c("vshort", "short", "medium", "long", "vlong")

# set defaults: if running interactively, adjust settings here
default_args = list(
  scenario    = "test",
  mcmc        = "vshort",
  rmd         = TRUE,
  sim         = TRUE,
  pp_check    = TRUE,
  lppd        = FALSE
)

# Build command line argument parser
# provides defaults and instructions
parser = arg_parser("Fit the GR-sslcm to Data", hide.opts = TRUE) |>
  
  add_argument("scenario", "Name of model run",
               type = "character", default = default_args$scenario) %>%
  
  add_argument("--mcmc", paste0("Overall run length, ranging from ~5mins to ~24hrs. Accepted options are: ", knitr::combine_words(accepted_mcmc, before = "'", and = " or ")),
               type = "character", default = default_args$mcmc) %>%
  
  add_argument("--sim", "Perform simulation as a validation?",
               type = "logical", default = default_args$sim, short = "sim") %>%
  
  add_argument("--rmd", "Render Rmarkdown output?",
               type = "logical", default = default_args$rmd, short = "rmd") %>%
  
  add_argument("--pp_check", "Include posterior predictive check calculations?",
               type = "logical", default = default_args$pp_check, short = "pp_check") %>%
  
  add_argument("--lppd", "Include log posterior predictive density calculations (for WAIC)?",
               type = "logical", default = default_args$lppd, short = "lppd") 

# handle arguments based on whether this is an interactive session or not
# if interactive, use all defaults
# code also discards irrelevant values
if (interactive()) {
  defaults = parser$defaults
  arg_names = parser$args
  arg_names = str_remove(arg_names, "--")
  names(defaults) = arg_names
  args = defaults[!(names(defaults) %in% c("", "help"))]
} else {
  args = parse_args(parser)
  args = args[!(names(args) %in% c("", "help"))]
}

# stop if MCMC is not accepted
if (!(args$mcmc %in% accepted_mcmc)) {
  stop ("value supplied to argument '--mcmc' must be one of: ",
        knitr::combine_words(accepted_mcmc, before = "'", and = " or "))
}

# build the scenario name
scenario = paste0(args$scenario, "_", args$mcmc, ifelse(args$sim, "_sim", ""))

# build and print the output file name
out_file = paste0("output_", scenario, ".rds")
my_cat("Output File: ", out_file)

##### STEP 1: PREPARE DATA FOR JAGS #####

# build JAGS data object
jags_data = create_jags_data_mult(c("CAT", "LOS", "MIN", "UGR"), first_y = ifelse(!args$sim, 1991, 2000), last_y = 2019)

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
if (!args$sim) {
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
  omega0_prior =          log(c(50, 100)),
  omega1_prior =              c(-0.5, 0.5),
  theta0_prior =          log(c(1, 2)),
  theta1_prior =              c(-0.5, 0.5)
))

# append all of this additional content to the data object
jags_data = append(jags_data, add_jags_data)

### ADD ON QUANTITIES FOR FORWARD SIMULATION ###

if (args$sim) {
  jags_data = append_values_for_sim(jags_data)
}

##### STEP 2: SPECIFY JAGS MODEL #####

# write the jags model source code to a text file
jags_source = "02-model/model-code.R"
jags_file = file.path("02-model", basename(tempfile(pattern = "model-", fileext = ".txt")))
write_model_code(jags_source, jags_file)

# toggle on/off the calculation of pp checks and lppd
toggle_data_diagnostics(args$lppd, args$pp_check, jags_file)

# toggle Rb_init for HOR fish
# (this is only for the sim vs. obs)
if (args$sim) {
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
  "Lphi_E_Pb_resid", "Lphi_Rb_Ra_resid",
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
  "Tau_Lphi_Rb_Ra"
)

# add these additional nodes if included in JAGS model
jags_params = c(jags_params, prior_params, Tau_params, paste0(Tau_params, "_pr"))
if (args$pp_check) jags_params = c(jags_params, pp_check_params)
if (args$lppd) jags_params = c(jags_params, lppd_params)

##### STEP 4: SELECT MCMC ATTRIBUTES #####

jags_dims = list(
  n_post  = switch(args$mcmc, "vshort" = 100, "short" = 2000, "medium" = 24000, "long" = 50000, "vlong" = 100000),
  n_burn  = switch(args$mcmc, "vshort" = 10,  "short" = 1000, "medium" = 20000, "long" = 30000, "vlong" = 50000),
  n_thin  = switch(args$mcmc, "vshort" = 1,   "short" = 3,    "medium" = 8,     "long" = 10,    "vlong" = 25),
  n_chain = switch(args$mcmc, "vshort" = 2,   "short" = 4,    "medium" = 4,     "long" = 4,     "vlong" = 4),
  n_adapt = switch(args$mcmc, "vshort" = 100, "short" = 100,  "medium" = 1000,  "long" = 3000,  "vlong" = 5000),
  parallel = TRUE
)

# calculate expected time to run and misc values for nice printing
hrs_per_10k = ifelse(args$sim, 0.9, 1.85)
iters = with(jags_dims, (n_burn + n_post) * ifelse(parallel, 1, n_chain))
pred_hrs = iters/1e4 * hrs_per_10k
units_c = ifelse(pred_hrs > 1, "hour(s)", "minute(s)")
units_m = ifelse(pred_hrs > 1, 1, 60)
units_r = ifelse(pred_hrs > 1, 2, 0)

##### STEP 5: GENERATE INITIAL VALUES #####

# set random seed so at least initial values are reproducible
set.seed(1234)

# generate initial values
jags_inits = lapply(1:jags_dims$n_chain, gen_initials, jags_data = jags_data)

##### STEP 6: CALL JAGS #####

# print summary of what is about to run
starttime_jags = Sys.time()
my_cat("MCMC Length Selected:", args$mcmc)
my_cat("Predicted MCMC Duration:", paste(round(pred_hrs * units_m, digits = units_r), units_c))
my_cat("Sim vs. Obs Comparison:", ifelse(args$sim, "Yes", "No"))
my_cat("Posterior Predictive Check:", ifelse(args$pp_check, "Yes", "No"))
my_cat("WAIC: ", ifelse(args$lppd, "Yes", "No"))
my_cat("Scenario Name: ", scenario)
my_cat("MCMC Started:", format(starttime_jags))
my_cat("MCMC Completion (predicted):", format(starttime_jags + 3600 * pred_hrs))
my_cat("MCMC Running:", "...please wait")

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

# print a stop message upon completion
stoptime_jags = Sys.time()
my_cat("MCMC Completion (actual):", format(stoptime_jags))
my_cat("    Elapsed:", format(round(stoptime_jags - starttime_jags, 1)))

##### STEP 7: HANDLE PROCESS COVARIANCE MATRICES #####

# print a start message
starttime = Sys.time()
my_cat("Tau Matrix Decomposition:", "...please wait")

# a function to get the posterior of the vcov components (sig and rhos) from a precision matrix
process_Tau = function(Tau_name) {
  cat("\r", paste(rep(" ", 100), collapse = ""), "\r        (", which(Tau_names == Tau_name), "/", length(Tau_names), ") ", Tau_name, sep = "")
  suppressMessages({
    vcov_decomp(
      post, param = paste0(Tau_name, "["), invert = TRUE,
      sigma_base_name = str_replace(Tau_name, "Tau", "sig"),
      rho_base_name = str_replace(Tau_name, "Tau", "rho")
    )
  })
}

# extract the names of all Tau nodes
Tau_names = match_params(post, "^Tau", type = "base_only")

# apply the above function to each one
Sig_components = lapply(Tau_names, function(x) {out = process_Tau(x); return(out)}); cat("\n")

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

# print a stop message when complete
stoptime = Sys.time()
my_cat("    Elapsed:", format(round(stoptime - starttime, 1)))

##### STEP 8: SAVE THE OUTPUT #####

starttime = Sys.time()
my_cat("Saving output file:", "...please wait")

# create the output directory if it doesn't already exist
if (!dir.exists(out_dir)) dir.create(out_dir)

# create the output object: stores data, posterior samples, and population name
out_obj = list(
  jags_model_code = readLines(jags_file),
  jags_data = jags_data,
  jags_inits = jags_inits,
  jags_dims = jags_dims,
  do_lppd = args$lppd,
  do_pp_check = args$pp_check,
  do_sim_vs_obs = args$sim,
  jags_time = c(starttime = format(starttime_jags), stoptime = format(stoptime_jags), elapsed = format(round(stoptime_jags - starttime_jags,2))),
  post = post,
  scenario = scenario
)

# save the file
saveRDS(out_obj, file.path(out_dir, out_file))

# delete the text file that contains the JAGS code
unlink(jags_file)

# print a stop message when complete
stoptime = Sys.time()
my_cat("    Elapsed:", format(round(stoptime - starttime, 1)))

# render the output plots if requested
if (args$rmd) {
  
  # print a start message
  starttime = Sys.time()
  my_cat("Rendering Rmd File: ", "output-plots.Rmd")
  
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
  if (interactive()) file.show(rmd_out_file)
  
  # print a message when complete
  stoptime = Sys.time()
  my_cat("    Elapsed:", format(round(stoptime - starttime, 1)))
  
  # set the working dir back
  setwd("../")
}

# render the simulation vs. observed time series plots if requested and applicable
if (args$rmd & args$sim) {
  
  # print a starting message
  starttime = Sys.time()
  my_cat("Rendering Rmd File: ", "compare-sim-to-obs.Rmd")
  
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
  if (interactive()) file.show(rmd_out_file)
  
  # print
  stoptime = Sys.time()
  my_cat("    Elapsed:", format(round(stoptime - starttime, 1)))
  
  # set the working dir back
  setwd("../")
}

##### STEP 10: ALL DONE #####

# stop the timer for all analyses
stoptime_all = Sys.time()

# print a completion message
my_cat("All Steps Complete:", "HOORAY!!")
my_cat("    Elapsed:", format(round(stoptime_all - starttime_all, 1)))
