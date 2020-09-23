# ::::::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT TO RUN THE JAGS MODEL FOR ONE POPULATION #
# ::::::::::::::::::::::::::::::::::::::::::::::: #

##### STEP 0: SET UP WORKSPACE #####

# load biological data (also loads packages)
source("00-data/prep-bio-data.R")

# load all necessary functions
invisible(sapply(list.files(path = "01-functions", pattern = "\\.R$", full.names = T), source))

# set the output directory: if it doesn't exist, a new directory will be created
out_dir = "02-model/model-output"

##### STEP 1: PREPARE DATA FOR JAGS #####

# select the population to fit
# pop = "UGR"
# pop = "CAT"
# pop = "MIN"  # don't do this one yet, the data prep steps (for age comp mostly) aren't correct for MIN
pop = "LOS"

# build JAGS data object
jags_data = create_jags_data_one(pop)
jags_data = append_no_na_indices(jags_data)

# some parameters we are assuming known (for now)
add_jags_data = list(
  mu_phi_M_O1 = 0.2,    # survival from arrival to estuary to next spring (become SWA1)
  mu_phi_O1_O2 = 0.8,   # survival from SWA1 to SWA2
  mu_phi_O2_O3 = 0.8,   # survival from SWA2 to SWA3
  mu_phi_Rb_Ra = 0.7   # survival upstream as adults in-river. mortality sources: sea lions, fishery, hydrosystem
)

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
  
  # hyperparameters: central tendency
  "mu_pi", "mu_phi_Pa_Mb", "mu_phi_Mb_Ma", "mu_phi_Ma_M",
  "mu_omega", "mu_psi_O1_Rb", "mu_psi_O2_Rb",
  
  # hyperparameters: inter-annual sd
  "sig_Lpi", "sig_Lphi_Pa_Mb", "sig_Lphi_Mb_Ma", "sig_Lphi_Ma_M",
  "sig_Lomega", "sig_Lpsi_O1_Rb", "sig_Lpsi_O2_Rb",
  
  # year-specific parameters
  "pi", "phi_Pa_Mb", "phi_Mb_Ma", "phi_Ma_M", "omega", 
  "psi_O1_Rb", "psi_O2_Rb",
  
  # derived survival terms
  "phi_Pb_Ma", "phi_Pa_Ma",
  
  # info for initialization
  "p_init", "mu_init_recruits", "sig_init_lrecruits", "init_recruits",
  
  # states
  "Pb", "Pa", "Mb", "Ma", "M", "O", "Rb",
  "Ra", "Sb", "Sa", "q", "Sa_tot", "Ra_tot"
)

##### STEP 4: SELECT MCMC ATTRIBUTES #####

jags_dims = list(
  n_post = 10000,
  n_burn = 5000,
  n_thin = 2, 
  n_chain = 2,
  n_adapt = 1000,
  parallel = F
)

##### STEP 5: GENERATE INITIAL VALUES #####

jags_inits = lapply(1:jags_dims$n_chain, gen_initials, jags_data)

##### STEP 6: CALL THE JAGS SAMPLER #####

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
  parallel = jags_dims$parallel
)

# delete the text file that contains the JAGS code
unlink(jags_file)

##### STEP 7: SAVE THE OUTPUT #####

if (!dir.exists(out_dir)) dir.create(out_dir)

# create the output object: stores data, posterior samples, and population name
out_obj = list(
  jags_data = jags_data,
  post = post,
  pop = pop
)

# create the output file name
out_file = paste0(pop, "-output.rds")

# save the file
saveRDS(out_obj, file.path(out_dir, out_file))
