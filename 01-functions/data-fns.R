# :::::::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT TO HOUSE FUNCTIONS FOR PREPARING RAW DATA #
# :::::::::::::::::::::::::::::::::::::::::::::::: #

##### OBTAIN LOGIT-SCALE STANDARD ERROR OF A PROPORTION #####
# some estimates have p_mean and CIs
# some estimates have p_mean and p_se
# either way, obtain CI's of proportion, then convert to a standard error on the logit scale
get_logit_se = function(p_mean, p_se, p_lwr, p_upr, alpha) {
  # turn se's into ci if that is what is available
  p_lwr = ifelse(is.na(p_se), p_lwr, p_mean + qnorm(alpha/2) * p_se)
  p_upr = ifelse(is.na(p_se), p_upr, p_mean + qnorm(1 - alpha/2) * p_se)
  
  # cap the CI
  p_lwr = ifelse(p_lwr <= 0, 0.001, p_lwr)
  p_upr = ifelse(p_upr >= 1, 0.999, p_upr)
  
  # convert CI into logit-normal standard error: M. Liermann's approximation
  (logit(p_upr) - logit(p_lwr))/(2 * qnorm(1 - alpha/2))
}



##### STANDARDIZE MEAN LENGTH DATA BASED ON MEDIAN CAPTURE DATE #####

# the median capture date varies widely within populations
# and in general, sampling occured later in the year in the earlier years
# b/c fish are growing over this whole period, this caused an apparent declining trend in mean length over time

# this function standardizes mean length based on the median capture date
# by fitting a growth function, obtaining the expected value if all years had the same median date, and adding the growth residual

# ARGUMENTS
# len_mean: numeric vector of observed mean lengths, elements represent years
# jday_med: numeric vector of median capture dates, elements represent years
# resid_type: character; either "mult" or "add"; do you wish to use additive or multiplicative residuals

standardize_mean_length = function(len_mean, jday_med, resid_type) {
  
  # step 1: fit a regression model; approximates average daily growth function
  fit = lm(len_mean ~ jday_med)
  
  # step 2: calculate residuals for each year
  if (resid_type != "mult" & resid_type != "add") stop ("resid_type must be one of 'mult' or 'add'")
  if (resid_type == "mult") r = len_mean/predict(fit, newdata = data.frame(jday_med))
  if (resid_type == "add") r = len_mean - fitted(fit)
  
  # step 3: calculate expected length on the average sampling day across years
  mu_len_mean = predict(fit, data.frame(jday_med = mean(jday_med, na.rm = TRUE)))
  
  # step 4: calculate the sampling-date-adjusted mean length each year by applying the
  # year-specific residuals to the standardized expected mean length
  if (resid_type == "mult") len_mean_adj = mu_len_mean * r
  if (resid_type == "add") len_mean_adj = mu_len_mean + r
  
  # return the sampling date-adjusted mean length vector
  return(len_mean_adj)
}

##### CREATE NAMES FOR A SPECIFIC COMPOSITION DATA SET #####
# o_names = c("NOR", "HOR")
# k_names = c(3, 4, 5)
# type = "weir";# type = "carc"; type = "rm"

create_comp_names = function(type, o_names, k_names) {
  # create combinations of origins and ages
  x = expand.grid(o = o_names, k = k_names)
  
  # sort them by origin
  x = x[order(x$o),]
  
  # combine into strings, along with the type: carc, weir, or rm
  x = apply(x, 1, function(x) paste(c(type, x), collapse = "_"))
  
  # build a list with the names for each origin type
  list(
    nat_names = unname(x[1:(length(k_names))]),
    hat_names = unname(x[(length(k_names) + 1):(length(o_names) * length(k_names))])
  )
}

##### CREATE DATA LIST FOR JAGS: ONE POPULATION #####

# formats the appropriate information from the bio_dat object
# to be used by a JAGS model, one population only

# pop: "CAT", "LOS", "MIN", or "UGR"
# first_y: the first return year with adult return data for any population
# last_y: the last return year with adult return data for any population

create_jags_data_one = function(pop, first_y = 1991, last_y = 2019) {
  
  ## ERROR HANDLER: make sure bio_dat object is available
  if (!exists("bio_dat")) {
    stop ("the 'bio_dat' main data.frame does not exist. Run source('00-data/prep-bio-data.R') to create it")
  }
  
  ## ERROR HANDLERS: ensure years and population are in data
  if (first_y %!in% bio_dat$brood_year) {
    stop ("supplied value of 'first_y' outside of range found in data")
  }
  if (last_y %!in% bio_dat$brood_year) {
    stop ("supplied value of 'last_y' outside of range found in data")
  }
  if (pop %!in% bio_dat$population) {
    stop ("supplied value of 'pop' is not a population found in data")
  }
  
  # extract all records/variables meeting the query specifications
  sub = subset(bio_dat, brood_year >= first_y & brood_year <= last_y & population == pop)
  all = subset(bio_dat, brood_year >= first_y & brood_year <= last_y & population == "ALL")
  
  # specify/calculate dimensional variables
  kmin = 3              # minimum age of return
  kmax = 5              # maximum age of return
  nk = kmax - kmin + 1  # number of ages of return
  ni = 2                # number of juvenile life history strategies
  no = 2                # number of origins
  nko = nk * no         # number of age/origin combinations
  ny = nrow(sub) + 1    # number of brood years tracked (add 1 to leave empty spot for year = 0; for AR(1) processes which require a previous year's residual)
  
  # assign these names
  y_names = (first_y - 1):last_y
  k_names = kmin:kmax
  i_names = paste0(c("fall", "spring"), "-mig")
  o_names = c("NOR", "HOR")
  ko_names = c(paste0(k_names, "-", o_names[1]), paste0(k_names, "-", o_names[2]))
  
  ### JUVENILE ABUNDANCE DATA ###
  # fall trap count
  Pa_obs = matrix(NA, ny, ni); dimnames(Pa_obs) = list(y_names, i_names)
  Pa_obs[y_names %in% sub$brood_year,i_names == "fall-mig"] = sub$fall_passage_est
  
  # fall trap sd(count)
  sig_Pa_obs = matrix(NA, ny, ni); dimnames(sig_Pa_obs) = list(y_names, i_names)
  sig_Pa_obs[y_names %in% sub$brood_year,i_names == "fall-mig"] = sub$fall_passage_log_se
  
  # spring trap count
  Mb_obs = matrix(NA, ny, ni); dimnames(Mb_obs) = list(y_names, i_names)
  Mb_obs[y_names %in% sub$brood_year,i_names == "spring-mig"] = sub$spring_passage_est
  
  # spring trap sd(count)
  sig_Mb_obs = matrix(NA, ny, ni); dimnames(sig_Mb_obs) = list(y_names, i_names)
  sig_Mb_obs[y_names %in% sub$brood_year,i_names == "spring-mig"] = sub$spring_passage_log_se
  
  # hatchery releases
  Mb_rel = matrix(NA, ny, ni); dimnames(Mb_rel) = list(y_names, i_names)
  Mb_rel[y_names %in% sub$brood_year,i_names == "spring-mig"] = sub$hatchery_smolt
  
  # uncertainty in hatchery releases
  sig_Mb_rel = matrix(NA, ny, ni); dimnames(sig_Mb_rel) = list(y_names, i_names)
  
  # combine natural and hatchery Mb
  Mb_obs = abind(Mb_obs, Mb_rel, along = 3)
  sig_Mb_obs = abind(sig_Mb_obs, sig_Mb_rel, along = 3)
  dimnames(Mb_obs)[[3]]  = dimnames(sig_Mb_obs)[[3]] = o_names
  
  ### JUVENILE LENGTH DATA ###
  # mean length at summer tagging
  L_Pb_obs = rep(NA, ny); names(L_Pb_obs) = y_names
  L_Pb_obs[y_names %in% sub$brood_year] = sub$length_mean_summer
  
  # se mean length at summer tagging: calculate lognormal se by calculating "cv" then applying transformation
  # this se is so small (sample size nearly always ~1000), may want to consider down-weighting these?
  sig_L_Pb_obs = rep(NA, ny); names(sig_L_Pb_obs) = y_names
  sig_L_Pb_obs[y_names %in% sub$brood_year] = cv2sig(sub$length_se_summer/sub$length_mean_summer)
  
  # mean length at spring tagging
  L_Mb_obs = rep(NA, ny); names(L_Mb_obs) = y_names
  L_Mb_obs[y_names %in% sub$brood_year] = sub$length_mean_spring
  
  # se mean length at summer tagging: calculate lognormal se by calculating "cv" then applying transformation
  # this se is so small (sample size nearly always ~1000), may want to consider down-weighting these?
  sig_L_Mb_obs = rep(NA, ny); names(sig_L_Mb_obs) = y_names
  sig_L_Mb_obs[y_names %in% sub$brood_year] = cv2sig(sub$length_se_spring/sub$length_mean_spring)
  
  ### JUVENILE SURVIVAL DATA ###
  # summer logit(surv) to LGD
  Lphi_obs_Pb_Ma = rep(NA, ny); names(Lphi_obs_Pb_Ma) = y_names
  Lphi_obs_Pb_Ma[y_names %in% sub$brood_year] = logit(sub$summer_surv_est)
  
  # sd summer logit(surv) to LGD
  sig_Lphi_obs_Pb_Ma = rep(NA, ny); names(sig_Lphi_obs_Pb_Ma) = y_names
  sig_Lphi_obs_Pb_Ma[y_names %in% sub$brood_year] = sub$summer_surv_logit_se
  
  # fall & winter logit(surv) to LGD
  Lphi_obs_Pa_Ma = matrix(NA, ny, ni); dimnames(Lphi_obs_Pa_Ma) = list(y_names, i_names)
  Lphi_obs_Pa_Ma[y_names %in% sub$brood_year,i_names == "fall-mig"] = logit(sub$fall_surv_est)
  Lphi_obs_Pa_Ma[y_names %in% sub$brood_year,i_names == "spring-mig"] = logit(sub$winter_surv_est)
  
  # sd fall & winter logit(surv) to LGD
  sig_Lphi_obs_Pa_Ma = matrix(NA, ny, ni); dimnames(sig_Lphi_obs_Pa_Ma) = list(y_names, i_names)
  sig_Lphi_obs_Pa_Ma[y_names %in% sub$brood_year,i_names == "fall-mig"] = sub$fall_surv_logit_se
  sig_Lphi_obs_Pa_Ma[y_names %in% sub$brood_year,i_names == "spring-mig"] = sub$winter_surv_logit_se
  
  # spring logit(surv) to LGD
  Lphi_obs_Mb_Ma = array(NA, dim = c(ny, ni, no)); dimnames(Lphi_obs_Mb_Ma) = list(y_names, i_names, o_names)
  Lphi_obs_Mb_Ma[y_names %in% sub$brood_year,i_names == "spring-mig",o_names[1]] = logit(sub$spring_surv_est)
  Lphi_obs_Mb_Ma[y_names %in% sub$brood_year,i_names == "spring-mig",o_names[2]] = logit(sub$hatchery_spring_surv_est)
  
  # sd spring logit(surv) to LGD
  sig_Lphi_obs_Mb_Ma = array(NA, dim = c(ny, ni, no)); dimnames(sig_Lphi_obs_Mb_Ma) = list(y_names, i_names, o_names)
  sig_Lphi_obs_Mb_Ma[y_names %in% sub$brood_year,i_names == "spring-mig",o_names[1]] = sub$spring_surv_logit_se
  sig_Lphi_obs_Mb_Ma[y_names %in% sub$brood_year,i_names == "spring-mig",o_names[2]] = sub$hatchery_spring_surv_logit_se
  
  # LGD to BON logit(surv)
  Lphi_obs_Ma_O0 = matrix(NA, ny, no); dimnames(Lphi_obs_Ma_O0) = list(y_names, o_names)
  Lphi_obs_Ma_O0[y_names %in% all$brood_year,o_names[1]] = logit(all$nat_hydro_est)
  Lphi_obs_Ma_O0[y_names %in% all$brood_year,o_names[2]] = logit(all$hat_hydro_est)
  
  # sd LGD to BON logit(surv)
  sig_Lphi_obs_Ma_O0 = matrix(NA, ny, no); dimnames(sig_Lphi_obs_Ma_O0) = list(y_names, o_names)
  sig_Lphi_obs_Ma_O0[y_names %in% all$brood_year,o_names[1]] = all$nat_hydro_logit_se
  sig_Lphi_obs_Ma_O0[y_names %in% all$brood_year,o_names[2]] = all$hat_hydro_logit_se
  
  ### ADULT AGE COMP: WEIR ###
  # obtain names of age comp variables
  weir_comp_names = create_comp_names("weir", o_names, k_names)
  
  # extract them by origin and coerce NA to zero
  nat_comp = sub[,weir_comp_names$nat_names]
  hat_comp = sub[,weir_comp_names$hat_names]
  nat_comp[is.na(nat_comp)] = 0
  hat_comp[is.na(hat_comp)] = 0
  
  # age frequencies: for fitting composition of returns
  x_Ra = matrix(NA, ny, no * nk); dimnames(x_Ra) = list(y_names, ko_names)
  x_Ra[y_names %in% sub$brood_year,] = as.matrix(cbind(nat_comp, hat_comp))
  nx_Ra = rowSums(x_Ra)
  
  ### ADULT AGE COMP: CARCASSES ###
  # obtain names of age comp variables
  carc_comp_names = create_comp_names("carc", o_names, k_names)
  
  # extract them by origin and coerce NA to zero
  nat_comp = sub[,carc_comp_names$nat_names]
  hat_comp = sub[,carc_comp_names$hat_names]
  nat_comp[is.na(nat_comp)] = 0
  hat_comp[is.na(hat_comp)] = 0
  
  # age frequencies: for fitting composition of returns
  x_Sa_prime = matrix(NA, ny, no * nk); dimnames(x_Sa_prime) = list(y_names, ko_names)
  x_Sa_prime[y_names %in% sub$brood_year,] = as.matrix(cbind(nat_comp, hat_comp))
  nx_Sa_prime = rowSums(x_Sa_prime)
  
  ### NUMBER OF RETURNING ADULTS REMOVED AT WEIR ###
  rm_comp_names = create_comp_names("rm", o_names, k_names)
  
  # extract compositions by type and coerce NAs to zero
  # this is the number of fish removed at weir each year by age/origin class
  rm_comp = sub[,c(rm_comp_names$nat_names, rm_comp_names$hat_names)]
  rm_comp[is.na(rm_comp)] = 0
  
  # place rm_comp in the correct location of B (name in model): same numbers just reformatted array structure used by model
  B = array(NA, dim = c(ny, nk, no)); dimnames(B) = list(y_names, k_names, o_names)
  B[y_names %in% sub$brood_year,,o_names[1]] = as.matrix(rm_comp[,paste("rm", o_names[1], k_names, sep = "_")])
  B[y_names %in% sub$brood_year,,o_names[2]] = as.matrix(rm_comp[,paste("rm", o_names[2], k_names, sep = "_")])

  ### ADULT HARVEST RATE BELOW BON ###
  # assume age-3 is harvested at a rate 25% of that for age-4 and age-5
  U = array(NA, dim = c(ny, nk, no)); dimnames(U) = list(y_names, k_names, o_names)
  U[y_names %in% all$brood_year,"3","NOR"] = all$HR_NOR * 0.25
  U[y_names %in% all$brood_year,c("4", "5"),"NOR"] = all$HR_NOR
  U[y_names %in% all$brood_year,"3","HOR"] = all$HR_HOR * 0.25
  U[y_names %in% all$brood_year,c("4", "5"),"HOR"] = all$HR_HOR
  
  ### ADULT SURVIVAL PAST SEA LIONS ###
  phi_SL = rep(NA, ny); names(phi_SL) = y_names
  phi_SL[y_names %in% sub$brood_year] = sub$surv_est_sea_lions
  
  ### COUNTS OF ADULT PIT TAG DETECTIONS AT BON ###
  x_BON = matrix(NA, ny, no); dimnames(x_BON) = list(y_names, o_names)
  x_BON[y_names %in% all$brood_year,o_names[1]] = all$NOR_BON_adults
  x_BON[y_names %in% all$brood_year,o_names[2]] = all$HOR_BON_adults
  
  ### COUNTS OF ADULT PIT TAG DETECTIONS AT LGR ###
  x_LGR = matrix(NA, ny, no); dimnames(x_LGR) = list(y_names, o_names)
  x_LGR[y_names %in% all$brood_year,o_names[1]] = all$NOR_LGR_adults
  x_LGR[y_names %in% all$brood_year,o_names[2]] = all$HOR_LGR_adults
  
  ### ADULT ABUNDANCE ###

  # total returning adults, nat + hat
  Ra_obs = rep(NA, ny); names(Ra_obs) = y_names
  Ra_obs[y_names %in% sub$brood_year] = sub$adults_at_weir
  
  # sd of total returning adults nat + hat
  sig_Ra_obs = rep(NA, ny); names(sig_Ra_obs) = y_names
  sig_Ra_obs[y_names %in% sub$brood_year] = sub$adult_log_se
  
  ### ADULT PRESPAWN DATA ###
  
  # number of carcasses sampled and found to have spawned successfully
  x_carcass_spawned = x_carcass_total = rep(NA, ny); names(x_carcass_spawned) = names(x_carcass_total) = y_names
  x_carcass_spawned[y_names %in% sub$brood_year] = sub$carcs_status_spawned
  x_carcass_total[y_names %in% sub$brood_year] = sub$carcs_samp_for_status
  
  ### INFORMATION ABOUT WHICH YEARS NEED STRAYING ###
  # the years in which strays will be needed
  if (pop != "MIN") {
    yrs = as.numeric(names(Mb_obs[,2,2]))
    first_brood_release = min(yrs[Mb_obs[,2,2] > 0], na.rm = T)
    first_adult_return = first_brood_release + kmax
    stray_yrs = 2:(which(yrs == (first_adult_return - 1)))
    not_stray_yrs = max(stray_yrs+1):ny
  } else {
    stray_yrs = 2:ny
    not_stray_yrs = numeric(0)
  }
  
  ### BUILD LIST TO RETURN ###
  
  out = list(
    
    ### DIMENSIONAL VARIABLES ###
    ny = ny,        # number of tracked brood years
    ny_obs = ny,    # number of tracked brood years
    nk = nk,        # number of ages of return
    nko = nko,    # number of age/origin classes of return
    ni = ni,        # number of life history strategies
    no = no,        # number of origin types
    kmin = kmin,    # minimum age of return
    kmax = kmax,    # maximum age of return
    
    ### JUVENILE ABUNDANCE ###
    # fall trap count
    Pa_obs = Pa_obs,
    sig_Pa_obs = sig_Pa_obs,
    
    # spring trap count
    Mb_obs = Mb_obs,
    sig_Mb_obs = sig_Mb_obs,
    
    ### JUVENILE LENGTH ###
    # mean length at summer tagging
    L_Pb_obs = L_Pb_obs,
    sig_L_Pb_obs = sig_L_Pb_obs,
    
    # mean length at spring tagging
    L_Mb_obs = L_Mb_obs,
    sig_L_Mb_obs = sig_L_Mb_obs,
    
    ### JUVENILE SURVIVAL ###
    # summer tagging to LGD
    Lphi_obs_Pb_Ma = Lphi_obs_Pb_Ma,
    sig_Lphi_obs_Pb_Ma = sig_Lphi_obs_Pb_Ma,
    
    # fall trap tagging to LGD
    Lphi_obs_Pa_Ma = Lphi_obs_Pa_Ma,
    sig_Lphi_obs_Pa_Ma = sig_Lphi_obs_Pa_Ma,
    
    # spring trap tagging to LGD
    # also, includes hatchery survival from release in spring to LGD
    Lphi_obs_Mb_Ma = Lphi_obs_Mb_Ma,
    sig_Lphi_obs_Mb_Ma = sig_Lphi_obs_Mb_Ma,
    
    # LGD to BON
    # includes both hatchery and natural origin
    Lphi_obs_Ma_O0 = Lphi_obs_Ma_O0,
    sig_Lphi_obs_Ma_O0 = sig_Lphi_obs_Ma_O0,
    
    ### ADULT SURVIVAL ###
    # adult survival past sea lions (not fitted; assumed known w/o error)
    phi_SL = phi_SL,
    
    # adult harvest rate below BON
    U = U,
    
    # counts of PIT tag detections at BON by origin
    x_BON = x_BON,
    
    # counts of PIT tag detections at LGR by origin
    x_LGR = x_LGR,
    
    ### ADULT ABUNDANCE ###
    # total adults arriving at "weir"
    Ra_obs = Ra_obs,
    sig_Ra_obs = sig_Ra_obs,
    
    # number removed at weir
    B = B,
    
    ### ADULT COMPOSITION ###
    # observed frequency of age/origin arriving at weir
    x_Ra = x_Ra,
    nx_Ra = nx_Ra,   # multinomial sample size
    
    # observed frequency of age/origin sampled as carcasses
    x_Sa_prime = x_Sa_prime,
    nx_Sa_prime = nx_Sa_prime,   # multinomial sample size
    
    # number of carcasses sampled for spawn status
    x_carcass_total = x_carcass_total,
    
    # number of carcasses found to have spawned successfully
    x_carcass_spawned = x_carcass_spawned,
    
    # weighted usable length
    wul = WUL[pop],
    
    # info about which years need straying
    stray_yrs = stray_yrs,
    not_stray_yrs = not_stray_yrs,
    n_stray_yrs = length(stray_yrs),
    n_not_stray_yrs = length(not_stray_yrs)
  )
  
  # return output
  return(out)
}

##### CREATE DATA LIST FOR JAGS: MULTIPLE POPULATIONS #####

# formats the appropriate information from the bio_dat object
# to be used by a JAGS model, multiple populations
# calls create_jags_data_one() on each supplied population
# and places the data in the right dimensions of each list element

# pops: a vector containing at least two of "CAT", "LOS", "MIN", "UGR"
# first_y: the first return year to model
# last_y: the last return year to model

create_jags_data_mult = function(pops, first_y = 1991, last_y = 2019) {
  # create the jags_data list for each population, store them as elements of a larger list
  main_list = lapply(as.list(pops), create_jags_data_one, first_y = first_y, last_y = last_y)
  
  # assign this list names for the population dimension
  names(main_list) = pops
  
  # extract the dimension variables from one of the populations
  dims_list = main_list[[1]][c("ny", "ny_obs", "nk", "nko", "ni", "no", "kmin", "kmax")]
  
  # add on nj to dimensions: number of populations
  dims_list = append(dims_list, list(nj = length(pops)))
  
  # for each data type, loop through populations extracting that data type
  # and abind it with the same data type from other populations
  obs_list = list(
    # fall trap count
    Pa_obs = abind(lapply(main_list, function(x) x$Pa_obs), along = 3),
    sig_Pa_obs = abind(lapply(main_list, function(x) x$sig_Pa_obs), along = 3),
    
    # spring trap count
    Mb_obs = abind(lapply(main_list, function(x) x$Mb_obs), along = 4),
    sig_Mb_obs = abind(lapply(main_list, function(x) x$sig_Mb_obs), along = 4),
    
    # mean length at summer tagging
    L_Pb_obs = abind(lapply(main_list, function(x) x$L_Pb_obs), along = 2),
    sig_L_Pb_obs = abind(lapply(main_list, function(x) x$sig_L_Pb_obs), along = 2),
    
    # mean length at spring tagging
    L_Mb_obs = abind(lapply(main_list, function(x) x$L_Mb_obs), along = 2),
    sig_L_Mb_obs = abind(lapply(main_list, function(x) x$sig_L_Mb_obs), along = 2),
    
    # summer tagging to LGD survival
    Lphi_obs_Pb_Ma = abind(lapply(main_list, function(x) x$Lphi_obs_Pb_Ma), along = 2),
    sig_Lphi_obs_Pb_Ma = abind(lapply(main_list, function(x) x$sig_Lphi_obs_Pb_Ma), along = 2),
    
    # fall trap tagging to LGD survival
    Lphi_obs_Pa_Ma = abind(lapply(main_list, function(x) x$Lphi_obs_Pa_Ma), along = 3),
    sig_Lphi_obs_Pa_Ma = abind(lapply(main_list, function(x) x$sig_Lphi_obs_Pa_Ma), along = 3),
    
    # spring trap tagging to LGD survival
    Lphi_obs_Mb_Ma = abind(lapply(main_list, function(x) x$Lphi_obs_Mb_Ma), along = 4),
    sig_Lphi_obs_Mb_Ma = abind(lapply(main_list, function(x) x$sig_Lphi_obs_Mb_Ma), along = 4),
    
    # hydropower survival
    Lphi_obs_Ma_O0 = main_list[[1]]$Lphi_obs_Ma_O0,
    sig_Lphi_obs_Ma_O0 = main_list[[1]]$sig_Lphi_obs_Ma_O0,
    
    # adult survival past sea lions (not fitted; assumed known w/o error)
    phi_SL = abind(lapply(main_list, function(x) x$phi_SL), along = 2),
    
    # adult PIT tag counts at BON
    x_BON = main_list[[1]]$x_BON,
    
    # adult PIT tag counts at LGR
    x_LGR = main_list[[1]]$x_LGR,
    
    # adult returns to tributary
    Ra_obs = abind(lapply(main_list, function(x) x$Ra_obs), along = 2),
    sig_Ra_obs = abind(lapply(main_list, function(x) x$sig_Ra_obs), along = 2),
    
    # age/origin composition of returns: weir
    x_Ra = abind(lapply(main_list, function(x) x$x_Ra), along = 3),
    nx_Ra = abind(lapply(main_list, function(x) x$nx_Ra), along = 2),
    
    # age/origin composition of returns: carcass
    x_Sa_prime = abind(lapply(main_list, function(x) x$x_Sa_prime), along = 3),
    nx_Sa_prime = abind(lapply(main_list, function(x) x$nx_Sa_prime), along = 2),
    
    # weir removals
    B = abind(lapply(main_list, function(x) x$B), along = 4),
    
    # number of carcasses sampled for spawn status
    x_carcass_total = abind(lapply(main_list, function(x) x$x_carcass_total), along = 2),
    
    # number of carcasses sampled with successful spawning status
    x_carcass_spawned = abind(lapply(main_list, function(x) x$x_carcass_spawned), along = 2),
    
    # weighted usable length
    wul = abind(lapply(main_list, function(x) x$wul), along = 1)
  )
  
  # add the harvest rate below BON (common to all pops)
  obs_list = append(obs_list, list(U = main_list[[1]]$U))
  
  # list containing information about years needing straying
  stray_yrs = lapply(main_list, function(x) x$stray_yrs); nmax_stray = max(sapply(stray_yrs, length))
  stray_yrs = abind(lapply(stray_yrs, function(x) c(x, rep(NA, nmax_stray - length(x)))), along = 2)
  n_stray_yrs = apply(stray_yrs, 2, function(x) sum(!is.na(x)))
  not_stray_yrs = lapply(main_list, function(x) x$not_stray_yrs); nmax_not_stray = max(sapply(not_stray_yrs, length))
  not_stray_yrs = abind(lapply(not_stray_yrs, function(x) c(x, rep(NA, nmax_not_stray - length(x)))), along = 2)
  n_not_stray_yrs = apply(not_stray_yrs, 2, function(x) sum(!is.na(x)))
  
  stray_yrs_list = list(
    stray_yrs = stray_yrs,
    n_stray_yrs = n_stray_yrs,
    not_stray_yrs = not_stray_yrs,
    n_not_stray_yrs = n_not_stray_yrs
  )
  
  # append dimension, observation, and stray years into lists together
  out_list = append(append(dims_list, obs_list), stray_yrs_list)
  
  # return the output
  return(out_list)
}

##### APPEND JAGS_DATA LIST WITH ELEMENTS THAT DO NOT CONTAIN NA VALUES #####

# adds information to a jags_data object that specifies which
# elements of data notes do not have NA values
# will allow us to include structural NAs to keep everything square
# while only fitting to non-NA values

append_no_na_indices = function(jags_data) {
  
  # obtain the elements in each dimension that do not have NA values
  # for each data type
  fit_list = with(jags_data, {
    list(
      fit_Pa = find_no_na_indices(Pa_obs),
      fit_Mb = find_no_na_indices(sig_Mb_obs),   # use sig for Mb b/c hatchery releases don't have sig, and shouldn't be fitted to
      fit_L_Pb = find_no_na_indices(L_Pb_obs),
      fit_L_Mb = find_no_na_indices(L_Mb_obs),
      fit_Lphi_Pb_Ma = find_no_na_indices(Lphi_obs_Pb_Ma),
      fit_Lphi_Pa_Ma = find_no_na_indices(Lphi_obs_Pa_Ma),
      fit_Lphi_Mb_Ma = find_no_na_indices(Lphi_obs_Mb_Ma),
      fit_Lphi_Ma_O0 = find_no_na_indices(Lphi_obs_Ma_O0),
      fit_x_LGR = find_no_na_indices(x_LGR),
      fit_Ra = find_no_na_indices(Ra_obs),
      fit_x_carcass_spawned = find_no_na_indices(x_carcass_spawned)
    )
  })
  
  # calculate the number of elements for each data type without NA values
  nfit_list = lapply(fit_list, function(x) xdim(x)[1])
  names(nfit_list) = paste0("n", names(nfit_list))
  
  # append these with the original jags_data object
  out = append(jags_data, append(fit_list, nfit_list))
  
  # return the output list: now ready for JAGS
  return(out)
}

append_values_for_sim = function(jags_data) {
  
  # last year in the observed time period
  last_obs_yr = max(as.numeric(rownames(jags_data$Pa_obs)))
  
  # how many years in the simulated period
  ny_sim = jags_data$ny
  
  # last year including the simulated period
  last_yr = ny_sim + last_obs_yr
  
  # bump up the total number of years to include the simulated period
  jags_data$ny = jags_data$ny + ny_sim
  
  ### KNOWN VALUES SOURCE 1: HATCHERY SMOLT RELEASES ###
  # extract the values in the observed period
  x = jags_data$Mb_obs
  
  # duplicate values for brood years not in data set
  take_yrs = as.character(last_obs_yr + c(-3,-2))
  fill_yrs = as.character(last_obs_yr + c(-1, 0))
  x[fill_yrs,"spring-mig","HOR",] = x[take_yrs,"spring-mig","HOR",]
  
  # duplicate whole object and change year names
  y = x
  dimnames(y)[[1]] = 1:ny_sim + last_obs_yr
  
  # insert a value in the old "year-0" position
  y[as.character(last_obs_yr + 1),"spring-mig","HOR",] = x[jags_data$ny_obs,"spring-mig","HOR",]
  
  # append with observed period to be passed to the model
  jags_data$Mb_obs = abind(x, y, along = 1)
  
  ### KNOWN VALUES SOURCE X: WEIR REMOVALS ###
  
  # extract the values from the observed period
  x = jags_data$B
  
  # insert a value in the old "year-0" position
  x[1,,,] = jags_data$B[jags_data$ny_obs,,,]
  
  # change year names
  dimnames(x)[[1]] = 1:ny_sim + last_obs_yr
  
  # append with observed period to be passed to the model
  jags_data$B = abind(jags_data$B, x, along = 1)
  
  ### KNOWN VALUES SOURCE X: STRAY YEARS ###
  
  jags_data$stray_yrs[,c("CAT","LOS","UGR")] = NA
  fill_yrs = (max(jags_data$stray_yrs[,"MIN"], na.rm = TRUE) + 1):jags_data$ny
  fill = matrix(NA, length(fill_yrs), jags_data$nj)
  colnames(fill) = colnames(jags_data$stray_yrs)
  fill[,"MIN"] = fill_yrs
  jags_data$stray_yrs = rbind(jags_data$stray_yrs[-(1:jags_data$kmax),], fill)
  jags_data$n_stray_yrs = apply(jags_data$stray_yrs, 2, function(x) sum(!is.na(x)))
  
  fill = matrix(2:6, 5, jags_data$nj)
  colnames(fill) = colnames(jags_data$not_stray_yrs)
  jags_data$not_stray_yrs = rbind(fill, jags_data$not_stray_yrs)
  fill = matrix(fill_yrs, length(fill_yrs), jags_data$nj)
  colnames(fill) = colnames(jags_data$not_stray_yrs)
  fill[,"MIN"] = NA
  jags_data$not_stray_yrs = rbind(jags_data$not_stray_yrs, fill)
  jags_data$n_not_stray_yrs = apply(jags_data$not_stray_yrs, 2, function(x) sum(!is.na(x)))
  
  ### KNOWN VALUES SOURCE X: SURVIVAL PAST SEA LIONS ###
  
  # homogenize and extract the values from the observed period
  jags_data$phi_SL["2000",] = jags_data$phi_SL["2001",]
  x = jags_data$phi_SL
  
  # insert a value in the old "year-0" position
  x[1,] = jags_data$phi_SL["2001",]
  
  # change year names
  dimnames(x)[[1]] = 1:ny_sim + last_obs_yr
  
  # append with observed period to be passed to the model
  jags_data$phi_SL = abind(jags_data$phi_SL, x, along = 1)
  
  ### KNOWN VALUES SOURCE X: HARVEST RATES ###
  
  # extract the values from the observed period
  x = jags_data$U
  
  # insert a value in the old "year-0" position
  x[1,,] = x[2,,]
  
  # change year names
  dimnames(x)[[1]] = 1:ny_sim + last_obs_yr
  
  # append with observed period to be passed to the model
  jags_data$U = abind(jags_data$U, x, along = 1)
  
  # return the updated jags_data object
  return(jags_data)
}
