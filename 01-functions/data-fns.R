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


##### CREATE NAMES FOR A SPECIFIC COMPOSITION DATA SET #####
# o_names = c("Nat", "Hat")
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
# first_y: the first return year to model
# last_y: the last return year to model

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
  ns = 2                # number of sexes of return
  ni = 2                # number of juvenile life history strategies
  no = 2                # number of origins
  nks = nk * ns         # number of age/sex classes
  nkso = nk * ns * no   # number of age/sex/origin classes
  nt = nrow(sub)        # number of return years tracked
  ny = nt + kmax        # number of brood years tracked
  
  # assign these names
  y_names = (first_y - kmax):last_y
  k_names = kmin:kmax
  i_names = paste0(c("fall", "spring"), "-mig")
  s_names = c("F", "M")
  o_names = c("Nat", "Hat")
  ks_names = c(paste0(s_names[1], k_names), paste0(s_names[2], k_names))
  kso_names = c(paste0(ks_names, "-", o_names[1]), paste0(ks_names, "-", o_names[2]))
  
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
  
  ### JUVENILE SURVIVAL DATA ###
  # summer logit(surv) to LGD
  Lphi_obs_Pb_Ma = rep(NA, ny); names(Lphi_obs_Pb_Ma) = y_names
  Lphi_obs_Pb_Ma[y_names %in% sub$brood_year] = logit(sub$summer_surv_est)
  
  # sd summer logit(surv) to LGD
  sig_Lphi_obs_Pb_Ma = rep(NA, ny); names(sig_Lphi_obs_Pb_Ma) = y_names
  sig_Lphi_obs_Pb_Ma[y_names %in% sub$brood_year] = sub$summer_surv_logit_se
  
  # fall logit(surv) to LGD
  Lphi_obs_Pa_Ma = matrix(NA, ny, ni); dimnames(Lphi_obs_Pa_Ma) = list(y_names, i_names)
  Lphi_obs_Pa_Ma[y_names %in% sub$brood_year,i_names == "fall-mig"] = logit(sub$fall_surv_est)
  
  # sd fall logit(surv) to LGD
  sig_Lphi_obs_Pa_Ma = matrix(NA, ny, ni); dimnames(sig_Lphi_obs_Pa_Ma) = list(y_names, i_names)
  sig_Lphi_obs_Pa_Ma[y_names %in% sub$brood_year,i_names == "fall-mig"] = sub$fall_surv_logit_se
  
  # spring logit(surv) to LGD
  Lphi_obs_Mb_Ma = array(NA, dim = c(ny, ni, no)); dimnames(Lphi_obs_Mb_Ma) = list(y_names, i_names, o_names)
  Lphi_obs_Mb_Ma[y_names %in% sub$brood_year,i_names == "spring-mig","Nat"] = logit(sub$spring_surv_est)
  Lphi_obs_Mb_Ma[y_names %in% sub$brood_year,i_names == "spring-mig","Hat"] = logit(sub$hatchery_spring_surv_est)
  
  # sd spring logit(surv) to LGD
  sig_Lphi_obs_Mb_Ma = array(NA, dim = c(ny, ni, no)); dimnames(sig_Lphi_obs_Mb_Ma) = list(y_names, i_names, o_names)
  sig_Lphi_obs_Mb_Ma[y_names %in% sub$brood_year,i_names == "spring-mig","Nat"] = sub$spring_surv_logit_se
  sig_Lphi_obs_Mb_Ma[y_names %in% sub$brood_year,i_names == "spring-mig","Hat"] = sub$hatchery_spring_surv_logit_se
  
  # LGD to BON logit(surv)
  Lphi_obs_Ma_O0 = matrix(NA, ny, no); dimnames(Lphi_obs_Ma_O0) = list(y_names, o_names)
  Lphi_obs_Ma_O0[y_names %in% all$brood_year,"Nat"] = logit(all$nat_hydro_est)
  Lphi_obs_Ma_O0[y_names %in% all$brood_year,"Hat"] = logit(all$hat_hydro_est)
  
  # sd LGD to BON logit(surv)
  sig_Lphi_obs_Ma_O0 = matrix(NA, ny, no); dimnames(sig_Lphi_obs_Ma_O0) = list(y_names, o_names)
  sig_Lphi_obs_Ma_O0[y_names %in% all$brood_year,"Nat"] = all$nat_hydro_logit_se
  sig_Lphi_obs_Ma_O0[y_names %in% all$brood_year,"Hat"] = all$hat_hydro_logit_se
  
  ### ADULT AGE COMP: WEIR ###
  # obtain names of age comp variables
  weir_comp_names = create_comp_names("weir", o_names, s_names, k_names)
  
  # extract them by origin and coerce NA to zero
  nat_comp = sub[,weir_comp_names$nat_names]
  hat_comp = sub[,weir_comp_names$hat_names]
  nat_comp[is.na(nat_comp)] = 0
  hat_comp[is.na(hat_comp)] = 0
  
  # age frequencies: for fitting composition of returns
  weir_x_obs = matrix(NA, ny, no * ns * nk); dimnames(weir_x_obs) = list(y_names, kso_names)
  weir_x_obs[y_names %in% sub$brood_year,] = as.matrix(cbind(nat_comp, hat_comp))
  weir_nx_obs = rowSums(weir_x_obs)
  
  ### ADULT AGE COMP: CARCASSES ###
  # obtain names of age comp variables
  carc_comp_names = create_comp_names("carc", o_names, s_names, k_names)
  
  # extract them by origin and coerce NA to zero
  nat_comp = sub[,carc_comp_names$nat_names]
  hat_comp = sub[,carc_comp_names$hat_names]
  nat_comp[is.na(nat_comp)] = 0
  hat_comp[is.na(hat_comp)] = 0
  
  # age frequencies: for fitting composition of returns
  carc_x_obs = matrix(NA, ny, no * ns * nk); dimnames(carc_x_obs) = list(y_names, kso_names)
  carc_x_obs[y_names %in% sub$brood_year,] = as.matrix(cbind(nat_comp, hat_comp))
  carc_nx_obs = rowSums(carc_x_obs)
  
  ### PROPORTION OF RETURNING ADULTS REMOVED AT WEIR ###
  # weir_comp_names = create_comp_names("weir", o_names, s_names, k_names)
  rm_comp_names = create_comp_names("rm", o_names, s_names, k_names)
  
  # extract compositions by type and coerce NAs to zero
  # this is the number of fish removed at weir each year by age/sex/origin class
  rm_comp = sub[,c(rm_comp_names$nat_names, rm_comp_names$hat_names)]
  rm_comp[is.na(rm_comp)] = 0
  
  # place rm_comp in the correct location of n_remove: same numbers just reformatted array structure used by model
  n_remove = array(NA, dim = c(ny, nk, ns, no)); dimnames(n_remove) = list(y_names, k_names, s_names, o_names)
  n_remove[y_names %in% sub$brood_year,,s_names[1],o_names[1]] = as.matrix(rm_comp[,paste("rm", o_names[1], s_names[1], k_names, sep = "_")])
  n_remove[y_names %in% sub$brood_year,,s_names[2],o_names[1]] = as.matrix(rm_comp[,paste("rm", o_names[1], s_names[2], k_names, sep = "_")])
  n_remove[y_names %in% sub$brood_year,,s_names[1],o_names[2]] = as.matrix(rm_comp[,paste("rm", o_names[2], s_names[1], k_names, sep = "_")])
  n_remove[y_names %in% sub$brood_year,,s_names[2],o_names[2]] = as.matrix(rm_comp[,paste("rm", o_names[2], s_names[2], k_names, sep = "_")])
  
  ### ADULT SURVIVAL PAST SEA LIONS ###
  phi_SL = rep(NA, ny); names(phi_SL) = y_names
  phi_SL[y_names %in% sub$brood_year] = sub$surv_est_sea_lions
  
  ### COUNTS OF ADULT PIT TAG DETECTIONS AT BON ###
  BON_adults = matrix(NA, ny, no); dimnames(BON_adults) = list(y_names, o_names)
  BON_adults[y_names %in% all$brood_year,o_names[1]] = all$NOR_BON_adults
  BON_adults[y_names %in% all$brood_year,o_names[2]] = all$HOR_BON_adults
  
  ### COUNTS OF ADULT PIT TAG DETECTIONS AT LGR ###
  LGR_adults = matrix(NA, ny, no); dimnames(LGR_adults) = list(y_names, o_names)
  LGR_adults[y_names %in% all$brood_year,o_names[1]] = all$NOR_LGR_adults
  LGR_adults[y_names %in% all$brood_year,o_names[2]] = all$HOR_LGR_adults
  
  ### ADULT ABUNDANCE ###

  # total returning adults, nat + hat
  Ra_obs = rep(NA, ny); names(Ra_obs) = y_names
  Ra_obs[y_names %in% sub$brood_year] = sub$adults_at_weir
  
  # sd of total returning adults nat + hat
  sig_Ra_obs = rep(NA, ny); names(sig_Ra_obs) = y_names
  sig_Ra_obs[y_names %in% sub$brood_year] = sub$adult_log_se
  
  ### ADULT PRESPAWN DATA ###
  
  # number of carcasses sampled and found to have spawned successfully
  carcs_spawned = carcs_sampled = rep(NA, ny); names(carcs_spawned) = names(carcs_sampled) = y_names
  carcs_spawned[y_names %in% sub$brood_year] = sub$carcs_status_spawned
  carcs_sampled[y_names %in% sub$brood_year] = sub$carcs_samp_for_status
  
  ### INFORMATION ABOUT WHICH YEARS NEED STRAYING ###
  # the years in which strays will be needed
  if (pop != "MIN") {
    yrs = as.numeric(names(Mb_obs[,2,2]))
    first_brood_release = min(yrs[Mb_obs[,2,2] > 0], na.rm = T)
    first_adult_return = first_brood_release + kmax
    stray_yrs = (kmax+1):(which(yrs == (first_adult_return - 1)))
    not_stray_yrs = max(stray_yrs+1):ny
  } else {
    stray_yrs = (kmax+1):ny
    not_stray_yrs = numeric(0)
  }
  
  ### BUILD LIST TO RETURN ###
  
  out = list(
    
    ### DIMENSIONAL VARIABLES ###
    ny = ny,        # number of tracked brood years
    nk = nk,        # number of ages of return
    ns = ns,        # number of sexes of return
    nks = nks,      # number of age/sex classes of return
    nkso = nkso,    # number of age/sex/origin classes of return
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
    
    # counts of PIT tag detections at BON by origin
    BON_adults = BON_adults,
    
    # counts of PIT tag detections at LGR by origin
    LGR_adults = LGR_adults,
    
    ### ADULT ABUNDANCE ###
    # total adults arriving at "weir"
    Ra_obs = Ra_obs,
    sig_Ra_obs = sig_Ra_obs,
    
    # number removed at weir
    n_remove = n_remove,
    
    ### ADULT COMPOSITION ###
    # observed frequency of age/sex/origin arriving at weir
    weir_x_obs = weir_x_obs,
    weir_nx_obs = weir_nx_obs,   # multinomial sample size
    
    # observed frequency of age/sex/origin sampled as carcasses
    carc_x_obs = carc_x_obs,
    carc_nx_obs = carc_nx_obs,   # multinomial sample size
    
    # number of carcasses sampled for spawn status
    carcs_sampled = carcs_sampled,
    
    # number of carcasses found to have spawned successfully
    carcs_spawned = carcs_spawned,
    
    # pool equivalent units
    peu = PEU[pop],
    
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
  dims_list = main_list[[1]][c("ny", "nk", "ns", "nks", "nkso", "ni", "no", "kmin", "kmax")]
  
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
    BON_adults = main_list[[1]]$BON_adults,
    
    # adult PIT tag counts at LGR
    LGR_adults = main_list[[1]]$LGR_adults,
    
    # adult returns to tributary
    Ra_obs = abind(lapply(main_list, function(x) x$Ra_obs), along = 2),
    sig_Ra_obs = abind(lapply(main_list, function(x) x$sig_Ra_obs), along = 2),
    
    # age/sex composition of returns: weir
    weir_x_obs = abind(lapply(main_list, function(x) x$weir_x_obs), along = 3),
    weir_nx_obs = abind(lapply(main_list, function(x) x$weir_nx_obs), along = 2),
    
    # age/sex composition of returns: carcass
    carc_x_obs = abind(lapply(main_list, function(x) x$carc_x_obs), along = 3),
    carc_nx_obs = abind(lapply(main_list, function(x) x$carc_nx_obs), along = 2),
    
    # weir removals removals
    n_remove = abind(lapply(main_list, function(x) x$n_remove), along = 5),
    
    # number of carcasses sampled for spawn status
    carcs_sampled = abind(lapply(main_list, function(x) x$carcs_sampled), along = 2),
    
    # number of carcasses sampled with successful spawning status
    carcs_spawned = abind(lapply(main_list, function(x) x$carcs_spawned), along = 2),
    
    # pool equivalent units
    peu = abind(lapply(main_list, function(x) x$peu), along = 1)
  )
  
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
      fit_Lphi_Pb_Ma = find_no_na_indices(Lphi_obs_Pb_Ma),
      fit_Lphi_Pa_Ma = find_no_na_indices(Lphi_obs_Pa_Ma),
      fit_Lphi_Mb_Ma = find_no_na_indices(Lphi_obs_Mb_Ma),
      fit_Lphi_Ma_O0 = find_no_na_indices(Lphi_obs_Ma_O0),
      fit_LGR_adults = find_no_na_indices(LGR_adults),
      fit_Ra = find_no_na_indices(Ra_obs),
      fit_spawned = find_no_na_indices(carcs_spawned)
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

