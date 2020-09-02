# :::::::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT TO HOUSE FUNCTIONS FOR PREPARING RAW DATA #
# :::::::::::::::::::::::::::::::::::::::::::::::: #

##### OBTAIN LOGIT-SCALE STANDARD ERROR OF A PROPORTION #####
# some estimates have p_mean and CIs
# some estimates have p_mean and p_se
# either way, obtain CI's of proportion, then convert to a standard error on the logit scale
get_logit_se = function(N, p_mean, p_se, p_lwr, p_upr) {
  # turn se's into ci if that is what is available
  p_lwr = ifelse(is.na(p_se), p_lwr, p_mean + qnorm(0.025) * sqrt((p_mean * (1 - p_mean))/N) - (0.5/N))
  p_upr = ifelse(is.na(p_se), p_upr, p_mean + qnorm(0.975) * sqrt((p_mean * (1 - p_mean))/N) + (0.5/N))
  
  # cap the CI
  p_lwr = ifelse(p_lwr <= 0, 0.001, p_lwr)
  p_upr = ifelse(p_upr >= 1, 0.999, p_upr)
  
  # convert CI into logit-normal standard error: M. Liermann's approximation
  (logit(p_upr) - logit(p_lwr))/(2 * qnorm(0.975))
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
  
  # specify/calculate dimensional variables
  kmin = 3              # minimum age of return
  kmax = 5              # maximum age of return
  nk = kmax - kmin + 1  # number of ages of return
  ns = 2                # number of sexes of return
  ni = 2                # number of juvenile life history strategies
  no = 2                # number of origins
  nt = nrow(sub)        # number of return years tracked
  ny = nt + kmax        # number of brood years tracked
  
  # assign these names
  y_names = (first_y - kmax):last_y
  k_names = kmin:kmax
  i_names = paste0(c("fall", "spring"), "-mig")
  s_names = c("F", "M")
  o_names = c("Nat", "Hat")
  ks_names = c(paste0(s_names[1], k_names), paste0(s_names[2], k_names))
  
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
  Lphi_obs_Mb_Ma = matrix(NA, ny, ni); dimnames(Lphi_obs_Mb_Ma) = list(y_names, i_names)
  Lphi_obs_Mb_Ma[y_names %in% sub$brood_year,i_names == "spring-mig"] = logit(sub$spring_surv_est)
  
  # sd spring logit(surv) to LGD
  sig_Lphi_obs_Mb_Ma = matrix(NA, ny, ni); dimnames(sig_Lphi_obs_Mb_Ma) = list(y_names, i_names)
  sig_Lphi_obs_Mb_Ma[y_names %in% sub$brood_year,i_names == "spring-mig"] = sub$spring_surv_logit_se
  
  ### ADULT AGE COMP: WEIR ###
  # obtain names of age comp variables
  x = expand.grid(o = o_names, s = s_names, k = k_names)
  x = x[order(x$o, x$s),]
  x = apply(x, 1, function(x) paste(c("weir", x), collapse = "_"))
  nat_names = x[1:(nk * ns)]
  hat_names = x[((nk * ns) + 1):(no * nk * ns)]
  
  # extract them by origin and coerce NA to zero
  nat_comp = sub[,nat_names]
  hat_comp = sub[,hat_names]
  nat_comp[is.na(nat_comp)] = 0
  hat_comp[is.na(hat_comp)] = 0
  
  # calculate percent of hatchery origin return
  # used to expand natural origin adult returns to obtain
  # total returns nat + hat. Could be done on an age-specific basis
  # too, but for now let's start here.
  p_HOR = rep(NA, ny); names(p_HOR) = y_names
  p_HOR[y_names %in% sub$brood_year] = rowSums(hat_comp)/(rowSums(nat_comp + hat_comp))
  p_HOR[is.na(p_HOR)] = 0
  
  # age frequencies: hat + nat; for fitting composition of returns
  x_obs = matrix(NA, ny, ns * nk); dimnames(x_obs) = list(y_names, ks_names)
  x_obs[y_names %in% sub$brood_year,] = as.matrix(hat_comp + nat_comp)
  nx_obs = rowSums(x_obs)
  
  ### ADULT ABUNDANCE ###
  # proportion of total adults arriving at weir that are removed for broodstock
  p_remove = rep(NA, ny); names(p_remove) = y_names
  p_remove[y_names %in% sub$brood_year] = (sub$adults_at_weir - sub$adults_above_weir)/sub$adults_at_weir
  
  # total returning adults, nat + hat
  Ra_obs = rep(NA, ny); names(Ra_obs) = y_names
  Ra_obs[y_names %in% sub$brood_year] = sub$adults_at_weir
  
  # sd of total returning adults nat + hat
  sig_Ra_obs = rep(NA, ny); names(sig_Ra_obs) = y_names
  sig_Ra_obs[y_names %in% sub$brood_year] = sub$adult_log_se
  
  ### BUILD LIST TO RETURN ###
  
  out = list(
    
    ### DIMENSIONAL VARIABLES ###
    ny = ny,        # number of tracked brood years
    nt = nt,        # number of tracked return years
    nk = nk,        # number of ages of return
    ns = ns,        # number of sexes of return
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
    Lphi_obs_Mb_Ma = Lphi_obs_Mb_Ma,
    sig_Lphi_obs_Mb_Ma = sig_Lphi_obs_Mb_Ma,
    
    ### ADULT ABUNDANCE ###
    # total adults arriving at "weir"
    Ra_obs = Ra_obs,
    sig_Ra_obs = sig_Ra_obs,
    
    # proportion removed for broodstock
    p_remove = p_remove,
    
    ### ADULT COMPOSITION ###
    # observed frequency of age/sex arriving at "weir"
    x_obs = x_obs,
    nx_obs = nx_obs,   # multinomial sample size
    
    # proportion of hatchery origin returns
    p_HOR = p_HOR
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
  dims_list = main_list[[1]][c("ny", "nt", "nk", "ns", "ni", "no", "kmin", "kmax")]
  
  # add on nj to dimensions: number of populations
  dims_list = append(dims_list, list(nj = length(pops)))
  
  # for each data type, loop through populations extracting that data type
  # and abind it with the same data type from other populations
  obs_list = list(
    # fall trap count
    Pa_obs = abind(lapply(main_list, function(x) x$Pa_obs), along = 3),
    sig_Pa_obs = abind(lapply(main_list, function(x) x$sig_Pa_obs), along = 3),
    
    # spring trap count
    Mb_obs = abind(lapply(main_list, function(x) x$Mb_obs), along = 3),
    sig_Mb_obs = abind(lapply(main_list, function(x) x$sig_Mb_obs), along = 3),
    
    # summer tagging to LGD survival
    Lphi_obs_Pb_Ma = abind(lapply(main_list, function(x) x$Lphi_obs_Pb_Ma), along = 2),
    sig_Lphi_obs_Pb_Ma = abind(lapply(main_list, function(x) x$sig_Lphi_obs_Pb_Ma), along = 2),
    
    # fall trap tagging to LGD survival
    Lphi_obs_Pa_Ma = abind(lapply(main_list, function(x) x$Lphi_obs_Pa_Ma), along = 3),
    sig_Lphi_obs_Pa_Ma = abind(lapply(main_list, function(x) x$sig_Lphi_obs_Pa_Ma), along = 3),
    
    # spring trap tagging to LGD survival
    Lphi_obs_Mb_Ma = abind(lapply(main_list, function(x) x$Lphi_obs_Mb_Ma), along = 3),
    sig_Lphi_obs_Mb_Ma = abind(lapply(main_list, function(x) x$sig_Lphi_obs_Mb_Ma), along = 3),
    
    # adult returns to tributary
    Ra_obs = abind(lapply(main_list, function(x) x$Ra_obs), along = 2),
    sig_Ra_obs = abind(lapply(main_list, function(x) x$sig_Ra_obs), along = 2),
    
    # age/sex composition of returns to tributary
    x_obs = abind(lapply(main_list, function(x) x$x_obs), along = 3),
    nx_obs = abind(lapply(main_list, function(x) x$nx_obs), along = 2),
    
    # broodstock removals
    p_remove = abind(lapply(main_list, function(x) x$p_remove), along = 2),
    
    # proportion of hatchery origin returns
    p_HOR = abind(lapply(main_list, function(x) x$p_HOR), along = 2)
  )
  
  # append dimension and observation lists together
  out_list = append(dims_list, obs_list)
  
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
      fit_Mb = find_no_na_indices(Mb_obs),
      fit_Lphi_Pb_Ma = find_no_na_indices(Lphi_obs_Pb_Ma),
      fit_Lphi_Pa_Ma = find_no_na_indices(Lphi_obs_Pa_Ma),
      fit_Lphi_Mb_Ma = find_no_na_indices(Lphi_obs_Mb_Ma),
      fit_Ra = find_no_na_indices(Ra_obs)
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

