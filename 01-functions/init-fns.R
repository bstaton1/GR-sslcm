# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT TO HOUSE FUNCTIONS FOR GENERATING INITIAL MCMC VALUES #
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #


# reconstruct total parr
get_Pb_obs = function(jags_data, fill_missing = FALSE, append_sim_yrs = FALSE) {
  
  # calculation:
  # (fall_trap_count * fall_toLGR_surv + spring_trap_count * spring_to_LGR_surv)/summer_toLGR_surv
  
  Pb = with(jags_data, {
    
    # extract survival from summer tagging to LGR
    summer_toLGR_surv = expit(Lphi_obs_Pb_Ma)
    
    # extract survival from fall tagging to LGR and count at trap
    fall_toLGR_surv = expit(Lphi_obs_Pa_Ma[,i_fall,])
    fall_at_trap = Pa_obs[,i_fall,]
    
    # extract survival from spring tagging to LGR and count at trap
    spring_toLGR_surv = expit(Lphi_obs_Mb_Ma[,i_spring,o_nor,])
    spring_at_trap = Mb_obs[1:ny_obs,i_spring,o_nor,]
    
    # calculate numbers reaching LGR
    fall_at_LGR = fall_at_trap * fall_toLGR_surv
    spring_at_LGR = spring_at_trap * spring_toLGR_surv
    tot_at_LGR = fall_at_LGR + spring_at_LGR
    
    # calculate initial parr abundance at summer tagging
    tot_summer = tot_at_LGR/summer_toLGR_surv
    
    # append years for simulation if requested and needed
    if (append_sim_yrs & (ny - ny_obs > 1)) {
      
      ny_sim = ny - ny_obs
      obs_yrs = as.numeric(rownames(Lphi_obs_Pb_Ma)[ny_obs])
      
      empty = matrix(NA, ny_sim, nj)
      rownames(empty) = obs_yrs + seq(1, ny_sim)
      colnames(empty) = colnames(tot_summer)
      tot_summer = rbind(tot_summer, empty)
    }
    
    tot_summer
  })
  
  # fill in missing values with the mean of all other years for that population
  if (fill_missing) {
    Pb_filled = apply(Pb, 2, function(x) {x[is.na(x)] = mean(x, na.rm = TRUE); x})
    dimnames(Pb_filled) = dimnames(Pb)
    Pb = Pb_filled
  }
  
  # first year should always be NA
  Pb[1,] = NA
  
  return(Pb)
}

# reconstruct total egg production
get_E_obs = function(jags_data, fill_missing = FALSE, append_sim_yrs = FALSE) {
  E = with(jags_data, {
    # get the total adult return
    total_return = Ra_obs
    
    # get the broodstock removed
    broodstock = B
    broodstock = lapply(1:nj, function(j) cbind(broodstock[1:ny_obs,,1,j], broodstock[1:ny_obs,,1,j]))
    broodstock = do.call(abind, append(broodstock, list(along = 3)))
    
    # get the pre-spawn survival
    prespawn_surv = x_carcass_spawned/x_carcass_total
    prespawn_surv[is.na(prespawn_surv)] = 0.75
    prespawn_surv[prespawn_surv < 0.5] = 0.5
    
    # get the age/origin compositon of the total adult return
    # assume 60% pHOS and 20%, 60%, 20% of age 3, 4, and 5
    p_hor_assume = 0.6
    p_age_assume = c(0.2, 0.6, 0.2)
    comp = c(p_age_assume * (1 - p_hor_assume), p_age_assume * p_hor_assume)
    
    # get the age/origin structured return
    age_origin_return = lapply(1:nj, function(j) {
      t(sapply(1:length(total_return[,j]), function(y) comp * total_return[y,j]))
    })
    age_origin_return = do.call(abind, append(age_origin_return, list(along = 3)))
    
    # remove brood_stock
    age_origin_above_weir = age_origin_return - broodstock
    age_origin_above_weir[age_origin_above_weir < 0] = 5
    
    # remove prespawn morts
    age_origin_survived = lapply(1:nj, function(j) {
      apply(age_origin_above_weir[,,j], 2, function(a) a * prespawn_surv[,j])
    })
    age_origin_survived = do.call(abind, append(age_origin_survived, list(along = 3)))
    
    # get number of females by age/origin
    p_female_age_origin = rbind(Omega, Omega)
    females_age_origin = lapply(1:nj, function(j) {
      t(apply(age_origin_survived[,,j], 1, function(y) y * p_female_age_origin[,j]))
    })
    females_age_origin = do.call(abind, append(females_age_origin, list(along = 3)))
    
    # get total egg output
    f_age_origin = c(f, f)
    eggs_age_origin = lapply(1:nj, function(j) {
      t(apply(females_age_origin[,,j], 1, function(y) y * f_age_origin))
    })
    eggs_age_origin = do.call(abind, append(eggs_age_origin, list(along = 3)))
    total_eggs = apply(eggs_age_origin, 3, rowSums)
    
    # add population names back in
    colnames(total_eggs) = colnames(total_return)
    
    # append years for simulation if requested and needed
    if (append_sim_yrs & (ny - ny_obs > 1)) {
      
      ny_sim = ny - ny_obs
      obs_yrs = as.numeric(rownames(Lphi_obs_Pb_Ma)[ny_obs])
      
      empty = matrix(NA, ny_sim, nj)
      rownames(empty) = obs_yrs + seq(1, ny_sim)
      colnames(empty) = colnames(total_eggs)
      total_eggs = rbind(total_eggs, empty)
    }
    
    total_eggs
  })
  
  if (fill_missing) {
    E_filled = apply(E, 2, function(x) {x[is.na(x)] = mean(x, na.rm = TRUE); x})
    dimnames(E_filled) = dimnames(E)
    E = E_filled
  }
  
  return(E)
}

# get egg to parr survival
get_phi_E_Pb_obs = function(jags_data, fill_missing = FALSE, append_sim_yrs = FALSE) {
  Pb = get_Pb_obs(jags_data, fill_missing = fill_missing, append_sim_yrs = append_sim_yrs)
  E = get_E_obs(jags_data, fill_missing = fill_missing, append_sim_yrs = append_sim_yrs)
  
  phi_E_Pb = Pb/E
  phi_E_Pb[phi_E_Pb >= 1] = 0.1
  
  phi_E_Pb
}

# get proportion of parr that become fall migrants (pi)
get_pi_obs = function(jags_data, fill_missing = FALSE, append_sim_yrs = FALSE) {
  
  pi = with(jags_data, {
    
    # reconstruct parr recruitment
    Pb = get_Pb_obs(jags_data, fill_missing = FALSE, append_sim_yrs = FALSE)
    
    # extract count of fall parr at trap
    Pa_fall = Pa_obs[,i_fall,]
    
    # divide to get fraction of all parr that become fall migrants
    pi = Pa_fall/Pb
    
    # append years for simulation if requested and needed
    if (append_sim_yrs & (ny - ny_obs > 1)) {
      
      ny_sim = ny - ny_obs
      obs_yrs = as.numeric(rownames(Lphi_obs_Pb_Ma)[ny_obs])
      
      empty = matrix(NA, ny_sim, nj)
      rownames(empty) = obs_yrs + seq(1, ny_sim)
      colnames(empty) = colnames(pi)
      pi = rbind(pi, empty)
    }
    pi
  })
  
  if (fill_missing) {
    pi_filled = apply(pi, 2, function(x) {x[is.na(x)] = mean(x, na.rm = TRUE); x})
    dimnames(pi_filled) = dimnames(pi)
    pi = pi_filled
  }
  
  return(pi)
}

# get overwinter survival
get_phi_Pa_Mb_obs = function(jags_data, fill_missing = FALSE, append_sim_yrs = FALSE) {
  
  phi_Pa_Mb = with(jags_data, {
    
    # extract spring migrant smolt
    Mb_spring = Mb_obs[1:ny_obs,i_spring,o_nor,]
    
    # extract smolt survival to LGR (assumed equal among LH types)
    phi_Mb_Ma = plogis(Lphi_obs_Mb_Ma[,i_spring,o_nor,])
    
    # reconstruct total parr
    Pb = get_Pb_obs(jags_data, fill_missing = fill_missing, append_sim_yrs = append_sim_yrs)
    
    # reconstruct proportion of total parr that migrate out as fall migrants
    pi = get_pi_obs(jags_data, fill_missing = fill_missing, append_sim_yrs = append_sim_yrs)
    
    # get parr that are fall migrants
    Pa_fall = Pb * pi
    
    # get parr that are spring migrants
    Pa_spring = Pb * (1 - pi)
    
    # get overwinter survival for spring migrants
    phi_Pa_Mb_spring = Mb_spring/Pa_spring[1:ny_obs,]
    
    # get fall migrant smolt at LGR
    Ma_fall = Pa_obs[,i_fall,] * plogis(Lphi_obs_Pa_Ma[,i_fall,])
    
    # get fall migrant smolt in trib
    Mb_fall = Ma_fall/phi_Mb_Ma
    
    # get overwinter survival for fall migrants
    phi_Pa_Mb_fall = Mb_fall/Pa_obs[,i_fall,]
    
    # combine into an array
    phi_Pa_Mb = abind(phi_Pa_Mb_fall, phi_Pa_Mb_spring, along = 3)
    dimnames(phi_Pa_Mb)[[3]] = c("fall-mig", "spring-mig")
    
    # append years for simulation if requested and needed
    if (append_sim_yrs & (ny - ny_obs > 1)) {
      
      ny_sim = ny - ny_obs
      obs_yrs = as.numeric(rownames(Lphi_obs_Pb_Ma)[ny_obs])
      
      empty = matrix(NA, ny_sim, nj)
      rownames(empty) = obs_yrs + seq(1, ny_sim)
      colnames(empty) = dimnames(phi_Pa_Mb)[[2]]
      empty = abind(empty, empty, along = 3)
      dimnames(empty)[[3]] = c("fall-mig", "spring-mig")
      phi_Pa_Mb = abind(phi_Pa_Mb, empty, along = 1)
    }
    phi_Pa_Mb
  })
  
  if (fill_missing) {
    phi_Pa_Mb_filled = lapply(1:jags_data$ni,  function(i) {
      apply(phi_Pa_Mb[,,i], 2, function(x) {x[is.na(x)] = mean(x, na.rm = TRUE); x})
    })
    phi_Pa_Mb = do.call(abind, append(phi_Pa_Mb_filled, list(along = 3)))
    dimnames(phi_Pa_Mb)[[3]] = c("fall-mig", "spring-mig")
  }
  
  # reorder dimensions so it is [year,LH_type,pop]
  phi_Pa_Mb = abind(phi_Pa_Mb[,1,], phi_Pa_Mb[,2,], phi_Pa_Mb[,3,], phi_Pa_Mb[,4,], along = 3)
  dimnames(phi_Pa_Mb)[[3]] = colnames(jags_data$Ra_obs)
  
  return(phi_Pa_Mb)
}

# get BH params
get_BH_params = function(jags_data) {
  
  # loop through populations and obtain estimates of alpha (max parr/egg), beta (max parr), and sigma (SD of logit(parr/egg))
  BH_params = t(sapply(1:jags_data$nj, function(j) {
    # reconstruct parr and eggs for this population
    parr = get_Pb_obs(jags_data)[,j]
    eggs = get_E_obs(jags_data)[,j]
    
    # get logit(phi_E_Pb) for this population
    logit_parr_per_egg = qlogis(parr/eggs)
    
    # fit the BH model for this population
    fit = nls(logit_parr_per_egg ~ qlogis(1/(1/plogis(logit_alpha) + eggs/exp(log_beta))),
              start = c(logit_alpha = qlogis(0.1), log_beta = log(1e5)))
    
    # return point estimates for this population
    c(alpha = unname(plogis(coef(fit)[1])), beta = unname(exp(coef(fit)[2])), sig_Lphi_E_Pb = summary(fit)$sigma)
  }))
  
  # assign rownames
  rownames(BH_params) = colnames(jags_data$Ra_obs)
  
  # return output
  return(BH_params)
}

gen_initials = function(CHAIN, jags_data) {
  
  with(jags_data, {
    # create random BH parameters
    BH_params = get_BH_params(jags_data)
    alpha_init = plogis(qlogis(BH_params[,"alpha"]) + rnorm(nj, 0, 0.1))
    lbeta_init = log(BH_params[,"beta"]) + rnorm(nj, 0, 0.1)
    sig_Lphi_E_Pb_init = exp(log(BH_params[,"sig_Lphi_E_Pb"]) + rnorm(nj, 0, 0.05))
    
    # create random egg-to-parr survival values
    phi_E_Pb = get_phi_E_Pb_obs(jags_data, TRUE, TRUE)
    Lphi_E_Pb_init = qlogis(phi_E_Pb) + matrix(rnorm(prod(dim(phi_E_Pb)), 0, 0.1), nrow(phi_E_Pb), ncol(phi_E_Pb))
    
    # create random parameters for capacity relationship
    lambda_init = runif(1, 8000, 12000)
    sig_lbeta_init = runif(1, 0.05, 0.15)
    
    # create random parameters for proportion of parr leaving in the fall
    mu_pi_init = plogis(qlogis(colMeans(get_pi_obs(jags_data), na.rm = TRUE)) + rnorm(nj, 0, 0.2))
    sig_Lpi_init = apply(qlogis(get_pi_obs(jags_data)), 2, sd, na.rm = TRUE) * rlnorm(nj, 0, 0.05)
    Lpi1_init = qlogis(get_pi_obs(jags_data, TRUE, TRUE)) + matrix(rnorm(ny * nj, 0, 0.4), ny, nj)
    Lpi1_init[1,] = NA
    mu_pi_init = rbind(mu_pi_init, rep(NA, nj))
    
    # create random parameters for omegas
    L_Pb_params = sapply(1:nj, function(j) {
      x = (get_E_obs(jags_data)/10000)[,j]/wul[j]
      y = L_Pb_obs[,j]
      fit = lm(log(y) ~ x)
      out = c(coef(fit), summary(fit)$sigma)
      names(out) = c("omega0", "omega1", "sigma")
      out
    })
    colnames(L_Pb_params) = colnames(jags_data$Ra_obs)
    omega0_init = L_Pb_params["omega0",] * rlnorm(nj, 0, 0.05)
    omega1_init = L_Pb_params["omega1",] * rlnorm(nj, 0, 0.05)
    sig_lL_Pb_init = L_Pb_params["sigma",] * rlnorm(nj, 0, 0.05)
    
    # create random parameters for gammas
    phi_Pa_Mb_params = lapply(1:nj, function(j) {
      x = L_Pb_obs[,j]
      y_fall = get_phi_Pa_Mb_obs(jags_data)[,i_fall,j]
      y_spring = get_phi_Pa_Mb_obs(jags_data)[,i_spring,j]
      fit_fall = lm(qlogis(y_fall) ~ x)
      fit_spring = lm(qlogis(y_spring) ~ x)
      
      out_fall = c(coef(fit_fall), summary(fit_fall)$sigma)
      out_spring = c(coef(fit_spring), summary(fit_spring)$sigma)
      out = rbind(out_fall, out_spring)
      rownames(out) = c("fall-mig", "spring-mig")
      colnames(out) = c("gamma0", "gamma1", "sigma")
      out
    })
    phi_Pa_Mb_params = do.call(abind, append(phi_Pa_Mb_params, list(along = 3)))
    dimnames(phi_Pa_Mb_params)[[3]] = colnames(Ra_obs)
    gamma0_init = phi_Pa_Mb_params[,"gamma0",] * matrix(rlnorm(ni * nj, 0, 0.05), ni, nj)
    gamma1_init = phi_Pa_Mb_params[,"gamma1",] * matrix(rlnorm(ni * nj, 0, 0.05), ni, nj)
    sig_Lphi_Pa_Mb_init = phi_Pa_Mb_params[,"sigma",] * matrix(rlnorm(ni * nj, 0, 0.05), ni, nj)
    
    # create random parameters for overwinter survival
    Lphi_Pa_Mb_init = qlogis(get_phi_Pa_Mb_obs(jags_data, TRUE, TRUE)) + array(rnorm(ny * ni * nj, 0, 0.2), dim = c(ny, ni, nj))
    Lphi_Pa_Mb_init[1,,] = NA
    
    # create random parameters for taus
    phi_Mb_Ma_params = sapply(1:nj, function(j) {
      x = L_Mb_obs[,j]
      y = Lphi_obs_Mb_Ma[,i_spring,o_nor,j]
      fit = lm(y ~ x)
      out = c(coef(fit), summary(fit)$sigma)
      names(out) = c("tau0", "tau1", "sigma")
      out
    })
    colnames(phi_Mb_Ma_params) = colnames(jags_data$Ra_obs)
    # tau0_init = phi_Mb_Ma_params["tau0",] * rlnorm(nj, 0, 0.05)
    # tau1_init = phi_Mb_Ma_params["tau1",] * rlnorm(nj, 0, 0.05)
    tau0_init = c(-6.31, -3.75, -2.68, -8.78) * rlnorm(nj, 0, 0.05)
    tau1_init = c(0.06, 0.04, 0.03, 0.09) * rlnorm(nj, 0, 0.05)
    
    sig_Lphi_Mb_Ma_init = phi_Mb_Ma_params["sigma",] * rlnorm(nj, 0, 0.05)
    sig_Lphi_Mb_Ma_init = rbind(NOR = sig_Lphi_Mb_Ma_init, HOR = rep(NA, nj))
    
    # create random parameters for trib to LGR migration survival
    # Lphi_obs_Mb_Ma
    
    # create random parameters for straying dynamics
    # p_G_init
    # G_random1_init
    # G_random2_init
    
    # create random parameters for carcass correction factors
    # mu_z_init
    # sig_z_init
    # z_init
    
    # create random parameters for NOR return abundance-at-age in years w/o juvenile process model
    # Rb_init
    
    # create random parameters for LGR to ocean migration
    mu_phi_Ma_O0_init = runif(no, 0.4, 0.6)
    sig_Lphi_Ma_O0_init = runif(no, 0.1, 0.5)
    Lphi_Ma_O0_init = sapply(1:no, function(o) rnorm(ny, qlogis(mu_phi_Ma_O0_init[o]), sig_Lphi_Ma_O0_init[o]))
    Lphi_Ma_O0_init[1,] = NA
    
    # create random parameters for first year ocean survival
    kappa_phi_O0_O1_init = runif(nj, 0.2, 0.5)
    mu_phi_O0_O1_init = matrix(c(runif(nj, 0.1, 0.2), rep(NA, nj)), no, nj, byrow = TRUE) 
    sig_Lphi_O0_O1_init = runif(nj, 0.1, 0.5)
    Lphi_O0_O1_resid_init = array(NA, dim = c(ny, no, nj))
    Lphi_O0_O1_resid_init[1,o_nor,] = rnorm(nj, 0, 0.05)
    Lphi_O0_O1_init = array(NA, dim = c(ny, no, nj))
    Lphi_O0_O1_init[2:ny,o_nor,] = rnorm((ny-1) * nj, qlogis(0.15), 0.2)
    
    # create random parameters for second/third year ocean survival
    mu_phi_O1_O2_init = array(runif(no * nj, 0.75, 0.85), dim = c(no, nj))
    mu_phi_O1_O2_init[o_hor,] = NA
    mu_phi_O2_O3_init = array(runif(no * nj, 0.75, 0.85), dim = c(no, nj))
    mu_phi_O2_O3_init[o_hor,] = NA
    
    # create random parameters for difference between NOR and HOR ocean survival
    delta_O0_O1_init = runif(nj, -1, -0.5)
    delta_O1_O2_init = runif(nj, -1, -0.5)
    delta_O2_O3_init = runif(nj, -1, -0.5)
    
    # create random parameters for age-3 maturation
    mu_psi_O1_init = array(runif(no * nj, 0.2, 0.3), dim = c(no, nj))
    sig_Lpsi_O1_init = array(runif(no * nj, 0.1, 0.5), dim = c(no, nj))
    Lpsi_O1_init = qlogis(array(runif(ny * no * nj, 0.1, 0.4), dim = c(ny, no, nj)))
    Lpsi_O1_init[1,,] = NA
    
    # create random parameters for age-4 maturation
    mu_psi_O2_init = array(runif(no * nj, 0.7, 0.8), dim = c(no, nj))
    sig_Lpsi_O2_init = array(runif(no * nj, 0.1, 0.5), dim = c(no, nj))
    Lpsi_O2_init = qlogis(array(runif(ny * no * nj, 0.6, 0.9), dim = c(ny, no, nj)))
    Lpsi_O2_init[1,,] = NA
    
    # create random parameters for BON to LGR migration
    mu_phi_Rb_Ra_init = runif(no, 0.75, 0.85)
    sig_Lphi_Rb_Ra_init = runif(no, 0.1, 0.5)
    Lphi_Rb_Ra_random_init = sapply(1:no, function(o) rnorm(ny, qlogis(mu_phi_Rb_Ra_init[o]), sig_Lphi_Rb_Ra_init[o]))
    Lphi_Rb_Ra_random_init[1,] = NA
    
    # create random parameters for prespawn survival
    mu_phi_Sb_Sa_init = runif(nj, 0.8, 0.95)
    sig_Lphi_Sb_Sa_init = runif(nj, 0.1, 0.5)
    Lphi_Sb_Sa_init = array(qlogis(runif(ny * nj, 0.7, 0.99)), dim = c(ny, nj))
    Lphi_Sb_Sa_init[1,] = NA
    
    # create random parameters for length at summer tagging
    lL_Pb_init = log(matrix(runif(nj * ny, 60, 90), ny, nj))
    lL_Pb_init[1,] = NA
    
    # create random parameters for growth between summer and spring tagging
    lgrowth_init = log(matrix(runif(nj * ny, 1, 1.5), ny, nj))
    lgrowth_init[1,] = NA
    
    # build the output as a list
    list(
      # BH parameters
      alpha = alpha_init,
      lbeta = lbeta_init,
      sig_Lphi_E_Pb = sig_Lphi_E_Pb_init,
      Lphi_E_Pb = Lphi_E_Pb_init,
      
      # habitat capacity parameters
      lambda = lambda_init,
      sig_lbeta = sig_lbeta_init,
      
      # summer length relationship
      omega0 = omega0_init,
      omega1 = omega1_init,
      sig_lL_Pb = sig_lL_Pb_init,
      lL_Pb = lL_Pb_init,
      
      # proportion of parr leaving in fall
      mu_pi = mu_pi_init,
      sig_Lpi = sig_Lpi_init,
      Lpi1 = Lpi1_init,
      
      # overwinter survival
      gamma0 = gamma0_init,
      gamma1 = gamma1_init,
      sig_Lphi_Pa_Mb = sig_Lphi_Pa_Mb_init,
      Lphi_Pa_Mb = Lphi_Pa_Mb_init,
      
      # summer to spring growth
      # mu_growth = mu_growth_init,
      # sig_lgrowth = sig_lgrowth_init,
      lgrowth = lgrowth_init,
      
      # migration to LGR
      # tau0 = tau0_init,
      # tau1 = tau1_init,
      sig_Lphi_Mb_Ma = sig_Lphi_Mb_Ma_init,
      # Lphi_Mb_Ma = Lphi_Mb_Ma_init,
      
      # LGR to ocean
      mu_phi_Ma_O0 = mu_phi_Ma_O0_init,
      sig_Lphi_Ma_O0 = sig_Lphi_Ma_O0_init,
      Lphi_Ma_O0 = Lphi_Ma_O0_init,
      
      # first year ocean survival
      kappa_phi_O0_O1 = kappa_phi_O0_O1_init,
      mu_phi_O0_O1 = mu_phi_O0_O1_init,
      sig_Lphi_O0_O1 = sig_Lphi_O0_O1_init,
      Lphi_O0_O1_resid = Lphi_O0_O1_resid_init,
      Lphi_O0_O1 = Lphi_O0_O1_init,
      
      # second year ocean survival
      # mu_phi_O1_O2 = mu_phi_O1_O2_init,
      
      # third year ocean survival
      # mu_phi_O2_O3 = mu_phi_O2_O3_init,
      
      # NOR:HOR ocean survival scalers
      delta_O0_O1 = delta_O0_O1_init,
      # delta_O1_O2 = delta_O1_O2_init,
      # delta_O2_O3 = delta_O2_O3_init,
      
      # age-3 maturity
      mu_psi_O1 = mu_psi_O1_init,
      sig_Lpsi_O1 = sig_Lpsi_O1_init,
      Lpsi_O1 = Lpsi_O1_init,
      
      # age-4 maturity
      mu_psi_O2 = mu_psi_O2_init,
      sig_Lpsi_O2 = sig_Lpsi_O2_init,
      Lpsi_O2 = Lpsi_O2_init,
      
      # LGR to BON migration
      mu_phi_Rb_Ra = mu_phi_Rb_Ra_init,
      sig_Lphi_Rb_Ra = sig_Lphi_Rb_Ra_init,
      Lphi_Rb_Ra_random = Lphi_Rb_Ra_random_init,
      
      # straying dynamics parameter
      # p_G = p_G_init,
      # G_random1 = G_random1_init,
      # G_random2 = G_random2_init,
      
      # prespawn survival
      mu_phi_Sb_Sa = mu_phi_Sb_Sa_init,
      sig_Lphi_Sb_Sa = sig_Lphi_Sb_Sa_init,
      Lphi_Sb_Sa = Lphi_Sb_Sa_init
    )
  })
}
