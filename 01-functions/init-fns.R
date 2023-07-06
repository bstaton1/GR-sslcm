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
    broodstock = lapply(1:nj, function(j) cbind(broodstock[1:ny_obs,,o_nor,j], broodstock[1:ny_obs,,o_hor,j]))
    broodstock = do.call(abind, append(broodstock, list(along = 3)))
    
    # get the harvest removed
    harvest = H
    harvest = lapply(1:nj, function(j) cbind(harvest[1:ny_obs,,o_nor,j], harvest[1:ny_obs,,o_hor,j]))
    harvest = do.call(abind, append(harvest, list(along = 3)))
    
    # get the pre-spawn survival
    prespawn_surv = colMeans(phi_Sb_Sa, na.rm = TRUE)
    
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
    
    # remove brood_stock & harvest
    age_origin_above_weir = age_origin_return - broodstock - harvest
    age_origin_above_weir[age_origin_above_weir < 0] = 5
    
    # remove prespawn morts
    age_origin_survived = lapply(1:nj, function(j) {
      apply(age_origin_above_weir[,,j], 2, function(a) a * prespawn_surv[j])
    })
    age_origin_survived = do.call(abind, append(age_origin_survived, list(along = 3)))
    
    # get number of females by age/origin
    p_female_age_origin = rbind(Omega, Omega)
    females_age_origin = lapply(1:nj, function(j) {
      t(apply(age_origin_survived[,,j], 1, function(y) y * p_female_age_origin[,j]))
    })
    females_age_origin = do.call(abind, append(females_age_origin, list(along = 3)))
    
    # get total egg output
    eggs_age_origin = lapply(1:nj, function(j) {
      f_age_origin = cbind(f[,,j], f[,,j])
      t(sapply(1:nrow(females_age_origin), function(y) females_age_origin[y,,j] * f_age_origin[y,]))
    })
    eggs_age_origin = do.call(abind, append(eggs_age_origin, list(along = 3)))
    total_eggs = apply(eggs_age_origin, 3, rowSums)
    
    # add population names back in
    colnames(total_eggs) = colnames(total_return)
    
    # add year names back in
    rownames(total_eggs) = rownames(total_return)
    
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
    
    # function to predict parr to egg survival based on BH params and egg abundance
    BH = function(eggs, alpha, beta) 1/(1/alpha + eggs/beta)
    
    # negative log likelihood function to minimize with optim
    nll = function(theta, parr, eggs) {
      # transform parameters from estimated scale to useful scale
      alpha = plogis(theta[1])
      beta = exp(theta[2])
      sigma = exp(theta[3])
      
      # produce predicted survival based on parameters
      preds = BH(eggs, alpha, beta)
      
      # produce observed values
      obs = parr/eggs
      
      # for CAT, one very large outcome prevented convergence
      # remove the largest value for all populations prior to fitting
      # remember, this is just for generating means for the initial value generator
      obs[which.max(obs)] = NA
      
      # return the negative log likelihood
      -sum(dnorm(qlogis(obs), qlogis(preds), sigma, log = TRUE), na.rm = TRUE)
    }
    
    # initial values for MLE
    inits = c(qlogis(0.2), log(max(parr, na.rm = TRUE)), log(0.5))
    
    # minimize nll to fit the model
    fit = optim(inits, nll, parr = parr, eggs = eggs, method = "L-BFGS-B",
                lower = c(-3, log(50000), log(0.1)), upper = c(3, log(1e6), log(1)))

    # return point estimates for this population
    c(alpha = plogis(fit$par[1]), beta = exp(fit$par[2]), sig_Lphi_E_Pb = exp(fit$par[3]))
  }))
  
  # assign rownames
  rownames(BH_params) = colnames(jags_data$Ra_obs)
  
  # return output
  return(BH_params)
}

# simulate a covariance matrix
# create a vcov matrix, supply the standard deviation vect (sqrt(diag(vcov))), and a type of correlation
sim_Sigma = function(sig, correlation) {
  
  # other functions
  cor_fns = list(
    none = function(n) {
      m = matrix(runif(n^2, 0, 0), n, n)
      m = (m * lower.tri(m)) + t(m * lower.tri(m)) 
      diag(m) = 1
      m
    },
    low = function(n) {
      m = matrix(runif(n^2, 0.01, 0.29), n, n)
      m = (m * lower.tri(m)) + t(m * lower.tri(m)) 
      diag(m) = 1
      m
    },
    moderate = function(n) {
      m = matrix(runif(n^2, 0.3, 0.49), n, n)
      m = (m * lower.tri(m)) + t(m * lower.tri(m)) 
      diag(m) = 1
      m
    },
    high = function(n) {
      m = matrix(runif(n^2, 0.5, 0.95), n, n)
      m = (m * lower.tri(m)) + t(m * lower.tri(m)) 
      diag(m) = 1
      m
    },
    near_perfect = function(n) {
      m = matrix(runif(n^2, 0.96, 0.99999), n, n)
      m = (m * lower.tri(m)) + t(m * lower.tri(m)) 
      diag(m) = 1
      m
    }
  )
  
  cor2cov = function(cor, sig) {
    n = length(sig)
    cov = matrix(NA, n, n)
    for (i in 1:n) {
      for (j in 1:n) {
        cov[i,j] = sig[i] * sig[j] * cor[i,j]
      }
    }
    cov
  }
  
  if (!(correlation %in% names(cor_fns))) {
    stop("invalid correlation type specified. The accepted values are: \n\n  ", 
         paste(paste(paste0("'", names(cor_fns), "'"), c(rep(", ", length(names(cor_fns)) - 2), ", or ", ""), sep = ""), collapse = "")
    )
  }
  
  # how many elements  
  n = length(sig)
  
  # continue creating matrices until one is found that's positive definite
  cor = cor_fns[[correlation]](length(sig))
  Sigma = cor2cov(cor, sig)
  while(!matrixcalc::is.positive.definite(Sigma)) {
    Sigma = cor2cov(cor_fns[[correlation]](length(sig)), sig)
  }
  Sigma
}

gen_initials = function(CHAIN, jags_data) {
  
  with(jags_data, {
    # create random BH parameters
    # BH_params = get_BH_params(jags_data)
    BH_params = matrix(c(0.3, 0.25, 0.2, 0.2, 1e5, 3e5, 5e5, 1e5, 0.4, 0.6, 0.6, 0.8), nrow = 4, ncol = 3)
    rownames(BH_params) = c("CAT", "LOS", "MIN", "UGR")
    colnames(BH_params) = c("alpha", "beta", "sig_Lphi_E_Pb")
    alpha_init = plogis(qlogis(BH_params[,"alpha"]) + rnorm(nj, 0, 0.1))
    lbeta_init = log(BH_params[,"beta"]) + rnorm(nj, 0, 0.1)
    sig_Lphi_E_Pb_init = exp(log(BH_params[,"sig_Lphi_E_Pb"]) + rnorm(nj, 0, 0.05))
    kappa_phi_E_Pb_init = runif(nj, 0.2, 0.5)

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
      fit = lm(log(y) ~ log(x))
      out = c(coef(fit), summary(fit)$sigma)
      names(out) = c("omega0", "omega1", "sigma")
      out
    })
    colnames(L_Pb_params) = colnames(jags_data$Ra_obs)
    omega0_init = L_Pb_params["omega0",] * rlnorm(nj, 0, 0.05)
    omega0_init[omega0_init < omega0_prior[1]] = log(exp(omega0_prior[1]) + 5)
    omega0_init[omega0_init > omega0_prior[2]] = log(exp(omega0_prior[2]) - 5)
    omega1_init = L_Pb_params["omega1",] * rlnorm(nj, 0, 0.05)
    sig_lL_Pb_init = L_Pb_params["sigma",] * rlnorm(nj, 0, 0.05)
    
    # create random parameters for gammas
    phi_Pa_Mb_params = lapply(1:nj, function(j) {
      x = (L_Pb_obs[,j] - L_Pb_center[j])/L_Pb_scale[j]
      y_fall = get_phi_Pa_Mb_obs(jags_data)[,i_fall,j]
      y_spring = get_phi_Pa_Mb_obs(jags_data)[,i_spring,j]
      y_mean = rowMeans(cbind(y_fall, y_spring), na.rm = TRUE)
      
      fit_fall = lm(qlogis(y_fall) ~ x)
      fit_spring = lm(qlogis(y_spring) ~ x)
      fit_mean = lm(qlogis(y_mean) ~ x)
      
      out_fall = c(coef(fit_fall), summary(fit_fall)$sigma)
      out_spring = c(coef(fit_spring), summary(fit_spring)$sigma)
      out_mean = c(coef(fit_mean), summary(fit_mean)$sigma)
      # out = rbind(out_fall, out_spring)
      out = rbind(out_mean, out_mean)
      rownames(out) = c("fall-mig", "spring-mig")
      colnames(out) = c("gamma0", "gamma1", "sigma")
      out
    })
    phi_Pa_Mb_params = do.call(abind, append(phi_Pa_Mb_params, list(along = 3)))
    dimnames(phi_Pa_Mb_params)[[3]] = colnames(Ra_obs)
    gamma0_init = phi_Pa_Mb_params[,"gamma0",] * matrix(rlnorm(ni * nj, 0, 0.05), ni, nj)
    gamma1_init = phi_Pa_Mb_params[,"gamma1",] * matrix(rlnorm(ni * nj, 0, 0.05), ni, nj)
    gamma1_init[i_spring,] = NA
    sig_Lphi_Pa_Mb_init = phi_Pa_Mb_params[,"sigma",] * matrix(rlnorm(ni * nj, 0, 0.05), ni, nj)
    
    # create random parameters for overwinter survival
    Lphi_Pa_Mb_init = qlogis(get_phi_Pa_Mb_obs(jags_data, TRUE, TRUE)) + array(rnorm(ny * ni * nj, 0, 0.2), dim = c(ny, ni, nj))
    Lphi_Pa_Mb_init[1,,] = NA
    
    # create random parameters for thetas
    lDelta_L_Pb_Mb_params = sapply(1:nj, function(j) {
      y1 = L_Pb_obs[,j]
      y2 = L_Mb_obs[,j]
      x = (L_Pb_obs[,j] - L_Pb_center[j])/L_Pb_scale[j]
      y = log(y2/y1)
      fit = lm(y ~ x)
      out = c(coef(fit), summary(fit)$sigma)
      names(out) = c("theta0", "theta1", "sig_lDelta_L_Pb_Mb")
      out
    })
    colnames(lDelta_L_Pb_Mb_params) = colnames(jags_data$Ra_obs)
    theta0_init = lDelta_L_Pb_Mb_params["theta0",] * rlnorm(nj, 0, 0.05)
    theta1_init = lDelta_L_Pb_Mb_params["theta1",] * rlnorm(nj, 0, 0.05)
    lDelta_L_Pb_Mb_init = array(NA, dim = c(ny, nj))
    lDelta_L_Pb_Mb_init[1:ny_obs,] = log(L_Mb_obs/L_Pb_obs)
    lDelta_L_Pb_Mb_init[is.na(lDelta_L_Pb_Mb_init)] = mean(lDelta_L_Pb_Mb_init, na.rm = TRUE)
    lDelta_L_Pb_Mb_init = lDelta_L_Pb_Mb_init + matrix(rnorm(ny * nj, 0, 0.03), ny, nj)
    lDelta_L_Pb_Mb_init[1,] = NA
    
    # create random parameters for taus
    phi_Mb_Ma_params = sapply(1:nj, function(j) {
      x = (L_Mb_obs[,j] - L_Mb_center[j])/L_Mb_scale[j]
      y = Lphi_obs_Mb_Ma[,i_spring,o_nor,j]
      fit = lm(y ~ x)
      out = c(coef(fit), summary(fit)$sigma)
      names(out) = c("tau0", "tau1", "sigma")
      out
    })
    colnames(phi_Mb_Ma_params) = colnames(jags_data$Ra_obs)
    tau0_init = phi_Mb_Ma_params["tau0",] * rlnorm(nj, 0, 0.05)
    tau1_init = phi_Mb_Ma_params["tau1",] * rlnorm(nj, 0, 0.05)
    # tau0_init = c(-6.31, -3.75, -2.68, -8.78) * rlnorm(nj, 0, 0.05)
    # tau1_init = c(0.06, 0.04, 0.03, 0.09) * rlnorm(nj, 0, 0.05)
    
    sig_Lphi_Mb_Ma_init = phi_Mb_Ma_params["sigma",] * rlnorm(nj, 0, 0.05)
    sig_Lphi_Mb_Ma_init = rbind(NOR = sig_Lphi_Mb_Ma_init, HOR = rep(NA, nj))
    
    # create random parameters for trib to LGR migration survival
    mu_phi_Mb_Ma_init = array(NA, dim = c(ni,no,nj))
    mu_phi_Mb_Ma_init[i_spring,o_hor,] = plogis(colMeans(Lphi_obs_Mb_Ma[,i_spring,o_hor,], na.rm = TRUE) + rnorm(nj, 0, 0.05))
    mu_phi_Mb_Ma_init[i_spring,o_hor,3] = 0.5
    
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
    # Lphi_O0_O1_resid_init = array(NA, dim = c(ny, no, nj))
    # Lphi_O0_O1_resid_init[1,o_nor,] = rnorm(nj, 0, 0.05)
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
    
    # create random parameters for length at summer tagging
    lL_Pb_init = log(matrix(runif(nj * ny, 60, 90), ny, nj))
    lL_Pb_init[1,] = NA
    
    Lphi_E_Pb_resid_init = matrix(NA, ny, nj)
    Lphi_E_Pb_resid_init[1,1:nj] = 0
    Lphi_O0_O1_resid_init = array(NA, dim = c(ny, no, nj))
    Lphi_O0_O1_resid_init[1,o_nor,1:nj] = 0
    
    
    # build the output as a list
    list(
      # BH parameters
      alpha = alpha_init,
      lbeta = lbeta_init,
      # sig_Lphi_E_Pb = sig_Lphi_E_Pb_init,
      kappa_phi_E_Pb = kappa_phi_E_Pb_init,
      Lphi_E_Pb = Lphi_E_Pb_init,
      # Lphi_E_Pb_resid = Lphi_E_Pb_resid_init,
      Tau_Lphi_E_Pb = solve(sim_Sigma(rep(0.3, nj), "low")),
      
      # habitat capacity parameters
      lambda = lambda_init,
      sig_lbeta = sig_lbeta_init,
      
      # summer length relationship
      omega0 = omega0_init,
      omega1 = omega1_init,
      # sig_lL_Pb = sig_lL_Pb_init,
      lL_Pb = lL_Pb_init,
      Tau_lL_Pb = solve(sim_Sigma(rep(0.05, nj), "low")),
      
      # proportion of parr leaving in fall
      mu_pi = mu_pi_init,
      # sig_Lpi = sig_Lpi_init,
      Lpi1 = Lpi1_init,
      Tau_Lpi = solve(sim_Sigma(rep(0.3, nj), "low")),
      
      # overwinter survival
      gamma0 = gamma0_init,
      gamma1 = gamma1_init,
      # sig_Lphi_Pa_Mb = sig_Lphi_Pa_Mb_init,
      Lphi_Pa_Mb = Lphi_Pa_Mb_init,
      Tau_Lphi_Pa_Mb = solve(sim_Sigma(rep(0.3, nj), "low")),
      
      # summer to spring growth
      theta0 = theta0_init,
      theta1 = theta1_init,
      lDelta_L_Pb_Mb = lDelta_L_Pb_Mb_init,
      Tau_lDelta_L_Pb_Mb = solve(sim_Sigma(rep(0.05, nj), "low")),
      
      # migration to LGR
      # tau0 = tau0_init,
      # tau1 = tau1_init,
      mu_phi_Mb_Ma = mu_phi_Mb_Ma_init,
      # sig_Lphi_Mb_Ma = sig_Lphi_Mb_Ma_init,
      # Lphi_Mb_Ma = Lphi_Mb_Ma_init,
      Tau_Lphi_Mb_Ma = solve(sim_Sigma(rep(0.3, nj), "low")),
      
      # LGR to ocean
      mu_phi_Ma_O0 = mu_phi_Ma_O0_init,
      sig_Lphi_Ma_O0 = sig_Lphi_Ma_O0_init,
      Lphi_Ma_O0 = Lphi_Ma_O0_init,
      Tau_Lphi_Ma_O0 = solve(sim_Sigma(rep(0.3, no), "low")),
      
      # first year ocean survival
      kappa_phi_O0_O1 = kappa_phi_O0_O1_init,
      mu_phi_O0_O1 = mu_phi_O0_O1_init,
      # sig_Lphi_O0_O1 = sig_Lphi_O0_O1_init,
      # Lphi_O0_O1_resid = Lphi_O0_O1_resid_init,
      # Lphi_O0_O1_resid = Lphi_O0_O1_resid_init,
      Lphi_O0_O1 = Lphi_O0_O1_init,
      Tau_Lphi_O0_O1 = solve(sim_Sigma(rep(0.15, nj), "low")),
      
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
      # sig_Lpsi_O1 = sig_Lpsi_O1_init,
      Lpsi_O1 = Lpsi_O1_init,
      Tau_Lpsi_O1 = solve(sim_Sigma(rep(0.15, nj), "low")),
      
      # age-4 maturity
      mu_psi_O2 = mu_psi_O2_init,
      # sig_Lpsi_O2 = sig_Lpsi_O2_init,
      Lpsi_O2 = Lpsi_O2_init,
      Tau_Lpsi_O2 = solve(sim_Sigma(rep(0.15, nj), "low")),
      
      # LGR to BON migration
      mu_phi_Rb_Ra = mu_phi_Rb_Ra_init,
      # sig_Lphi_Rb_Ra = sig_Lphi_Rb_Ra_init,
      Lphi_Rb_Ra_random = Lphi_Rb_Ra_random_init,
      Tau_Lphi_Rb_Ra = solve(sim_Sigma(rep(0.15, no), "low")),
      
      # straying dynamics parameter
      # p_G = p_G_init,
      # G_random1 = G_random1_init,
      # G_random2 = G_random2_init,
      
      # prespawn survival
      mu_phi_Sb_Sa = mu_phi_Sb_Sa_init
    )
  })
}
