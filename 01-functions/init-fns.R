# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT TO HOUSE FUNCTIONS FOR GENERATING INITIAL MCMC VALUES #
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #

##### RECONSTRUCT PARR AND SPAWNER "OBSERVATIONS" AND FIT BH MODEL #####

# function to fit a basic spawner to parr BH function
fit_basic_BH = function(jags_data, plot = FALSE) {
  out = with(jags_data, {
    # get "observed" smolt at LGD
    Ma_obs = Pa_obs[,"fall-mig"] * expit(Lphi_obs_Pa_Ma[,"fall-mig"]) + 
      Mb_obs[,"spring-mig",1] * expit(Lphi_obs_Mb_Ma[,"spring-mig",1])
    
    # get "observed" summer parr
    Pb_obs = Ma_obs/expit(Lphi_obs_Pb_Ma)
    
    # get "observed" pre-spawn mortality
    mu_phi_Sb_Sa = mean(carcs_spawned/carcs_sampled, na.rm = T)
    
    # get "observed" spawners
    Sa_obs = Ra_obs * mu_phi_Sb_Sa
    
    # fit a basic BH model (assume lognormal errors)
    fit = nls(log(Pb_obs) ~ log(BH(Sa_obs, exp(log_alpha), exp(log_beta))),
              start = c(log_alpha = log(300), log_beta = log(300000)))
    
    # extract relevant information and return it
    list(BH_dat = data.frame(brood_year = as.numeric(names(Pb_obs)),
                             Pb_obs = Pb_obs, Sa_obs = Sa_obs),
         BH_fit = fit)
  })
  
  # draw a basic plot if requested
  if (plot) {
    newx = seq(0, max(out$BH_dat$Sa_obs, na.rm = T) * 1.2, length = 30)
    pred = exp(predict(out$BH_fit, data.frame(Sa_obs = newx)))
    plot(out$BH_dat$Pb_obs ~ out$BH_dat$Sa_obs, pch = 16, col = "blue",
         ylim = c(0, max(out$BH_dat$Pb_obs, na.rm = T)),
         xlim = c(0, max(out$BH_dat$Sa_obs, na.rm = T)) * 1.2,
         xlab = "Spawners", ylab = "Parr")
    lines(pred ~ newx, lwd = 2)
  }
  
  # return output
  return(out)
}

# function to fit a basic egg to parr BH function
fit_fecund_BH = function(jags_data, pop, plot = FALSE) {
  out = with(jags_data, {
    # get "observed" smolt at LGD
    Ma_obs = Pa_obs[,"fall-mig",pop] * expit(Lphi_obs_Pa_Ma[,"fall-mig",pop]) + 
      Mb_obs[,"spring-mig",1,pop] * expit(Lphi_obs_Mb_Ma[,"spring-mig",1,pop])
    
    # get "observed" summer parr
    Pb_obs = Ma_obs/expit(Lphi_obs_Pb_Ma[,pop])
    
    # get "observed" pre-spawn mortality
    mu_phi_Sb_Sa = mean(carcs_spawned[,pop]/carcs_sampled[,pop], na.rm = TRUE)
    
    # get "observed" spawners
    Sa_obs = Ra_obs[,pop] * mu_phi_Sb_Sa
    
    # set approximate composition (proportion of all spawners that are females of age 3, 4, or 5)
    q = c(0, 0.4, 0.05)
    
    # get approximate egg production
    f_tot = sapply(Sa_obs, function(x) {
      sum(x * q * f)
    })
    
    # fit a basic BH model (assume lognormal errors)
    fit = nls(log(Pb_obs) ~ log(BH(f_tot, exp(log_alpha), exp(log_beta))),
              start = c(log_alpha = log(0.1), log_beta = log(300000)))
    
    # extract relevant information and return it
    list(BH_dat = data.frame(brood_year = as.numeric(names(Pb_obs)),
                             Pb_obs = Pb_obs, Sa_obs = Sa_obs, f_tot = f_tot),
         BH_fit = fit)
  })
  
  # draw a basic plot if requested
  if (plot) {
    newx = seq(0, max(out$BH_dat$f_tot, na.rm = T) * 1.2, length = 30)
    pred = exp(predict(out$BH_fit, data.frame(f_tot = newx)))
    plot(out$BH_dat$Pb_obs ~ out$BH_dat$f_tot, pch = 16, col = "blue",
         ylim = c(0, max(out$BH_dat$Pb_obs, na.rm = T)),
         xlim = c(0, max(out$BH_dat$f_tot, na.rm = T)) * 1.2,
         xlab = "Total Egg Production", ylab = "Parr")
    lines(pred ~ newx, lwd = 2)
  }
  
  # return output
  return(out)
}

##### GENERATE SOME INITIAL VALUES FOR MCMC #####

gen_initials = function(c, jags_data) {
  # fit a simple BH model to each population
  BH_ests = lapply(1:jags_data$nj, function(j) {
    summary(fit_fecund_BH(jags_data, j)$BH_fit)$coef
  })
  
  Pb_ests = sapply(1:jags_data$nj, function(j) {
    Pb_ests = fit_fecund_BH(jags_data, j)$BH_dat$Pb_obs
    mn_Pb = mean(Pb_ests, na.rm = TRUE)
    fill_yrs = (jags_data$kmax+1):jags_data$ny
    Pb_ests[fill_yrs][is.na(Pb_ests[fill_yrs])] = mn_Pb
    Pb_ests
  })

  # generate random initial values for all years where strays can be present
  n_stray_tot = sapply(1:jags_data$nj, function(j) {
    stray_yrs = as.numeric(na.omit(jags_data$stray_yrs[,j]))
    comp = t(apply(jags_data$carc_x_obs[stray_yrs,,j], 1, function(x) x/sum(x)))
    p_nat = rowSums(comp[,1:jags_data$nk])
    p_nat[is.na(p_nat)] = 1
    n_strays = jags_data$Ra_obs[stray_yrs,j] * (1 - p_nat)
    out = rep(NA, jags_data$ny)
    out[stray_yrs] = runif(jags_data$n_stray_yrs[j], n_strays * 0.8, n_strays + mean(n_strays) * 1.2)
    
    out
  })
  
  mu_list = with(jags_data, {

    # hydropower survival
    mu_phi_Ma_O0 = runif(no, 0.4, 0.6)
    
    # migration survival from trib to LGR
    mu_phi_Mb_Ma = array(NA, dim = c(ni, no, nj))
    mu_phi_Mb_Ma[i_spring,1:no,1:nj] = runif(no*nj, 0.4, 0.6)
    
    # first year ocean survival
    mu_phi_O0_O1 = matrix(c(runif(nj, 0.05, 0.15), rep(NA, nj)), no, nj, byrow = TRUE)
    
    # second year ocean survival
    mu_phi_O1_O2 = matrix(c(runif(nj, 0.6, 0.8), rep(NA, nj)), no, nj, byrow = TRUE)
    
    # pre-spawn survival
    mu_phi_Sb_Sa = runif(nj, 0.8, 1)
    
    # LH apportionment
    mu_pi = array(NA, dim = c(ni, nj))
    mu_pi[i_fall,] = runif(nj, 0.2, 0.4)
    
    # maturity: at total age 3
    mu_psi_O1_Rb = array(NA, dim = c(no,nj))
    mu_psi_O1_Rb[1:no,1:nj] = runif(no*nj, 0.2, 0.3)

    # maturity: at total age 4
    mu_psi_O2_Rb = array(NA, dim = c(no,nj))
    mu_psi_O2_Rb[1:no,1:nj] = runif(no*nj, 0.7, 0.9)

    list(
      mu_phi_Ma_O0 = mu_phi_Ma_O0,
      mu_phi_Mb_Ma = mu_phi_Mb_Ma, 
      mu_phi_O0_O1 = mu_phi_O0_O1,
      mu_phi_O1_O2 = mu_phi_O1_O2,
      mu_phi_Sb_Sa = mu_phi_Sb_Sa,
      mu_pi = mu_pi,
      mu_psi_O1_Rb = mu_psi_O1_Rb,
      mu_psi_O2_Rb = mu_psi_O2_Rb
    )
  })
  
  # average adult recruits
  mean_Ra = colMeans(jags_data$Ra_obs, na.rm = TRUE)
  mu_init_recruits = runif(jags_data$nj, mean_Ra * 0.8, mean_Ra * 1.2)
  
  # generate random initial values for one chain
  base_list = list(
    alpha = sapply(BH_ests, function(ests) exp(rnorm(1, ests["log_alpha","Estimate"], ests["log_alpha","Std. Error"]))),
    log_beta = sapply(BH_ests, function(ests) rnorm(1, ests["log_beta","Estimate"], ests["log_beta","Std. Error"])),
    Pb = Pb_ests * matrix(exp(rnorm(jags_data$ny * jags_data$nj, 0, 0.1)), jags_data$ny, jags_data$nj),
    sigma_Pb = runif(jags_data$nj, 0.3, 0.4),
    n_stray_tot = n_stray_tot, 
    mu_init_recruits = mu_init_recruits
  )
  
  append(base_list, mu_list)
}
