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
fit_fecund_BH = function(jags_data, plot = FALSE) {
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
    
    # set approximate [age,sex] composition
    q = matrix(c(0, 0.4, 0.05, 0.2, 0.3, 0.05), nrow = 3, ncol = 2)
    
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
  # fit a simple BH model
  BH_fit = fit_fecund_BH(jags_data, plot = FALSE)$BH_fit
  
  # extract the parameter summaries
  BH_ests = summary(BH_fit)
  
  # generate random initial values for one chain
  list(
    alpha = exp(rnorm(1, BH_ests$coef["log_alpha","Estimate"], BH_ests$coef["log_alpha","Std. Error"])),
    log_beta = rnorm(1, BH_ests$coef["log_beta","Estimate"], BH_ests$coef["log_beta","Std. Error"]),
    sigma_Pb = runif(1, BH_ests$sigma - 0.1, BH_ests$sigma + 0.1)
  )
}
