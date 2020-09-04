jags_model_code = function() {
  
  ### PRIORS: RECRUITMENT FUNCTION ###
  # spawner to aggregate parr recruitment function
  log_alpha ~ dnorm(0, 0.001)             # log productivity
  alpha <- exp(log_alpha)
  log_beta ~ dnorm(0, 0.001) %_% T(,15)    # log capacity. bound to prevent nonsensically large draws
  beta <- exp(log_beta)
  sigma_Pb ~ dunif(0, 5)
  
  ### PRIORS: FRESHWATER PARAMETERS ###
  # aggregate parr to LH-specific parr
  mu_pi[1] ~ dbeta(1, 1)
  sig_Lpi ~ dunif(0, 5)
  mu_pi[2] <- 1 - mu_pi[1]
  
  # overwinter survival: LH-specific
  for (i in 1:ni) {
    mu_phi_Pa_Mb[i] ~ dbeta(1, 1)  
    sig_Lphi_Pa_Mb[i] ~ dunif(0, 5)
  }

  # movement survival (trib to LGD): estimate for spring migrants and assume the same value for fall migrants
  mu_phi_Mb_Ma[2] ~ dbeta(1, 1)
  sig_Lphi_Mb_Ma[2] ~ dunif(0, 5)
  mu_phi_Mb_Ma[1] <- mu_phi_Mb_Ma[2]
  sig_Lphi_Mb_Ma[1] <- sig_Lphi_Mb_Ma[2]
  
  # movement survival (LGD to estuary): same for both LH types
  mu_phi_Ma_M ~ dbeta(1, 1)
  sig_Lphi_Ma_M ~ dunif(0, 5)
  
  ### PRIORS: MATURATION ###
  # probability of returning as female ([1]) or male ([2])
  mu_omega[1] ~ dbeta(1, 1)
  mu_omega[2] <- 1 - mu_omega[1]
  sig_Lomega ~ dunif(0, 5)
  
  # sex-specific maturation probabilities
  for (s in 1:ns) {
    # pr(return at SWA1)
    mu_psi_O1_Rb[s] ~ dbeta(1, 1)     
    sig_Lpsi_O1_Rb[s] ~ dunif(0, 5)
    
    # pr(return at SWA2|not returned at SWA1)
    mu_psi_O2_Rb[s] ~ dbeta(1, 1)
    sig_Lpsi_O2_Rb[s] ~ dunif(0, 5)
    
    # pr(return at SWA3|not returned at SWA1 or SWA2)
    mu_psi_O3_Rb[s] <- 1
    sig_Lpsi_O3_Rb[s] <- 0
  }
  
}
