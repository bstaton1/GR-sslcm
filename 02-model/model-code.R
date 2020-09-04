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
  
  ### PRIORS: BROOD-YEAR-SPECIFIC PARAMETERS ###
  for (y in (kmax+1):ny) {
    # aggregate parr to LH-specific parr
    Lpi1[y] ~ dnorm(logit(mu_pi[1]), 1/sig_Lpi^2)
    pi[y,1] <- ilogit(Lpi1[y])
    pi[y,2] <- 1 - pi[y,1]
    
    # overwinter survival
    for (i in 1:ni) {
      Lphi_Pa_Mb[y,i] ~ dnorm(logit(mu_phi_Pa_Mb[i]), 1/sig_Lphi_Pa_Mb[i]^2)
      phi_Pa_Mb[y,i] <- ilogit(Lphi_Pa_Mb[y,i])
    }
    
    # movement survival: trib to LGD
    # assume equal between LH types
    Lphi_Mb_Ma[y,2] ~ dnorm(logit(mu_phi_Mb_Ma[2]), 1/sig_Lphi_Mb_Ma[2]^2)
    phi_Mb_Ma[y,2] <- ilogit(Lphi_Mb_Ma[y,2])
    phi_Mb_Ma[y,1] <- phi_Mb_Ma[y,2]
    
    # movement survival: LGD to estuary
    Lphi_Ma_M[y] ~ dnorm(logit(mu_phi_Ma_M), 1/sig_Lphi_Ma_M^2)
    phi_Ma_M[y] <- ilogit(Lphi_Ma_M[y])
    
    # probability of returning as female ([1]) or male ([2])
    Lomega1[y] ~ dnorm(logit(mu_omega[1]), 1/sig_Lomega^2)
    omega[y,1] <- ilogit(Lomega1[y])
    omega[y,2] <- 1 - omega[y,1]
    
    # sex-specific maturation probabilities
    for (s in 1:ns) {
      # pr(return at SWA1)
      Lpsi_O1_Rb[y,s] ~ dnorm(logit(mu_psi_O1_Rb[s]), 1/sig_Lpsi_O1_Rb[s]^2)
      psi_O1_Rb[y,s] <- ilogit(Lpsi_O1_Rb[y,s])
      
      # pr(return at SWA2|not returned at SWA1)
      Lpsi_O2_Rb[y,s] ~ dnorm(logit(mu_psi_O2_Rb[s]), 1/sig_Lpsi_O2_Rb[s]^2)
      psi_O2_Rb[y,s] <- ilogit(Lpsi_O2_Rb[y,s])
      
      # pr(return at SWA3|not returned at SWA1 or SWA2)
      psi_O3_Rb[y,s] <- 1
    }
    # ocean survival
    phi_M_O1[y] <- mu_phi_M_O1
    phi_O1_O2[y] <- mu_phi_O1_O2
    phi_O2_O3[y] <- mu_phi_O2_O3
    
    # upstream adult survival
    phi_Rb_Ra[y] <- mu_phi_Rb_Ra
    
    # pre-spawn survival
    phi_Sb_Sa[y] <- mu_phi_Sb_Sa
  }
  
  ### PROCESS MODEL: INITIALIZATION ###
  # need to do something special with first kmax brood years prior to data collection
  # to populate adult states for reproduction/observation in the first years with data
  # THIS APPROACH: 
  #  Estimate the mean and year-specific TOTAL ADULT RECRUITS (arriving to mouth of river) 
  #  for the first kmax brood years. Then, based on the hyperparams of return-by-sex, return-by-age, and ocean survival
  #  (this is the p_init_prime and p_init stuff below)
  #  derive the probability of returning at age for each sex that will apply equally to all kmax of the first brood years
  #  Then place these adult recruits in the appropriate return year by age/sex
  
  # obtain expected probability of adult recruits returning at age for each sex
  for (s in 1:ns) {
    p_init_prime[1,s] <- mu_phi_M_O1 * mu_psi_O1_Rb[s]
    p_init_prime[2,s] <- mu_phi_M_O1 * (1 - mu_psi_O1_Rb[s]) * mu_phi_O1_O2 * mu_psi_O2_Rb[s]
    p_init_prime[3,s] <- mu_phi_M_O1 * (1 - mu_psi_O1_Rb[s]) * mu_phi_O1_O2 * (1 - mu_psi_O2_Rb[s]) * mu_phi_O2_O3 * mu_psi_O3_Rb[s]
    for (k in 1:nk) {
      p_init[k,s] <- p_init_prime[k,s]/sum(p_init_prime[1:nk,s])
    }
  }
  
  # hyperparameters for initialization
  mu_init_recruits ~ dunif(0, max_init_recruits)
  sig_init_lrecruits ~ dunif(0,10)
  
  for (y in 1:kmax) {
    # random, constant-mean adult recruits for first kmax brood years
    init_recruits[y] ~ dlnorm(log(mu_init_recruits), 1/sig_init_lrecruits^2) %_% T(,max_init_recruits)
    for (s in 1:ns) {
      # apportion them to each sex
      init_recruits_sex[y,s] <- init_recruits[y] * mu_omega[s]
      for (k in 1:nk) {
        # apportion them to return year, age, and sex
        Rb[y+kmin+k-1,k,s] <- init_recruits_sex[y,s] * p_init[k,s] 
      }
    }
  }
  
  ### PROCESS MODEL: COMPLETE LIFE CYCLE FOR REMAINING BROOD YEARS ###
  for (y in (kmax+1):ny) {
    
    # reproductive link: total summer parr
    Pb_pred[y] <- Sa_tot[y]/(1/alpha + Sa_tot[y]/beta)
    Pb[y] ~ dlnorm(log(Pb_pred[y]), 1/sigma_Pb^2)
    
    for (i in 1:ni) {
      # apportion to LH-strategy: parr after fall migration
      Pa[y,i] <- Pb[y] * pi[y,i]

      # survive over winter: smolt before spring migration
      Mb[y,i] <- Pa[y,i] * phi_Pa_Mb[y,i]

      # move to LGD: smolt after spring migration, at top of LGD
      Ma[y,i] <- Mb[y,i] * phi_Mb_Ma[y,i]
      
      # derived survival for fitting: fall trap to LGD
      phi_Pa_Ma[y,i] <- Ma[y,i]/Pa[y,i]
    }
    
    # derived survival for fitting: summer tagging to LGD
    phi_Pb_Ma[y] <- sum(Ma[y,1:ni])/Pb[y]
    
    # sex-specific processes
    for (s in 1:ns) {
      # move to estuary and assign to sex
      M[y,s] <- sum(Ma[y,1:ni]) * phi_Ma_M[y] * omega[y,s]

      # move juveniles through ocean ages and survivals
      O[y,1,s] <- M[y,s] * phi_M_O1[y] # survive first winter at sea. now SWA1, TA3
      O[y,2,s] <- O[y,1,s] * (1 - psi_O1_Rb[y,s]) * phi_O1_O2[y] # don't mature at SWA1 and survive second winter at sea. now SWA2, TA4
      O[y,3,s] <- O[y,2,s] * (1 - psi_O2_Rb[y,s]) * phi_O2_O3[y] # don't mature at SWA2 and survive third winter at sea. now SWA3, TA5
    
      # mature and return to river in appropriate year at age/sex
      Rb[y+kmin+1-1,1,s] <- O[y,1,s] * psi_O1_Rb[y,s]
      Rb[y+kmin+2-1,2,s] <- O[y,2,s] * psi_O2_Rb[y,s]
      Rb[y+kmin+3-1,3,s] <- O[y,3,s] * psi_O3_Rb[y,s]
    }
    
    # adult in-river processes
    # y now represents brood year these fish are returning in
    for (k in 1:nk) {
      for (s in 1:ns) {
        # survive upstream migration
        Ra[y,k,s,1] <- Rb[y,k,s] * phi_Rb_Ra[y]
        
        # expansion to obtain returning hatchery adults
        # natural origin fish only accounted for up until this point
        # could (should) make p_HOR[y] sex-specific - but problems when p_HOR[y,s] == 1?
        Ra[y,k,s,2] <- (Ra[y,k,s,1]/(1 - p_HOR[y])) - Ra[y,k,s,1]
        
        for (o in 1:no) {
          # remove fish for brood stock
          # currently assumed to apply equally to hatchery and wild origin returns
          # and to be age/sex non-selective
          Sb[y,k,s,o] <- Ra[y,k,s,o] * (1 - p_remove[y])
          
          # survive pre-spawn mortality
          # same assumptions as p_remove
          Sa[y,k,s,o] <- Sb[y,k,s,o] * phi_Sb_Sa[y]
        }
      }
    }
    
    # total adults returned to trib
    Ra_tot[y] <- sum(Ra[y,1:nk,1:ns,1:no])
    
    # total spawning adults: used as spawning stock
    # assumes all spawners contribute equally to progeny
    Sa_tot[y] <- sum(Sa[y,1:nk,1:ns,1:no])
    
    # reformat for age comp calculation/fitting
    for (k in 1:nk) {
      # female return by age, natural + hatchery
      Ra_2d[y,k] <- Ra[y,k,1,1] + Ra[y,k,1,2]
      # male return by age, natural + hatchery
      Ra_2d[y,nk+k] <- Ra[y,k,2,1] + Ra[y,k,2,2]
    }
    
    # calculate age/sex composition
    for (ks in 1:nks) {
      q[y,ks] <- Ra_2d[y,ks]/sum(Ra_2d[y,1:nks])
    }
  }
  
}
