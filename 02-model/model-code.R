jags_model_code = function() {
  
  ### PRIORS: RECRUITMENT FUNCTION ###
  # spawner to aggregate parr recruitment function
  log_alpha ~ dnorm(0, 0.001) %_% T(,9.2) # log productivity. bound to prevent nonsensically large draws
  alpha <- exp(log_alpha)
  log_beta ~ dnorm(0, 0.001) %_% T(,15)   # log capacity. bound to prevent nonsensically large draws
  beta <- exp(log_beta)
  sigma_Pb ~ dunif(0, 5)
  
  ### PRIORS: FRESHWATER PARAMETERS ###
  # aggregate parr to LH-specific parr
  mu_pi[i_fall] ~ dbeta(1, 1)
  sig_Lpi ~ dunif(0, 5)
  mu_pi[i_spring] <- 1 - mu_pi[i_fall]
  
  # overwinter survival parameters: LH-specific
  # uses logistic relationship to introduce density-depdendence
  # gamma0: LH-specific intercepts
  # gamma1: LH-common slopes
  for (i in 1:ni) {
    gamma0[i] ~ dnorm(0, 1e-3)
    sig_Lphi_Pa_Mb[i] ~ dunif(0, 5)
  }
  gamma1[i_fall] ~ dnorm(0, 1e-3)
  gamma1[i_spring] <- gamma1[i_fall]

  # natural origin movement survival (trib to LGD): estimate for spring migrants and assume the same value for fall migrants
  mu_phi_Mb_Ma[i_spring,o_nat] ~ dbeta(1, 1)
  sig_Lphi_Mb_Ma[i_spring,o_nat] ~ dunif(0, 5)
  mu_phi_Mb_Ma[i_fall,o_nat] <- mu_phi_Mb_Ma[i_spring,o_nat]
  sig_Lphi_Mb_Ma[i_fall,o_nat] <- sig_Lphi_Mb_Ma[i_spring,o_nat]
  
  # hatchery origin movement survival (trib to LGD): have spring migrants only
  mu_phi_Mb_Ma[i_spring,o_hat] ~ dbeta(1, 1)
  sig_Lphi_Mb_Ma[i_spring,o_hat] ~ dunif(0, 5)
  
  # movement survival (LGD to estuary): same for both LH types, different for origin types
  for (o in 1:no) {
    mu_phi_Ma_M[o] ~ dbeta(1, 1)
    sig_Lphi_Ma_M[o] ~ dunif(0, 5)
  }
  
  # pre-spawn survival (after brood-stock removal to successful spawning)
  mu_phi_Sb_Sa ~ dbeta(1, 1)
  sig_Lphi_Sb_Sa ~ dunif(0, 5)
  
  ### PRIORS: SEX-ASSIGNMENT ###
  for (o in 1:no) {
    mu_omega[s_female,o] ~ dbeta(1, 1)
    mu_omega[s_male,o] <- 1 - mu_omega[s_female,o]
    sig_Lomega[o] ~ dunif(0, 5)
  }
  
  ### PRIORS: MATURATION ###
  # sex/origin-specific maturation probabilities
  for (o in 1:no) {
    for (s in 1:ns) {
      # pr(return at SWA1)
      mu_psi_O1_Rb[s,o] ~ dbeta(1, 1)     
      sig_Lpsi_O1_Rb[s,o] ~ dunif(0, 5)
      
      # pr(return at SWA2|not returned at SWA1)
      mu_psi_O2_Rb[s,o] ~ dbeta(1, 1)
      sig_Lpsi_O2_Rb[s,o] ~ dunif(0, 5)
      
      # pr(return at SWA3|not returned at SWA1 or SWA2)
      mu_psi_O3_Rb[s,o] <- 1
      sig_Lpsi_O3_Rb[s,o] <- 0
    }
  }
  
  ### PRIORS: OCEAN SURVIVAL ###
  # natural origin
  mu_phi_M_O1[o_nat] ~ dbeta(1, 1)             # first winter at sea: to become SWA1
  mu_phi_O1_O2[o_nat] ~ dbeta(1, 1)            # second winter at sea: to become SWA2
  mu_phi_O2_O3[o_nat] <- mu_phi_O1_O2[o_nat]   # third winter at sea: to become SW3
  sig_Lphi_M_O1[o_nat] ~ dunif(0, 5)
  sig_Lphi_O1_O2[o_nat] ~ dunif(0, 5)
  sig_Lphi_O2_O3[o_nat] <- sig_Lphi_O1_O2[o_nat]
  
  # hatchery origin: use a scaler that adjusts natural origin survival to get hatchery survival
  O_phi_scaler_nat_hat ~ dt(0, 1/1.566^2, 7.763)
  logit(mu_phi_M_O1[o_hat]) <- logit(mu_phi_M_O1[o_nat]) + O_phi_scaler_nat_hat
  logit(mu_phi_O1_O2[o_hat]) <- logit(mu_phi_O1_O2[o_nat]) + O_phi_scaler_nat_hat
  logit(mu_phi_O2_O3[o_hat]) <- logit(mu_phi_O1_O2[o_nat]) + O_phi_scaler_nat_hat
  
  ### PRIORS: AR(1) COEFFICIENTS
  kappa_Pb ~ dunif(-0.99,0.99)                # total summer parr recruitment
  kappa_pi ~ dunif(-0.99,0.99)                # proportion of summer parr that are fall migrants
  kappa_phi_Mb_Ma[o_nat] ~ dunif(-0.99,0.99)  # movement survival in spring to LGR (natural origin)
  kappa_phi_Mb_Ma[o_hat] ~ dunif(-0.99,0.99)  # movement survival in spring to LGR (hatchery origin)
  kappa_phi_Ma_M[o_nat] ~ dunif(-0.99,0.99)   # movement survival from LGR thru BON (natural origin)
  kappa_phi_Ma_M[o_hat] ~ dunif(-0.99,0.99)   # movement survival from LGR thru BON (hatchery origin)
  kappa_phi_M_O1 ~ dunif(-0.99,0.99)          # overwinter survival for first winter at sea
  
  ### PRIORS: YEAR-0 RESIDUALS FOR ALL TERMS THAT USE AR(1) PROCESS
  lPb_resid[kmax] ~ dnorm(0, (1/sigma_Pb^2) * (1 - kappa_Pb^2))
  Lpi_resid[kmax] ~ dnorm(0, (1/sig_Lpi^2) * (1 - kappa_pi^2))
  Lphi_Mb_Ma_resid[kmax,o_nat] ~ dnorm(0, (1/sig_Lphi_Mb_Ma[i_spring,o_nat]^2) * (1 - kappa_phi_Mb_Ma[o_nat]^2))
  Lphi_Mb_Ma_resid[kmax,o_hat] ~ dnorm(0, (1/sig_Lphi_Mb_Ma[i_spring,o_hat]^2) * (1 - kappa_phi_Mb_Ma[o_hat]^2))
  Lphi_Ma_M_resid[kmax,o_nat] ~ dnorm(0, (1/sig_Lphi_Ma_M[o_nat]^2) * (1 - kappa_phi_Ma_M[o_nat]^2)) 
  Lphi_Ma_M_resid[kmax,o_hat] ~ dnorm(0, (1/sig_Lphi_Ma_M[o_hat]^2) * (1 - kappa_phi_Ma_M[o_hat]^2))   
  Lphi_M_O1_resid[kmax] ~ dnorm(0, (1/sig_Lphi_M_O1[o_nat]^2) * (1 - kappa_phi_M_O1^2))
  
  ### PRIORS: BROOD-YEAR-SPECIFIC PARAMETERS ###
  for (y in (kmax+1):ny) {
    # aggregate parr to LH-specific parr
    Lpi1[y] ~ dnorm(logit(mu_pi[i_fall]) + Lpi_resid[y-1] * kappa_pi, 1/sig_Lpi^2)
    pi[y,i_fall] <- ilogit(Lpi1[y])
    pi[y,i_spring] <- 1 - pi[y,i_fall]
    
    # overwinter survival: density-dependent
    for (i in 1:ni) {
      Lphi_Pa_Mb[y,i] ~ dnorm(gamma0[i] + gamma1[i] * (Pa[y,i]/peu), 1/sig_Lphi_Pa_Mb[i]^2) 
      phi_Pa_Mb[y,i] <- ilogit(Lphi_Pa_Mb[y,i])
    }
    
    # natural origin movement survival: trib to LGD
    # assume equal between LH types
    Lphi_Mb_Ma[y,i_spring,o_nat] ~ dnorm(logit(mu_phi_Mb_Ma[i_spring,o_nat]) + Lphi_Mb_Ma_resid[y-1,o_nat] * kappa_phi_Mb_Ma[o_nat], 1/sig_Lphi_Mb_Ma[i_spring,o_nat]^2)
    phi_Mb_Ma[y,i_spring,o_nat] <- ilogit(Lphi_Mb_Ma[y,i_spring,o_nat])
    phi_Mb_Ma[y,i_fall,o_nat] <- phi_Mb_Ma[y,i_spring,o_nat]
    
    # hatchery origin movement survival: trib to LGD
    # spring migrants only
    Lphi_Mb_Ma[y,i_spring,o_hat] ~ dnorm(logit(mu_phi_Mb_Ma[i_spring,o_hat]) + Lphi_Mb_Ma_resid[y-1,o_hat] * kappa_phi_Mb_Ma[o_hat], 1/sig_Lphi_Mb_Ma[i_spring,o_hat]^2)
    phi_Mb_Ma[y,i_spring,o_hat] <- ilogit(Lphi_Mb_Ma[y,i_spring,o_hat])
    
    # movement survival: LGD to estuary
    # separate for each origin type
    for (o in 1:no) {
      Lphi_Ma_M[y,o] ~ dnorm(logit(mu_phi_Ma_M[o]) + Lphi_Ma_M_resid[y-1,o] * kappa_phi_Ma_M[o], 1/sig_Lphi_Ma_M[o]^2)
      phi_Ma_M[y,o] <- ilogit(Lphi_Ma_M[y,o])
    }
    
    # probability of returning as female or male by origin
    for (o in 1:no) {
      Lomega1[y,o] ~ dnorm(logit(mu_omega[s_female,o]), 1/sig_Lomega[o]^2)
      omega[y,s_female,o] <- ilogit(Lomega1[y,o])
      omega[y,s_male,o] <- 1 - omega[y,s_female,o]
    }
    
    # sex/origin-specific maturation probabilities
    for (s in 1:ns) {
      for (o in 1:no) {
        # pr(return at SWA1)
        Lpsi_O1_Rb[y,s,o] ~ dnorm(logit(mu_psi_O1_Rb[s,o]), 1/sig_Lpsi_O1_Rb[s,o]^2)
        psi_O1_Rb[y,s,o] <- ilogit(Lpsi_O1_Rb[y,s,o])
        
        # pr(return at SWA2|not returned at SWA1)
        Lpsi_O2_Rb[y,s,o] ~ dnorm(logit(mu_psi_O2_Rb[s,o]), 1/sig_Lpsi_O2_Rb[s,o]^2)
        psi_O2_Rb[y,s,o] <- ilogit(Lpsi_O2_Rb[y,s,o])
        
        # pr(return at SWA3|not returned at SWA1 or SWA2)
        psi_O3_Rb[y,s,o] <- 1
      }
    }
    
    # natural origin ocean survival SWA0 -> SWA1
    Lphi_M_O1[y,o_nat] ~ dnorm(logit(mu_phi_M_O1[o_nat]) + Lphi_M_O1_resid[y-1] * kappa_phi_M_O1, 1/sig_Lphi_M_O1[o_nat]^2)
    phi_M_O1[y,o_nat] <- ilogit(Lphi_M_O1[y,o_nat])
    
    # natural origin ocean survival SWA1 -> SWA2
    Lphi_O1_O2[y,o_nat] ~ dnorm(logit(mu_phi_O1_O2[o_nat]), 1/sig_Lphi_O1_O2[o_nat]^2)
    phi_O1_O2[y,o_nat] <- ilogit(Lphi_O1_O2[y,o_nat])
    
    # natural origin ocean survival SWA2 -> SWA3
    phi_O2_O3[y-1,o_nat] <- phi_O1_O2[y,o_nat]
    
    # hatchery origin ocean survival: use scaler
    logit(phi_M_O1[y,o_hat]) <- logit(phi_M_O1[y,o_nat]) + O_phi_scaler_nat_hat
    logit(phi_O1_O2[y,o_hat]) <- logit(phi_O1_O2[y,o_nat]) + O_phi_scaler_nat_hat
    logit(phi_O2_O3[y-1,o_hat]) <- logit(phi_O2_O3[y-1,o_nat]) + O_phi_scaler_nat_hat
    
    # pre-spawn survival
    Lphi_Sb_Sa[y] ~ dnorm(logit(mu_phi_Sb_Sa), 1/sig_Lphi_Sb_Sa^2)
    phi_Sb_Sa[y] <- ilogit(Lphi_Sb_Sa[y])
  }
  
  # populate the last year/age of ocean survival with the mean across years
  phi_O2_O3[ny,o_nat] <- mu_phi_O2_O3[o_nat]
  phi_O2_O3[ny,o_hat] <- mu_phi_O2_O3[o_hat]
  
  ### PRIORS: HATCHERY STRAYS RETURNING IN YEARS WITH NO ASSOCIATED SMOLT RELEASE ###
  # the number of hatchery strays: only estimate in years where no other mechanism for generating hatchery fish
  for (i in 1:n_stray_yrs) {
    n_stray_tot[stray_yrs[i]] ~ dunif(0, 5000)
  }
  for (i in 1:n_not_stray_yrs) {
    n_stray_tot[not_stray_yrs[i]] <- 0
  }
  
  # when n_stray_tot > 0, what age/sex are they? (assume all hatchery fish)
  stray_comp_2d[1:nks] ~ ddirich(rep(1, nks))
  
  # put stray_comp_2d into array format for looping in process model
  for (k in 1:nk) {
    for (s in 1:ns) {
      # assume no strays for natural fish
      stray_comp[k,s,o_nat] <- 0
      
      # place stray_comp_2d in the right index locations
      stray_comp[k,s,o_hat] <- stray_comp_2d[k + nk * (s-1)]
    }
  }
  
  ### PRIORS: OBSERVATION MODEL ###
  
  # carcass vs. weir composition correction factor coefficients
  # for (i in 1:3) {
  #   z[i] ~ dunif(-10,10)
  # }
  
  # if fitting MIN, use this instead of the priors above
  # params not estimable for MIN alone
  z[1] <- -0.76
  z[2] <- -0.01
  z[3] <- -0.28
  
  # if fitting MIN, use this instead of the priors above
  # params not estimable for MIN alone
  # z[1] <- -0.76
  # z[2] <- -0.01
  # z[3] <- -0.28
  
  # calculate correction factor
  for (k in 1:3) {
    for (s in 1:2) {
      log(carc_adj[k,s]) <- z[1] * age3[k] + z[2] * age5[k] + z[3] * male[s]
    }
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
  # natural origin fish only
  for (s in 1:ns) {
    p_init_prime[1,s] <- mu_phi_M_O1[o_nat] * mu_psi_O1_Rb[s,o_nat]
    p_init_prime[2,s] <- mu_phi_M_O1[o_nat] * (1 - mu_psi_O1_Rb[s,o_nat]) * mu_phi_O1_O2[o_nat] * mu_psi_O2_Rb[s,o_nat]
    p_init_prime[3,s] <- mu_phi_M_O1[o_nat] * (1 - mu_psi_O1_Rb[s,o_nat]) * mu_phi_O1_O2[o_nat] * (1 - mu_psi_O2_Rb[s,o_nat]) * mu_phi_O2_O3[o_nat] * mu_psi_O3_Rb[s,o_nat]
    for (k in 1:nk) {
      p_init[k,s] <- p_init_prime[k,s]/sum(p_init_prime[1:nk,s])
    }
  }
  
  # hyperparameters for initialization
  mu_init_recruits ~ dunif(0, max_init_recruits)
  sig_init_lrecruits ~ dunif(0,10)
  
  # for natural origin
  for (y in 1:kmax) {
    # random, constant-mean adult recruits for first kmax brood years
    init_recruits[y] ~ dlnorm(log(mu_init_recruits), 1/sig_init_lrecruits^2) %_% T(,max_init_recruits)
    for (s in 1:ns) {
      # apportion them to each sex
      init_recruits_sex[y,s] <- init_recruits[y] * mu_omega[s,o_nat]
      for (k in 1:nk) {
        # apportion them to return year, age, and sex
        Rb[y+kmin+k-1,k,s,o_nat] <- init_recruits_sex[y,s] * p_init[k,s] 
      }
    }
  }
  
  # for hatchery origin
  for (y in 1:kmax) {
    for (s in 1:ns) {
      for (k in 1:nk) {
        # there were no hatchery operations in these years, but still need to populate with a number
        Rb[y+kmin+k-1,k,s,o_hat] <- 0 
      }
    }
  }
  
  ### PROCESS MODEL: COMPLETE LIFE CYCLE FOR REMAINING BROOD YEARS ###
  for (y in (kmax+1):ny) {
    
    # reproductive link: total summer parr
    Pb_pred[y] <- Sa_tot[y]/(1/alpha + Sa_tot[y]/beta)
    Pb[y] ~ dlnorm(log(Pb_pred[y]) + lPb_resid[y-1] * kappa_Pb, 1/sigma_Pb^2)
    
    # natural origin tributary-to-LGD dynamics
    for (i in 1:ni) {
      # apportion to LH-strategy: parr after fall migration
      Pa[y,i] <- Pb[y] * pi[y,i]

      # survive over winter: smolt before spring migration
      Mb[y,i,o_nat] <- Pa[y,i] * phi_Pa_Mb[y,i]

      # move to LGD: smolt after spring migration, at top of LGD
      Ma[y,i,o_nat] <- Mb[y,i,o_nat] * phi_Mb_Ma[y,i,o_nat]
      
      # derived survival for fitting: fall trap to LGD
      phi_Pa_Ma[y,i] <- Ma[y,i,o_nat]/Pa[y,i]
    }
    
    # derived survival for fitting: summer tagging to LGD
    phi_Pb_Ma[y] <- sum(Ma[y,1:ni,o_nat])/Pb[y]
    
    # put hatchery smolts in tributary
    Mb[y,i_spring,o_hat] <- Mb_obs[y,i_spring,o_hat]
    
    # move hatchery smolts from tributary to LGD
    Ma[y,i_spring,o_hat] <- Mb[y,i_spring,o_hat] * phi_Mb_Ma[y,i_spring,o_hat]
    
    # create zeros for fall migrant hatchery fish at LGD
    # needed because we sum over this dimension below
    Ma[y,i_fall,o_hat] <- 0
    
    # sex/origin-specific processes
    for (s in 1:ns) {
      for (o in 1:no) {
        # move origin-specific smolts from LGD to estuary and assign to sex
        M[y,s,o] <- sum(Ma[y,1:ni,o]) * phi_Ma_M[y,o] * omega[y,s,o]
        
        # move juveniles through ocean ages and survivals
        O[y,1,s,o] <- M[y,s,o] * phi_M_O1[y,o] # survive first winter at sea. now SWA1, TA3
        O[y,2,s,o] <- O[y,1,s,o] * (1 - psi_O1_Rb[y,s,o]) * phi_O1_O2[y,o] # don't mature at SWA1 and survive second winter at sea. now SWA2, TA4
        O[y,3,s,o] <- O[y,2,s,o] * (1 - psi_O2_Rb[y,s,o]) * phi_O2_O3[y,o] # don't mature at SWA2 and survive third winter at sea. now SWA3, TA5
        
        # mature and return to river in appropriate year at age/sex
        Rb[y+kmin+1-1,1,s,o] <- O[y,1,s,o] * psi_O1_Rb[y,s,o]
        Rb[y+kmin+2-1,2,s,o] <- O[y,2,s,o] * psi_O2_Rb[y,s,o]
        Rb[y+kmin+3-1,3,s,o] <- O[y,3,s,o] * psi_O3_Rb[y,s,o]
      }
    }
    
    # adult in-river processes
    # y now represents brood year these fish are returning in
    for (k in 1:nk) {
      for (s in 1:ns) {
        for (o in 1:no) {
          
          # survive upstream migration (survive sea lions * survive fisheries * survive through dams and to tributary) and add strays
          Ra[y,k,s,o] <- Rb[y,k,s,o] * phi_SL[y] * (1 - Ub[y,k,o]) * (1 - Ua[y,k,o]) * phi_D1_D4 * phi_D4_D5 * phi_D5_D8 * phi_D8_Ra + n_stray_tot[y] * stray_comp[k,s,o]
          
          # remove fish for brood stock
          Sb[y,k,s,o] <- max(Ra[y,k,s,o] - n_remove[y,k,s,o], 0)

          # survive pre-spawn mortality
          # currently assumed to apply equally to hatchery and wild origin returns
          # and to be age/sex non-selective
          Sa[y,k,s,o] <- Sb[y,k,s,o] * phi_Sb_Sa[y]
          
          # calculate "adjusted carcasses": accounts for sampling bias relative to weir
          Sa_adj[y,k,s,o] <- Sa[y,k,s,o] * carc_adj[k,s]
        }
      }
    }
    
    # total adults returned to trib
    Ra_tot[y] <- sum(Ra[y,1:nk,1:ns,1:no])
    
    # total spawning adults: used as spawning stock
    # assumes all spawners contribute equally to progeny
    Sa_tot[y] <- sum(Sa[y,1:nk,1:ns,1:no])
    
    # reformat returner by age/sex/origin for fitting to weir comp data
    Ra_2d[y,1:nk] <- Ra[y,1:nk,s_female,o_nat]            # nat. females, all ages
    Ra_2d[y,(nk+1):(2*nk)] <- Ra[y,1:nk,s_male,o_nat]     # nat. males, all ages
    Ra_2d[y,(2*nk+1):(3*nk)] <- Ra[y,1:nk,s_female,o_hat] # hat. females, all ages
    Ra_2d[y,(3*nk+1):(4*nk)] <- Ra[y,1:nk,s_male,o_hat]   # hat. males, all ages
    
    # reformat "adjusted carcasses" by age/sex/origin for fitting to carcass comp data
    Sa_adj_2d[y,1:nk] <- Sa_adj[y,1:nk,s_female,o_nat]            # nat. females, all ages
    Sa_adj_2d[y,(nk+1):(2*nk)] <- Sa_adj[y,1:nk,s_male,o_nat]     # nat. males, all ages
    Sa_adj_2d[y,(2*nk+1):(3*nk)] <- Sa_adj[y,1:nk,s_female,o_hat] # hat. females, all ages
    Sa_adj_2d[y,(3*nk+1):(4*nk)] <- Sa_adj[y,1:nk,s_male,o_hat]   # hat. males, all ages

    # calculate age/sex compositions
    for (kso in 1:nkso) {
      q_Ra[y,kso] <- Ra_2d[y,kso]/sum(Ra_2d[y,1:nkso])
      q_Sa_adj[y,kso] <- Sa_adj_2d[y,kso]/sum(Sa_adj_2d[y,1:nkso])
    }
  }

  ### OBSERVATION MODEL ###
  
  # age/sex composition
  for (y in (kmax+1):ny) {
    # at weir
    weir_x_obs[y,1:nkso] ~ dmulti(q_Ra[y,1:nkso], weir_nx_obs[y])
    
    # for carcasses
    carc_x_obs[y,1:nkso] ~ dmulti(q_Sa_adj[y,1:nkso], carc_nx_obs[y])
  }
  
  # the "fit" objects allow using square objects but skipping over missing obs in likelihoods
  # fall trap abundance
  for (d in 1:nfit_Pa) {
    Pa_obs[fit_Pa[d,1],fit_Pa[d,2]] ~ dlnorm(log(Pa[fit_Pa[d,1],fit_Pa[d,2]]), 1/sig_Pa_obs[fit_Pa[d,1],fit_Pa[d,2]]^2)
  }
  
  # spring trap abundance
  for (d in 1:nfit_Mb) {
    Mb_obs[fit_Mb[d,1],fit_Mb[d,2],fit_Mb[d,3]] ~ dlnorm(log(Mb[fit_Mb[d,1],fit_Mb[d,2],fit_Mb[d,3]]), 1/sig_Mb_obs[fit_Mb[d,1],fit_Mb[d,2],fit_Mb[d,3]]^2)
  }
  
  # adult abundance
  for (d in 1:nfit_Ra) {
    Ra_obs[fit_Ra[d]] ~ dlnorm(log(Ra_tot[fit_Ra[d]]), 1/sig_Ra_obs[fit_Ra[d]]^2)
  }
  
  # summer tagging to LGD
  for (d in 1:nfit_Lphi_Pb_Ma) {
    Lphi_obs_Pb_Ma[fit_Lphi_Pb_Ma[d]] ~ dnorm(logit(phi_Pb_Ma[fit_Lphi_Pb_Ma[d]]), 1/sig_Lphi_obs_Pb_Ma[fit_Lphi_Pb_Ma[d]]^2)
  }
  
  # fall tagging to LGD
  for (d in 1:nfit_Lphi_Pa_Ma) {
    Lphi_obs_Pa_Ma[fit_Lphi_Pa_Ma[d,1],fit_Lphi_Pa_Ma[d,2]] ~ dnorm(logit(phi_Pa_Ma[fit_Lphi_Pa_Ma[d,1],fit_Lphi_Pa_Ma[d,2]]), 1/sig_Lphi_obs_Pa_Ma[fit_Lphi_Pa_Ma[d,1],fit_Lphi_Pa_Ma[d,2]]^2)
  }
  
  # spring tagging/smolt releases to LGD
  for (d in 1:nfit_Lphi_Mb_Ma) {
    Lphi_obs_Mb_Ma[fit_Lphi_Mb_Ma[d,1],fit_Lphi_Mb_Ma[d,2],fit_Lphi_Mb_Ma[d,3]] ~ dnorm(logit(phi_Mb_Ma[fit_Lphi_Mb_Ma[d,1],fit_Lphi_Mb_Ma[d,2],fit_Lphi_Mb_Ma[d,3]]), 1/sig_Lphi_obs_Mb_Ma[fit_Lphi_Mb_Ma[d,1],fit_Lphi_Mb_Ma[d,2],fit_Lphi_Mb_Ma[d,3]]^2)
  }
  
  # hydrosystem survival
  for (d in 1:nfit_Lphi_Ma_M) {
    Lphi_obs_Ma_M[fit_Lphi_Ma_M[d,1],fit_Lphi_Ma_M[d,2]] ~ dnorm(logit(phi_Ma_M[fit_Lphi_Ma_M[d,1],fit_Lphi_Ma_M[d,2]]), 1/sig_Lphi_obs_Ma_M[fit_Lphi_Ma_M[d,1],fit_Lphi_Ma_M[d,2]]^2)
  }
  
  # pre-spawn survival
  for (d in 1:nfit_spawned) {
    carcs_spawned[fit_spawned[d]] ~ dbin(phi_Sb_Sa[fit_spawned[d]], carcs_sampled[fit_spawned[d]])
  }
  
  ### CALCULATE ALL PROCESS MODEL RESIDUALS ###
  
  for (y in (kmax+1):ny) {
    
    # total summer parr recruitment
    lPb_resid[y] <- log(Pb[y]) - log(Pb_pred[y])
    
    # proportion of summer parr that become fall migrants
    Lpi_resid[y] <- logit(pi[y,i_fall]) - logit(mu_pi[i_fall])
    
    # LH-specific overwinter survival
    for (i in 1:ni) {
      Lphi_Pa_Mb_resid[y,i] <- Lphi_Pa_Mb[y,i] - (gamma0[i] + gamma1[i] * (Pa[y,i]/peu))
    }
    
    # origin-specific quantities
    for (o in 1:no) {
      # movement survival to LGR
      Lphi_Mb_Ma_resid[y,o] <- Lphi_Mb_Ma[y,i_spring,o] - logit(mu_phi_Mb_Ma[i_spring,o])
      
      # movement survival from LGR thru BON
      Lphi_Ma_M_resid[y,o] <- Lphi_Ma_M[y,o] - logit(mu_phi_Ma_M[o])
      
      # probability of returning as female
      Lomega_resid[y,o] <- Lomega1[y,o] - logit(mu_omega[s_female,o])
      
      # sex-specific quantities: maturity
      for (s in 1:ns) {
        # pr(return at SWA1)
        Lpsi_O1_Rb_resid[y,s,o] <- Lpsi_O1_Rb[y,s,o] - logit(mu_psi_O1_Rb[s,o])
        
        # pr(return at SWA2|not returned at SWA1)
        Lpsi_O2_Rb_resid[y,s,o] <- Lpsi_O2_Rb[y,s,o] - logit(mu_psi_O2_Rb[s,o])
      }
    }
    
    # natural origin ocean survival SWA0 -> SWA1
    Lphi_M_O1_resid[y] <- Lphi_M_O1[y,o_nat] - logit(mu_phi_M_O1[o_nat])
    
    # natural origin ocean survival SWA1 -> SWA2
    Lphi_O1_O2_resid[y] <- Lphi_O1_O2[y,o_nat] - logit(mu_phi_O1_O2[o_nat])
    
    # pre-spawn survival
    Lphi_Sb_Sa_resid[y] <- Lphi_Sb_Sa[y] - logit(mu_phi_Sb_Sa)
  }
}
