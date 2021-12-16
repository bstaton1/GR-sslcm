jags_model_code = function() {
  
  ### PRIORS: capacity vs. weighted usable habitat length relationship ###
  mu_beta_per_wul ~ dnorm(0, 1e-8) %_% T(0,)
  sig_lbeta ~ dunif(0, 5)
  
  for (j in 1:nj) {
    ### PRIORS: RECRUITMENT FUNCTION ###
    # total egg production to aggregate parr recruitment function
    alpha[j] ~ dbeta(1, 1)
    sig_Pb[j] ~ dunif(0, 5)
    lbeta[j] ~ dnorm(log(mu_beta_per_wul * wul[j]), 1/sig_lbeta^2) %_% T(,15) # log capacity. bound to prevent nonsensically large draws
    beta[j] <- exp(lbeta[j])
    beta_per_wul[j] <- beta[j]/wul[j]

    ### PRIORS: FRESHWATER PARAMETERS ###
    # aggregate parr to LH-specific parr
    mu_pi[i_fall,j] ~ dbeta(1, 1)
    sig_Lpi[j] ~ dunif(0, 5)
    mu_pi[i_spring,j] <- 1 - mu_pi[i_fall,j]
    
    # overwinter survival parameters: LH-specific
    # uses logistic relationship to introduce density-dependence
    # gamma0: LH-specific intercepts
    # gamma1: LH-common slopes
    for (i in 1:ni) {
      gamma0[i,j] ~ dnorm(0, 1e-3)
      sig_Lphi_Pa_Mb[i,j] ~ dunif(0, 5)
    }
    gamma1[i_fall,j] ~ dnorm(0, 1e-3)
    gamma1[i_spring,j] <- gamma1[i_fall,j]
    
    # natural origin movement survival (trib to LGD): estimate for spring migrants and assume the same value for fall migrants
    mu_phi_Mb_Ma[i_spring,o_nat,j] ~ dbeta(1, 1)
    sig_Lphi_Mb_Ma[o_nat,j] ~ dunif(0, 5)
    mu_phi_Mb_Ma[i_fall,o_nat,j] <- mu_phi_Mb_Ma[i_spring,o_nat,j]

    # hatchery origin movement survival (trib to LGD): have spring migrants only
    mu_phi_Mb_Ma[i_spring,o_hat,j] ~ dbeta(1, 1)
    sig_Lphi_Mb_Ma[o_hat,j] ~ dunif(0, 5)
    
    # pre-spawn survival (after brood-stock removal to successful spawning)
    mu_phi_Sb_Sa[j] ~ dbeta(1, 1)
    sig_Lphi_Sb_Sa[j] ~ dunif(0, 5)
    
    ### PRIORS: MATURATION ###
    # age/origin-specific maturation probabilities
    for (o in 1:no) {
      
      # pr(return at SWA1)
      mu_psi_O1_Rb[o,j] ~ dbeta(1, 1)
      sig_Lpsi_O1_Rb[o,j] ~ dunif(0, 5)
      
      # pr(return at SWA2|not returned at SWA1)
      mu_psi_O2_Rb[o,j] ~ dbeta(1, 1)
      sig_Lpsi_O2_Rb[o,j] ~ dunif(0, 5)
    }
    
    ### PRIORS: OCEAN SURVIVAL ###
    # mean survival by ocean year transition for natural origin
    mu_phi_O0_O1[o_nat,j] ~ dbeta(1,1)               # first winter at sea: to become SWA1
    mu_phi_O1_O2[o_nat,j] ~ dbeta(60,15)             # second winter at sea: to become SWA2
    mu_phi_O2_O3[o_nat,j] <- mu_phi_O1_O2[o_nat,j]   # third winter at sea: to become SW3
    
    # log odds ratio between natural and hatchery origin
    O_phi_scaler_nat_hat[j] ~ dt(0, 1/1.566^2, 7.763)
    
    # AR(1) coefficient for first year ocean survival
    kappa_phi_O0_O1[j] ~ dunif(-0.99,0.99)
    
    # standard deviation of white noise residuals
    sig_Lphi_O0_O1[j] ~ dunif(0, 5)
    
    # standard deviation of total year 1 residual (i.e., includes autocorrelation)
    sig_Lphi_O0_O1_init[j] <- sqrt(sig_Lphi_O0_O1[j]^2/(1 - kappa_phi_O0_O1[j]^2))
    
    # mean survival by ocean year transition for hatchery origin
    logit(mu_phi_O0_O1[o_hat,j]) <- logit(mu_phi_O0_O1[o_nat,j]) + O_phi_scaler_nat_hat[j]
    logit(mu_phi_O1_O2[o_hat,j]) <- logit(mu_phi_O1_O2[o_nat,j]) + O_phi_scaler_nat_hat[j]
    logit(mu_phi_O2_O3[o_hat,j]) <- logit(mu_phi_O2_O3[o_nat,j]) + O_phi_scaler_nat_hat[j]
    
    ### PRIORS: HATCHERY STRAYS RETURNING IN YEARS WITH NO ASSOCIATED SMOLT RELEASE ###
    # the number of hatchery strays: only estimate in years where no other mechanism for generating hatchery fish
    for (i in 1:n_stray_yrs[j]) {
      n_stray_tot[stray_yrs[i,j],j] ~ dunif(0, 5000)
    }
    for (i in 1:n_not_stray_yrs[j]) {
      n_stray_tot[not_stray_yrs[i,j],j] <- 0
    }
    
    # when n_stray_tot > 0, what age are they? (assume all hatchery fish)
    stray_comp_2d[1:nk,j] ~ ddirich(rep(1, nk))
    
    # put stray_comp_2d into array format for looping in process model
    for (k in 1:nk) {
      # assume no strays for natural fish
      stray_comp[k,o_nat,j] <- 0
      
      # place stray_comp_2d in the right index locations
      stray_comp[k,o_hat,j] <- stray_comp_2d[k,j]
    }
  }
  
  # survival between summer tagging/LH apportionment and winter tagging for spring migrants
  # parameter is needed to allow inclusion of winter tagging data, this is not a key pop dyn parameter
  phi_Pb_Pa[i_fall,1:nj] <- rep(1, nj)  # doesn't apply to fall migrants, but loops over i dimension
  phi_Pb_Pa[i_spring,1] ~ dunif(0, 1)
  phi_Pb_Pa[i_spring,2] ~ dunif(0, 1)
  phi_Pb_Pa[i_spring,3] <- 1            # don't estimate for MIN -- no winter tagging data
  phi_Pb_Pa[i_spring,4] ~ dunif(0, 1)
  
  ### PRIORS: PARAMETERS THAT ARE CONSTANT ACROSS POPULATIONS ###
  
  # juvenile movement survival (LGR to estuary [ocean-age 0]): same for both LH types, different for origin types
  for (o in 1:no) {
    mu_phi_Ma_O0[o] ~ dbeta(1, 1)
    sig_Lphi_Ma_O0[o] ~ dunif(0, 5)
  }
  
  # correlation parameters
  # see toggle_rho_estimation() function; 01-functions/util-fns.R
  rho_lPb <- 0
  rho_Lpi <- 0
  rho_Lphi_Pa_Mb[i_fall] <- 0
  rho_Lphi_Pa_Mb[i_spring] <- 0
  rho_Lphi_Mb_Ma[o_nat] <- 0
  rho_Lphi_Mb_Ma[o_hat] <- 0
  rho_Lphi_Ma_O0 <- 0
  rho_Lphi_O0_O1 <- 0
  rho_Lpsi_O1_Rb[o_nat] <- 0
  rho_Lpsi_O1_Rb[o_hat] <- 0
  rho_Lpsi_O2_Rb[o_nat] <- 0
  rho_Lpsi_O2_Rb[o_hat] <- 0
  rho_Lphi_Rb_Ra <- 0
  rho_Lphi_Sb_Sa <- 0
  
  # construct all population covariance matrices
  for (i in 1:nj) {
    for (j in 1:nj) {
      # covariance matrix of parr recruitment
      Sig_lPb[i,j] <- sig_Pb[i] * sig_Pb[j] * ifelse(i == j, 1, rho_lPb)
      
      # covariance matrix of LH apportionment
      Sig_Lpi[i,j] <- sig_Lpi[i] * sig_Lpi[j] * ifelse(i == j, 1, rho_Lpi)
      
      # covariance matrices of LH-specific overwinter survival (g used for LH index)
      for (g in 1:ni) {
        Sig_Lphi_Pa_Mb[i,j,g] <- sig_Lphi_Pa_Mb[g,i] * sig_Lphi_Pa_Mb[g,j] * ifelse(i == j, 1, rho_Lphi_Pa_Mb[g])
      }
      
      for (o in 1:no) {
        # covariance matrices of movement survival to LGR
        Sig_Lphi_Mb_Ma[i,j,o] <- sig_Lphi_Mb_Ma[o,i] * sig_Lphi_Mb_Ma[o,j] * ifelse(i == j, 1, rho_Lphi_Mb_Ma[o])
        
        # covariance matrices of pr(mature at SWA1)
        Sig_Lpsi_O1_Rb[i,j,o] <- sig_Lpsi_O1_Rb[o,i] * sig_Lpsi_O1_Rb[o,j] * ifelse(i == j, 1, rho_Lpsi_O1_Rb[o])
        
        # covariance matrices of pr(mature at SWA1)
        Sig_Lpsi_O2_Rb[i,j,o] <- sig_Lpsi_O2_Rb[o,i] * sig_Lpsi_O2_Rb[o,j] * ifelse(i == j, 1, rho_Lpsi_O2_Rb[o])
      }
      
      # covariance matrix of white noise portion of yr1 ocean survival
      Sig_Lphi_O0_O1[i,j] <- sig_Lphi_O0_O1[i] * sig_Lphi_O0_O1[j] * ifelse(i == j, 1, rho_Lphi_O0_O1)
      
      # covariance matrix of total process noise in yr1 ocean survival
      Sig_Lphi_O0_O1_init[i,j] <- sig_Lphi_O0_O1_init[i] * sig_Lphi_O0_O1_init[j] * ifelse(i == j, 1, rho_Lphi_O0_O1)
      
      # covariance matrix of pre-spawn survival
      Sig_Lphi_Sb_Sa[i,j] <- sig_Lphi_Sb_Sa[i] * sig_Lphi_Sb_Sa[j] * ifelse(i == j, 1, rho_Lphi_Sb_Sa)
    }
  }
  
  # construct all origin covariance matrices
  for (i in 1:no) {
    for (j in 1:no) {
      # movement survival (juveniles LGR to BON)
      Sig_Lphi_Ma_O0[i,j] <- sig_Lphi_Ma_O0[i] * sig_Lphi_Ma_O0[j] * ifelse(i == j, 1, rho_Lphi_Ma_O0)
      
      # movement survival (adults BON to LGR)
      Sig_Lphi_Rb_Ra[i,j] <- sig_Lphi_Rb_Ra[i] * sig_Lphi_Rb_Ra[j] * ifelse(i == j, 1, rho_Lphi_Rb_Ra)
    }
  }
  
  # year 0 residuals for yr1 ocean survival (needed for AR(1) process)
  Lphi_O0_O1_resid[kmax,o_nat,1:nj] ~ dmnorm.vcov(rep(0, nj), Sig_Lphi_O0_O1_init[1:nj,1:nj])

  # migration survival adults from BON to LGR
  for (o in 1:no) {
    mu_phi_Rb_Ra[o] ~ dbeta(1, 1)
    sig_Lphi_Rb_Ra[o] ~ dunif(0, 5)
  }
  
  # carcass vs. weir composition correction factor coefficients
  # age/population-specific coefficients of correction factor
  # only for populations with both weir and carcass data (j_z)
  for (i in 1:2) {
    # hyper-params
    mu_z[i] ~ dunif(-10,10)
    sig_z[i] ~ dunif(0, 5)
    
    # pop-specific effects
    for (j in 1:nj_z) {
      z[i,j_z[j]] ~ dnorm(mu_z[i], 1/sig_z[i]^2)
    }
  }
  
  # calculate correction factor for populations with both weir and carcass data
  for (j in 1:nj_z) {
    for (k in 1:3) {
      log(carc_adj[k,j_z[j]]) <- z[1,j_z[j]] * age3[k] + z[2,j_z[j]] * age5[k]
    }
  }
  
  # calculate correction factor for population(s) without both weir and carcass data (Minam)
  carc_adj[1,3] <- (carc_adj[1,1] + carc_adj[1,2] + carc_adj[1,4])/3
  carc_adj[2,3] <- (carc_adj[2,1] + carc_adj[2,2] + carc_adj[2,4])/3
  carc_adj[3,3] <- (carc_adj[3,1] + carc_adj[3,2] + carc_adj[3,4])/3
  
  ### PRIORS: BROOD-YEAR-SPECIFIC PARAMETERS ###
  for (y in (kmax+1):ny) {
    
    # parr recruitment
    lPb[y,1:nj] ~ dmnorm.vcov(log(Pb_mean[y,1:nj]), Sig_lPb[1:nj,1:nj])

    # LH apportionment
    Lpi1[y,1:nj] ~ dmnorm.vcov(logit(mu_pi[i_fall,1:nj]), Sig_Lpi[1:nj,1:nj])

    # overwinter survival by LH type: density dependent
    Lphi_Pa_Mb[y,i_fall,1:nj] ~ dmnorm.vcov(gamma0[i_fall,1:nj] + gamma1[i_fall,1:nj] * (Pa[y,i_fall,1:nj]/wul[1:nj]), Sig_Lphi_Pa_Mb[1:nj,1:nj,i_fall])
    Lphi_Pa_Mb[y,i_spring,1:nj] ~ dmnorm.vcov(gamma0[i_spring,1:nj] + gamma1[i_spring,1:nj] * (Pa[y,i_spring,1:nj]/wul[1:nj]), Sig_Lphi_Pa_Mb[1:nj,1:nj,i_spring])

    for (o in 1:no) {
      # migration survival from trib to LGR by origin
      Lphi_Mb_Ma[y,i_spring,o,1:nj] ~ dmnorm.vcov(logit(mu_phi_Mb_Ma[i_spring,o,1:nj]), Sig_Lphi_Mb_Ma[1:nj,1:nj,o])

      # pr(return at SWA1)
      Lpsi_O1_Rb[y,o,1:nj] ~ dmnorm.vcov(logit(mu_psi_O1_Rb[o,1:nj]), Sig_Lpsi_O1_Rb[1:nj,1:nj,o])

      # pr(return at SWA2)
      Lpsi_O2_Rb[y,o,1:nj] ~ dmnorm.vcov(logit(mu_psi_O2_Rb[o,1:nj]), Sig_Lpsi_O2_Rb[1:nj,1:nj,o])
    }
    
    # yr1 NOR ocean survival: includes AR(1) process
    Lphi_O0_O1[y,o_nat,1:nj] ~ dmnorm.vcov(logit(mu_phi_O0_O1[o_nat,1:nj]) + Lphi_O0_O1_resid[y-1,o_nat,1:nj] * kappa_phi_O0_O1[1:nj], Sig_Lphi_O0_O1[1:nj,1:nj])

    # yr2/yr3 NOR ocean survival: time constant
    Lphi_O1_O2[y,o_nat,1:nj] <- logit(mu_phi_O1_O2[o_nat,1:nj])
    Lphi_O2_O3[y,o_nat,1:nj] <- logit(mu_phi_O2_O3[o_nat,1:nj])
    
    # yr1/yr2/yr3 HOR ocean survival: same as NOR but adjusted by a time-constant log odds ratio
    Lphi_O0_O1[y,o_hat,1:nj] <- Lphi_O0_O1[y,o_nat,1:nj] + O_phi_scaler_nat_hat[1:nj]
    Lphi_O1_O2[y,o_hat,1:nj] <- Lphi_O1_O2[y,o_nat,1:nj] + O_phi_scaler_nat_hat[1:nj]
    Lphi_O2_O3[y,o_hat,1:nj] <- Lphi_O2_O3[y,o_nat,1:nj] + O_phi_scaler_nat_hat[1:nj]
    
    # pre-spawn survival
    Lphi_Sb_Sa[y,1:nj] ~ dmnorm.vcov(logit(mu_phi_Sb_Sa[1:nj]), Sig_Lphi_Sb_Sa[1:nj,1:nj])

    # movement survival (juveniles LGR to BON)
    Lphi_Ma_O0[y,1:no] ~ dmnorm.vcov(logit(mu_phi_Ma_O0[1:no]), Sig_Lphi_Ma_O0[1:no,1:no])

    # movement survival (adults BON to LGR)
    Lphi_Rb_Ra_random[y,1:no] ~ dmnorm.vcov(logit(mu_phi_Rb_Ra[1:no]), Sig_Lphi_Rb_Ra[1:no,1:no])

    # apply inverse logit-transformation to all of these brood year-specific parameters
    # ilogit() and exp() cannot be applied to vectors unfortunately; must loop over population dimension
    for (j in 1:nj) {
      
      # transform parr recruitment
      Pb[y,j] <- exp(lPb[y,j])
      
      # transform LH apportionment
      pi[y,i_fall,j] <- ilogit(Lpi1[y,j])
      pi[y,i_spring,j] <- 1 - pi[y,i_fall,j]
      
      # transform overwinter survival
      phi_Pa_Mb[y,i_fall,j] <- ilogit(Lphi_Pa_Mb[y,i_fall,j])
      phi_Pa_Mb[y,i_spring,j] <- ilogit(Lphi_Pa_Mb[y,i_spring,j])

      for (o in 1:no) {
        # transform movement survival trib to LGR
        phi_Mb_Ma[y,i_spring,o,j] <- ilogit(Lphi_Mb_Ma[y,i_spring,o,j])
        
        # transform pr(return at SWA1)
        psi_O1_Rb[y,o,j] <- ilogit(Lpsi_O1_Rb[y,o,j])
        
        # transform pr(return at SWA2|not returned at SWA1)
        psi_O2_Rb[y,o,j] <- ilogit(Lpsi_O2_Rb[y,o,j])
        
        # transform ocean survival terms
        phi_O0_O1[y,o,j] <- ilogit(Lphi_O0_O1[y,o,j])
        phi_O1_O2[y,o,j] <- ilogit(Lphi_O1_O2[y,o,j])
        phi_O2_O3[y,o,j] <- ilogit(Lphi_O2_O3[y,o,j])
      }
      
      # assume movement survival trib to LGR for NOR fish is equal between LH types
      phi_Mb_Ma[y,i_fall,o_nat,j] <- phi_Mb_Ma[y,i_spring,o_nat,j]
      
      # pre-spawn survival
      phi_Sb_Sa[y,j] <- ilogit(Lphi_Sb_Sa[y,j])
    }
    
    # transform quantities common to all pops but different by origin
    for (o in 1:no) {
      # transform movement survival (juveniles LGR to BON)
      phi_Ma_O0[y,o] <- ilogit(Lphi_Ma_O0[y,o])
      
      # transform movement survival (adults BON to LGR)
      Lphi_Rb_Ra[y,o] <- ifelse(y < first_LGR_adults, logit(mu_phi_Rb_Ra[o]), Lphi_Rb_Ra_random[y,o])
      phi_Rb_Ra[y,o] <- ilogit(Lphi_Rb_Ra[y,o])
    }
  }
  
  ### PROCESS MODEL: INITIALIZATION ###
  # need to do something special with first kmax brood years prior to data collection
  # to populate adult states for reproduction/observation in the first years with data
  # THIS APPROACH:
  #  Estimate the mean and year-specific TOTAL ADULT RECRUITS (arriving to mouth of river)
  #  for the first kmax brood years. Then, based on the hyperparams of return-by-age and ocean survival
  #  (this is the p_init_prime and p_init stuff below)
  #  derive the probability of returning at age will apply equally to all kmax of the first brood years
  #  Then place these adult recruits in the appropriate return year by age
  
  # obtain expected probability of adult recruits returning at age
  # natural origin fish only
  for (j in 1:nj) {
    for (k in 1:nk) {
      for (y in (kmax+1):(kmax+kmin+k-1)) {
        # Rb[y,k,o_nat,j] <- Rb_init[y,k,j] 
        Rb[y,k,o_nat,j] ~ dunif(0, max_Rb_init[k])
        Rb[y,k,o_hat,j] <- 0
      }
    }
  }
  
  ### PROCESS MODEL: COMPLETE LIFE CYCLE FOR REMAINING BROOD YEARS ###
  for (j in 1:nj) {
    for (y in (kmax+1):ny) {
      # reproductive link: total summer parr
      # Pb[y,j] <- f_tot[y,j]/(1/alpha[j] + f_tot[y,j]/beta[j]) * exp(lPb_resid[y,j])
      Pb_mean[y,j] <- f_tot[y,j]/(1/alpha[j] + f_tot[y,j]/beta[j])
      # Pb[y,j] ~ dlnorm(log(Pb_pred[y,j]), 1/sig_Pb[j]^2)
      
      # natural origin tributary-to-LGD dynamics
      for (i in 1:ni) {
        # apportion to LH-strategy: parr after fall migration
        Pa[y,i,j] <- Pb[y,j] * pi[y,i,j]
        
        # survive over winter: smolt before spring migration
        Mb[y,i,o_nat,j] <- Pa[y,i,j] * phi_Pa_Mb[y,i,j]
        
        # move to LGD: smolt after spring migration, at top of LGD
        Ma[y,i,o_nat,j] <- Mb[y,i,o_nat,j] * phi_Mb_Ma[y,i,o_nat,j]
        
        # derived survival for fitting: fall trap to LGD
        phi_Pa_Ma[y,i,j] <- Ma[y,i,o_nat,j]/max(Pa[y,i,j] * phi_Pb_Pa[i,j], Ma[y,i,o_nat,j])
      }
      
      # flag that tells us how often max constraint is violated
      # will remove this eventually. the max constraint is needed for MCMC during early tuning (crashes w/o it)
      # but none of the converged samples have this occur
      bad_flag[y,j] <- ifelse(Pa[y,i_spring,j] * phi_Pb_Pa[i_spring,j] < Ma[y,i_spring,o_nat,j], 1, 0)
      
      # derived survival for fitting: summer tagging to LGD
      phi_Pb_Ma[y,j] <- sum(Ma[y,1:ni,o_nat,j])/Pb[y,j]
      
      # put hatchery smolts in tributary
      Mb[y,i_spring,o_hat,j] <- Mb_obs[y,i_spring,o_hat,j]
      
      # move hatchery smolts from tributary to LGD
      Ma[y,i_spring,o_hat,j] <- Mb[y,i_spring,o_hat,j] * phi_Mb_Ma[y,i_spring,o_hat,j]
      
      # create zeros for fall migrant hatchery fish at LGD
      # needed because we sum over this dimension below
      Ma[y,i_fall,o_hat,j] <- 0
      
      # origin-specific processes
      for (o in 1:no) {
        # move origin-specific smolts from LGD to estuary (ocean age 0)
        O0[y,o,j] <- sum(Ma[y,1:ni,o,j]) * phi_Ma_O0[y,o]
        
        # move juveniles through ocean ages and survivals
        O[y,1,o,j] <- O0[y,o,j] * phi_O0_O1[y,o,j] # survive first winter at sea. now SWA1, TA3
        O[y,2,o,j] <- O[y,1,o,j] * (1 - psi_O1_Rb[y,o,j]) * phi_O1_O2[y,o,j] # don't mature at SWA1 and survive second winter at sea. now SWA2, TA4
        O[y,3,o,j] <- O[y,2,o,j] * (1 - psi_O2_Rb[y,o,j]) * phi_O2_O3[y,o,j] # don't mature at SWA2 and survive third winter at sea. now SWA3, TA5
        
        # mature and return to river in appropriate year at age
        Rb[y+kmin+1-1,1,o,j] <- O[y,1,o,j] * psi_O1_Rb[y,o,j]
        Rb[y+kmin+2-1,2,o,j] <- O[y,2,o,j] * psi_O2_Rb[y,o,j]
        Rb[y+kmin+3-1,3,o,j] <- O[y,3,o,j]
        
        # adult in-river processes: all age/origin specific
        # for these adult stages, y represents the brood year fish returned in
        for (k in 1:nk) {
          # returning adults making it to BON (survive sea lions * survive fisheries downstream of BON)
          Rb_BON[y,k,o,j] <- Rb[y,k,o,j] * phi_SL[y,j] * (1 - Ub[y,k,o])
          
          # survive upstream migration from BON to LGR and add strays
          Ra_LGR[y,k,o,j] <- Rb_BON[y,k,o,j] * phi_Rb_Ra[y,o]
          
          # add strays
          Ra[y,k,o,j] <- Ra_LGR[y,k,o,j] + n_stray_tot[y,j] * stray_comp[k,o,j]
          
          # remove fish at weir: use max() to ensure that fewer fish were removed than existed
          Sb[y,k,o,j] <- max(Ra[y,k,o,j] - n_remove[y,k,o,j], 1)
          
          # survive pre-spawn mortality
          Sa[y,k,o,j] <- Sb[y,k,o,j] * phi_Sb_Sa[y,j]
          
          # calculate egg production
          eggs[y,k,o,j] <- Sa[y,k,o,j] * p_female[k,j] * f[k]
          
          # calculate "adjusted carcasses": accounts for sampling bias relative to weir
          Sa_adj[y,k,o,j] <- Sa[y,k,o,j] * carc_adj[k,j]
        }
      }
      
      # total adults returned to trib
      Ra_tot[y,j] <- sum(Ra[y,1:nk,1:no,j])
      
      # total spawning adults
      Sa_tot[y,j] <- sum(Sa[y,1:nk,1:no,j])
      
      # total egg production: used as spawning stock
      f_tot[y,j] <- sum(eggs[y,1:nk,1:no,j])
      
      # reformat returns by age/origin for fitting to weir comp data (basically cbind two array slices)
      Ra_2d[y,1:nk,j] <- Ra[y,1:nk,o_nat,j]            # nat., all ages
      Ra_2d[y,(nk+1):(2*nk),j] <- Ra[y,1:nk,o_hat,j]   # hat., all ages

      # reformat "adjusted carcasses" by age/origin for fitting to carcass comp data (basically cbind two array slices)
      Sa_adj_2d[y,1:nk,j] <- Sa_adj[y,1:nk,o_nat,j]            # nat., all ages
      Sa_adj_2d[y,(nk+1):(2*nk),j] <- Sa_adj[y,1:nk,o_hat,j]   # hat., all ages

      # calculate age compositions
      for (ko in 1:nko) {
        q_Ra[y,ko,j] <- Ra_2d[y,ko,j]/sum(Ra_2d[y,1:nko,j])
        q_Sa_adj[y,ko,j] <- Sa_adj_2d[y,ko,j]/sum(Sa_adj_2d[y,1:nko,j])
      }
      
      # calculate misc derived quantities
      Pb_per_Sa_tot[y,j] <- Pb[y,j]/Sa_tot[y,j]                  # parr per spawner
      Pb_per_f_tot[y,j] <- Pb[y,j]/f_tot[y,j]                    # parr per egg
      Mb_per_Sa_tot[y,j] <- sum(Mb[y,1:ni,o_nat,j])/Sa_tot[y,j]  # smolt per spawner
    }
    
    # BON -> BON survival -- can't be calculated for all brood years in model
    for (y in (kmax+1):(ny-kmax)) {
      for (o in 1:no) {
        # calculate year/origin/pop specific survival
        phi_O0_Rb_BON[y,o,j] <- (
          Rb_BON[y+kmin+1-1,1,o,j] + 
            Rb_BON[y+kmin+2-1,2,o,j] + 
            Rb_BON[y+kmin+3-1,3,o,j])/
          ifelse(O0[y,o,j] == 0, 1, O0[y,o,j])
      }
    }
    
    # spawners per spawner -- can't be calculated for all brood years in model
    for (y in (kmax+1):(ny-kmax)) {
      Sa_tot_per_Sa_tot[y,j] <- (
        sum(Sa[y+kmin+1-1,1,1:no,j]) +       # age 3 adults produced by spawners in brood year y
          sum(Sa[y+kmin+2-1,2,1:no,j]) +     # age 4 adults produced by spawners in brood year y
          sum(Sa[y+kmin+3-1,3,1:no,j]))/     # age 5 adults produced by spawners in brood year y
        Sa_tot[y,j]                               # total spawners in brood year y
    }
    
    # adults per smolt -- can't be calculated for all brood years in model
    for (y in (kmax+1):(ny-kmax)) {
      for (o in 1:no) {
        Ra_per_Ma[y,o,j] <- (
          Ra_LGR[y+kmin+1-1,1,o,j] +       # age 3 adults back at trib that were once smolts from brood year y
            Ra_LGR[y+kmin+2-1,2,o,j] +     # age 4 adults back at trib that were once smolts from brood year y
            Ra_LGR[y+kmin+3-1,3,o,j])/     # age 5 adults back at trib that were once smolts from brood year y
          ifelse(sum(Ma[y,1:ni,o,j]) == 0, 1, sum(Ma[y,1:ni,o,j]))                    # total smolts at LGR from brood year y
      }
    }
  }
  
  ### OBSERVATION MODEL ###
  
  # the content under the "data likelihood" header is the only portion needed to fit the model
  # the other content is for evaluating model fit and parsimony
  # objects labeled "new" pertain to simulated data, for evaluating model fit to the observed data compared to "perfect data"
  # objects labeled "dev" represent measure of the deviation of the data from the model
  # objects labeled "lppd" represent log posterior predictive density; used in calculating WAIC
  
  # the long form index notation (i.e., for(d in 1:n_fit){}) enables skipping over missing data while keeping all object dimensions identical
  
  # age/origin composition
  for (j in 1:nj) {
    for (y in (kmax+1):ny_obs) {
      # AT WEIR
      # data likelihood
      weir_x_obs[y,1:nko,j] ~ dmulti(q_Ra[y,1:nko,j], weir_nx_obs[y,j])
      
      # simulate new data
      weir_x_obs_new[y,1:nko,j] ~ dmulti(q_Ra[y,1:nko,j], weir_nx_obs[y,j])
      
      # calculate expected count add small number to avoid division by zero
      expected_weir_x_obs[y,1:nko,j] <- q_Ra[y,1:nko,j] * weir_nx_obs[y,j] + 1e-6
      
      # calculate fit statistic: chi-squared statistic
      weir_x_obs_dev[y,j] <- sum(((weir_x_obs[y,1:nko,j] - expected_weir_x_obs[y,1:nko,j])^2)/expected_weir_x_obs[y,1:nko,j])
      weir_x_obs_new_dev[y,j] <- sum(((weir_x_obs_new[y,1:nko,j] - expected_weir_x_obs[y,1:nko,j])^2)/expected_weir_x_obs[y,1:nko,j])
      
      # calculate log posterior predictive density
      weir_x_obs_lppd[y,j] <- logdensity.multi(weir_x_obs[y,1:nko,j], q_Ra[y,1:nko,j], weir_nx_obs[y,j])
      
      # FOR CARCASSES
      # data likelihood
      carc_x_obs[y,1:nko,j] ~ dmulti(q_Sa_adj[y,1:nko,j], carc_nx_obs[y,j])
      
      # simulate new data
      carc_x_obs_new[y,1:nko,j] ~ dmulti(q_Ra[y,1:nko,j], carc_nx_obs[y,j])
      
      # calculate expected count add small number to avoid division by zero
      expected_carc_x_obs[y,1:nko,j] <- q_Ra[y,1:nko,j] * carc_nx_obs[y,j] + 1e-6
      
      # calculate fit statistic: chi-squared statistic
      carc_x_obs_dev[y,j] <- sum(((carc_x_obs[y,1:nko,j] - expected_carc_x_obs[y,1:nko,j])^2)/expected_carc_x_obs[y,1:nko,j])
      carc_x_obs_new_dev[y,j] <- sum(((carc_x_obs_new[y,1:nko,j] - expected_carc_x_obs[y,1:nko,j])^2)/expected_carc_x_obs[y,1:nko,j])
      
      # calculate log posterior predictive density
      carc_x_obs_lppd[y,j] <- logdensity.multi(carc_x_obs[y,1:nko,j], q_Ra[y,1:nko,j], carc_nx_obs[y,j])
    }
  }
  
  # the "fit" objects allow using square objects but skipping over missing obs in likelihoods
  # fall trap abundance
  for (d in 1:nfit_Pa) {
    # data likelihood
    Pa_obs[fit_Pa[d,1],fit_Pa[d,2],fit_Pa[d,3]] ~ dlnorm(log(Pa[fit_Pa[d,1],fit_Pa[d,2],fit_Pa[d,3]]), 1/sig_Pa_obs[fit_Pa[d,1],fit_Pa[d,2],fit_Pa[d,3]]^2)
    
    # simulate new data
    Pa_obs_new[fit_Pa[d,1],fit_Pa[d,2],fit_Pa[d,3]] ~ dlnorm(log(Pa[fit_Pa[d,1],fit_Pa[d,2],fit_Pa[d,3]]), 1/sig_Pa_obs[fit_Pa[d,1],fit_Pa[d,2],fit_Pa[d,3]]^2)
    
    # calculate fit statistic: squared log deviates
    Pa_obs_dev[fit_Pa[d,1],fit_Pa[d,2],fit_Pa[d,3]] <- (log(Pa_obs[fit_Pa[d,1],fit_Pa[d,2],fit_Pa[d,3]]) - log(Pa[fit_Pa[d,1],fit_Pa[d,2],fit_Pa[d,3]]))^2
    Pa_obs_new_dev[fit_Pa[d,1],fit_Pa[d,2],fit_Pa[d,3]] <- (log(Pa_obs_new[fit_Pa[d,1],fit_Pa[d,2],fit_Pa[d,3]]) - log(Pa[fit_Pa[d,1],fit_Pa[d,2],fit_Pa[d,3]]))^2
    
    # calculate log posterior predictive density
    Pa_obs_lppd[fit_Pa[d,1],fit_Pa[d,2],fit_Pa[d,3]] <- logdensity.lnorm(Pa_obs[fit_Pa[d,1],fit_Pa[d,2],fit_Pa[d,3]], log(Pa[fit_Pa[d,1],fit_Pa[d,2],fit_Pa[d,3]]), 1/sig_Pa_obs[fit_Pa[d,1],fit_Pa[d,2],fit_Pa[d,3]]^2)
  }
  
  # spring trap abundance
  for (d in 1:nfit_Mb) {
    # data likelihood
    Mb_obs[fit_Mb[d,1],fit_Mb[d,2],fit_Mb[d,3],fit_Mb[d,4]] ~ dlnorm(log(Mb[fit_Mb[d,1],fit_Mb[d,2],fit_Mb[d,3],fit_Mb[d,4]]), 1/sig_Mb_obs[fit_Mb[d,1],fit_Mb[d,2],fit_Mb[d,3],fit_Mb[d,4]]^2)
    
    # simulate new data
    Mb_obs_new[fit_Mb[d,1],fit_Mb[d,2],fit_Mb[d,3],fit_Mb[d,4]] ~ dlnorm(log(Mb[fit_Mb[d,1],fit_Mb[d,2],fit_Mb[d,3],fit_Mb[d,4]]), 1/sig_Mb_obs[fit_Mb[d,1],fit_Mb[d,2],fit_Mb[d,3],fit_Mb[d,4]]^2)
    
    # calculate fit statistic: squared log deviates
    Mb_obs_dev[fit_Mb[d,1],fit_Mb[d,2],fit_Mb[d,3],fit_Mb[d,4]] <- (log(Mb_obs[fit_Mb[d,1],fit_Mb[d,2],fit_Mb[d,3],fit_Mb[d,4]]) - log(Mb[fit_Mb[d,1],fit_Mb[d,2],fit_Mb[d,3],fit_Mb[d,4]]))^2
    Mb_obs_new_dev[fit_Mb[d,1],fit_Mb[d,2],fit_Mb[d,3],fit_Mb[d,4]] <- (log(Mb_obs_new[fit_Mb[d,1],fit_Mb[d,2],fit_Mb[d,3],fit_Mb[d,4]]) - log(Mb[fit_Mb[d,1],fit_Mb[d,2],fit_Mb[d,3],fit_Mb[d,4]]))^2
    
    # calculate log posterior predictive density
    Mb_obs_lppd[fit_Mb[d,1],fit_Mb[d,2],fit_Mb[d,3],fit_Mb[d,4]] <- logdensity.lnorm(Mb_obs[fit_Mb[d,1],fit_Mb[d,2],fit_Mb[d,3],fit_Mb[d,4]], log(Mb[fit_Mb[d,1],fit_Mb[d,2],fit_Mb[d,3],fit_Mb[d,4]]), 1/sig_Mb_obs[fit_Mb[d,1],fit_Mb[d,2],fit_Mb[d,3],fit_Mb[d,4]]^2)
  }
  
  # adult abundance
  for (d in 1:nfit_Ra) {
    # data likelihood
    Ra_obs[fit_Ra[d,1],fit_Ra[d,2]] ~ dlnorm(log(Ra_tot[fit_Ra[d,1],fit_Ra[d,2]]), 1/sig_Ra_obs[fit_Ra[d,1],fit_Ra[d,2]]^2)
    
    # simulate new data
    Ra_obs_new[fit_Ra[d,1],fit_Ra[d,2]] ~ dlnorm(log(Ra_tot[fit_Ra[d,1],fit_Ra[d,2]]), 1/sig_Ra_obs[fit_Ra[d,1],fit_Ra[d,2]]^2)
    
    # calculate fit statistic: squared log deviates
    Ra_obs_dev[fit_Ra[d,1],fit_Ra[d,2]] <- (log(Ra_obs[fit_Ra[d,1],fit_Ra[d,2]]) - log(Ra_tot[fit_Ra[d,1],fit_Ra[d,2]]))^2
    Ra_obs_new_dev[fit_Ra[d,1],fit_Ra[d,2]] <- (log(Ra_obs_new[fit_Ra[d,1],fit_Ra[d,2]]) - log(Ra_tot[fit_Ra[d,1],fit_Ra[d,2]]))^2
    
    # calculate log posterior predictive density
    Ra_obs_lppd[fit_Ra[d,1],fit_Ra[d,2]] <- logdensity.lnorm(Ra_obs[fit_Ra[d,1],fit_Ra[d,2]], log(Ra_tot[fit_Ra[d,1],fit_Ra[d,2]]), 1/sig_Ra_obs[fit_Ra[d,1],fit_Ra[d,2]]^2)
  }
  
  # summer tagging to LGD
  for (d in 1:nfit_Lphi_Pb_Ma) {
    # data likelihood
    Lphi_obs_Pb_Ma[fit_Lphi_Pb_Ma[d,1],fit_Lphi_Pb_Ma[d,2]] ~ dnorm(logit(phi_Pb_Ma[fit_Lphi_Pb_Ma[d,1],fit_Lphi_Pb_Ma[d,2]]), 1/sig_Lphi_obs_Pb_Ma[fit_Lphi_Pb_Ma[d,1],fit_Lphi_Pb_Ma[d,2]]^2)
    
    # simulate new data
    Lphi_obs_new_Pb_Ma[fit_Lphi_Pb_Ma[d,1],fit_Lphi_Pb_Ma[d,2]] ~ dnorm(logit(phi_Pb_Ma[fit_Lphi_Pb_Ma[d,1],fit_Lphi_Pb_Ma[d,2]]), 1/sig_Lphi_obs_Pb_Ma[fit_Lphi_Pb_Ma[d,1],fit_Lphi_Pb_Ma[d,2]]^2)
    
    # calculate fit statistic: squared logit deviates
    Lphi_obs_Pb_Ma_dev[fit_Lphi_Pb_Ma[d,1],fit_Lphi_Pb_Ma[d,2]] <- (Lphi_obs_Pb_Ma[fit_Lphi_Pb_Ma[d,1],fit_Lphi_Pb_Ma[d,2]] - logit(phi_Pb_Ma[fit_Lphi_Pb_Ma[d,1],fit_Lphi_Pb_Ma[d,2]]))^2
    Lphi_obs_new_Pb_Ma_dev[fit_Lphi_Pb_Ma[d,1],fit_Lphi_Pb_Ma[d,2]] <- (Lphi_obs_new_Pb_Ma[fit_Lphi_Pb_Ma[d,1],fit_Lphi_Pb_Ma[d,2]] - logit(phi_Pb_Ma[fit_Lphi_Pb_Ma[d,1],fit_Lphi_Pb_Ma[d,2]]))^2
    
    # calculate log posterior predictive density
    Lphi_obs_Pb_Ma_lppd[fit_Lphi_Pb_Ma[d,1],fit_Lphi_Pb_Ma[d,2]] <- logdensity.norm(Lphi_obs_Pb_Ma[fit_Lphi_Pb_Ma[d,1],fit_Lphi_Pb_Ma[d,2]], logit(phi_Pb_Ma[fit_Lphi_Pb_Ma[d,1],fit_Lphi_Pb_Ma[d,2]]), 1/sig_Lphi_obs_Pb_Ma[fit_Lphi_Pb_Ma[d,1],fit_Lphi_Pb_Ma[d,2]]^2)
  }
  
  # fall tagging to LGD
  for (d in 1:nfit_Lphi_Pa_Ma) {
    # data likelihood
    Lphi_obs_Pa_Ma[fit_Lphi_Pa_Ma[d,1],fit_Lphi_Pa_Ma[d,2],fit_Lphi_Pa_Ma[d,3]] ~ dnorm(logit(phi_Pa_Ma[fit_Lphi_Pa_Ma[d,1],fit_Lphi_Pa_Ma[d,2],fit_Lphi_Pa_Ma[d,3]]), 1/sig_Lphi_obs_Pa_Ma[fit_Lphi_Pa_Ma[d,1],fit_Lphi_Pa_Ma[d,2],fit_Lphi_Pa_Ma[d,3]]^2)
    
    # simulate new data
    Lphi_obs_new_Pa_Ma[fit_Lphi_Pa_Ma[d,1],fit_Lphi_Pa_Ma[d,2],fit_Lphi_Pa_Ma[d,3]] ~ dnorm(logit(phi_Pa_Ma[fit_Lphi_Pa_Ma[d,1],fit_Lphi_Pa_Ma[d,2],fit_Lphi_Pa_Ma[d,3]]), 1/sig_Lphi_obs_Pa_Ma[fit_Lphi_Pa_Ma[d,1],fit_Lphi_Pa_Ma[d,2],fit_Lphi_Pa_Ma[d,3]]^2)
    
    # calculate fit statistic: squared logit deviates
    Lphi_obs_Pa_Ma_dev[fit_Lphi_Pa_Ma[d,1],fit_Lphi_Pa_Ma[d,2],fit_Lphi_Pa_Ma[d,3]] <- (Lphi_obs_Pa_Ma[fit_Lphi_Pa_Ma[d,1],fit_Lphi_Pa_Ma[d,2],fit_Lphi_Pa_Ma[d,3]] - logit(phi_Pa_Ma[fit_Lphi_Pa_Ma[d,1],fit_Lphi_Pa_Ma[d,2],fit_Lphi_Pa_Ma[d,3]]))^2
    Lphi_obs_new_Pa_Ma_dev[fit_Lphi_Pa_Ma[d,1],fit_Lphi_Pa_Ma[d,2],fit_Lphi_Pa_Ma[d,3]] <- (Lphi_obs_new_Pa_Ma[fit_Lphi_Pa_Ma[d,1],fit_Lphi_Pa_Ma[d,2],fit_Lphi_Pa_Ma[d,3]] - logit(phi_Pa_Ma[fit_Lphi_Pa_Ma[d,1],fit_Lphi_Pa_Ma[d,2],fit_Lphi_Pa_Ma[d,3]]))^2
    
    # calculate log posterior predictive density
    Lphi_obs_Pa_Ma_lppd[fit_Lphi_Pa_Ma[d,1],fit_Lphi_Pa_Ma[d,2],fit_Lphi_Pa_Ma[d,3]] <- logdensity.norm(Lphi_obs_Pa_Ma[fit_Lphi_Pa_Ma[d,1],fit_Lphi_Pa_Ma[d,2],fit_Lphi_Pa_Ma[d,3]], logit(phi_Pa_Ma[fit_Lphi_Pa_Ma[d,1],fit_Lphi_Pa_Ma[d,2],fit_Lphi_Pa_Ma[d,3]]), 1/sig_Lphi_obs_Pa_Ma[fit_Lphi_Pa_Ma[d,1],fit_Lphi_Pa_Ma[d,2],fit_Lphi_Pa_Ma[d,3]]^2)
  }
  
  # spring tagging/smolt releases to LGD
  for (d in 1:nfit_Lphi_Mb_Ma) {
    # data likelihood
    Lphi_obs_Mb_Ma[fit_Lphi_Mb_Ma[d,1],fit_Lphi_Mb_Ma[d,2],fit_Lphi_Mb_Ma[d,3],fit_Lphi_Mb_Ma[d,4]] ~ dnorm(logit(phi_Mb_Ma[fit_Lphi_Mb_Ma[d,1],fit_Lphi_Mb_Ma[d,2],fit_Lphi_Mb_Ma[d,3],fit_Lphi_Mb_Ma[d,4]]), 1/sig_Lphi_obs_Mb_Ma[fit_Lphi_Mb_Ma[d,1],fit_Lphi_Mb_Ma[d,2],fit_Lphi_Mb_Ma[d,3],fit_Lphi_Mb_Ma[d,4]]^2)
    
    # simulate new data
    Lphi_obs_new_Mb_Ma[fit_Lphi_Mb_Ma[d,1],fit_Lphi_Mb_Ma[d,2],fit_Lphi_Mb_Ma[d,3],fit_Lphi_Mb_Ma[d,4]] ~ dnorm(logit(phi_Mb_Ma[fit_Lphi_Mb_Ma[d,1],fit_Lphi_Mb_Ma[d,2],fit_Lphi_Mb_Ma[d,3],fit_Lphi_Mb_Ma[d,4]]), 1/sig_Lphi_obs_Mb_Ma[fit_Lphi_Mb_Ma[d,1],fit_Lphi_Mb_Ma[d,2],fit_Lphi_Mb_Ma[d,3],fit_Lphi_Mb_Ma[d,4]]^2)
    
    # calculate fit statistic: squared logit deviates
    Lphi_obs_Mb_Ma_dev[fit_Lphi_Mb_Ma[d,1],fit_Lphi_Mb_Ma[d,2],fit_Lphi_Mb_Ma[d,3],fit_Lphi_Mb_Ma[d,4]] <- (Lphi_obs_Mb_Ma[fit_Lphi_Mb_Ma[d,1],fit_Lphi_Mb_Ma[d,2],fit_Lphi_Mb_Ma[d,3],fit_Lphi_Mb_Ma[d,4]] - logit(phi_Mb_Ma[fit_Lphi_Mb_Ma[d,1],fit_Lphi_Mb_Ma[d,2],fit_Lphi_Mb_Ma[d,3],fit_Lphi_Mb_Ma[d,4]]))^2
    Lphi_obs_new_Mb_Ma_dev[fit_Lphi_Mb_Ma[d,1],fit_Lphi_Mb_Ma[d,2],fit_Lphi_Mb_Ma[d,3],fit_Lphi_Mb_Ma[d,4]] <- (Lphi_obs_new_Mb_Ma[fit_Lphi_Mb_Ma[d,1],fit_Lphi_Mb_Ma[d,2],fit_Lphi_Mb_Ma[d,3],fit_Lphi_Mb_Ma[d,4]] - logit(phi_Mb_Ma[fit_Lphi_Mb_Ma[d,1],fit_Lphi_Mb_Ma[d,2],fit_Lphi_Mb_Ma[d,3],fit_Lphi_Mb_Ma[d,4]]))^2
    
    # calculate log posterior predictive density
    Lphi_obs_Mb_Ma_lppd[fit_Lphi_Mb_Ma[d,1],fit_Lphi_Mb_Ma[d,2],fit_Lphi_Mb_Ma[d,3],fit_Lphi_Mb_Ma[d,4]] <- logdensity.norm(Lphi_obs_Mb_Ma[fit_Lphi_Mb_Ma[d,1],fit_Lphi_Mb_Ma[d,2],fit_Lphi_Mb_Ma[d,3],fit_Lphi_Mb_Ma[d,4]], logit(phi_Mb_Ma[fit_Lphi_Mb_Ma[d,1],fit_Lphi_Mb_Ma[d,2],fit_Lphi_Mb_Ma[d,3],fit_Lphi_Mb_Ma[d,4]]), 1/sig_Lphi_obs_Mb_Ma[fit_Lphi_Mb_Ma[d,1],fit_Lphi_Mb_Ma[d,2],fit_Lphi_Mb_Ma[d,3],fit_Lphi_Mb_Ma[d,4]]^2)
  }
  
  # hydrosystem survival
  for (d in 1:nfit_Lphi_Ma_O0) {
    # data likelihood
    Lphi_obs_Ma_O0[fit_Lphi_Ma_O0[d,1],fit_Lphi_Ma_O0[d,2]] ~ dnorm(logit(phi_Ma_O0[fit_Lphi_Ma_O0[d,1],fit_Lphi_Ma_O0[d,2]]), 1/sig_Lphi_obs_Ma_O0[fit_Lphi_Ma_O0[d,1],fit_Lphi_Ma_O0[d,2]]^2)
    
    # simulate new data
    Lphi_obs_new_Ma_O0[fit_Lphi_Ma_O0[d,1],fit_Lphi_Ma_O0[d,2]] ~ dnorm(logit(phi_Ma_O0[fit_Lphi_Ma_O0[d,1],fit_Lphi_Ma_O0[d,2]]), 1/sig_Lphi_obs_Ma_O0[fit_Lphi_Ma_O0[d,1],fit_Lphi_Ma_O0[d,2]]^2)
    
    # calculate fit statistic: squared logit deviates
    Lphi_obs_Ma_O0_dev[fit_Lphi_Ma_O0[d,1],fit_Lphi_Ma_O0[d,2]] <- (Lphi_obs_Ma_O0[fit_Lphi_Ma_O0[d,1],fit_Lphi_Ma_O0[d,2]] - logit(phi_Ma_O0[fit_Lphi_Ma_O0[d,1],fit_Lphi_Ma_O0[d,2]]))^2
    Lphi_obs_new_Ma_O0_dev[fit_Lphi_Ma_O0[d,1],fit_Lphi_Ma_O0[d,2]] <- (Lphi_obs_new_Ma_O0[fit_Lphi_Ma_O0[d,1],fit_Lphi_Ma_O0[d,2]] - logit(phi_Ma_O0[fit_Lphi_Ma_O0[d,1],fit_Lphi_Ma_O0[d,2]]))^2
    
    # calculate log posterior predictive density
    Lphi_obs_Ma_O0_lppd[fit_Lphi_Ma_O0[d,1],fit_Lphi_Ma_O0[d,2]] <- logdensity.norm(Lphi_obs_Ma_O0[fit_Lphi_Ma_O0[d,1],fit_Lphi_Ma_O0[d,2]], logit(phi_Ma_O0[fit_Lphi_Ma_O0[d,1],fit_Lphi_Ma_O0[d,2]]), 1/sig_Lphi_obs_Ma_O0[fit_Lphi_Ma_O0[d,1],fit_Lphi_Ma_O0[d,2]]^2)
  }
  
  # movement survival from BON to LGR
  for (d in 1:nfit_LGR_adults) {
    # data likelihood
    LGR_adults[fit_LGR_adults[d,1],fit_LGR_adults[d,2]] ~ dbin(phi_Rb_Ra[fit_LGR_adults[d,1],fit_LGR_adults[d,2]], BON_adults[fit_LGR_adults[d,1],fit_LGR_adults[d,2]])
    
    # simulate new data
    LGR_adults_new[fit_LGR_adults[d,1],fit_LGR_adults[d,2]] ~ dbin(phi_Rb_Ra[fit_LGR_adults[d,1],fit_LGR_adults[d,2]], BON_adults[fit_LGR_adults[d,1],fit_LGR_adults[d,2]])
    
    # calculate expected count
    expected_LGR_adults[fit_LGR_adults[d,1],fit_LGR_adults[d,2]] <- phi_Rb_Ra[fit_LGR_adults[d,1],fit_LGR_adults[d,2]] * BON_adults[fit_LGR_adults[d,1],fit_LGR_adults[d,2]]
    
    # calculate fit statistic: chi-squared statistic
    LGR_adults_dev[fit_LGR_adults[d,1],fit_LGR_adults[d,2]] <- ((LGR_adults[fit_LGR_adults[d,1],fit_LGR_adults[d,2]] - expected_LGR_adults[fit_LGR_adults[d,1],fit_LGR_adults[d,2]])^2)/expected_LGR_adults[fit_LGR_adults[d,1],fit_LGR_adults[d,2]]
    LGR_adults_new_dev[fit_LGR_adults[d,1],fit_LGR_adults[d,2]] <- ((LGR_adults_new[fit_LGR_adults[d,1],fit_LGR_adults[d,2]] - expected_LGR_adults[fit_LGR_adults[d,1],fit_LGR_adults[d,2]])^2)/expected_LGR_adults[fit_LGR_adults[d,1],fit_LGR_adults[d,2]]
    
    # calculate log posterior predictive density
    LGR_adults_lppd[fit_LGR_adults[d,1],fit_LGR_adults[d,2]] <- logdensity.bin(LGR_adults[fit_LGR_adults[d,1],fit_LGR_adults[d,2]], phi_Rb_Ra[fit_LGR_adults[d,1],fit_LGR_adults[d,2]], BON_adults[fit_LGR_adults[d,1],fit_LGR_adults[d,2]])
  }
  
  # pre-spawn survival
  for (d in 1:nfit_spawned) {
    # data likelihood
    carcs_spawned[fit_spawned[d,1],fit_spawned[d,2]] ~ dbin(phi_Sb_Sa[fit_spawned[d,1],fit_spawned[d,2]], carcs_sampled[fit_spawned[d,1],fit_spawned[d,2]])
    
    # simulate new data
    carcs_spawned_new[fit_spawned[d,1],fit_spawned[d,2]] ~ dbin(phi_Sb_Sa[fit_spawned[d,1],fit_spawned[d,2]], carcs_sampled[fit_spawned[d,1],fit_spawned[d,2]])
    
    # calculate expected count
    expected_carcs_spawned[fit_spawned[d,1],fit_spawned[d,2]] <- phi_Sb_Sa[fit_spawned[d,1],fit_spawned[d,2]] * carcs_sampled[fit_spawned[d,1],fit_spawned[d,2]]
    
    # calculate fit statistic: chi-squared statistic
    carcs_spawned_dev[fit_spawned[d,1],fit_spawned[d,2]] <- ((carcs_spawned[fit_spawned[d,1],fit_spawned[d,2]] - expected_carcs_spawned[fit_spawned[d,1],fit_spawned[d,2]])^2)/expected_carcs_spawned[fit_spawned[d,1],fit_spawned[d,2]]
    carcs_spawned_new_dev[fit_spawned[d,1],fit_spawned[d,2]] <- ((carcs_spawned_new[fit_spawned[d,1],fit_spawned[d,2]] - expected_carcs_spawned[fit_spawned[d,1],fit_spawned[d,2]])^2)/expected_carcs_spawned[fit_spawned[d,1],fit_spawned[d,2]]
    
    # calculate log posterior predictive density
    carcs_spawned_lppd[fit_spawned[d,1],fit_spawned[d,2]] <- logdensity.bin(carcs_spawned[fit_spawned[d,1],fit_spawned[d,2]], phi_Sb_Sa[fit_spawned[d,1],fit_spawned[d,2]], carcs_sampled[fit_spawned[d,1],fit_spawned[d,2]])
  }
  
  ### CALCULATE ALL PROCESS MODEL RESIDUALS ###
  
  for (y in (kmax+1):ny) {
    for (j in 1:nj) {
      
      # total summer parr recruitment
      lPb_resid[y,j] <- log(Pb[y,j]) - log(Pb_mean[y,j])
      
      # proportion of summer parr that become fall migrants
      Lpi_resid[y,j] <- logit(pi[y,i_fall,j]) - logit(mu_pi[i_fall,j])
      
      # LH-specific overwinter survival
      for (i in 1:ni) {
        Lphi_Pa_Mb_resid[y,i,j] <- Lphi_Pa_Mb[y,i,j] - (gamma0[i,j] + gamma1[i,j] * (Pa[y,i,j]/wul[j]))
      }
      
      # origin-specific quantities
      for (o in 1:no) {
        # movement survival to LGR
        Lphi_Mb_Ma_resid[y,o,j] <- Lphi_Mb_Ma[y,i_spring,o,j] - logit(mu_phi_Mb_Ma[i_spring,o,j])
        
        # pr(return at SWA1)
        Lpsi_O1_Rb_resid[y,o,j] <- Lpsi_O1_Rb[y,o,j] - logit(mu_psi_O1_Rb[o,j])
        
        # pr(return at SWA2|not returned at SWA1)
        Lpsi_O2_Rb_resid[y,o,j] <- Lpsi_O2_Rb[y,o,j] - logit(mu_psi_O2_Rb[o,j])
      }
      
      # pre-spawn survival
      Lphi_Sb_Sa_resid[y,j] <- Lphi_Sb_Sa[y,j] - logit(mu_phi_Sb_Sa[j])
    }
    
    for (o in 1:no) {
      # movement survival from LGR thru BON
      Lphi_Ma_O0_resid[y,o] <- Lphi_Ma_O0[y,o] - logit(mu_phi_Ma_O0[o])
      
      # movement survival from BON thru LGR
      Lphi_Rb_Ra_resid[y,o] <- Lphi_Rb_Ra[y,o] - logit(mu_phi_Rb_Ra[o])
      
      for (j in 1:nj) {
        # ocean survival
        Lphi_O0_O1_resid[y,o,j] <- Lphi_O0_O1[y,o,j] - logit(mu_phi_O0_O1[o,j])
        Lphi_O1_O2_resid[y,o,j] <- Lphi_O1_O2[y,o,j] - logit(mu_phi_O1_O2[o,j])
        Lphi_O2_O3_resid[y,o,j] <- Lphi_O2_O3[y,o,j] - logit(mu_phi_O2_O3[o,j])
      }
    }
  }
}
