jags_model_code = function() {
  
  for (j in 1:nj) {
    ### PRIORS: RECRUITMENT FUNCTION ###
    # total egg production to aggregate parr recruitment function
    alpha[j] ~ dbeta(1, 1)
    log_beta[j] ~ dnorm(0, 0.001) %_% T(,15)   # log capacity. bound to prevent nonsensically large draws
    beta[j] <- exp(log_beta[j])
    sigma_Pb[j] ~ dunif(0, 5)
    beta_per_peu[j] <- beta[j]/peu[j]
    
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
    sig_Lphi_Mb_Ma[i_spring,o_nat,j] ~ dunif(0, 5)
    mu_phi_Mb_Ma[i_fall,o_nat,j] <- mu_phi_Mb_Ma[i_spring,o_nat,j]
    sig_Lphi_Mb_Ma[i_fall,o_nat,j] <- sig_Lphi_Mb_Ma[i_spring,o_nat,j]
    
    # hatchery origin movement survival (trib to LGD): have spring migrants only
    mu_phi_Mb_Ma[i_spring,o_hat,j] ~ dbeta(1, 1)
    sig_Lphi_Mb_Ma[i_spring,o_hat,j] ~ dunif(0, 5)
    
    # pre-spawn survival (after brood-stock removal to successful spawning)
    mu_phi_Sb_Sa[j] ~ dbeta(1, 1)
    sig_Lphi_Sb_Sa[j] ~ dunif(0, 5)
    
    ### PRIORS: SEX-ASSIGNMENT ###
    for (o in 1:no) {
      mu_omega[s_female,o,j] ~ dbeta(1, 1)
      mu_omega[s_male,o,j] <- 1 - mu_omega[s_female,o,j]
      sig_Lomega[o,j] ~ dunif(0, 5)
    }
    
    ### PRIORS: MATURATION ###
    # sex/origin-specific maturation probabilities
    for (o in 1:no) {
      for (s in 1:ns) {
        # pr(return at SWA1)
        mu_psi_O1_Rb[s,o,j] ~ dbeta(1, 1)
        sig_Lpsi_O1_Rb[s,o,j] ~ dunif(0, 5)
        
        # pr(return at SWA2|not returned at SWA1)
        mu_psi_O2_Rb[s,o,j] ~ dbeta(1, 1)
        sig_Lpsi_O2_Rb[s,o,j] ~ dunif(0, 5)
        
        # pr(return at SWA3|not returned at SWA1 or SWA2)
        mu_psi_O3_Rb[s,o,j] <- 1
        sig_Lpsi_O3_Rb[s,o,j] <- 0
      }
    }
    
    ### PRIORS: HATCHERY STRAYS RETURNING IN YEARS WITH NO ASSOCIATED SMOLT RELEASE ###
    # the number of hatchery strays: only estimate in years where no other mechanism for generating hatchery fish
    for (i in 1:n_stray_yrs[j]) {
      n_stray_tot[stray_yrs[i,j],j] ~ dunif(0, 5000)
    }
    for (i in 1:n_not_stray_yrs[j]) {
      n_stray_tot[not_stray_yrs[i,j],j] <- 0
    }
    
    # when n_stray_tot > 0, what age/sex are they? (assume all hatchery fish)
    stray_comp_2d[1:nks,j] ~ ddirich(rep(1, nks))
    
    # put stray_comp_2d into array format for looping in process model
    for (k in 1:nk) {
      for (s in 1:ns) {
        # assume no strays for natural fish
        stray_comp[k,s,o_nat,j] <- 0
        
        # place stray_comp_2d in the right index locations
        stray_comp[k,s,o_hat,j] <- stray_comp_2d[k+nk*(s-1),j]
      }
    }
  }
  
  ### PRIORS: PARAMETERS THAT ARE CONSTANT ACROSS POPULATIONS ###
  
  # juvenile movement survival (LGD to estuary [ocean-age 0]): same for both LH types, different for origin types
  for (o in 1:no) {
    mu_phi_Ma_O0[o] ~ dbeta(1, 1)
    sig_Lphi_Ma_O0[o] ~ dunif(0, 5)
  }
  
  # ocean survival: natural origin
  mu_phi_O0_O1[o_nat] ~ dbeta(1, 1)            # first winter at sea: to become SWA1
  mu_phi_O1_O2[o_nat] ~ dbeta(1, 1)            # second winter at sea: to become SWA2
  mu_phi_O2_O3[o_nat] <- mu_phi_O1_O2[o_nat]   # third winter at sea: to become SW3
  sig_Lphi_O0_O1[o_nat] ~ dunif(0, 5)
  sig_Lphi_O1_O2[o_nat] ~ dunif(0, 5)
  sig_Lphi_O2_O3[o_nat] <- sig_Lphi_O1_O2[o_nat]
  
  # ocean survival: hatchery origin; use a scaler that adjusts natural origin survival to get hatchery survival
  O_phi_scaler_nat_hat ~ dt(0, 1/1.566^2, 7.763)
  logit(mu_phi_O0_O1[o_hat]) <- logit(mu_phi_O0_O1[o_nat]) + O_phi_scaler_nat_hat
  logit(mu_phi_O1_O2[o_hat]) <- logit(mu_phi_O1_O2[o_nat]) + O_phi_scaler_nat_hat
  logit(mu_phi_O2_O3[o_hat]) <- logit(mu_phi_O1_O2[o_nat]) + O_phi_scaler_nat_hat
  
  # AR(1) coefficient and year 0 residual
  kappa_phi_O0_O1 ~ dunif(-0.99,0.99)
  Lphi_O0_O1_resid[kmax] ~ dnorm(0, (1/sig_Lphi_O0_O1[o_nat]^2) * (1 - kappa_phi_O0_O1^2))
  
  # carcass vs. weir composition correction factor coefficients
  for (i in 1:3) {
    z[i] ~ dunif(-10,10)
  }
  
  # calculate correction factor
  for (k in 1:3) {
    for (s in 1:2) {
      log(carc_adj[k,s]) <- z[1] * age3[k] + z[2] * age5[k] + z[3] * male[s]
    }
  }
  
  ### PRIORS: BROOD-YEAR-SPECIFIC PARAMETERS ###
  for (y in (kmax+1):ny) {
    for (j in 1:nj) {
      # aggregate parr to LH-specific parr
      Lpi1[y,j] ~ dnorm(logit(mu_pi[i_fall,j]), 1/sig_Lpi[j]^2)
      pi[y,i_fall,j] <- ilogit(Lpi1[y,j])
      pi[y,i_spring,j] <- 1 - pi[y,i_fall,j]
      
      # overwinter survival: density-dependent
      for (i in 1:ni) {
        Lphi_Pa_Mb[y,i,j] ~ dnorm(gamma0[i,j] + gamma1[i,j] * (Pa[y,i,j]/peu[j]), 1/sig_Lphi_Pa_Mb[i,j]^2)
        phi_Pa_Mb[y,i,j] <- ilogit(Lphi_Pa_Mb[y,i,j])
      }
      
      # natural origin movement survival: trib to LGD
      # assume equal between LH types
      Lphi_Mb_Ma[y,i_spring,o_nat,j] ~ dnorm(logit(mu_phi_Mb_Ma[i_spring,o_nat,j]), 1/sig_Lphi_Mb_Ma[i_spring,o_nat,j]^2)
      phi_Mb_Ma[y,i_spring,o_nat,j] <- ilogit(Lphi_Mb_Ma[y,i_spring,o_nat,j])
      phi_Mb_Ma[y,i_fall,o_nat,j] <- phi_Mb_Ma[y,i_spring,o_nat,j]
      
      # hatchery origin movement survival: trib to LGD
      # spring migrants only
      Lphi_Mb_Ma[y,i_spring,o_hat,j] ~ dnorm(logit(mu_phi_Mb_Ma[i_spring,o_hat,j]), 1/sig_Lphi_Mb_Ma[i_spring,o_hat,j]^2)
      phi_Mb_Ma[y,i_spring,o_hat,j] <- ilogit(Lphi_Mb_Ma[y,i_spring,o_hat,j])
      
      for (o in 1:no) {
        # probability of returning as female or male by origin
        Lomega1[y,o,j] ~ dnorm(logit(mu_omega[s_female,o,j]), 1/sig_Lomega[o,j]^2)
        omega[y,s_female,o,j] <- ilogit(Lomega1[y,o,j])
        omega[y,s_male,o,j] <- 1 - omega[y,s_female,o,j]
        
        for (s in 1:ns) {
          # pr(return at SWA1)
          Lpsi_O1_Rb[y,s,o,j] ~ dnorm(logit(mu_psi_O1_Rb[s,o,j]), 1/sig_Lpsi_O1_Rb[s,o,j]^2)
          psi_O1_Rb[y,s,o,j] <- ilogit(Lpsi_O1_Rb[y,s,o,j])
          
          # pr(return at SWA2|not returned at SWA1)
          Lpsi_O2_Rb[y,s,o,j] ~ dnorm(logit(mu_psi_O2_Rb[s,o,j]), 1/sig_Lpsi_O2_Rb[s,o,j]^2)
          psi_O2_Rb[y,s,o,j] <- ilogit(Lpsi_O2_Rb[y,s,o,j])
          
          # pr(return at SWA3|not returned at SWA1 or SWA2)
          psi_O3_Rb[y,s,o,j] <- 1
        }
      }
      
      # pre-spawn survival
      Lphi_Sb_Sa[y,j] ~ dnorm(logit(mu_phi_Sb_Sa[j]), 1/sig_Lphi_Sb_Sa[j]^2)
      phi_Sb_Sa[y,j] <- ilogit(Lphi_Sb_Sa[y,j])
    }
    
    # movement survival: LGD to estuary
    # separate for each origin type
    for (o in 1:no) {
      Lphi_Ma_O0[y,o] ~ dnorm(logit(mu_phi_Ma_O0[o]), 1/sig_Lphi_Ma_O0[o]^2)
      phi_Ma_O0[y,o] <- ilogit(Lphi_Ma_O0[y,o])
    }
    
    # natural origin ocean survival SWA0 -> SWA1 (uses AR(1) on process noise)
    Lphi_O0_O1[y,o_nat] ~ dnorm(logit(mu_phi_O0_O1[o_nat]) + Lphi_O0_O1_resid[y-1] * kappa_phi_O0_O1, 1/sig_Lphi_O0_O1[o_nat]^2)
    phi_O0_O1[y,o_nat] <- ilogit(Lphi_O0_O1[y,o_nat])
    
    # natural origin ocean survival SWA1 -> SWA2
    Lphi_O1_O2[y,o_nat] ~ dnorm(logit(mu_phi_O1_O2[o_nat]), 1/sig_Lphi_O1_O2[o_nat]^2)
    phi_O1_O2[y,o_nat] <- ilogit(Lphi_O1_O2[y,o_nat])
    
    # natural origin ocean survival SWA2 -> SWA3
    phi_O2_O3[y-1,o_nat] <- phi_O1_O2[y,o_nat]
    
    # hatchery origin ocean survival: use scaler
    logit(phi_O0_O1[y,o_hat]) <- logit(phi_O0_O1[y,o_nat]) + O_phi_scaler_nat_hat
    logit(phi_O1_O2[y,o_hat]) <- logit(phi_O1_O2[y,o_nat]) + O_phi_scaler_nat_hat
    logit(phi_O2_O3[y-1,o_hat]) <- logit(phi_O2_O3[y-1,o_nat]) + O_phi_scaler_nat_hat
  }
  
  # populate the last year/age of ocean survival with the mean across years
  phi_O2_O3[ny,o_nat] <- mu_phi_O2_O3[o_nat]
  phi_O2_O3[ny,o_hat] <- mu_phi_O2_O3[o_hat]
  
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
  for (j in 1:nj) {
    for (s in 1:ns) {
      p_init_prime[1,s,j] <- mu_phi_O0_O1[o_nat] * mu_psi_O1_Rb[s,o_nat,j]
      p_init_prime[2,s,j] <- mu_phi_O0_O1[o_nat] * (1 - mu_psi_O1_Rb[s,o_nat,j]) * mu_phi_O1_O2[o_nat] * mu_psi_O2_Rb[s,o_nat,j]
      p_init_prime[3,s,j] <- mu_phi_O0_O1[o_nat] * (1 - mu_psi_O1_Rb[s,o_nat,j]) * mu_phi_O1_O2[o_nat] * (1 - mu_psi_O2_Rb[s,o_nat,j]) * mu_phi_O2_O3[o_nat] * mu_psi_O3_Rb[s,o_nat,j]
      for (k in 1:nk) {
        p_init[k,s,j] <- p_init_prime[k,s,j]/sum(p_init_prime[1:nk,s,j])
      }
    }
    
    # hyperparameters for initialization
    mu_init_recruits[j] ~ dunif(0, max_init_recruits[j])
    sig_init_lrecruits[j] ~ dunif(0,10)
    
    # for natural origin
    for (y in 1:kmax) {
      # random, constant-mean adult recruits for first kmax brood years
      init_recruits[y,j] ~ dlnorm(log(mu_init_recruits[j]), 1/sig_init_lrecruits[j]^2) %_% T(,max_init_recruits[j])
      for (s in 1:ns) {
        # apportion them to each sex
        init_recruits_sex[y,s,j] <- init_recruits[y,j] * mu_omega[s,o_nat,j]
        for (k in 1:nk) {
          # apportion them to return year, age, and sex
          Rb[y+kmin+k-1,k,s,o_nat,j] <- init_recruits_sex[y,s,j] * p_init[k,s,j]
          
          # for hatchery origin - no operations these years, but still need to populate with a number
          Rb[y+kmin+k-1,k,s,o_hat,j] <- 0
        }
      }
    }
  }
  
  ### PROCESS MODEL: COMPLETE LIFE CYCLE FOR REMAINING BROOD YEARS ###
  for (j in 1:nj) {
    for (y in (kmax+1):ny) {
      # reproductive link: total summer parr
      Pb_pred[y,j] <- f_tot[y,j]/(1/alpha[j] + f_tot[y,j]/beta[j])
      Pb[y,j] ~ dlnorm(log(Pb_pred[y,j]), 1/sigma_Pb[j]^2)
      
      # natural origin tributary-to-LGD dynamics
      for (i in 1:ni) {
        # apportion to LH-strategy: parr after fall migration
        Pa[y,i,j] <- Pb[y,j] * pi[y,i,j]
        
        # survive over winter: smolt before spring migration
        Mb[y,i,o_nat,j] <- Pa[y,i,j] * phi_Pa_Mb[y,i,j]
        
        # move to LGD: smolt after spring migration, at top of LGD
        Ma[y,i,o_nat,j] <- Mb[y,i,o_nat,j] * phi_Mb_Ma[y,i,o_nat,j]
        
        # derived survival for fitting: fall trap to LGD
        phi_Pa_Ma[y,i,j] <- Ma[y,i,o_nat,j]/Pa[y,i,j]
      }
      
      # derived survival for fitting: summer tagging to LGD
      phi_Pb_Ma[y,j] <- sum(Ma[y,1:ni,o_nat,j])/Pb[y,j]
      
      # put hatchery smolts in tributary
      Mb[y,i_spring,o_hat,j] <- Mb_obs[y,i_spring,o_hat,j]
      
      # move hatchery smolts from tributary to LGD
      Ma[y,i_spring,o_hat,j] <- Mb[y,i_spring,o_hat,j] * phi_Mb_Ma[y,i_spring,o_hat,j]
      
      # create zeros for fall migrant hatchery fish at LGD
      # needed because we sum over this dimension below
      Ma[y,i_fall,o_hat,j] <- 0
      
      # sex/origin-specific processes
      for (s in 1:ns) {
        for (o in 1:no) {
          # move origin-specific smolts from LGD to estuary (ocean age 0) and assign to sex
          O0[y,s,o,j] <- sum(Ma[y,1:ni,o,j]) * phi_Ma_O0[y,o] * omega[y,s,o,j]
          
          # move juveniles through ocean ages and survivals
          O[y,1,s,o,j] <- O0[y,s,o,j] * phi_O0_O1[y,o] # survive first winter at sea. now SWA1, TA3
          O[y,2,s,o,j] <- O[y,1,s,o,j] * (1 - psi_O1_Rb[y,s,o,j]) * phi_O1_O2[y,o] # don't mature at SWA1 and survive second winter at sea. now SWA2, TA4
          O[y,3,s,o,j] <- O[y,2,s,o,j] * (1 - psi_O2_Rb[y,s,o,j]) * phi_O2_O3[y,o] # don't mature at SWA2 and survive third winter at sea. now SWA3, TA5
          
          # mature and return to river in appropriate year at age/sex
          Rb[y+kmin+1-1,1,s,o,j] <- O[y,1,s,o,j] * psi_O1_Rb[y,s,o,j]
          Rb[y+kmin+2-1,2,s,o,j] <- O[y,2,s,o,j] * psi_O2_Rb[y,s,o,j]
          Rb[y+kmin+3-1,3,s,o,j] <- O[y,3,s,o,j] * psi_O3_Rb[y,s,o,j]
          
          # adult in-river processes: all age/origin/sex specific
          # for these adult stages, y represents the brood year fish returned in
          for (k in 1:nk) {
            # survive upstream migration (survive sea lions * survive fisheries * survive through dams and to tributary) and add strays
            Ra[y,k,s,o,j] <- Rb[y,k,s,o,j] * phi_SL[y,j] * (1 - Ub[y,k,o]) * (1 - Ua[y,k,o]) * phi_D1_D4 * phi_D4_D5 * phi_D5_D8 * phi_D8_Ra + n_stray_tot[y,j] * stray_comp[k,s,o,j]
            
            # remove fish at weir: use max() to ensure that fewer fish were removed than existed
            Sb[y,k,s,o,j] <- max(Ra[y,k,s,o,j] - n_remove[y,k,s,o,j], 1)
            
            # survive pre-spawn mortality
            Sa[y,k,s,o,j] <- Sb[y,k,s,o,j] * phi_Sb_Sa[y,j]
            
            # calculate egg production
            eggs[y,k,s,o,j] <- Sa[y,k,s,o,j] * f[k,s]
            
            # calculate "adjusted carcasses": accounts for sampling bias relative to weir
            Sa_adj[y,k,s,o,j] <- Sa[y,k,s,o,j] * carc_adj[k,s]
          }
        }
      }
      
      # total adults returned to trib
      Ra_tot[y,j] <- sum(Ra[y,1:nk,1:ns,1:no,j])
      
      # total spawning adults
      Sa_tot[y,j] <- sum(Sa[y,1:nk,1:ns,1:no,j])
      
      # total egg production: used as spawning stock
      f_tot[y,j] <- sum(eggs[y,1:nk,1:ns,1:no,j])
      
      # reformat returns by age/sex/origin for fitting to weir comp data
      Ra_2d[y,1:nk,j] <- Ra[y,1:nk,s_female,o_nat,j]            # nat. females, all ages
      Ra_2d[y,(nk+1):(2*nk),j] <- Ra[y,1:nk,s_male,o_nat,j]     # nat. males, all ages
      Ra_2d[y,(2*nk+1):(3*nk),j] <- Ra[y,1:nk,s_female,o_hat,j] # hat. females, all ages
      Ra_2d[y,(3*nk+1):(4*nk),j] <- Ra[y,1:nk,s_male,o_hat,j]   # hat. males, all ages
      
      # reformat "adjusted carcasses" by age/sex/origin for fitting to carcass comp data
      Sa_adj_2d[y,1:nk,j] <- Sa_adj[y,1:nk,s_female,o_nat,j]            # nat. females, all ages
      Sa_adj_2d[y,(nk+1):(2*nk),j] <- Sa_adj[y,1:nk,s_male,o_nat,j]     # nat. males, all ages
      Sa_adj_2d[y,(2*nk+1):(3*nk),j] <- Sa_adj[y,1:nk,s_female,o_hat,j] # hat. females, all ages
      Sa_adj_2d[y,(3*nk+1):(4*nk),j] <- Sa_adj[y,1:nk,s_male,o_hat,j]   # hat. males, all ages
      
      # calculate age/sex compositions
      for (kso in 1:nkso) {
        q_Ra[y,kso,j] <- Ra_2d[y,kso,j]/sum(Ra_2d[y,1:nkso,j])
        q_Sa_adj[y,kso,j] <- Sa_adj_2d[y,kso,j]/sum(Sa_adj_2d[y,1:nkso,j])
      }
      
      # calculate misc derived quantities
      Pb_per_Sa_tot[y,j] <- Pb[y,j]/Sa_tot[y,j]                  # parr per spawner
      Pb_per_f_tot[y,j] <- Pb[y,j]/f_tot[y,j]                    # parr per egg
      Mb_per_Sa_tot[y,j] <- sum(Mb[y,1:ns,o_nat,j])/Sa_tot[y,j]  # smolt per spawner
    }
    
    # spawners per spawner -- can't be calculated for all brood years in model
    for (y in (kmax+1):(ny-kmax)) {
      Sa_tot_per_Sa_tot[y,j] <- (
        sum(Sa[y+kmin+1-1,1,1:ns,1:no,j]) +       # age 3 adults produced by spawners in brood year y
          sum(Sa[y+kmin+2-1,2,1:ns,1:no,j]) +     # age 4 adults produced by spawners in brood year y
          sum(Sa[y+kmin+3-1,3,1:ns,1:no,j]))/     # age 5 adults produced by spawners in brood year y
        Sa_tot[y,j]                               # total spawners in brood year y
    }
    
    # adults per smolt -- can't be calculated for all brood years in model
    for (y in (kmax+1):(ny-kmax)) {
      for (o in 1:no) {
        Ra_per_Ma[y,o,j] <- (
          sum(Ra[y+kmin+1-1,1,1:ns,o,j]) +       # age 3 adults back at trib that were once smolts from brood year y
            sum(Ra[y+kmin+2-1,2,1:ns,o,j]) +     # age 4 adults back at trib that were once smolts from brood year y
            sum(Ra[y+kmin+3-1,3,1:ns,o,j]))/     # age 5 adults back at trib that were once smolts from brood year y
          ifelse(sum(Ma[y,1:ni,o,j]) == 0, 1, sum(Ma[y,1:ni,o,j]))                    # total smolts at LGR from brood year y
      }
    }
  }
  
  ### OBSERVATION MODEL ###
  
  # age/sex composition
  for (j in 1:nj) {
    for (y in (kmax+1):ny) {
      # at weir
      weir_x_obs[y,1:nkso,j] ~ dmulti(q_Ra[y,1:nkso,j], weir_nx_obs[y,j])
      
      # for carcasses
      carc_x_obs[y,1:nkso,j] ~ dmulti(q_Sa_adj[y,1:nkso,j], carc_nx_obs[y,j])
    }
  }
  
  # the "fit" objects allow using square objects but skipping over missing obs in likelihoods
  # fall trap abundance
  for (d in 1:nfit_Pa) {
    Pa_obs[fit_Pa[d,1],fit_Pa[d,2],fit_Pa[d,3]] ~ dlnorm(log(Pa[fit_Pa[d,1],fit_Pa[d,2],fit_Pa[d,3]]), 1/sig_Pa_obs[fit_Pa[d,1],fit_Pa[d,2],fit_Pa[d,3]]^2)
  }
  
  # spring trap abundance
  for (d in 1:nfit_Mb) {
    Mb_obs[fit_Mb[d,1],fit_Mb[d,2],fit_Mb[d,3],fit_Mb[d,4]] ~ dlnorm(log(Mb[fit_Mb[d,1],fit_Mb[d,2],fit_Mb[d,3],fit_Mb[d,4]]), 1/sig_Mb_obs[fit_Mb[d,1],fit_Mb[d,2],fit_Mb[d,3],fit_Mb[d,4]]^2)
  }
  
  # adult abundance
  for (d in 1:nfit_Ra) {
    Ra_obs[fit_Ra[d,1],fit_Ra[d,2]] ~ dlnorm(log(Ra_tot[fit_Ra[d,1],fit_Ra[d,2]]), 1/sig_Ra_obs[fit_Ra[d,1],fit_Ra[d,2]]^2)
  }
  
  # summer tagging to LGD
  for (d in 1:nfit_Lphi_Pb_Ma) {
    Lphi_obs_Pb_Ma[fit_Lphi_Pb_Ma[d,1],fit_Lphi_Pb_Ma[d,2]] ~ dnorm(logit(phi_Pb_Ma[fit_Lphi_Pb_Ma[d,1],fit_Lphi_Pb_Ma[d,2]]), 1/sig_Lphi_obs_Pb_Ma[fit_Lphi_Pb_Ma[d,1],fit_Lphi_Pb_Ma[d,2]]^2)
  }
  
  # fall tagging to LGD
  for (d in 1:nfit_Lphi_Pa_Ma) {
    Lphi_obs_Pa_Ma[fit_Lphi_Pa_Ma[d,1],fit_Lphi_Pa_Ma[d,2],fit_Lphi_Pa_Ma[d,3]] ~ dnorm(logit(phi_Pa_Ma[fit_Lphi_Pa_Ma[d,1],fit_Lphi_Pa_Ma[d,2],fit_Lphi_Pa_Ma[d,3]]), 1/sig_Lphi_obs_Pa_Ma[fit_Lphi_Pa_Ma[d,1],fit_Lphi_Pa_Ma[d,2],fit_Lphi_Pa_Ma[d,3]]^2)
  }
  
  # spring tagging/smolt releases to LGD
  for (d in 1:nfit_Lphi_Mb_Ma) {
    Lphi_obs_Mb_Ma[fit_Lphi_Mb_Ma[d,1],fit_Lphi_Mb_Ma[d,2],fit_Lphi_Mb_Ma[d,3],fit_Lphi_Mb_Ma[d,4]] ~ dnorm(logit(phi_Mb_Ma[fit_Lphi_Mb_Ma[d,1],fit_Lphi_Mb_Ma[d,2],fit_Lphi_Mb_Ma[d,3],fit_Lphi_Mb_Ma[d,4]]), 1/sig_Lphi_obs_Mb_Ma[fit_Lphi_Mb_Ma[d,1],fit_Lphi_Mb_Ma[d,2],fit_Lphi_Mb_Ma[d,3],fit_Lphi_Mb_Ma[d,4]]^2)
  }
  
  # hydrosystem survival
  for (d in 1:nfit_Lphi_Ma_O0) {
    Lphi_obs_Ma_O0[fit_Lphi_Ma_O0[d,1],fit_Lphi_Ma_O0[d,2]] ~ dnorm(logit(phi_Ma_O0[fit_Lphi_Ma_O0[d,1],fit_Lphi_Ma_O0[d,2]]), 1/sig_Lphi_obs_Ma_O0[fit_Lphi_Ma_O0[d,1],fit_Lphi_Ma_O0[d,2]]^2)
  }
  
  # pre-spawn survival
  for (d in 1:nfit_spawned) {
    carcs_spawned[fit_spawned[d,1],fit_spawned[d,2]] ~ dbin(phi_Sb_Sa[fit_spawned[d,1],fit_spawned[d,2]], carcs_sampled[fit_spawned[d,1],fit_spawned[d,2]])
  }
  
  ### CALCULATE ALL PROCESS MODEL RESIDUALS ###
  
  for (y in (kmax+1):ny) {
    for (j in 1:nj) {
      
      # total summer parr recruitment
      lPb_resid[y,j] <- log(Pb[y,j]) - log(Pb_pred[y,j])
      
      # proportion of summer parr that become fall migrants
      Lpi_resid[y,j] <- logit(pi[y,i_fall,j]) - logit(mu_pi[i_fall,j])
      
      # LH-specific overwinter survival
      for (i in 1:ni) {
        Lphi_Pa_Mb_resid[y,i,j] <- Lphi_Pa_Mb[y,i,j] - (gamma0[i,j] + gamma1[i,j] * (Pa[y,i,j]/peu[j]))
      }
      
      # origin-specific quantities
      for (o in 1:no) {
        # movement survival to LGR
        Lphi_Mb_Ma_resid[y,o,j] <- Lphi_Mb_Ma[y,i_spring,o,j] - logit(mu_phi_Mb_Ma[i_spring,o,j])
        
        # probability of returning as female
        Lomega_resid[y,o,j] <- Lomega1[y,o,j] - logit(mu_omega[s_female,o,j])
        
        # sex-specific quantities: maturity
        for (s in 1:ns) {
          # pr(return at SWA1)
          Lpsi_O1_Rb_resid[y,s,o,j] <- Lpsi_O1_Rb[y,s,o,j] - logit(mu_psi_O1_Rb[s,o,j])
          
          # pr(return at SWA2|not returned at SWA1)
          Lpsi_O2_Rb_resid[y,s,o,j] <- Lpsi_O2_Rb[y,s,o,j] - logit(mu_psi_O2_Rb[s,o,j])
        }
      }
      
      # pre-spawn survival
      Lphi_Sb_Sa_resid[y,j] <- Lphi_Sb_Sa[y,j] - logit(mu_phi_Sb_Sa[j])
    }
    
    for (o in 1:no) {
      # movement survival from LGR thru BON
      Lphi_Ma_O0_resid[y,o] <- Lphi_Ma_O0[y,o] - logit(mu_phi_Ma_O0[o])
    }
    
    # natural origin ocean survival SWA0 -> SWA1
    Lphi_O0_O1_resid[y] <- Lphi_O0_O1[y,o_nat] - logit(mu_phi_O0_O1[o_nat])
    
    # natural origin ocean survival SWA1 -> SWA2
    Lphi_O1_O2_resid[y] <- Lphi_O1_O2[y,o_nat] - logit(mu_phi_O1_O2[o_nat])
  }
}
