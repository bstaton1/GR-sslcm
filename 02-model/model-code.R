jags_model_code = function() {
  
  ### PRIORS: capacity vs. weighted usable habitat length relationship ###
  lambda ~ dnorm(0, 1e-8) %_% T(0,)
  sig_lbeta ~ dunif(0, 5)
  
  for (j in 1:nj) {
    ### PRIORS: RECRUITMENT FUNCTION ###
    # total egg production to aggregate parr recruitment function
    alpha[j] ~ dbeta(1, 1)
    lbeta[j] ~ dnorm(log(lambda * wul[j]), 1/sig_lbeta^2) %_% T(,15) # log capacity. bound to prevent nonsensically large draws
    beta[j] <- exp(lbeta[j])
    sig_Lphi_E_Pb[j] ~ dunif(0, 5)
    lambda_pop[j] <- beta[j]/wul[j]  # derived pop-specific lambda, includes process noise in regression

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
    mu_phi_Mb_Ma[i_spring,o_nor,j] ~ dbeta(1, 1)
    sig_Lphi_Mb_Ma[o_nor,j] ~ dunif(0, 5)
    mu_phi_Mb_Ma[i_fall,o_nor,j] <- mu_phi_Mb_Ma[i_spring,o_nor,j]

    # hatchery origin movement survival (trib to LGD): have spring migrants only
    mu_phi_Mb_Ma[i_spring,o_hor,j] ~ dbeta(1, 1)
    sig_Lphi_Mb_Ma[o_hor,j] ~ dunif(0, 5)
    
    # pre-spawn survival (after brood-stock removal to successful spawning)
    mu_phi_Sb_Sa[j] ~ dbeta(1, 1)
    sig_Lphi_Sb_Sa[j] ~ dunif(0, 5)
    
    ### PRIORS: MATURATION ###
    # age/origin-specific maturation probabilities
    for (o in 1:no) {
      
      # pr(return at SWA1)
      mu_psi_O1[o,j] ~ dbeta(1, 1)
      sig_Lpsi_O1[o,j] ~ dunif(0, 5)
      
      # pr(return at SWA2|not returned at SWA1)
      mu_psi_O2[o,j] ~ dbeta(1, 1)
      sig_Lpsi_O2[o,j] ~ dunif(0, 5)
    }
    
    ### PRIORS: OCEAN SURVIVAL ###
    # mean survival by ocean year transition for natural origin
    mu_phi_O0_O1[o_nor,j] ~ dbeta(1,1)               # first winter at sea: to become SWA1
    mu_phi_O1_O2[o_nor,j] ~ dbeta(60,15)             # second winter at sea: to become SWA2
    mu_phi_O2_O3[o_nor,j] <- mu_phi_O1_O2[o_nor,j]   # third winter at sea: to become SW3
    
    # log odds ratio between natural and hatchery origin
    delta[j] ~ dt(0, 1/1.566^2, 7.763)
    
    # AR(1) coefficient for first year ocean survival
    kappa_phi_O0_O1[j] ~ dunif(-0.99,0.99)
    
    # standard deviation of white noise residuals
    sig_Lphi_O0_O1[j] ~ dunif(0, 5)
    
    # standard deviation of total year 1 residual (i.e., includes autocorrelation)
    sig_Lphi_O0_O1_init[j] <- sqrt(sig_Lphi_O0_O1[j]^2/(1 - kappa_phi_O0_O1[j]^2))
    
    # mean survival by ocean year transition for hatchery origin
    logit(mu_phi_O0_O1[o_hor,j]) <- logit(mu_phi_O0_O1[o_nor,j]) + delta[j]
    logit(mu_phi_O1_O2[o_hor,j]) <- logit(mu_phi_O1_O2[o_nor,j]) + delta[j]
    logit(mu_phi_O2_O3[o_hor,j]) <- logit(mu_phi_O2_O3[o_nor,j]) + delta[j]
    
    ### PRIORS: HATCHERY STRAYS RETURNING IN YEARS WITH NO ASSOCIATED SMOLT RELEASE ###
    
    # age composition of strays: only estimate for hatchery origin
    p_G[1:nk,o_nor,j] <- rep(0, nk)
    p_G[1:nk,o_hor,j] ~ ddirich(rep(1, nk))
    
    # the number of hatchery strays: only estimate in years where no other mechanism for generating hatchery fish
    for (i in 1:n_stray_yrs[j]) {
      G[stray_yrs[i,j],o_nor,j] <- 0
      G_random1[stray_yrs[i,j],o_hor,j] ~ dunif(0, 500)  # prior if an observed year
      G_random2[stray_yrs[i,j],o_hor,j] ~ dunif(50, 150) # prior if a simulated year, only applies for MIN
      G[stray_yrs[i,j],o_hor,j] <- ifelse(stray_yrs[i,j] <= ny_obs, G_random1[stray_yrs[i,j],o_hor,j], G_random2[stray_yrs[i,j],o_hor,j])
    }
    
    # force zero strays in years where they aren't needed
    for (o in 1:no) {
      for (i in 1:n_not_stray_yrs[j]) {
        G[not_stray_yrs[i,j],o,j] <- 0
      }
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
  rho_Lphi_E_Pb <- 0
  rho_Lpi <- 0
  rho_Lphi_Pa_Mb[i_fall] <- 0
  rho_Lphi_Pa_Mb[i_spring] <- 0
  rho_Lphi_Mb_Ma[o_nor] <- 0
  rho_Lphi_Mb_Ma[o_hor] <- 0
  rho_Lphi_Ma_O0 <- 0
  rho_Lphi_O0_O1 <- 0
  rho_Lpsi_O1[o_nor] <- 0
  rho_Lpsi_O1[o_hor] <- 0
  rho_Lpsi_O2[o_nor] <- 0
  rho_Lpsi_O2[o_hor] <- 0
  rho_Lphi_Rb_Ra <- 0
  rho_Lphi_Sb_Sa <- 0
  
  # construct all population covariance matrices
  for (i in 1:nj) {
    for (j in 1:nj) {
      # covariance matrix of parr recruitment
      Sig_Lphi_E_Pb[i,j] <- sig_Lphi_E_Pb[i] * sig_Lphi_E_Pb[j] * ifelse(i == j, 1, rho_Lphi_E_Pb)
      
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
        Sig_Lpsi_O1[i,j,o] <- sig_Lpsi_O1[o,i] * sig_Lpsi_O1[o,j] * ifelse(i == j, 1, rho_Lpsi_O1[o])
        
        # covariance matrices of pr(mature at SWA1)
        Sig_Lpsi_O2[i,j,o] <- sig_Lpsi_O2[o,i] * sig_Lpsi_O2[o,j] * ifelse(i == j, 1, rho_Lpsi_O2[o])
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
  Lphi_O0_O1_resid[1,o_nor,1:nj] ~ dmnorm.vcov(rep(0, nj), Sig_Lphi_O0_O1_init[1:nj,1:nj])

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
      log(zeta[k,j_z[j]]) <- z[1,j_z[j]] * age3[k] + z[2,j_z[j]] * age5[k]
    }
  }
  
  # calculate correction factor for population(s) without both weir and carcass data (Minam)
  zeta[1,3] <- (zeta[1,1] + zeta[1,2] + zeta[1,4])/3
  zeta[2,3] <- (zeta[2,1] + zeta[2,2] + zeta[2,4])/3
  zeta[3,3] <- (zeta[3,1] + zeta[3,2] + zeta[3,4])/3
  
  ### PRIORS: BROOD-YEAR-SPECIFIC PARAMETERS ###
  for (y in 2:ny) {
    
    # egg to parr survival: Beverton-Holt relationship
    Lphi_E_Pb[y,1:nj] ~ dmnorm.vcov(logit(1/(1/alpha[1:nj] + E[y,1:nj]/beta[1:nj])), Sig_Lphi_E_Pb[1:nj,1:nj])

    # LH apportionment
    Lpi1[y,1:nj] ~ dmnorm.vcov(logit(mu_pi[i_fall,1:nj]), Sig_Lpi[1:nj,1:nj])

    # overwinter survival by LH type: density dependent
    Lphi_Pa_Mb[y,i_fall,1:nj] ~ dmnorm.vcov(gamma0[i_fall,1:nj] + gamma1[i_fall,1:nj] * (Pa[y,i_fall,1:nj]/wul[1:nj]), Sig_Lphi_Pa_Mb[1:nj,1:nj,i_fall])
    Lphi_Pa_Mb[y,i_spring,1:nj] ~ dmnorm.vcov(gamma0[i_spring,1:nj] + gamma1[i_spring,1:nj] * (Pa[y,i_spring,1:nj]/wul[1:nj]), Sig_Lphi_Pa_Mb[1:nj,1:nj,i_spring])

    for (o in 1:no) {
      # migration survival from trib to LGR by origin
      Lphi_Mb_Ma[y,i_spring,o,1:nj] ~ dmnorm.vcov(logit(mu_phi_Mb_Ma[i_spring,o,1:nj]), Sig_Lphi_Mb_Ma[1:nj,1:nj,o])

      # pr(return at SWA1)
      Lpsi_O1[y,o,1:nj] ~ dmnorm.vcov(logit(mu_psi_O1[o,1:nj]), Sig_Lpsi_O1[1:nj,1:nj,o])

      # pr(return at SWA2)
      Lpsi_O2[y,o,1:nj] ~ dmnorm.vcov(logit(mu_psi_O2[o,1:nj]), Sig_Lpsi_O2[1:nj,1:nj,o])
    }
    
    # yr1 NOR ocean survival: includes AR(1) process
    Lphi_O0_O1[y,o_nor,1:nj] ~ dmnorm.vcov(logit(mu_phi_O0_O1[o_nor,1:nj]) + Lphi_O0_O1_resid[y-1,o_nor,1:nj] * kappa_phi_O0_O1[1:nj], Sig_Lphi_O0_O1[1:nj,1:nj])

    # yr2/yr3 NOR ocean survival: time constant
    Lphi_O1_O2[y,o_nor,1:nj] <- logit(mu_phi_O1_O2[o_nor,1:nj])
    Lphi_O2_O3[y,o_nor,1:nj] <- logit(mu_phi_O2_O3[o_nor,1:nj])
    
    # yr1/yr2/yr3 HOR ocean survival: same as NOR but adjusted by a time-constant log odds ratio
    Lphi_O0_O1[y,o_hor,1:nj] <- Lphi_O0_O1[y,o_nor,1:nj] + delta[1:nj]
    Lphi_O1_O2[y,o_hor,1:nj] <- Lphi_O1_O2[y,o_nor,1:nj] + delta[1:nj]
    Lphi_O2_O3[y,o_hor,1:nj] <- Lphi_O2_O3[y,o_nor,1:nj] + delta[1:nj]
    
    # pre-spawn survival
    Lphi_Sb_Sa[y,1:nj] ~ dmnorm.vcov(logit(mu_phi_Sb_Sa[1:nj]), Sig_Lphi_Sb_Sa[1:nj,1:nj])

    # movement survival (juveniles LGR to BON)
    Lphi_Ma_O0[y,1:no] ~ dmnorm.vcov(logit(mu_phi_Ma_O0[1:no]), Sig_Lphi_Ma_O0[1:no,1:no])

    # movement survival (adults BON to LGR)
    Lphi_Rb_Ra_random[y,1:no] ~ dmnorm.vcov(logit(mu_phi_Rb_Ra[1:no]), Sig_Lphi_Rb_Ra[1:no,1:no])

    # apply inverse logit-transformation to all of these brood year-specific parameters
    # ilogit() and exp() cannot be applied to vectors unfortunately; must loop over population dimension
    for (j in 1:nj) {
      
      # transform egg to parr survival
      phi_E_Pb[y,j] <- ilogit(Lphi_E_Pb[y,j])

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
        psi_O1[y,o,j] <- ilogit(Lpsi_O1[y,o,j])
        
        # transform pr(return at SWA2|not returned at SWA1)
        psi_O2[y,o,j] <- ilogit(Lpsi_O2[y,o,j])
        
        # transform ocean survival terms
        phi_O0_O1[y,o,j] <- ilogit(Lphi_O0_O1[y,o,j])
        phi_O1_O2[y,o,j] <- ilogit(Lphi_O1_O2[y,o,j])
        phi_O2_O3[y,o,j] <- ilogit(Lphi_O2_O3[y,o,j])
      }
      
      # assume movement survival trib to LGR for NOR fish is equal between LH types
      phi_Mb_Ma[y,i_fall,o_nor,j] <- phi_Mb_Ma[y,i_spring,o_nor,j]
      
      # pre-spawn survival
      phi_Sb_Sa[y,j] <- ilogit(Lphi_Sb_Sa[y,j])
    }
    
    # transform quantities common to all pops but different by origin
    for (o in 1:no) {
      # transform movement survival (juveniles LGR to BON)
      phi_Ma_O0[y,o] <- ilogit(Lphi_Ma_O0[y,o])
      
      # transform movement survival (adults BON to LGR)
      Lphi_Rb_Ra[y,o] <- ifelse(y < first_x_LGR, logit(mu_phi_Rb_Ra[o]), Lphi_Rb_Ra_random[y,o])
      phi_Rb_Ra[y,o] <- ilogit(Lphi_Rb_Ra[y,o])
    }
  }
  
  ### PROCESS MODEL: INITIALIZATION ###
  # need to do something special with first kmax brood years prior to data collection
  # to populate adult states for reproduction/observation in the first years with data
  # THIS APPROACH: estimate year/age-specific NOR return to Columbia R. mouth
  # there are 12 year/age-specific returns that are not linked to pop dyn in the model
  # this approach estimates (with a fairly restrictive prior) the cells with "??"
  # NA cells are a placeholder year for y = 0 residuals in AR(1) processes
  # do NOR here only, HOR returns in these years handled by "straying" model
  #
  # ---------------------------
  #   year    age3    age4 age5
  # ---------------------------
  #   1990      NA      NA   NA
  #   1991      ??      ??   ??
  #   1992      ??      ??   ??
  #   1993      ??      ??   ??
  #   1994 modeled      ??   ??
  #   1995 modeled modeled   ??
  # ---------------------------

  for (j in 1:nj) {
    for (k in 1:nk) {
      for (y in 2:(kmin+k)) {
        Rb[y,k,o_nor,j] ~ dunif(0, max_Rb_init[k])
        Rb[y,k,o_hor,j] <- 0
      }
    }
  }
  
  ### PROCESS MODEL: COMPLETE LIFE CYCLE FOR REMAINING BROOD YEARS ###
  for (j in 1:nj) {
    for (y in 2:ny) {
      # reproductive link: expected egg-to-parr survival is a BH function and parr recruitment is egg production times a latent survival
      Pb[y,j] <- E[y,j] * phi_E_Pb[y,j]

      # natural origin tributary-to-LGD dynamics
      for (i in 1:ni) {
        # apportion to LH-strategy: parr after fall migration
        Pa[y,i,j] <- Pb[y,j] * pi[y,i,j]
        
        # survive over winter: smolt before spring migration
        Mb[y,i,o_nor,j] <- Pa[y,i,j] * phi_Pa_Mb[y,i,j]
        
        # move to LGD: smolt after spring migration, at top of LGD
        Ma[y,i,o_nor,j] <- Mb[y,i,o_nor,j] * phi_Mb_Ma[y,i,o_nor,j]

        # derived survival for fitting: fall trap to LGD
        phi_Pa_Ma[y,i,j] <- Ma[y,i,o_nor,j]/max(Pa[y,i,j] * phi_Pb_Pa[i,j], Ma[y,i,o_nor,j])
      }
      
      # flag that tells us how often max constraint is violated
      # will remove this eventually. the max constraint is needed for MCMC during early tuning (crashes w/o it)
      # but none of the converged samples have this occur
      bad_flag[y,j] <- ifelse(Pa[y,i_spring,j] * phi_Pb_Pa[i_spring,j] < Ma[y,i_spring,o_nor,j], 1, 0)
      
      # derived survival for fitting: summer tagging to LGD
      phi_Pb_Ma[y,j] <- sum(Ma[y,1:ni,o_nor,j])/Pb[y,j]
      
      # put hatchery smolts in tributary
      Mb[y,i_spring,o_hor,j] <- Mb_obs[y,i_spring,o_hor,j]
      
      # move hatchery smolts from tributary to LGD
      Ma[y,i_spring,o_hor,j] <- Mb[y,i_spring,o_hor,j] * phi_Mb_Ma[y,i_spring,o_hor,j]
      
      # create zeros for fall migrant hatchery fish at LGD
      # needed because we sum over this dimension below
      Ma[y,i_fall,o_hor,j] <- 0
      
      # origin-specific processes
      for (o in 1:no) {
        # move origin-specific smolts from LGD to estuary (ocean age 0)
        O0[y,o,j] <- sum(Ma[y,1:ni,o,j]) * phi_Ma_O0[y,o]
        
        # move juveniles through ocean ages and survivals
        O[y,1,o,j] <- O0[y,o,j] * phi_O0_O1[y,o,j] # survive first winter at sea. now SWA1, TA3
        O[y,2,o,j] <- O[y,1,o,j] * (1 - psi_O1[y,o,j]) * phi_O1_O2[y,o,j] # don't mature at SWA1 and survive second winter at sea. now SWA2, TA4
        O[y,3,o,j] <- O[y,2,o,j] * (1 - psi_O2[y,o,j]) * phi_O2_O3[y,o,j] # don't mature at SWA2 and survive third winter at sea. now SWA3, TA5
        
        # mature and return to river in appropriate year at age
        Rb[y+kmin+1-1,1,o,j] <- O[y,1,o,j] * psi_O1[y,o,j]
        Rb[y+kmin+2-1,2,o,j] <- O[y,2,o,j] * psi_O2[y,o,j]
        Rb[y+kmin+3-1,3,o,j] <- O[y,3,o,j]
        
        # adult in-river processes: all age/origin specific
        # for these adult stages, y represents the brood year fish returned in
        for (k in 1:nk) {
          # returning adults making it to BON (survive sea lions * survive fisheries downstream of BON)
          Rb_BON[y,k,o,j] <- Rb[y,k,o,j] * phi_SL[y,j] * (1 - U[y,k,o])
          
          # survive upstream migration from BON to LGR
          Ra_LGR[y,k,o,j] <- Rb_BON[y,k,o,j] * phi_Rb_Ra[y,o]
          
          # add strays
          Ra[y,k,o,j] <- Ra_LGR[y,k,o,j] + (G[y,o,j] * p_G[k,o,j])

          # remove fish at weir: use max() to ensure that fewer fish were removed than existed
          Sb[y,k,o,j] <- max(Ra[y,k,o,j] - B[y,k,o,j], 1)
          
          # survive pre-spawn mortality
          Sa[y,k,o,j] <- Sb[y,k,o,j] * phi_Sb_Sa[y,j]
          
          # calculate egg production (eggs separated by age/origin)
          E_sep[y,k,o,j] <- Sa[y,k,o,j] * Omega[k,j] * f[k]
          
          # calculate "adjusted carcasses": accounts for sampling bias relative to weir
          Sa_prime[y,k,o,j] <- Sa[y,k,o,j] * zeta[k,j]
        }
      }
      
      # total adults returned to trib
      Ra_tot[y,j] <- sum(Ra[y,1:nk,1:no,j])
      
      # total spawning adults
      Sa_tot[y,j] <- sum(Sa[y,1:nk,1:no,j])
      
      # total egg production summed across ages/origins: used as spawning stock index
      E[y,j] <- sum(E_sep[y,1:nk,1:no,j])
      
      # reformat returns by age/origin for fitting to weir comp data (basically cbind two array slices)
      Ra_2d[y,1:nk,j] <- Ra[y,1:nk,o_nor,j]            # nat., all ages
      Ra_2d[y,(nk+1):(2*nk),j] <- Ra[y,1:nk,o_hor,j]   # hat., all ages

      # reformat "adjusted carcasses" by age/origin for fitting to carcass comp data (basically cbind two array slices)
      Sa_prime_2d[y,1:nk,j] <- Sa_prime[y,1:nk,o_nor,j]            # nat., all ages
      Sa_prime_2d[y,(nk+1):(2*nk),j] <- Sa_prime[y,1:nk,o_hor,j]   # hat., all ages

      # calculate age compositions
      for (ko in 1:nko) {
        p_Ra[y,ko,j] <- Ra_2d[y,ko,j]/sum(Ra_2d[y,1:nko,j])
        p_Sa_prime[y,ko,j] <- Sa_prime_2d[y,ko,j]/sum(Sa_prime_2d[y,1:nko,j])
      }
      
      # calculate misc derived quantities
      Pb_per_Sa_tot[y,j] <- Pb[y,j]/Sa_tot[y,j]                  # parr per spawner
      Pb_per_E[y,j] <- Pb[y,j]/E[y,j]                    # parr per egg
      Mb_per_Sa_tot[y,j] <- sum(Mb[y,1:ni,o_nor,j])/Sa_tot[y,j]  # smolt per spawner
    }
    
    # BON -> BON survival -- can't be calculated for all brood years in model
    for (y in 2:(ny-kmax)) {
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
    for (y in 2:(ny-kmax)) {
      Sa_tot_per_Sa_tot[y,j] <- (
        sum(Sa[y+kmin+1-1,1,1:no,j]) +       # age 3 adults produced by spawners in brood year y
          sum(Sa[y+kmin+2-1,2,1:no,j]) +     # age 4 adults produced by spawners in brood year y
          sum(Sa[y+kmin+3-1,3,1:no,j]))/     # age 5 adults produced by spawners in brood year y
        Sa_tot[y,j]                          # total spawners in brood year y
    }
    
    # adults per smolt -- can't be calculated for all brood years in model
    for (y in 2:(ny-kmax)) {
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
    for (y in 2:ny_obs) {
      # AT WEIR
      # data likelihood
      x_Ra[y,1:nko,j] ~ dmulti(p_Ra[y,1:nko,j], nx_Ra[y,j])
      
      # simulate new data
      x_Ra_new[y,1:nko,j] ~ dmulti(p_Ra[y,1:nko,j], nx_Ra[y,j])
      
      # calculate expected count add small number to avoid division by zero
      expected_x_Ra[y,1:nko,j] <- p_Ra[y,1:nko,j] * nx_Ra[y,j] + 1e-6
      
      # calculate fit statistic: chi-squared statistic
      x_Ra_dev[y,j] <- sum(((x_Ra[y,1:nko,j] - expected_x_Ra[y,1:nko,j])^2)/expected_x_Ra[y,1:nko,j])
      x_Ra_new_dev[y,j] <- sum(((x_Ra_new[y,1:nko,j] - expected_x_Ra[y,1:nko,j])^2)/expected_x_Ra[y,1:nko,j])
      
      # calculate log posterior predictive density
      x_Ra_lppd[y,j] <- logdensity.multi(x_Ra[y,1:nko,j], p_Ra[y,1:nko,j], nx_Ra[y,j])
      
      # FOR CARCASSES
      # data likelihood
      x_Sa_prime[y,1:nko,j] ~ dmulti(p_Sa_prime[y,1:nko,j], nx_Sa_prime[y,j])
      
      # simulate new data
      x_Sa_prime_new[y,1:nko,j] ~ dmulti(p_Ra[y,1:nko,j], nx_Sa_prime[y,j])
      
      # calculate expected count add small number to avoid division by zero
      expected_x_Sa_prime[y,1:nko,j] <- p_Ra[y,1:nko,j] * nx_Sa_prime[y,j] + 1e-6
      
      # calculate fit statistic: chi-squared statistic
      x_Sa_prime_dev[y,j] <- sum(((x_Sa_prime[y,1:nko,j] - expected_x_Sa_prime[y,1:nko,j])^2)/expected_x_Sa_prime[y,1:nko,j])
      x_Sa_prime_new_dev[y,j] <- sum(((x_Sa_prime_new[y,1:nko,j] - expected_x_Sa_prime[y,1:nko,j])^2)/expected_x_Sa_prime[y,1:nko,j])
      
      # calculate log posterior predictive density
      x_Sa_prime_lppd[y,j] <- logdensity.multi(x_Sa_prime[y,1:nko,j], p_Ra[y,1:nko,j], nx_Sa_prime[y,j])
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
  for (d in 1:nfit_x_LGR) {
    # data likelihood
    x_LGR[fit_x_LGR[d,1],fit_x_LGR[d,2]] ~ dbin(phi_Rb_Ra[fit_x_LGR[d,1],fit_x_LGR[d,2]], x_BON[fit_x_LGR[d,1],fit_x_LGR[d,2]])
    
    # simulate new data
    x_LGR_new[fit_x_LGR[d,1],fit_x_LGR[d,2]] ~ dbin(phi_Rb_Ra[fit_x_LGR[d,1],fit_x_LGR[d,2]], x_BON[fit_x_LGR[d,1],fit_x_LGR[d,2]])
    
    # calculate expected count
    expected_x_LGR[fit_x_LGR[d,1],fit_x_LGR[d,2]] <- phi_Rb_Ra[fit_x_LGR[d,1],fit_x_LGR[d,2]] * x_BON[fit_x_LGR[d,1],fit_x_LGR[d,2]]
    
    # calculate fit statistic: chi-squared statistic
    x_LGR_dev[fit_x_LGR[d,1],fit_x_LGR[d,2]] <- ((x_LGR[fit_x_LGR[d,1],fit_x_LGR[d,2]] - expected_x_LGR[fit_x_LGR[d,1],fit_x_LGR[d,2]])^2)/expected_x_LGR[fit_x_LGR[d,1],fit_x_LGR[d,2]]
    x_LGR_new_dev[fit_x_LGR[d,1],fit_x_LGR[d,2]] <- ((x_LGR_new[fit_x_LGR[d,1],fit_x_LGR[d,2]] - expected_x_LGR[fit_x_LGR[d,1],fit_x_LGR[d,2]])^2)/expected_x_LGR[fit_x_LGR[d,1],fit_x_LGR[d,2]]
    
    # calculate log posterior predictive density
    x_LGR_lppd[fit_x_LGR[d,1],fit_x_LGR[d,2]] <- logdensity.bin(x_LGR[fit_x_LGR[d,1],fit_x_LGR[d,2]], phi_Rb_Ra[fit_x_LGR[d,1],fit_x_LGR[d,2]], x_BON[fit_x_LGR[d,1],fit_x_LGR[d,2]])
  }
  
  # pre-spawn survival
  for (d in 1:nfit_x_carcass_spawned) {
    # data likelihood
    x_carcass_spawned[fit_x_carcass_spawned[d,1],fit_x_carcass_spawned[d,2]] ~ dbin(phi_Sb_Sa[fit_x_carcass_spawned[d,1],fit_x_carcass_spawned[d,2]], x_carcass_total[fit_x_carcass_spawned[d,1],fit_x_carcass_spawned[d,2]])
    
    # simulate new data
    x_carcass_spawned_new[fit_x_carcass_spawned[d,1],fit_x_carcass_spawned[d,2]] ~ dbin(phi_Sb_Sa[fit_x_carcass_spawned[d,1],fit_x_carcass_spawned[d,2]], x_carcass_total[fit_x_carcass_spawned[d,1],fit_x_carcass_spawned[d,2]])
    
    # calculate expected count
    expected_x_carcass_spawned[fit_x_carcass_spawned[d,1],fit_x_carcass_spawned[d,2]] <- phi_Sb_Sa[fit_x_carcass_spawned[d,1],fit_x_carcass_spawned[d,2]] * x_carcass_total[fit_x_carcass_spawned[d,1],fit_x_carcass_spawned[d,2]]
    
    # calculate fit statistic: chi-squared statistic
    x_carcass_spawned_dev[fit_x_carcass_spawned[d,1],fit_x_carcass_spawned[d,2]] <- ((x_carcass_spawned[fit_x_carcass_spawned[d,1],fit_x_carcass_spawned[d,2]] - expected_x_carcass_spawned[fit_x_carcass_spawned[d,1],fit_x_carcass_spawned[d,2]])^2)/expected_x_carcass_spawned[fit_x_carcass_spawned[d,1],fit_x_carcass_spawned[d,2]]
    x_carcass_spawned_new_dev[fit_x_carcass_spawned[d,1],fit_x_carcass_spawned[d,2]] <- ((x_carcass_spawned_new[fit_x_carcass_spawned[d,1],fit_x_carcass_spawned[d,2]] - expected_x_carcass_spawned[fit_x_carcass_spawned[d,1],fit_x_carcass_spawned[d,2]])^2)/expected_x_carcass_spawned[fit_x_carcass_spawned[d,1],fit_x_carcass_spawned[d,2]]
    
    # calculate log posterior predictive density
    x_carcass_spawned_lppd[fit_x_carcass_spawned[d,1],fit_x_carcass_spawned[d,2]] <- logdensity.bin(x_carcass_spawned[fit_x_carcass_spawned[d,1],fit_x_carcass_spawned[d,2]], phi_Sb_Sa[fit_x_carcass_spawned[d,1],fit_x_carcass_spawned[d,2]], x_carcass_total[fit_x_carcass_spawned[d,1],fit_x_carcass_spawned[d,2]])
  }
  
  ### CALCULATE ALL PROCESS MODEL RESIDUALS ###
  
  for (y in 2:ny) {
    for (j in 1:nj) {
      
      # total summer parr recruitment
      Lphi_E_Pb_resid[y,j] <- Lphi_E_Pb[y,j] - logit(1/(1/alpha[j] + E[y,j]/beta[j]))

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
        Lpsi_O1_resid[y,o,j] <- Lpsi_O1[y,o,j] - logit(mu_psi_O1[o,j])
        
        # pr(return at SWA2|not returned at SWA1)
        Lpsi_O2_resid[y,o,j] <- Lpsi_O2[y,o,j] - logit(mu_psi_O2[o,j])
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
