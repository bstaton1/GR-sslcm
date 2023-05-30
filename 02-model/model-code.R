jags_model_code = function() {
  
  ### PRIORS: capacity vs. weighted usable habitat length relationship ###
  lambda ~ dnorm(0, 1e-8) %_% T(0,)
  lambda_pop[1:nj] <- beta[1:nj]/wul[1:nj]
  sig_lbeta ~ dunif(0, 5)
  
  ### PRIORS: Precision Matrices for Process Noise
  
  # egg to parr survival
  Tau_Lphi_E_Pb[1:nj,1:nj] ~ dscaled.wishart(rep(Tau_Lphi_E_Pb_prior[1], nj), Tau_Lphi_E_Pb_prior[2])
  Tau_Lphi_E_Pb_pr[1:nj,1:nj] ~ dscaled.wishart(rep(Tau_Lphi_E_Pb_prior[1], nj), Tau_Lphi_E_Pb_prior[2])
  
  # parr mean length
  Tau_lL_Pb[1:nj,1:nj] ~ dscaled.wishart(rep(Tau_lL_Pb_prior[1], nj), Tau_lL_Pb_prior[2])
  Tau_lL_Pb_pr[1:nj,1:nj] ~ dscaled.wishart(rep(Tau_lL_Pb_prior[1], nj), Tau_lL_Pb_prior[2])

  # proportion of parr that are fall migrants
  Tau_Lpi[1:nj,1:nj] ~ dscaled.wishart(rep(Tau_Lpi_prior[1], nj), Tau_Lpi_prior[2])
  Tau_Lpi_pr[1:nj,1:nj] ~ dscaled.wishart(rep(Tau_Lpi_prior[1], nj), Tau_Lpi_prior[2])

  # overwinter survival
  Tau_Lphi_Pa_Mb[1:nj,1:nj] ~ dscaled.wishart(rep(Tau_Lphi_Pa_Mb_prior[1], nj), Tau_Lphi_Pa_Mb_prior[2])
  Tau_Lphi_Pa_Mb_pr[1:nj,1:nj] ~ dscaled.wishart(rep(Tau_Lphi_Pa_Mb_prior[1], nj), Tau_Lphi_Pa_Mb_prior[2])

  # overwinter change in mean size
  Tau_lDelta_L_Pb_Mb[1:nj,1:nj] ~ dscaled.wishart(rep(Tau_lDelta_L_Pb_Mb_prior[1], nj), Tau_lDelta_L_Pb_Mb_prior[2])
  Tau_lDelta_L_Pb_Mb_pr[1:nj,1:nj] ~ dscaled.wishart(rep(Tau_lDelta_L_Pb_Mb_prior[1], nj), Tau_lDelta_L_Pb_Mb_prior[2])

  # migration survival to LGR
  Tau_Lphi_Mb_Ma[1:nj,1:nj] ~ dscaled.wishart(rep(Tau_Lphi_Mb_Ma_prior[1], nj), Tau_Lphi_Mb_Ma_prior[2])
  Tau_Lphi_Mb_Ma_pr[1:nj,1:nj] ~ dscaled.wishart(rep(Tau_Lphi_Mb_Ma_prior[1], nj), Tau_Lphi_Mb_Ma_prior[2])

  # migration survival from LGR to ocean
  Tau_Lphi_Ma_O0[1:no,1:no] ~ dscaled.wishart(rep(Tau_Lphi_Ma_O0_prior[1], no), Tau_Lphi_Ma_O0_prior[2])
  Tau_Lphi_Ma_O0_pr[1:no,1:no] ~ dscaled.wishart(rep(Tau_Lphi_Ma_O0_prior[1], no), Tau_Lphi_Ma_O0_prior[2])
  
  # Yr1 ocean survival
  Tau_Lphi_O0_O1[1:nj,1:nj] ~ dscaled.wishart(rep(Tau_Lphi_O0_O1_prior[1], nj), Tau_Lphi_O0_O1_prior[2])
  Tau_Lphi_O0_O1_pr[1:nj,1:nj] ~ dscaled.wishart(rep(Tau_Lphi_O0_O1_prior[1], nj), Tau_Lphi_O0_O1_prior[2])
  
  # Pr(mature after 1 year at sea)
  Tau_Lpsi_O1[1:nj,1:nj] ~ dscaled.wishart(rep(Tau_Lpsi_O1_prior[1], nj), Tau_Lpsi_O1_prior[2])
  Tau_Lpsi_O1_pr[1:nj,1:nj] ~ dscaled.wishart(rep(Tau_Lpsi_O1_prior[1], nj), Tau_Lpsi_O1_prior[2])
  
  # Pr(mature after 2 years at sea)
  Tau_Lpsi_O2[1:nj,1:nj] ~ dscaled.wishart(rep(Tau_Lpsi_O2_prior[1], nj), Tau_Lpsi_O2_prior[2])
  Tau_Lpsi_O2_pr[1:nj,1:nj] ~ dscaled.wishart(rep(Tau_Lpsi_O2_prior[1], nj), Tau_Lpsi_O2_prior[2])
  
  # migration survival BON to LGR
  Tau_Lphi_Rb_Ra[1:no,1:no] ~ dscaled.wishart(rep(Tau_Lphi_Rb_Ra_prior[1], no), Tau_Lphi_Rb_Ra_prior[2])
  Tau_Lphi_Rb_Ra_pr[1:no,1:no] ~ dscaled.wishart(rep(Tau_Lphi_Rb_Ra_prior[1], no), Tau_Lphi_Rb_Ra_prior[2])
  
  for (j in 1:nj) {
    ### PRIORS: RECRUITMENT FUNCTION ###
    # total egg production to aggregate parr recruitment function
    # max egg-to-parr survival
    alpha[j] ~ dbeta(alpha_prior[1], alpha_prior[2])
    
    # max parr recruitment
    lbeta[j] ~ dnorm(log(lambda * wul[j]), 1/sig_lbeta^2) %_% T(,15) # log capacity. bound to prevent nonsensically large draws
    beta[j] <- exp(lbeta[j])
    
    # AR(1) coefficient for egg-to-parr survival
    kappa_phi_E_Pb[j] ~ dunif(-0.99,0.99)
    
    ### PRIORS: DENSITY-DEPENDENT GROWTH ###
    # length at end of summer vs. eggs relationship: density-dependent spring/summer growth
    omega0[j] ~ dunif(omega0_prior[1], omega0_prior[2])  # intercept
    omega1[j] ~ dunif(omega1_prior[1], omega1_prior[2])  # slope

    ### PRIORS: FRESHWATER PARAMETERS ###
    # aggregate parr to LH-specific parr
    mu_pi[i_fall,j] ~ dbeta(1, 1)
    mu_pi[i_spring,j] <- 1 - mu_pi[i_fall,j]
    
    # overwinter survival parameters: LH-specific
    # uses logistic relationship to introduce size-dependence
    # gamma0: LH-specific intercepts
    # gamma1: LH-common slopes
    for (i in 1:ni) {
      gamma0[i,j] ~ dt(0, 1/1.566^2, 7.763)
    }
    gamma1[i_fall,j] ~ dt(0, 1/1.566^2, 7.763)
    gamma1[i_spring,j] <- gamma1[i_fall,j]

    # coefficients for obtaining summer mean length to spring mean length growth factor
    theta0[j] ~ dunif(theta0_prior[1], theta0_prior[2])
    theta1[j] ~ dunif(theta1_prior[1], theta1_prior[2])
    
    # size-dependent NOR migration to LGR survival
    # assumed equal among migration types
    tau0[j] ~ dt(0, 1/1.566^2, 7.763)
    tau1[j] ~ dt(0, 1/1.566^2, 7.763)
    
    # hatchery origin movement survival (trib to LGD): have spring migrants only
    mu_phi_Mb_Ma[i_spring,o_hor,j] ~ dbeta(1, 1)
    
    # pre-spawn survival (after brood-stock removal to successful spawning)
    mu_phi_Sb_Sa[j] ~ dbeta(1, 1)
    
    ### PRIORS: MATURATION ###
    # age/origin-specific maturation probabilities
    for (o in 1:no) {
      # pr(return at SWA1)
      mu_psi_O1[o,j] ~ dbeta(mu_psi_O1_prior[1], mu_psi_O1_prior[2])

      # pr(return at SWA2|not returned at SWA1)
      mu_psi_O2[o,j] ~ dbeta(mu_psi_O2_prior[1], mu_psi_O2_prior[2])
    }
    
    ### PRIORS: OCEAN SURVIVAL ###
    # mean survival by ocean year transition for natural origin
    mu_phi_O0_O1[o_nor,j] ~ dbeta(mu_phi_O0_O1_prior[1], mu_phi_O0_O1_prior[2]) # first winter at sea: to become SWA1
    mu_phi_O1_O2[o_nor,j] ~ dbeta(mu_phi_O1_O2_prior[1], mu_phi_O1_O2_prior[2]) # second winter at sea: to become SWA2
    mu_phi_O2_O3[o_nor,j] ~ dbeta(mu_phi_O2_O3_prior[1], mu_phi_O2_O3_prior[2]) # third winter at sea: to become SW3

    # log odds ratio between natural and hatchery origin
    delta_O0_O1[j] ~ dt(0, 1/1.566^2, 7.763)
    delta_O1_O2[j] <- delta_O0_O1[j]
    delta_O2_O3[j] <- delta_O0_O1[j]
    
    # AR(1) coefficient for first year ocean survival
    kappa_phi_O0_O1[j] ~ dunif(-0.99,0.99)
    
    # mean survival by ocean year transition for hatchery origin
    logit(mu_phi_O0_O1[o_hor,j]) <- logit(mu_phi_O0_O1[o_nor,j]) + delta_O0_O1[j]
    logit(mu_phi_O1_O2[o_hor,j]) <- logit(mu_phi_O1_O2[o_nor,j]) + delta_O1_O2[j]
    logit(mu_phi_O2_O3[o_hor,j]) <- logit(mu_phi_O2_O3[o_nor,j]) + delta_O2_O3[j]
    
    ### PRIORS: HATCHERY STRAYS RETURNING IN YEARS WITH NO ASSOCIATED SMOLT RELEASE ###
    
    # age composition of strays: only estimate for hatchery origin
    p_G[1:nk,o_nor,j] <- rep(0, nk)
    p_G[1:nk,o_hor,j] ~ ddirich(rep(1, nk))
    
    # the number of hatchery strays: only estimate in years where no other mechanism for generating hatchery fish
    for (i in 1:n_stray_yrs[j]) {
      G[stray_yrs[i,j],o_nor,j] <- 0
      G_random1[stray_yrs[i,j],o_hor,j] ~ dunif(0, 500)  # prior if an observed year
      G_random2[stray_yrs[i,j],o_hor,j] ~ dunif(20,110) # prior if a simulated year, only applies for MIN
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
  }

  # migration survival adults from BON to LGR
  for (o in 1:no) {
    mu_phi_Rb_Ra[o] ~ dbeta(1, 1)
    # sig_Lphi_Rb_Ra[o] ~ dunif(0, 1)
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
    phi_E_Pb_dot[y,1:nj] <- 1/(1/alpha[1:nj] + E[y,1:nj]/beta[1:nj])
    Lphi_E_Pb[y,1:nj] ~ dmnorm(logit(phi_E_Pb_dot[y,1:nj]) + Lphi_E_Pb_resid[y-1,1:nj] * kappa_phi_E_Pb[1:nj], Tau_Lphi_E_Pb[1:nj,1:nj])

    # density-dependent length at end of summer
    lL_Pb_dot[y,1:nj] <- omega0[1:nj] + omega1[1:nj] * log((E[y,1:nj]/E_scale)/wul[1:nj])
    lL_Pb[y,1:nj] ~ dmnorm(lL_Pb_dot[y,1:nj], Tau_lL_Pb[1:nj,1:nj])
    
    # LH apportionment
    Lpi1[y,1:nj] ~ dmnorm(logit(mu_pi[i_fall,1:nj]), Tau_Lpi[1:nj,1:nj])

    # overwinter survival by LH type: function of size
    for (i in 1:ni) {
      Lphi_Pa_Mb_dot[y,i,1:nj] <- gamma0[i,1:nj] + gamma1[i,1:nj] * L_Pb_star[y,1:nj]
      Lphi_Pa_Mb[y,i,1:nj] ~ dmnorm(Lphi_Pa_Mb_dot[y,i,1:nj], Tau_Lphi_Pa_Mb[1:nj,1:nj])
    }

    # summer to spring growth factor
    lDelta_L_Pb_Mb_dot[y,1:nj] <- theta0[1:nj] + theta1[1:nj] * L_Pb_star[y,1:nj]
    lDelta_L_Pb_Mb[y,1:nj] ~ dmnorm(lDelta_L_Pb_Mb_dot[y,1:nj], Tau_lDelta_L_Pb_Mb[1:nj,1:nj])
    
    # size-based migration survival from trib to LGR: NOR
    Lphi_Mb_Ma_dot[y,i_spring,o_nor,1:nj] <- tau0[1:nj] + tau1[1:nj] * L_Mb_star[y,1:nj]
    Lphi_Mb_Ma[y,i_spring,o_nor,1:nj] ~ dmnorm(Lphi_Mb_Ma_dot[y,i_spring,o_nor,1:nj], Tau_Lphi_Mb_Ma[1:nj,1:nj])
    
    # non-size-based migration survival from trib to LGR: HOR
    Lphi_Mb_Ma[y,i_spring,o_hor,1:nj] ~ dmnorm(logit(mu_phi_Mb_Ma[i_spring,o_hor,1:nj]), Tau_Lphi_Mb_Ma[1:nj,1:nj])
    
    for (o in 1:no) {
      # pr(return at SWA1)
      Lpsi_O1[y,o,1:nj] ~ dmnorm(logit(mu_psi_O1[o,1:nj]), Tau_Lpsi_O1[1:nj,1:nj])

      # pr(return at SWA2)
      Lpsi_O2[y,o,1:nj] ~ dmnorm(logit(mu_psi_O2[o,1:nj]), Tau_Lpsi_O2[1:nj,1:nj])
    }
    
    # yr1 NOR ocean survival: includes AR(1) process
    Lphi_O0_O1_dot2[y,o_nor,1:nj] <- logit(mu_phi_O0_O1[o_nor,1:nj]) + Lphi_O0_O1_resid[y-1,o_nor,1:nj] * kappa_phi_O0_O1[1:nj]
    Lphi_O0_O1[y,o_nor,1:nj] ~ dmnorm(Lphi_O0_O1_dot2[y,o_nor,1:nj], Tau_Lphi_O0_O1[1:nj,1:nj])

    # yr2/yr3 NOR ocean survival: time constant
    Lphi_O1_O2[y,o_nor,1:nj] <- logit(mu_phi_O1_O2[o_nor,1:nj])
    Lphi_O2_O3[y,o_nor,1:nj] <- logit(mu_phi_O2_O3[o_nor,1:nj])
    
    # yr1/yr2/yr3 HOR ocean survival: same as NOR but adjusted by a time-constant log odds ratio
    Lphi_O0_O1[y,o_hor,1:nj] <- Lphi_O0_O1[y,o_nor,1:nj] + delta_O0_O1[1:nj]
    Lphi_O1_O2[y,o_hor,1:nj] <- Lphi_O1_O2[y,o_nor,1:nj] + delta_O1_O2[1:nj]
    Lphi_O2_O3[y,o_hor,1:nj] <- Lphi_O2_O3[y,o_nor,1:nj] + delta_O2_O3[1:nj]
    
    # pre-spawn survival
    Lphi_Sb_Sa[y,1:nj] <- logit(mu_phi_Sb_Sa[1:nj])

    # movement survival (juveniles LGR to BON)
    Lphi_Ma_O0[y,1:no] ~ dmnorm(logit(mu_phi_Ma_O0[1:no]), Tau_Lphi_Ma_O0[1:no,1:no])

    # movement survival (adults BON to LGR)
    Lphi_Rb_Ra_random[y,1:no] ~ dmnorm(logit(mu_phi_Rb_Ra[1:no]), Tau_Lphi_Rb_Ra[1:no,1:no])

    # apply inverse logit-transformation to all of these brood year-specific parameters
    # ilogit() and exp() cannot be applied to vectors unfortunately; must loop over population dimension
    for (j in 1:nj) {
      
      # transform egg to parr survival
      phi_E_Pb[y,j] <- ilogit(Lphi_E_Pb[y,j])
      
      # transform summer length and obtain scaled and centered version
      L_Pb[y,j] <- exp(lL_Pb[y,j])
      L_Pb_dot[y,j] <- exp(lL_Pb_dot[y,j])
      L_Pb_star[y,j] <- (L_Pb[y,j] - L_Pb_center[j])/L_Pb_scale[j]

      # transform LH apportionment
      pi[y,i_fall,j] <- ilogit(Lpi1[y,j])
      pi[y,i_spring,j] <- 1 - pi[y,i_fall,j]
      
      # transform overwinter survival
      for (i in 1:ni) {
        phi_Pa_Mb_dot[y,i,j] <- ilogit(Lphi_Pa_Mb_dot[y,i,j])
        phi_Pa_Mb[y,i,j] <- ilogit(Lphi_Pa_Mb[y,i,j])
      }

      # transform summer to spring growth factor
      Delta_L_Pb_Mb_dot[y,j] <- exp(lDelta_L_Pb_Mb_dot[y,j])
      Delta_L_Pb_Mb[y,j] <- exp(lDelta_L_Pb_Mb[y,j])
      
      # obtain spring length and scaled and centered version
      L_Mb[y,j] <- L_Pb[y,j] * Delta_L_Pb_Mb[y,j]
      L_Mb_star[y,j] <- (L_Mb[y,j] - L_Mb_center[j])/L_Mb_scale[j]

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
      
      phi_Mb_Ma_dot[y,i_spring,o_nor,j] <- ilogit(Lphi_Mb_Ma_dot[y,i_spring,o_nor,j])
      phi_O0_O1_dot2[y,o_nor,j] <- ilogit(Lphi_O0_O1_dot2[y,o_nor,j])
      
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
          # returning adults making it to BON (survive fisheries downstream of BON)
          Rb_BON[y,k,o,j] <- Rb[y,k,o,j] * (1 - U[y,k,o])
          
          # survive upstream migration from BON to LGR
          Ra_LGR[y,k,o,j] <- Rb_BON[y,k,o,j] * phi_Rb_Ra[y,o]
          
          # add strays
          Ra[y,k,o,j] <- Ra_LGR[y,k,o,j] + (G[y,o,j] * p_G[k,o,j])

          # remove fish at weir: use max() to ensure that fewer fish were removed than existed
          Sb[y,k,o,j] <- max(Ra[y,k,o,j] - B[y,k,o,j] - H[y,k,o,j], 1)
          
          # survive pre-spawn mortality
          Sa[y,k,o,j] <- Sb[y,k,o,j] * phi_Sb_Sa[y,j]
          
          # calculate egg production (eggs separated by age/origin)
          E_sep[y,k,o,j] <- Sa[y,k,o,j] * Omega[k,j] * f[y,k,j]
          
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
        sum(Sa[y+kmin+1-1,1,o_nor,j]) +       # age 3 adults produced by spawners in brood year y
          sum(Sa[y+kmin+2-1,2,o_nor,j]) +     # age 4 adults produced by spawners in brood year y
          sum(Sa[y+kmin+3-1,3,o_nor,j]))/     # age 5 adults produced by spawners in brood year y
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
    
    # eggs per spawner
    for (y in 2:ny) {
      E_per_Sa[y,j] <- E[y,j]/Sa_tot[y,j]
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
      expected_x_Ra[y,1:nko,j] <- p_Ra[y,1:nko,j] * nx_Ra[y,j]

      # calculate fit statistic: Freeman-Tukey statistic
      x_Ra_dev[y,j] <- sum((sqrt(x_Ra[y,1:nko,j]) - sqrt(expected_x_Ra[y,1:nko,j]))^2)
      x_Ra_new_dev[y,j] <- sum((sqrt(x_Ra_new[y,1:nko,j]) - sqrt(expected_x_Ra[y,1:nko,j]))^2)
      
      # calculate log posterior predictive density
      x_Ra_lppd[y,j] <- logdensity.multi(x_Ra[y,1:nko,j], p_Ra[y,1:nko,j], nx_Ra[y,j])
      
      # calculate observation model residual
      # x_Ra_obs_resid[y,1:nko,j] <- sqrt(x_Ra[y,1:nko,j]) - sqrt(expected_x_Ra[y,1:nko,j])
      for (ko in 1:nko) {
        x_Ra_obs_qresid[y,ko,j] <- pbin(x_Ra[y,ko,j], p_Ra[y,ko,j], nx_Ra[y,j])
      }
      
      # FOR CARCASSES
      # data likelihood
      x_Sa_prime[y,1:nko,j] ~ dmulti(p_Sa_prime[y,1:nko,j], nx_Sa_prime[y,j])
      
      # simulate new data
      x_Sa_prime_new[y,1:nko,j] ~ dmulti(p_Sa_prime[y,1:nko,j], nx_Sa_prime[y,j])
      
      # calculate expected count
      expected_x_Sa_prime[y,1:nko,j] <- p_Sa_prime[y,1:nko,j] * nx_Sa_prime[y,j]

      # calculate fit statistic: Freeman-Tukey statistic
      x_Sa_prime_dev[y,j] <- sum((sqrt(x_Sa_prime[y,1:nko,j]) - sqrt(expected_x_Sa_prime[y,1:nko,j]))^2)
      x_Sa_prime_new_dev[y,j] <- sum((sqrt(x_Sa_prime_new[y,1:nko,j]) - sqrt(expected_x_Sa_prime[y,1:nko,j]))^2)
      
      # calculate log posterior predictive density
      x_Sa_prime_lppd[y,j] <- logdensity.multi(x_Sa_prime[y,1:nko,j], p_Sa_prime[y,1:nko,j], nx_Sa_prime[y,j])
      
      # calculate observation model residual
      # x_Sa_prime_obs_resid[y,1:nko,j] <- sqrt(x_Sa_prime[y,1:nko,j]) - sqrt(expected_x_Sa_prime[y,1:nko,j])
      for (ko in 1:nko) {
        x_Sa_prime_obs_qresid[y,ko,j] <- pbin(x_Sa_prime[y,ko,j], p_Sa_prime[y,ko,j], nx_Sa_prime[y,j])
      }
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
    
    # calculate observation model residual
    lPa_obs_qresid[fit_Pa[d,1],fit_Pa[d,2],fit_Pa[d,3]] <- plnorm(Pa_obs[fit_Pa[d,1],fit_Pa[d,2],fit_Pa[d,3]], log(Pa[fit_Pa[d,1],fit_Pa[d,2],fit_Pa[d,3]]), 1/sig_Pa_obs[fit_Pa[d,1],fit_Pa[d,2],fit_Pa[d,3]]^2)
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
    
    # calculate observation model residual
    lMb_obs_qresid[fit_Mb[d,1],fit_Mb[d,2],fit_Mb[d,3],fit_Mb[d,4]] <- plnorm(Mb_obs[fit_Mb[d,1],fit_Mb[d,2],fit_Mb[d,3],fit_Mb[d,4]], log(Mb[fit_Mb[d,1],fit_Mb[d,2],fit_Mb[d,3],fit_Mb[d,4]]), 1/sig_Mb_obs[fit_Mb[d,1],fit_Mb[d,2],fit_Mb[d,3],fit_Mb[d,4]]^2)
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
    
    # calculate observation model residual
    lRa_obs_qresid[fit_Ra[d,1],fit_Ra[d,2]] <- plnorm(Ra_obs[fit_Ra[d,1],fit_Ra[d,2]], log(Ra_tot[fit_Ra[d,1],fit_Ra[d,2]]), 1/sig_Ra_obs[fit_Ra[d,1],fit_Ra[d,2]]^2)
  }
  
  # mean length at summer tagging 
  for (d in 1:nfit_L_Pb) {
    # data likelihood
    L_Pb_obs[fit_L_Pb[d,1],fit_L_Pb[d,2]] ~ dlnorm(log(L_Pb[fit_L_Pb[d,1],fit_L_Pb[d,2]]), 1/sig_L_Pb_obs[fit_L_Pb[d,1],fit_L_Pb[d,2]]^2)
    
    # simulate new data
    L_Pb_obs_new[fit_L_Pb[d,1],fit_L_Pb[d,2]] ~ dlnorm(log(L_Pb[fit_L_Pb[d,1],fit_L_Pb[d,2]]), 1/sig_L_Pb_obs[fit_L_Pb[d,1],fit_L_Pb[d,2]]^2)
    
    # calculate fit statistic: squared log deviates
    L_Pb_obs_dev[fit_L_Pb[d,1],fit_L_Pb[d,2]] <- (log(L_Pb_obs[fit_L_Pb[d,1],fit_L_Pb[d,2]]) - log(L_Pb[fit_L_Pb[d,1],fit_L_Pb[d,2]]))^2
    L_Pb_obs_new_dev[fit_L_Pb[d,1],fit_L_Pb[d,2]] <- (log(L_Pb_obs_new[fit_L_Pb[d,1],fit_L_Pb[d,2]]) - log(L_Pb[fit_L_Pb[d,1],fit_L_Pb[d,2]]))^2
    
    # calculate log posterior predictive density
    L_Pb_obs_lppd[fit_L_Pb[d,1],fit_L_Pb[d,2]] <- logdensity.lnorm(L_Pb_obs[fit_L_Pb[d,1],fit_L_Pb[d,2]], log(L_Pb[fit_L_Pb[d,1],fit_L_Pb[d,2]]), 1/sig_L_Pb_obs[fit_L_Pb[d,1],fit_L_Pb[d,2]]^2)
    
    # calculate observation model residual
    lL_Pb_obs_qresid[fit_L_Pb[d,1],fit_L_Pb[d,2]] <- plnorm(L_Pb_obs[fit_L_Pb[d,1],fit_L_Pb[d,2]], log(L_Pb[fit_L_Pb[d,1],fit_L_Pb[d,2]]), 1/sig_L_Pb_obs[fit_L_Pb[d,1],fit_L_Pb[d,2]]^2)
  }
  
  # mean length at spring tagging 
  for (d in 1:nfit_L_Mb) {
    # data likelihood
    L_Mb_obs[fit_L_Mb[d,1],fit_L_Mb[d,2]] ~ dlnorm(log(L_Mb[fit_L_Mb[d,1],fit_L_Mb[d,2]]), 1/sig_L_Mb_obs[fit_L_Mb[d,1],fit_L_Mb[d,2]]^2)
    
    # simulate new data
    L_Mb_obs_new[fit_L_Mb[d,1],fit_L_Mb[d,2]] ~ dlnorm(log(L_Mb[fit_L_Mb[d,1],fit_L_Mb[d,2]]), 1/sig_L_Mb_obs[fit_L_Mb[d,1],fit_L_Mb[d,2]]^2)
    
    # calculate fit statistic: squared log deviates
    L_Mb_obs_dev[fit_L_Mb[d,1],fit_L_Mb[d,2]] <- (log(L_Mb_obs[fit_L_Mb[d,1],fit_L_Mb[d,2]]) - log(L_Mb[fit_L_Mb[d,1],fit_L_Mb[d,2]]))^2
    L_Mb_obs_new_dev[fit_L_Mb[d,1],fit_L_Mb[d,2]] <- (log(L_Mb_obs_new[fit_L_Mb[d,1],fit_L_Mb[d,2]]) - log(L_Mb[fit_L_Mb[d,1],fit_L_Mb[d,2]]))^2
    
    # calculate log posterior predictive density
    L_Mb_obs_lppd[fit_L_Mb[d,1],fit_L_Mb[d,2]] <- logdensity.lnorm(L_Mb_obs[fit_L_Mb[d,1],fit_L_Mb[d,2]], log(L_Mb[fit_L_Mb[d,1],fit_L_Mb[d,2]]), 1/sig_L_Mb_obs[fit_L_Mb[d,1],fit_L_Mb[d,2]]^2)
    
    # calculate observation model residual
    lL_Mb_obs_qresid[fit_L_Mb[d,1],fit_L_Mb[d,2]] <- plnorm(L_Mb_obs[fit_L_Mb[d,1],fit_L_Mb[d,2]], log(L_Mb[fit_L_Mb[d,1],fit_L_Mb[d,2]]), 1/sig_L_Mb_obs[fit_L_Mb[d,1],fit_L_Mb[d,2]]^2)
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
    
    # calculate observation model residual
    Lphi_obs_Pb_Ma_qresid[fit_Lphi_Pb_Ma[d,1],fit_Lphi_Pb_Ma[d,2]] <- pnorm(Lphi_obs_Pb_Ma[fit_Lphi_Pb_Ma[d,1],fit_Lphi_Pb_Ma[d,2]], logit(phi_Pb_Ma[fit_Lphi_Pb_Ma[d,1],fit_Lphi_Pb_Ma[d,2]]), 1/sig_Lphi_obs_Pb_Ma[fit_Lphi_Pb_Ma[d,1],fit_Lphi_Pb_Ma[d,2]]^2)
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
    
    # calculate observation model residual
    Lphi_obs_Pa_Ma_qresid[fit_Lphi_Pa_Ma[d,1],fit_Lphi_Pa_Ma[d,2],fit_Lphi_Pa_Ma[d,3]] <- pnorm(Lphi_obs_Pa_Ma[fit_Lphi_Pa_Ma[d,1],fit_Lphi_Pa_Ma[d,2],fit_Lphi_Pa_Ma[d,3]], logit(phi_Pa_Ma[fit_Lphi_Pa_Ma[d,1],fit_Lphi_Pa_Ma[d,2],fit_Lphi_Pa_Ma[d,3]]), 1/sig_Lphi_obs_Pa_Ma[fit_Lphi_Pa_Ma[d,1],fit_Lphi_Pa_Ma[d,2],fit_Lphi_Pa_Ma[d,3]]^2)
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
    
    # calculate observation model residual
    Lphi_obs_Mb_Ma_qresid[fit_Lphi_Mb_Ma[d,1],fit_Lphi_Mb_Ma[d,2],fit_Lphi_Mb_Ma[d,3],fit_Lphi_Mb_Ma[d,4]] <- pnorm(Lphi_obs_Mb_Ma[fit_Lphi_Mb_Ma[d,1],fit_Lphi_Mb_Ma[d,2],fit_Lphi_Mb_Ma[d,3],fit_Lphi_Mb_Ma[d,4]], logit(phi_Mb_Ma[fit_Lphi_Mb_Ma[d,1],fit_Lphi_Mb_Ma[d,2],fit_Lphi_Mb_Ma[d,3],fit_Lphi_Mb_Ma[d,4]]), 1/sig_Lphi_obs_Mb_Ma[fit_Lphi_Mb_Ma[d,1],fit_Lphi_Mb_Ma[d,2],fit_Lphi_Mb_Ma[d,3],fit_Lphi_Mb_Ma[d,4]]^2)
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
    
    # calculate observation model residual
    Lphi_obs_Ma_O0_qresid[fit_Lphi_Ma_O0[d,1],fit_Lphi_Ma_O0[d,2]] <- pnorm(Lphi_obs_Ma_O0[fit_Lphi_Ma_O0[d,1],fit_Lphi_Ma_O0[d,2]], logit(phi_Ma_O0[fit_Lphi_Ma_O0[d,1],fit_Lphi_Ma_O0[d,2]]), 1/sig_Lphi_obs_Ma_O0[fit_Lphi_Ma_O0[d,1],fit_Lphi_Ma_O0[d,2]]^2)
  }
  
  # movement survival from BON to LGR
  for (d in 1:nfit_x_LGR) {
    # data likelihood
    x_LGR[fit_x_LGR[d,1],fit_x_LGR[d,2]] ~ dbin(phi_Rb_Ra[fit_x_LGR[d,1],fit_x_LGR[d,2]], x_BON[fit_x_LGR[d,1],fit_x_LGR[d,2]])
    
    # simulate new data
    x_LGR_new[fit_x_LGR[d,1],fit_x_LGR[d,2]] ~ dbin(phi_Rb_Ra[fit_x_LGR[d,1],fit_x_LGR[d,2]], x_BON[fit_x_LGR[d,1],fit_x_LGR[d,2]])
    
    # calculate expected count
    expected_x_LGR[fit_x_LGR[d,1],fit_x_LGR[d,2]] <- phi_Rb_Ra[fit_x_LGR[d,1],fit_x_LGR[d,2]] * x_BON[fit_x_LGR[d,1],fit_x_LGR[d,2]]

    # calculate fit statistic: Freeman-Tukey statistic
    x_LGR_dev[fit_x_LGR[d,1],fit_x_LGR[d,2]] <- (sqrt(x_LGR[fit_x_LGR[d,1],fit_x_LGR[d,2]]) - sqrt(expected_x_LGR[fit_x_LGR[d,1],fit_x_LGR[d,2]]))^2
    x_LGR_new_dev[fit_x_LGR[d,1],fit_x_LGR[d,2]] <- (sqrt(x_LGR_new[fit_x_LGR[d,1],fit_x_LGR[d,2]]) - sqrt(expected_x_LGR[fit_x_LGR[d,1],fit_x_LGR[d,2]]))^2

    # calculate log posterior predictive density
    x_LGR_lppd[fit_x_LGR[d,1],fit_x_LGR[d,2]] <- logdensity.bin(x_LGR[fit_x_LGR[d,1],fit_x_LGR[d,2]], phi_Rb_Ra[fit_x_LGR[d,1],fit_x_LGR[d,2]], x_BON[fit_x_LGR[d,1],fit_x_LGR[d,2]])
    
    # calculate observation model residuals
    x_LGR_obs_qresid[fit_x_LGR[d,1],fit_x_LGR[d,2]] <- pbin(x_LGR[fit_x_LGR[d,1],fit_x_LGR[d,2]], phi_Rb_Ra[fit_x_LGR[d,1],fit_x_LGR[d,2]], x_BON[fit_x_LGR[d,1],fit_x_LGR[d,2]])
  }
  
  # pre-spawn survival
  for (d in 1:nfit_x_carcass_spawned) {
    # data likelihood
    x_carcass_spawned[fit_x_carcass_spawned[d,1],fit_x_carcass_spawned[d,2]] ~ dbin(phi_Sb_Sa[fit_x_carcass_spawned[d,1],fit_x_carcass_spawned[d,2]], x_carcass_total[fit_x_carcass_spawned[d,1],fit_x_carcass_spawned[d,2]])
    
    # simulate new data
    x_carcass_spawned_new[fit_x_carcass_spawned[d,1],fit_x_carcass_spawned[d,2]] ~ dbin(phi_Sb_Sa[fit_x_carcass_spawned[d,1],fit_x_carcass_spawned[d,2]], x_carcass_total[fit_x_carcass_spawned[d,1],fit_x_carcass_spawned[d,2]])
    
    # calculate expected count
    expected_x_carcass_spawned[fit_x_carcass_spawned[d,1],fit_x_carcass_spawned[d,2]] <- phi_Sb_Sa[fit_x_carcass_spawned[d,1],fit_x_carcass_spawned[d,2]] * x_carcass_total[fit_x_carcass_spawned[d,1],fit_x_carcass_spawned[d,2]]

    # calculate fit statistic: Freeman-Tukey statistic
    x_carcass_spawned_dev[fit_x_carcass_spawned[d,1],fit_x_carcass_spawned[d,2]] <- (sqrt(x_carcass_spawned[fit_x_carcass_spawned[d,1],fit_x_carcass_spawned[d,2]]) - sqrt(expected_x_carcass_spawned[fit_x_carcass_spawned[d,1],fit_x_carcass_spawned[d,2]]))^2
    x_carcass_spawned_new_dev[fit_x_carcass_spawned[d,1],fit_x_carcass_spawned[d,2]] <- (sqrt(x_carcass_spawned_new[fit_x_carcass_spawned[d,1],fit_x_carcass_spawned[d,2]]) - sqrt(expected_x_carcass_spawned[fit_x_carcass_spawned[d,1],fit_x_carcass_spawned[d,2]]))^2

    # calculate log posterior predictive density
    x_carcass_spawned_lppd[fit_x_carcass_spawned[d,1],fit_x_carcass_spawned[d,2]] <- logdensity.bin(x_carcass_spawned[fit_x_carcass_spawned[d,1],fit_x_carcass_spawned[d,2]], phi_Sb_Sa[fit_x_carcass_spawned[d,1],fit_x_carcass_spawned[d,2]], x_carcass_total[fit_x_carcass_spawned[d,1],fit_x_carcass_spawned[d,2]])
    
    # calculate observation model residuals
    x_carcass_spawned_obs_qresid[fit_x_carcass_spawned[d,1],fit_x_carcass_spawned[d,2]] <- pbin(x_carcass_spawned[fit_x_carcass_spawned[d,1],fit_x_carcass_spawned[d,2]], phi_Sb_Sa[fit_x_carcass_spawned[d,1],fit_x_carcass_spawned[d,2]], x_carcass_total[fit_x_carcass_spawned[d,1],fit_x_carcass_spawned[d,2]])
  }
  
  ### CALCULATE ALL PROCESS MODEL RESIDUALS ###
  
  # year 0 residuals for process that include AR(1)
  Lphi_E_Pb_resid[1,1:nj] <- rep(0, nj)
  Lphi_O0_O1_resid[1,o_nor,1:nj] <- rep(0, nj)
  
  # all processes vary annually
  for (y in 2:ny) {
    
    # processes that vary by population
    for (j in 1:nj) {
      
      # total summer parr recruitment
      Lphi_E_Pb_resid[y,j] <- Lphi_E_Pb[y,j] - logit(phi_E_Pb_dot[y,j])
      logit(phi_E_Pb_dot2[y,j]) <- logit(phi_E_Pb_dot[y,j]) + kappa_phi_E_Pb[j] * Lphi_E_Pb_resid[y-1,j]
      
      # summer mean length
      lL_Pb_resid[y,j] <- lL_Pb[y,j] - lL_Pb_dot[y,j]
      
      # proportion of summer parr that become fall migrants
      Lpi_resid[y,j] <- logit(pi[y,i_fall,j]) - logit(mu_pi[i_fall,j])
      
      # LH-specific overwinter survival
      for (i in 1:ni) {
        Lphi_Pa_Mb_resid[y,i,j] <- Lphi_Pa_Mb[y,i,j] - Lphi_Pa_Mb_dot[y,i,j]
      }
      
      # summer to spring growth rate
      lDelta_L_Pb_Mb_resid[y,j] <- lDelta_L_Pb_Mb[y,j] - lDelta_L_Pb_Mb_dot[y,j]
      
      # movement survival to LGR: NOR
      Lphi_Mb_Ma_resid[y,o_nor,j] <- Lphi_Mb_Ma[y,i_spring,o_nor,j] - Lphi_Mb_Ma_dot[y,i_spring,o_nor,j]
      
      # movement survival to LGR: HOR
      Lphi_Mb_Ma_resid[y,o_hor,j] <- Lphi_Mb_Ma[y,i_spring,o_hor,j] - logit(mu_phi_Mb_Ma[i_spring,o_hor,j])
      
      # origin-specific quantities
      for (o in 1:no) {
        # pr(return at SWA1)
        Lpsi_O1_resid[y,o,j] <- Lpsi_O1[y,o,j] - logit(mu_psi_O1[o,j])
        
        # pr(return at SWA2|not returned at SWA1)
        Lpsi_O2_resid[y,o,j] <- Lpsi_O2[y,o,j] - logit(mu_psi_O2[o,j])
        
        # ocean survival
        Lphi_O0_O1_resid[y,o,j] <- Lphi_O0_O1[y,o,j] - logit(mu_phi_O0_O1[o,j])
      }
    }
    
    # processes that are constant across populations but vary by origin
    for (o in 1:no) {
      # movement survival from LGR thru BON
      Lphi_Ma_O0_resid[y,o] <- Lphi_Ma_O0[y,o] - logit(mu_phi_Ma_O0[o])
      
      # movement survival from BON thru LGR
      Lphi_Rb_Ra_resid[y,o] <- Lphi_Rb_Ra[y,o] - logit(mu_phi_Rb_Ra[o])
    }
  }
  
  # for parameters that may have non-vague priors, sample the prior without a link to likelihood
  # for comparison of posterior vs. prior density
  alpha_pr ~ dbeta(alpha_prior[1], alpha_prior[2])
  mu_psi_O1_pr ~ dbeta(mu_psi_O1_prior[1], mu_psi_O1_prior[2])
  mu_psi_O2_pr ~ dbeta(mu_psi_O2_prior[1], mu_psi_O2_prior[2])
  mu_phi_O0_O1_pr ~ dbeta(mu_phi_O0_O1_prior[1], mu_phi_O0_O1_prior[2])
  mu_phi_O1_O2_pr ~ dbeta(mu_phi_O1_O2_prior[1], mu_phi_O1_O2_prior[2])
  mu_phi_O2_O3_pr ~ dbeta(mu_phi_O2_O3_prior[1], mu_phi_O2_O3_prior[2])
}
