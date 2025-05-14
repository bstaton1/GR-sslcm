# scenario name
sslcm_scenario = "final_vlong"

# set the proportion of MCMC samples to retain
keep_percent = 1

# load the model output
source("load-model-info.R")
id_mat = postpack:::id_mat(post)

# basic BH parameter summaries
post_summ(post, "^alpha[", digits = 2)  # theoretical max egg-to-parr survival
post_summ(post, "^beta[", digits = 0)   # theoretical max parr rearing capacity

# increase in parr rearing capacity per 1 unit increase in WUL
post_summ(post, "^lambda$", digits = 0)

# How many parr produced per spawner at 50, 100 spawners?
E_per_Sa = post_subset(post, "E_per_Sa", matrix = TRUE)
mean_E_per_Sa = apply(E_per_Sa, 1, function(x) colMeans(array_format(x), na.rm = TRUE))
mean_E_per_Sa_samps = t(mean_E_per_Sa)
colMeans(mean_E_per_Sa_samps)   # time series average eggs-per-spawner
eggs_50 = mean_E_per_Sa_samps * 50
eggs_100 = mean_E_per_Sa_samps * 100
alpha_samps = post_subset(post, "^alpha[", matrix = TRUE)
beta_samps = post_subset(post, "^beta[", matrix = TRUE)
BH_eggs = function(eggs, alpha, beta) {
  (1/((1/alpha) + (eggs/beta))) * eggs
}
parr_per_spawner_50 = do.call(cbind, lapply(1:4, function(j) {
  sapply(1:nrow(eggs_50), function(i) BH_eggs(eggs_50[i,j], alpha_samps[i,j], beta_samps[i,j])/50)
})); colnames(parr_per_spawner_50) = paste0("parr_per_spawner_50[", 1:4, "]")
parr_per_spawner_100 = do.call(cbind, lapply(1:4, function(j) {
  sapply(1:nrow(eggs_100), function(i) BH_eggs(eggs_100[i,j], alpha_samps[i,j], beta_samps[i,j])/100)
})); colnames(parr_per_spawner_100) = paste0("parr_per_spawner_100[", 1:4, "]")
parr_per_spawner_post = post_convert(cbind(id_mat, parr_per_spawner_50, parr_per_spawner_100))
post_summ(parr_per_spawner_post, "50", digits = 0)
post_summ(parr_per_spawner_post, "100", digits = 0)

# percent change moving from 100 to 50 spawners
mat100 = post_subset(parr_per_spawner_post, "100", matrix = TRUE)
mat50 = post_subset(parr_per_spawner_post, "50", matrix = TRUE)
pdiff = (mat50 - mat100)/mat100
colnames(pdiff) = paste0("pdiff[", 1:4, "]")
pdiff_post = post_convert(cbind(id_mat, pdiff))
post_summ(pdiff_post, ".", digits = 2) * 100

# autoregressive coefficients for BH relationship
kappa_samps = post_subset(post, "kappa_phi_E_Pb", matrix = TRUE)
kappa_mean_samps = post_convert(cbind(id_mat, mean_kappa = rowMeans(kappa_samps)))
post_summ(kappa_mean_samps, ".", digits = 2)

# mean probability of parr becoming fall migrants
post_summ(post, "mu_pi[1,", digits = 2) * 100

# overwinter survival differences between spring and fall migrants
# reported approximate values based on Figure 4

# statements about fraction of years that over-winter survival for fall migrants overlapped that of spring migrants
fall_mig = post_summ(post, "phi_Pa_Mb[.+,1,.]")
spring_mig = post_summ(post, "phi_Pa_Mb[.+,2,.]")
fall_mean = array_format(fall_mig["mean",])[,1,]
fall_lwr = array_format(fall_mig["2.5%",])[,1,]
fall_upr = array_format(fall_mig["97.5%",])[,1,]
spring_mean = array_format(spring_mig["mean",])[,2,]
spring_lwr = array_format(spring_mig["2.5%",])[,2,]
spring_upr = array_format(spring_mig["97.5%",])[,2,]
dimnames(fall_mean) = dimnames(fall_lwr) = dimnames(fall_upr) = 
  dimnames(spring_mean) = dimnames(spring_lwr) = dimnames(spring_upr) = 
  list(dimnames(jags_data$Pa_obs)[[1]], dimnames(jags_data$Pa_obs)[[3]])

fall_mean_is_lower = fall_mean < spring_mean
CRL_overlaps = fall_lwr < spring_upr
fall_mean_is_lower
round(colMeans(CRL_overlaps, na.rm = TRUE), 2) * 100
msdown::percentize(sapply(1:4, function(j) cor(fall_mean[-1,j], spring_mean[-1,j])))

# average NOR and HOR survival migration to LGR
# reported approximate values from figure 4

# mean survival from LGR to ocean
post_summ(post, "mu_phi_Ma_O0", digits = 2)

# first year ocean survival for NOR fish
post_summ(post, "mu_phi_O0_O1[1,", digits = 2)

# "rare" ocean survival years reported approximate values from figure 4

# odds ratio to obtain HOR first year ocean survival
post_summ(post_convert(cbind(id_mat, exp(post_subset(post, "delta_O0_O1", matrix = TRUE)[,-3]))), ".", digits = 2)

# average proportion of age-3/4 maturation reported as approximate values from figure 6

# survival rate along upstream migration reported as approximate values from figure 4

# harvest rates below BON
jags_data$U[,,"NOR"]  # typically <5% for NOR
jags_data$U[,,"HOR"]  # typically 10-15% for HOR

# tributary harvest summaries as fraction of total return
x = apply(jags_data$H[,,"NOR",] + jags_data$H[,,"HOR",], 3, rowSums)/jags_data$Ra_obs
mean(x[,"LOS"] < 0.1, na.rm = TRUE)  # 95% of years <10%
max(x[,"LOS"], na.rm = TRUE)         # max percentage
rownames(x)[which.max(x[,"LOS"])]    # year with

# brood stock removal summaries
x = apply(jags_data$B[,,"NOR",] + jags_data$B[,,"HOR",], 3, rowSums)/jags_data$Ra_obs
colMeans(x[all_yrs >= 2000,], na.rm = TRUE)

# prespawn mortality summaries
ps_mort = 1 - jags_data$phi_Sb_Sa[-1,]
apply(ps_mort, 2, range)
apply(ps_mort, 2, mean)

# overall survival from estuary to spawning grounds
# taken from Supplement B, section 9.5

# change in odds of overwinter survival for 1 SD increase in mean length
jags_data$L_Pb_scale  # 1SD mean length
egamma_post = post_convert(cbind(id_mat, exp(post_subset(post, "gamma1[1,", matrix = TRUE))))
post_summ(egamma_post, ".", digits = 2)

# change in odds of outmigration survival for 1 SD increase in mean length
jags_data$L_Mb_scale
x = post_subset(post, "tau1[", matrix = TRUE)
x = cbind(x, rowMeans(x)); colnames(x) = paste0("tau1[", 1:5, "]")
etau_post = post_convert(cbind(id_mat, exp(x)))
post_summ(etau_post, ".", digits = 2)

# example of out-migration survival at two sizes
j = 1
LMb_start = (85 - jags_data$L_Mb_center[j])/jags_data$L_Mb_scale[j]
LMb_end = (95 - jags_data$L_Mb_center[j])/jags_data$L_Mb_scale[j]
tau_post = post_subset(post, sub_index("tau.[pop]", pop = j), matrix = TRUE)
p_start = plogis(tau_post[,1] + tau_post[,2] * LMb_start)
p_end = plogis(tau_post[,1] + tau_post[,2] * LMb_end)
p_change_post = post_convert(cbind(id_mat, "phi_start_L" = p_start, "phi_end_L" = p_end))
post_summ(p_change_post, ".", digits = 2)

# correlation in process noise in egg-to-parr survival
# given from Supplement B, section 10.2.2

## DISCUSSION STATEMENTS ABOUT FALL MIGRANTS BEING MORE ABUNDANT AS LGR SMOLT THAN AS REARING PARR -----

# probability of being fall migrant, at parr stage
pi = post_summ(post, "^pi[")["mean",] |> 
  array_format()
pi_Pb = pi[,i_fall,]

# abundance of smolt reaching LGR
Ma = post_summ(post, "^Ma[")["mean",] |> 
  array_format()
Ma_fall = Ma[,i_fall,o_nor,]     # extract NOR fall migrants
Ma_spring = Ma[,i_spring,o_nor,] # extract NOR spring migrants

# calculate proportion of fall migrants at LGR
pi_Ma = Ma_fall/(Ma_fall+Ma_spring)

# compare across year means
colMeans(pi_Pb, na.rm = TRUE)  # proportion fall migrants as parr
colMeans(pi_Ma, na.rm = TRUE)  # proportion fall migrants as LGR smolt

## DISCUSSION STATEMENTS ABOUT HOR FISH BEING LESS COMMON AS ADULTS THAN AS SMOLTS -----

# abundance of smolt reaching LGR
Ma = post_summ(post, "^Ma[")["mean",] |> 
  array_format()
Ma_nor = Ma[,i_fall,o_nor,] + Ma[,i_spring,o_nor,]    # extract NOR fall migrants
Ma_hor = Ma[,i_spring,o_hor,] # extract NOR spring migrants

# proportion of NOR fish upon reaching LGR as smolts
p_NOR_Ma = Ma_nor/(Ma_nor + Ma_hor)

# abundance of adults returning to river
Ra = post_summ(post, "^Ra[")["mean",] |> 
  array_format()

# sum across ages
Ra_tot_nor = apply(Ra[,,o_nor,], 3, rowSums)
Ra_tot_hor = apply(Ra[,,o_hor,], 3, rowSums)

# proportion of NOR fish upon return to river
p_NOR_Ra = Ra_tot_nor/(Ra_tot_nor + Ra_tot_hor)

# years to keep (basically all years HOR fish present for all pops with data)
keep = !all_yrs %in% c(1990:2003, 2021:2022)

# calculate average proportions for each
round(colMeans(1 - p_NOR_Ma[keep,-j_min]), 2)  # proportion NOR as LGR smolts
round(colMeans(1 - p_NOR_Ra[keep,-j_min]), 2)  # proportion NOR as adults

## CALCULATE WITHIN-POPULATION CORRELATIONS (specific results not reported in-text) -----

resid_params = match_params(post, "resid", type = "base_only")
resid_params[!stringr::str_detect(resid_params, "obs")]

within_pop_resid_corr = function(j) {
  
  resid_params = c(
    "Egg-to-Parr Surv." = "Lphi_E_Pb_resid[year,pop]",
    "Parr Growth" = "lL_Pb_resid[year,pop]",
    "Prop. Fall Mig." = "Lpi_resid[year,pop]",
    "Overwinter Surv. (Fall)" = "Lphi_Pa_Mb_resid[year,1,pop]",
    "Overwinter Surv. (Spring)" = "Lphi_Pa_Mb_resid[year,2,pop]",
    "Overwinter Growth" = "lDelta_L_Pb_Mb_resid[.+,pop]",
    "Outmigration Surv. (NOR)" = "Lphi_Mb_Ma_resid[.+,1,pop]",
    "Outmigration Surv. (HOR)" = "Lphi_Mb_Ma_resid[.+,2,pop]",
    "Juv. Hydro. Surv. (NOR)" = "Lphi_Ma_O0_resid[.+,1]",
    "Juv. Hydro. Surv. (HOR)" = "Lphi_Ma_O0_resid[.+,2]",
    "Yr1 Ocean Surv." = "Lphi_O0_O1_resid[year,1,pop]", 
    "Prop. Age-3 (NOR)" = "Lpsi_O1_resid[year,1,pop]",
    "Prop. Age-3 (HOR)" = "Lpsi_O1_resid[year,2,pop]",
    "Prop. Age-4 (NOR)" = "Lpsi_O2_resid[year,1,pop]",
    "Prop. Age-4 (HOR)" = "Lpsi_O2_resid[year,2,pop]",
    "Adult Hydro. Surv. (NOR)" = "Lphi_Rb_Ra_resid[year,1]",
    "Adult Hydro. Surv. (HOR)" = "Lphi_Rb_Ra_resid[year,2]"
  )
  
  # dimensions
  n_terms = length(resid_params)
  n_samps = post_dim(post, types = "saved")
  
  # get posterior samples of each residual term
  resid_posts = lapply(resid_params, function(param) {
    # subset out all posterior samples of all years for this residual for this population
    post_sub = post_subset(post, sub_index(param, pop = j, year = ".+"), matrix = TRUE)
    
    # remove the first year if for egg-to-parr & year-1 ocean survival residuals
    # has extra year because AR(1) process
    if (param %in% c("Lphi_E_Pb_resid[year,pop]", "Lphi_O0_O1_resid[year,1,pop]")) post_sub = post_sub[,-1]
    
    # return
    return(post_sub)
  })
  
  # loop through MCMC iters, calculating correlation term between all residual terms
  # store in 3D array [term,term,mcmc_iter]
  rho_post = array(NA, c(n_terms, n_terms, n_samps))
  for (iter in 1:n_samps) {
    cat("\r", iter)
    for (row in 1:n_terms) {
      for (col in 1:n_terms) {
        rho_post[row,col,iter] = cor(resid_posts[[row]][iter,], resid_posts[[col]][iter,])
      }
    }
  }
  
  # summarize posteriors of correlation matrix
  rho_mean = rho_lwr = rho_upr = rho_sig = matrix(NA, n_terms, n_terms)
  for (row in 1:n_terms) {
    for (col in 1:n_terms) {
      rho_mean[row,col] = mean(rho_post[row,col,])
      rho_lwr[row,col] = quantile(rho_post[row,col,], 0.025)
      rho_upr[row,col] = quantile(rho_post[row,col,], 0.975)
      rho_sig[row,col] = !(rho_lwr[row,col] < 0 & rho_upr[row,col] > 0)
    }
  }
  
  # assign dimension names to these summary objects
  dimnames(rho_mean) = dimnames(rho_lwr) = dimnames(rho_upr) = dimnames(rho_sig) = 
    list(names(resid_params), names(resid_params))
  
  
  # create a matrix with NAs in all correlation matrix cells we don't want to keep
  # upper triangle (duplicated) and diagonal (all 1s)
  NA_matrix = matrix(1, n_terms, n_terms)
  diag(NA_matrix) = NA
  NA_matrix[upper.tri(NA_matrix)] = NA
  dimnames(NA_matrix) = dimnames(rho_mean)
  
  # model already explicitly models these correlations, can exclude here
  NA_matrix["Juv. Hydro. Surv. (HOR)","Juv. Hydro. Surv. (NOR)"] = NA
  NA_matrix["Adult Hydro. Surv. (HOR)","Adult Hydro. Surv. (NOR)"] = NA
  
  # format the summary objects
  rho_mean = round(rho_mean, 2) * NA_matrix
  rho_lwr = round(rho_lwr, 2) * NA_matrix
  rho_upr = round(rho_upr, 2) * NA_matrix
  rho_sig = rho_sig * NA_matrix
  
  # build row and column IDs: ensures correct labeling when converted to long form
  row_name = matrix(rep(names(resid_params), each = n_terms), nrow = n_terms, ncol = n_terms, byrow = TRUE)
  col_name = matrix(rep(names(resid_params), each = n_terms), nrow = n_terms, ncol = n_terms, byrow = FALSE)
  
  # convert summaries to long form data.frame
  rho_summary = data.frame(
    pop = pops[j],
    term1 = as.character(row_name),
    term2 = as.character(col_name),
    mean = as.numeric(rho_mean),
    lwr = as.numeric(rho_lwr),
    upr = as.numeric(rho_upr),
    sig = as.logical(rho_sig)
  )
  
  # remove all non-reported rows (see NA_matrix above)
  rho_summary[-which(is.na(rho_summary$sig)),]
}

out = do.call(rbind, lapply(1:4, within_pop_resid_corr))

out[order(abs(out$mean), decreasing = TRUE),] |> 
  subset(sig)
