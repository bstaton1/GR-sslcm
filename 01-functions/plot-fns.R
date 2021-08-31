# :::::::::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT TO HOUSE FUNCTIONS FOR PLOTTING INFORMATION #
# :::::::::::::::::::::::::::::::::::::::::::::::::: #

##### FUNCTIONS TO GET UNCERTAINTY IN DATA FOR PLOTTING #####

## FOR A DATASET OF MULTINOMIAL RANDOM VARIABLES
# x is the matrix of counts (rows = years; columns = group), i is the group id
get_obs_ests_multinomial = function(x, i) {
  
  # calculate the proportion point estimate
  p_ests = t(apply(x, 1, function(y) y/sum(y)))
  p_ests[p_ests == "NaN"] = NA
  
  # calculate number of trials
  N_tot = rowSums(x)
  
  # calculate standard error of proportion point estimate
  se_p_ests = apply((p_ests * (1 - p_ests)), 2, function(x) sqrt(x/N_tot))
  
  # calculate confidence intervals; impose boundaries
  lwr = p_ests - 1.96 * se_p_ests; lwr[lwr < 0] = 0.001
  upr = p_ests + 1.96 * se_p_ests; upr[upr > 1] = 0.999
  
  out = cbind(
    mean = p_ests[,i],
    lwr95 = lwr[,i],
    upr95 = upr[,i]
  )
  return(out)
}

## FOR A DATASET OF LOGIT-NORMAL RANDOM VARIABLES
# Lmean is the logit-scale point estimate; Lsig is the logit-normal standard error
get_obs_ests_logit_normal = function(Lmean, Lsig) {
  out = cbind(
    mean = expit(Lmean),
    lwr95 = expit(qnorm(0.025, Lmean, Lsig)),
    upr95 = expit(qnorm(0.975, Lmean, Lsig))
  )
  return(out)
}

## FOR A DATASET OF LOGNORMAL RANDOM VARIABLES
# lmean is the log-scale point estiamte; lsig is the log-normal standard error
get_obs_ests_log_normal = function(lmean, lsig) {
  out = cbind(
    mean = exp(lmean),
    lwr95 = qlnorm(0.025, lmean, lsig),
    upr95 = qlnorm(0.975, lmean, lsig)
  )
  return(out)
}

##### PLOT A TIME SERIES OF ESTIMATES, WITH DATA IF AVAILABLE #####

# est: output from postpack::post_summ()
# obs: optional vector of observed data to plot on top of estimates
# main: plot title
# xaxis: plot the xaxis?
# yaxis_side: 2 or 4 or NULL, where should the y-axis be placed?
# set_par: set some par() options internally?

plot_tseries = function(est, obs = NULL, main = NULL, xaxis = T, yaxis_side = 2, set_par = T, ylim = NULL, yrs = NULL) {
  
  # extract year labels
  if (is.null(yrs)) {
    yrs = as.numeric(rownames(obs))
  }
  
  # set y limits
  if (is.null(ylim)) {
    ylim = range(cbind(t(est[c("2.5%", "97.5%"),]), obs), na.rm = T)
  }
  
  # set the graphics device parameters
  if (set_par) par(mar = c(1.5,1.5,1.5,0.5), tcl = -0.15, mgp = c(2,0.35,0))
  
  # create empty plot with correct dimensions/labels
  plot(1,1, type = "n", 
       xlim = range(yrs), ylim = ylim,
       xlab = "", las = 1, ylab = "", main = main, xaxt = "n", yaxt = "n")
  
  # draw band representing 95% credible region
  polygon(x = c(yrs, rev(yrs)), y = c(est["2.5%",], rev(est["97.5%",])), border = NA, col = alpha("salmon", 0.5))
  lines(est["2.5%",] ~ yrs, col = "red", lty = 2)
  lines(est["97.5%",] ~ yrs, col = "red", lty = 2)
  
  # draw posterior median
  lines(est["50%",] ~ yrs, col = "red", lwd = 2)
  
  # draw data if provided
  if (!is.null(obs)) {
    segments(yrs, obs[,"lwr95"], yrs, obs[,"upr95"], col = "blue")
    points(obs[,"mean"] ~ yrs, pch = 21, col = "blue", bg = alpha("skyblue2", 0.5), cex = 1.2)
  } 
  
  # draw axes/labels
  if (xaxis) {
    at_yrs = c(seq(min(yrs), max(yrs), 4), max(yrs))
    axis(side = 1, at = at_yrs, labels = paste0("", substr(at_yrs, 3, 4)))
  } 
  if (!is.null(yaxis_side)) axis(side = yaxis_side)
}
