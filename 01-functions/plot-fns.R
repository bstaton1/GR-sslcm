# :::::::::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT TO HOUSE FUNCTIONS FOR PLOTTING INFORMATION #
# :::::::::::::::::::::::::::::::::::::::::::::::::: #

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
    yrs = as.numeric(names(obs))
  }
  
  # set y limits
  if (is.null(ylim)) {
    ylim = range(rbind(est[c("2.5%", "97.5%"),], obs), na.rm = T)
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
  if (!is.null(obs)) points(obs ~ yrs, pch = 21, col = "blue", bg = alpha("skyblue2", 0.5), cex = 1.2)
  
  # draw axes/labels
  if (xaxis) {
    at_yrs = c(seq(min(yrs), max(yrs), 4), max(yrs))
    axis(side = 1, at = at_yrs, labels = paste0("", substr(at_yrs, 3, 4)))
  } 
  if (!is.null(yaxis_side)) axis(side = yaxis_side)
}
