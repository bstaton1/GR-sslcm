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

##### PLOTTING INFRASTRUCTURE FUNCTIONS #####

# create a nice axis limit
make_lim = function(..., buffer = c(0, 0.1)) {
  lim = range(..., na.rm = TRUE)
  lim_diff = diff(lim)
  lim = lim + lim_diff * buffer
  lim
}

# include a label in the corner of a plot
panel_label = function(label_text, label_loc = "topleft", cex = 1,
                       font = 2, x_inp = 0.015, y_inp = 0.05, col = par("col.axis"), ...) {
  
  accepted = c("topleft", "topright", "bottomleft", "bottomright")
  
  usr = par("usr"); xdiff = diff(usr[1:2]); ydiff = diff(usr[3:4])
  
  is_left = stringr::str_detect(label_loc, "left")
  is_top = stringr::str_detect(label_loc, "top")
  xsign = ifelse(is_left, -1, 1)
  ysign = ifelse(is_top, -1, 1)
  xi = ifelse(is_left, 1, 2)
  yi = ifelse(is_top, 4, 3)
  pos = ifelse(is_left, 4, 2)
  
  text(x = usr[xi] + xsign * xdiff * x_inp,
       y = usr[yi] + ysign * ydiff * y_inp,
       labels = label_text, cex = cex, pos = pos, font = font, col = col, ...)
}

# draw axes labels
axis_labels = function(xlab = NULL, ylab = NULL, xline = 0.4, yline = 0.4, outer = TRUE, ...) {
  mtext(side = 1, line = xline, text = xlab, outer = outer, ...)
  mtext(side = 2, line = yline, text = ylab, outer = outer, ...)
} 

# default par()
mypar = function(mar = c(1,1,0.5,0.5),
                 mfrow = c(2,2), 
                 oma = c(1.5,1.5,0,0), mgp = c(200,0.05,0), 
                 tcl = 0, lend = "square",
                 ljoin = "mitre",
                 col.axis = "grey40", cex.axis = 0.85, ...) {
  
  do.call(par, c(as.list(environment()), list(...)))
}

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
    ylim = range(list(est[c("2.5%", "97.5%"),], obs), na.rm = T)
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
    obs_yrs = as.numeric(rownames(obs))
    segments(obs_yrs, obs[,"lwr95"], obs_yrs, obs[,"upr95"], col = alpha("blue", 0.75))
    points(obs[,"mean"] ~ obs_yrs, pch = 21, col = alpha("blue", 0.75), bg = alpha("skyblue2", 0.5), cex = 1.2)
  } 
  
  # draw axes/labels
  if (xaxis) {
    at_yrs = c(seq(min(yrs), max(yrs), 4), max(yrs))
    axis(side = 1, at = yrs, labels = FALSE)
    axis(side = 1, at = at_yrs, labels = paste0("", substr(at_yrs, 3, 4)), tcl = -0.3)
  } 
  if (!is.null(yaxis_side)) axis(side = yaxis_side)
}

add_tseries = function(est, yrs, col = "royalblue") {
  polygon(x = c(all_yrs[ts_yrs], rev(all_yrs[ts_yrs])),
          y = c(est["2.5%",], rev(est["97.5%",])),
          col = scales::alpha(col, 0.15), border = NA)
  lines(est["mean",] ~ all_yrs[ts_yrs], col = col, lwd = 2)
  lines(est["2.5%",] ~ all_yrs[ts_yrs], col = col, lty = 2)
  lines(est["97.5%",] ~ all_yrs[ts_yrs], col = col, lty = 2)
}

##### FUNCTIONS FOR COMPARING POSTERIORS AMONG 2 OR MORE MODELS #####

## EXAMPLE USAGE OF THESE FUNCTIONS ##
# store two or more mcmc.lists as list elements, call it post_list

# legend = list(a = "topleft", b = NULL, c = NULL, d = NULL)
# main = c("CAT", "LOS", "MIN", "UGR")
# 
# # COMPARE ALPHA ESTIMATES FOR EACH POPULATION ACROSS MODELS
# compare_param(post_list, "^alpha[pop]", pop = 1:4, main = "Max Egg-to-Parr Survival", legend = "topleft")
# 
# # COMPARE ADULT RETURNS TO EACH POPULATION ACROSS MODELS, WITH DATA+OBS VARIABILITY SHOWN
# par(mfrow = c(2,2), oma = c(2,2,0,0))
# sapply(1:4, function(j) {
#   compare_tseries(post_list, "^Ra_tot[.+,pop]", yrs = 1991:2019, 
#                   obs = get_obs_ests_log_normal(log(jags_data$Ra_obs[as.character(1991:2019),j]), jags_data$sig_Ra_obs[as.character(1991:2019),j]),
#                   pop = j, origin = 1, main = main[j], legend = legend[[j]], y_scale = 1000)
# })
# mtext(side = 1, outer = TRUE, line = 0.5, "Return Year")
# mtext(side = 2, outer = TRUE, line = 0.5, "Total Adults (000s)")
# 
# # COMPARE 1ST YEAR NOR OCEAN SURVIVAL ACROSS MODELS, NO DATA TO SHOW
# par(mfrow = c(2,2), oma = c(2,2,0,0))
# sapply(1:4, function(j) {
#   compare_tseries(post_list, "^phi_O0_O1[.+,origin,pop]", yrs = 1991:2019, pop = j, origin = 1, main = main[j], legend = legend[[j]])
# })
# mtext(side = 1, outer = TRUE, line = 0.5, "Brood Year")
# mtext(side = 2, outer = TRUE, line = 0.5, "Yr1 Ocean Survival (NOR)")
# 
# # COMPARE THE BH RELATIONSHIP FOR EACH POPULATION ACROSS MODELS
# par(mfrow = c(2,2), oma = c(2,2,0,0))
# sapply(1:4, function(j) {
#   compare_BH(post_list, pop = j, main = main[j], legend = legend[[j]])
# })
# mtext(side = 1, outer = TRUE, line = 0.5, "Total Egg Production (000 000s)")
# mtext(side = 2, outer = TRUE, line = 0.5, "Total Parr Recruitment (000s)")

# function to create a plot comparing posterior summmaries of the same quantity among models
compare_param = function(post_list, param, ylab = "", main = "", legend = "topleft", y_scale = 1, pt_est = "mean", cols = NULL, ylim = NULL, ...) {
  
  # get identifier for each element that will be summarized
  param_ID = dim_IDs(param, ...)
  
  # get posterior summaries from all models
  ests = lapply(post_list, function(post) post_summ(post, sub_index(param, ...), probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))
  
  # convert summaries into an array format rather than a list format (models are third dimension)
  ests = do.call(abind, append(ests, list(along = 3)))
  dimnames(ests)[[2]] = param_ID
  
  # calculate the y-axis limits if not supplied
  if (is.null(ylim)) {
    y_range = range(ests[stringr::str_detect(dimnames(ests)[[1]], "%"),,])
    y_diff = diff(y_range)
    ylim = y_range + c(-1,1) * 0.05 * y_diff
  }
  
  # define colors if not supplied
  if (is.null(cols)) {
    all_cols = c("salmon", "skyblue2", "orange", "forestgreen")
    cols = all_cols[1:length(post_list)]
  } 
  
  # graphics parameters
  par(lend = "square", mgp = c(2,0.35,0), tcl = -0.15, mar = c(2,2,2,1))
  
  # empty plot
  mp = barplot(t(ests[pt_est,,]), beside = TRUE, col = "white", border = NA, ylim = ylim, yaxt = "n", ylab = ylab, main = main)
  usr = par("usr")
  segments(usr[1], usr[3], usr[2], usr[3], xpd = TRUE)
  segments(usr[1], usr[3], usr[1], usr[4], xpd = TRUE)
  
  # draw uncertainty intervals
  segments(mp, t(ests["2.5%",,]), mp, t(ests["97.5%",,]), col = cols)
  segments(mp, t(ests["25%",,]), mp, t(ests["75%",,]), lwd = 6, col = cols)
  
  # draw point estimates
  points(x = mp, y = t(ests[pt_est,,]), cex = 1.5, pch = 3, col = cols)
  
  # draw axes
  at_y = axisTicks(usr[3:4], log = FALSE)
  axis(side = 1, at = colSums(mp)/nrow(mp), labels = FALSE)
  axis(side = 2, at = at_y, labels = at_y/y_scale, las = 2)
  
  if (!is.null(legend)) {
    legend(legend, legend = names(post_list), pch = 15, col = cols, pt.cex = 2, bty = "n")
  }
}

# function to create a plot comparing a posterior summaries of the same time series among models
compare_tseries = function(post_list, param, yrs, obs = NULL, ylab = "", xlab = "", main = "", legend = "topleft", y_scale = 1, pt_est = "mean", cols = NULL, ylim = NULL,  ...) {
  
  # get posterior summaries from all models
  ests = lapply(post_list, function(post) post_summ(post, sub_index(param, ...), probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))
  
  # convert summaries into an array format rather than a list format (models are third dimension)
  ests = do.call(abind, append(ests, list(along = 3)))
  
  # calculate the y-axis limits
  if (is.null(ylim)) {
    y_range = range(ests[stringr::str_detect(dimnames(ests)[[1]], "%"),,], obs)
    y_diff = diff(y_range)
    ylim = y_range + c(-1,1) * 0.05 * y_diff
    ylim[1] = 0
  }
  
  # define colors if not supplied
  if (is.null(cols)) {
    all_cols = c("salmon", "skyblue2", "orange", "forestgreen")
    cols = all_cols[1:length(post_list)]
  } 
  tcols = scales::alpha(cols, 0.3)
  
  # graphics parameters
  par(mar = c(2,2,2,1), lend = "square", mgp = c(2,0.35,0), tcl = -0.15)
  
  plot(1,1, type = "n", xlim = range(yrs), ylim = ylim, main = main, ylab = ylab, xaxt = "n", yaxt = "n", xlab = xlab)
  usr = par("usr")
  junk = sapply(1:dim(ests)[3], function(i) {
    polygon(x = c(yrs, rev(yrs)), y = c(ests["2.5%",,i], rev(ests["97.5%",,i])), col = tcols[i], border = NA)
    lines(ests[pt_est,,i] ~ yrs, col = cols[i], lwd = 2)
    lines(ests["2.5%",,i] ~ yrs, col = cols[i], lty = 2)
    lines(ests["97.5%",,i] ~ yrs, col = cols[i], lty = 2)
  })
  
  # draw data if provided
  if (!is.null(obs)) {
    segments(yrs, obs[,"lwr95"], yrs, obs[,"upr95"], col = alpha("grey25", 0.25))
    points(obs[,"mean"] ~ yrs, pch = 21, col = alpha("grey25", 0.5), bg = alpha("grey25", 0.25), cex = 1.2)
  } 
  
  # draw axes
  at_x = axisTicks(usr[1:2], log = FALSE)
  at_y = axisTicks(usr[3:4], log = FALSE)
  axis(side = 1, at = at_x, labels = at_x)
  axis(side = 2, at = at_y, labels = at_y/y_scale, las = 2)
  
  if (!is.null(legend)) {
    if (is.null(obs)) {
      legend(legend, legend = names(post_list), pch = 15, col = cols, pt.cex = 2, bty = "n")
    } else {
      legend(legend, legend = c(names(post_list), "Data"), pch = 15, col = c(cols, "grey"), pt.cex = 2, bty = "n")
    }
  }
}

# function to compare BH relationships for one population across multiple models
compare_BH = function(post_list, pop, main = "", legend = "topleft", cols = NULL) {
  
  # function to get a sequence of egg abundances to predict at
  # calculates a 30 element vector counting up from zero to the maximum 97.5% quantile of eggs by a given population across models
  get_pred_eggs = function(j) {
    max_eggs = max(unlist(lapply(post_list, function(post) max(post_summ(post, sub_index("f_tot[.+,pop]", pop = j))["97.5%",]))))
    pred_eggs = seq(0, max_eggs, length = 30)
    return(pred_eggs)
  }
  
  # function to posterior summaries of predicted parr at each egg abundance
  get_pred_parr = function(post, j) {
    alpha_post = post_subset(post, sub_index("alpha[pop]", pop = j), TRUE)
    beta_post = post_subset(post, sub_index("beta[pop]", pop = j), TRUE)
    pred_eggs = get_pred_eggs(j)
    pred_parr = sapply(1:nrow(alpha_post), function(i) BH(pred_eggs, alpha_post[i], beta_post[i]))
    apply(pred_parr, 1, function(x) c(mean = mean(x), quantile(x, c(0.025, 0.975))))
  }
  
  # get the sequence of egg abundances to predict parr abundances at
  pred_eggs = get_pred_eggs(pop)
  
  # get posterior of predicted parr abundances
  pred_parr = lapply(post_list, get_pred_parr, j = pop)
  
  # define colors if not supplied
  if (is.null(cols)) {
    all_cols = c("salmon", "skyblue2", "orange", "forestgreen")
    cols = all_cols[1:length(post_list)]
  } 
  tcols = scales::alpha(cols, 0.3)
  
  # graphics parameters
  par(mar = c(2,2,2,1), lend = "square", mgp = c(2,0.35,0), tcl = -0.15)
  
  plot(1,1, type = "n", xlim = range(pred_eggs), ylim = c(0, max(pred_parr[[1]], pred_parr[[2]])),
       main = main, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  usr = par("usr")
  junk = sapply(1:length(post_list), function(i) {
    polygon(x = c(pred_eggs, rev(pred_eggs)), y = c(pred_parr[[i]]["2.5%",], rev(pred_parr[[i]]["97.5%",])), col = tcols[i], border = NA)
    lines(pred_parr[[i]]["mean",] ~ pred_eggs, col = cols[i], lwd = 2)
    lines(pred_parr[[i]]["2.5%",] ~ pred_eggs, col = cols[i], lty = 2)
    lines(pred_parr[[i]]["97.5%",] ~ pred_eggs, col = cols[i], lty = 2)
  })
  
  # draw axes
  at_x = axisTicks(usr[1:2], log = FALSE)
  at_y = axisTicks(usr[3:4], log = FALSE)
  axis(side = 1, at = at_x, labels = at_x/1e6)
  axis(side = 2, at = at_y, labels = at_y/1000, las = 2)
  
  if (!is.null(legend)) {
    legend(legend, legend = names(post_list), pch = 15, col = cols, pt.cex = 2, bty = "n")
  }
}

##### GENERALIZED PLOTTING FUNCTIONS #####

# given a parameter name, return the axis limits that would include the 95% CRI
make_lim_param = function(post, param) {
  
  # which rows of post_summ() output should be used in calculating axis limits
  use_in_lim = c(FALSE, FALSE, FALSE, TRUE, TRUE)
  
  # get the range of relevant posterior summaries
  sub_index(sub_index(param, year = ts_yrs)) |>
    post_summ(post, param = _) |>
    subset(use_in_lim) |>
    range()
}
