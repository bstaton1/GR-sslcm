
# load model info if not already in workspace
if (!exists("post")) source(file.path(this_dir, "load-model-info.R"))

##### EGG-TO-PARR SURVIVAL VS. EGG DENSITY #####

phi_E_Pa_VS_E_per_WUL = function(j, xlim = NULL, ylim = NULL) {
  
  # extract posterior summaries of total egg production by brood year
  E = post_summ(post, sub_index("^E[.+,pop]", pop = j))
  
  # extract posterior summaries of egg to parr survival rates
  phi_E_Pb = post_summ(post, sub_index("^phi_E_Pb[.+,pop]", pop = j))
  
  # extract posterior samples of the BH parameters
  bh_params = post_subset(post, sub_index(c("alpha[pop]", "^beta[pop]$"), pop = j), matrix = TRUE)
  
  # construct vector of egg abundances to predict at
  E_pred = seq(0.1, max(E[1,]), length = 30)
  
  # obtain predicted curves for each posterior sample
  phi_E_Pb_pred = t(sapply(1:post_dim(post, "saved"), function(i) 1/((1/bh_params[i,1]) + E_pred/bh_params[i,2])))
  
  # obtain posterior summary of predicted curve
  colnames(phi_E_Pb_pred) = paste0("phi_E_Pb_pred[", 1:length(E_pred), "]")
  phi_E_Pb_pred = post_convert(cbind(postpack:::id_mat(post), phi_E_Pb_pred))
  phi_E_Pb_pred = post_summ(phi_E_Pb_pred, "phi_E_Pb_pred")
  
  # set xaxis limits
  if (is.null(xlim)) xlim = c(0, max(E[1,]))
  
  # set yaxis limits
  if (is.null(ylim)) ylim = range(0, phi_E_Pb["mean",], phi_E_Pb_pred["2.5%",], phi_E_Pb_pred["97.5%",])
  
  # empty plot with correct dimensions
  plot(phi_E_Pb["mean",] ~ E["mean",],
       ylim = ylim,
       xlim = xlim, type = "n",
       xaxt = "n")
  
  # panel label
  panel_label(paste0("(", paste0("a", j), ")"), "topright", y_inp = 0.075)
  
  # draw prettier x-axis
  x_ticks = axisTicks(par("usr")[1:2], log = F)
  axis(side = 1, at = x_ticks, labels = x_ticks/1e6)
  
  # add credible region for predicted curve
  polygon(c(E_pred, rev(E_pred)), c(phi_E_Pb_pred["2.5%",], rev(phi_E_Pb_pred["97.5%",])), col = tranp_col, border = NA)
  
  # add posterior mean predicted curve
  lines(phi_E_Pb_pred["mean",] ~ E_pred, lwd = 2, col = solid_col2)
  
  # add posterior mean of realized pairs
  points(phi_E_Pb["mean",] ~ E["mean",], pch = 21, col = solid_col, bg = tranp_col, cex = pt_cex)
  
  # draw boundary box
  box()
}

##### PARR RECRUITMENT VS. EGG PRODUCTION #####

Pb_VS_E = function(j, xlim = NULL, ylim = NULL) {
  
  # extract posterior summaries of total egg production by brood year
  E = post_summ(post, sub_index("^E[.+,pop]", pop = j))
  
  # extract posterior summaries of total parr production by brood year
  Pb = post_summ(post, sub_index("^Pb[.+,pop]", pop = j))
  
  # extract posterior samples of the BH parameters
  bh_params = post_subset(post, sub_index(c("alpha[pop]", "^beta[pop]$"), pop = j), matrix = TRUE)
  
  # construct vector of egg abundances to predict at
  E_pred = seq(0.1, max(E[1,]), length = 30)
  
  # obtain predicted curves for each posterior sample
  Pb_pred = t(sapply(1:post_dim(post, "saved"), function(i) E_pred/((1/bh_params[i,1]) + E_pred/bh_params[i,2])))
  
  # obtain posterior summary of predicted curve
  colnames(Pb_pred) = paste0("Pb_pred[", 1:length(E_pred), "]")
  Pb_pred = post_convert(cbind(postpack:::id_mat(post), Pb_pred))
  Pb_pred = post_summ(Pb_pred, "Pb_pred")
  
  # set xaxis limits
  if (is.null(xlim)) xlim = c(0, max(E[1,]))
  
  # set yaxis limit
  if (is.null(ylim)) ylim = c(0, max(Pb[1,]))
  
  # empty plot with correct dimensions
  plot(Pb["mean",] ~ E["mean",],
       ylim = ylim,
       xlim = xlim, type = "n",
       xaxt = "n",
       yaxt = "n")
  
  # panel label
  panel_label(paste0("(", paste0("b", j), ")"), "topright", y_inp = 0.075)
  
  # draw prettier x-axis
  x_ticks = axisTicks(par("usr")[1:2], log = F)
  axis(side = 1, at = x_ticks, labels = x_ticks/1e6)
  
  # draw prettier y-axis
  y_ticks = axisTicks(par("usr")[3:4], log = F)
  axis(side = 2, at = y_ticks, labels = y_ticks/1e3)
  
  # add credible region for predicted curve
  polygon(c(E_pred, rev(E_pred)), c(Pb_pred["2.5%",], rev(Pb_pred["97.5%",])), col = tranp_col, border = NA)
  
  # add posterior mean predicted curve
  lines(Pb_pred["mean",] ~ E_pred, lwd = 2, col = solid_col2)
  
  # add posterior mean of realized pairs
  points(Pb["mean",] ~ E["mean",], pch = 21, col = solid_col, bg = tranp_col, cex = pt_cex)
  
  # draw boundary box
  box()
}

##### PARR SIZE VS. EGG DENSITY #####

L_Pb_VS_E_per_WUL = function(j, xlim = NULL, ylim = NULL) {
  # obtain posterior samples of total egg abundance scaled to WUL
  E = post_subset(post, sub_index("^E[.+,pop]", pop = j), matrix = TRUE)
  E_scaled = (E/jags_data$E_scale)/jags_data$wul[j]
  colnames(E_scaled) = gsub("E", "E_scaled", x = colnames(E_scaled))
  post_E = post_convert(cbind(postpack:::id_mat(post), E_scaled))
  
  # summarize scaled egg abundance
  E_scaled_mean = post_summ(post_E, "E_scaled")["mean",]
  
  # create a vector of scaled egg abundances to predict at
  E_scaled_seq = seq(min(E_scaled_mean), max(E_scaled_mean), length = 30)
  
  # extract posteriors of coefficients
  omega0 = post_subset(post, sub_index("omega0[pop]", pop = j), matrix = TRUE)
  omega1 = post_subset(post, sub_index("omega1[pop]", pop = j), matrix = TRUE)
  
  # extract posterior summaries of summer length outcomes
  L_Pb_mean = post_summ(post, sub_index("^L_Pb[.+,pop]", pop = j))["mean",]
  
  # function to create predicted survival curves for one posterior sample
  pred_fn = function(i) {
    pred_L_Pb = exp(omega0[i] + omega1[i] * log(E_scaled_seq))
    names(pred_L_Pb) = paste0("pred_L_Pb[", 1:30, "1]")
    pred_L_Pb
  }
  
  # calculate predicted length
  pred_L_Pb = t(sapply(1:post_dim(post, "saved"), pred_fn))
  post_pred_L_Pb = post_convert(cbind(postpack:::id_mat(post), pred_L_Pb))
  
  # extract summarized predicted length
  pred_L_Pb_mean = post_summ(post_pred_L_Pb, ".")["mean",]
  pred_L_Pb_lwr = post_summ(post_pred_L_Pb, ".")["2.5%",]
  pred_L_Pb_upr = post_summ(post_pred_L_Pb, ".")["97.5%",]
  
  E_scaled_seq = (E_scaled_seq * jags_data$wul[j] * jags_data$E_scale)/1e6
  E_scaled_mean = (E_scaled_mean * jags_data$wul[j] * jags_data$E_scale)/1e6
  
  # set xaxis limit
  if (is.null(xlim)) xlim = range(0, E_scaled_seq)
  
  # set yaxis limit
  if (is.null(ylim)) ylim = range(L_Pb_mean, pred_L_Pb_lwr, pred_L_Pb_upr, na.rm = TRUE)
  
  # empty plot with correct dimensions
  plot(1,1, type = "n", 
       ylim = ylim,
       xlim = xlim
  )
  
  # panel label
  panel_label(paste0("(", paste0("c", j), ")"), "topright", y_inp = 0.075)
  
  # add credible region for predicted curve
  polygon(c(E_scaled_seq, rev(E_scaled_seq)), c(pred_L_Pb_lwr, rev(pred_L_Pb_upr)), border = NA, col = tranp_col)
  
  # add posterior mean predicted curve
  lines(pred_L_Pb_mean ~ E_scaled_seq, col = solid_col2, lwd = 2)
  
  # add posterior mean of realized pairs
  points(L_Pb_mean ~ E_scaled_mean, pch = 21, col = solid_col, bg = tranp_col, cex = pt_cex)
  
  # draw boundary box
  box()
}

##### OVERWINTER SURVIVAL VS. PARR SIZE #####

phi_Pa_Mb_VS_L_Pb = function(j, i, xlim = NULL, ylim = NULL) {
  
  pch = c(21, 23)
  lty = c(1, 2)
  
  # extract posterior summaries of summer length outcomes
  L_Pb_mean = post_summ(post, sub_index("^L_Pb[year,pop]", year = ".+", pop = j))["mean",]
  
  # obtain scaled/centered versions
  L_Pb_star_mean = (L_Pb_mean - jags_data$L_Pb_center[j])/jags_data$L_Pb_scale[j]
  
  # summarize posterior of overwinter survival outcomes
  phi_Pa_Mb_mean = array_format(post_summ(post, sub_index("^phi_Pa_Mb[.+,LH_type,pop]", pop = j, LH_type = i))["mean",])[-1,i,j]
  
  # create vector of lengths to predict at
  L_Pb_seq = seq(min(L_Pb_mean), max(L_Pb_mean), length = 30)
  L_Pb_star_seq = seq(min(L_Pb_star_mean), max(L_Pb_star_mean), length = 30)
  
  # extract posteriors of coefficients
  gamma0 = post_subset(post, sub_index("gamma0[LH_type,pop]", pop = j, LH_type = i), matrix = TRUE)
  gamma1 = post_subset(post, sub_index("gamma1[LH_type,pop]", pop = j, LH_type = i), matrix = TRUE)
  
  # function to create predicted survival curve for one posterior sample
  pred_fn = function(i) {
    pred_phi_Pa_Mb = plogis(gamma0[i] + gamma1[i] * L_Pb_star_seq)
    names(pred_phi_Pa_Mb) = paste0("pred_phi_Pa_Mb[", 1:30, "]")
    pred_phi_Pa_Mb
  }
  
  # calculate predicted survivals and add to posterior samples
  pred_phi_Pa_Mb = t(sapply(1:post_dim(post, "saved"), pred_fn))
  post_pred_phi_Pa_Mb = post_convert(cbind(postpack:::id_mat(post), pred_phi_Pa_Mb))
  
  # summarize predicted survival
  pred_phi_Pa_Mb_mean = post_summ(post_pred_phi_Pa_Mb, ".")["mean",]
  pred_phi_Pa_Mb_lwr = post_summ(post_pred_phi_Pa_Mb, ".")["2.5%",]
  pred_phi_Pa_Mb_upr = post_summ(post_pred_phi_Pa_Mb, ".")["97.5%",]
  
  # set xaxis limits
  if (is.null(xlim)) xlim = range(L_Pb_seq)
  
  # set yaxis limits
  if (is.null(ylim)) ylim = range(0, pred_phi_Pa_Mb_lwr, pred_phi_Pa_Mb_upr, phi_Pa_Mb_mean)
  
  # empty plot with correct dimensions
  if (i == 1) {
    plot(1,1, type = "n", 
         ylim = ylim,
         xlim = xlim,
         xlab = "", ylab = ""
    )
  }
  
  # panel label
  panel_label(paste0("(", paste0("d", j), ")"), "topright", y_inp = 0.075)
  
  # add credible region for predicted curve
  polygon(c(L_Pb_seq, rev(L_Pb_seq)), c(pred_phi_Pa_Mb_lwr, rev(pred_phi_Pa_Mb_upr)), border = NA, col = tranp_col)
  
  # add posterior mean predicted curve
  lines(pred_phi_Pa_Mb_mean ~ L_Pb_seq, col = solid_col2, lwd = 2, lty = lty[i])
  
  # add posterior mean of realized pairs
  points(phi_Pa_Mb_mean ~ L_Pb_mean, pch = pch[i], col = solid_col, bg = tranp_col, cex = pt_cex)
  
  legend_pop = 2
  if (j == legend_pop & i == 1) {
    legend("topleft", legend = c("Fall", "Spring"), title = "Mig. Type", pch = pch, pt.bg = tranp_col, col = solid_col, pt.cex = pt_cex, cex = 0.8, bty = "n")
  }
  
  # draw boundary box
  box()
}

##### SMOLT SIZE VS PARR SIZE #####

L_Mb_VS_L_Pb = function(j, xlim = NULL, ylim = NULL) {
  
  # extract posterior summaries of summer length outcomes
  L_Pb_mean = post_summ(post, sub_index("^L_Pb[year,pop]", year = ".+", pop = j))["mean",]
  
  # extract posterior summaries of spring length outcomes
  L_Mb_mean = post_summ(post, sub_index("^L_Mb[year,pop]", year = ".+", pop = j))["mean",]
  
  # create a vector of scaled summer lengths to predict at
  L_Pb_seq = seq(min(L_Pb_mean), max(L_Pb_mean), length = 30)
  L_Pb_star_seq = (L_Pb_seq - jags_data$L_Pb_center[j])/jags_data$L_Pb_scale[j]
  
  # extract posteriors of coefficients
  theta0 = post_subset(post, sub_index("theta0[pop]", pop = j), matrix = TRUE)
  theta1 = post_subset(post, sub_index("theta1[pop]", pop = j), matrix = TRUE)
  
  # function to create predicted curves for one posterior sample
  pred_fn = function(i) {
    pred_L_Mb = exp(theta0[i] + theta1[i] * L_Pb_star_seq) * L_Pb_seq
    names(pred_L_Mb) = paste0("pred_L_Mb[", 1:30, "1]")
    pred_L_Mb
  }
  
  # calculate predicted growth values and add to posterior samples
  pred_L_Mb = t(sapply(1:post_dim(post, "saved"), pred_fn))
  post_pred_L_Mb = post_convert(cbind(postpack:::id_mat(post), pred_L_Mb))
  
  # summarize predicted growth
  pred_L_Mb_mean = post_summ(post_pred_L_Mb, ".")["mean",]
  pred_L_Mb_lwr = post_summ(post_pred_L_Mb, ".")["2.5%",]
  pred_L_Mb_upr = post_summ(post_pred_L_Mb, ".")["97.5%",]
  
  # set xaxis limits
  if (is.null(xlim)) xlim = range(L_Pb_seq)
  
  # set yaxis limits
  if (is.null(ylim)) ylim = range(pred_L_Mb_lwr, pred_L_Mb_upr, L_Mb_mean, na.rm = TRUE)
  
  # empty plot with the correct dimensions
  plot(1,1, type = "n", 
       ylim = ylim,
       xlim = xlim, 
       xlab = "", ylab = ""
  )
  
  # panel label
  panel_label(paste0("(", paste0("e", j), ")"), "topright", y_inp = 0.075)
  
  # add credible region for predicted curve
  polygon(c(L_Pb_seq, rev(L_Pb_seq)), c(pred_L_Mb_lwr, rev(pred_L_Mb_upr)), border = NA, col = tranp_col)
  
  # add posterior mean of predicted curve
  lines(pred_L_Mb_mean ~ L_Pb_seq, col = solid_col2, lwd = 2)
  
  # add posterior mean of realized pairs
  points(L_Mb_mean ~ L_Pb_mean, pch = 21, col = solid_col, bg = tranp_col, cex = pt_cex)
  
  # draw boundary box
  box()
}

##### SURVIVAL TO LGR VS. SMOLT SIZE #####

phi_Mb_Ma_VS_L_Mb = function(j, xlim = NULL, ylim = NULL) {
  
  # extract posterior summaries of spring length outcomes
  L_Mb_mean = post_summ(post, sub_index("^L_Mb[year,pop]", year = ".+", pop = j))["mean",]
  
  # obtain scaled/centered versions
  L_Mb_star_mean = (L_Mb_mean - jags_data$L_Mb_center[j])/jags_data$L_Mb_scale[j]
  
  # summarize posterior of migration to LGR survival outcomes
  phi_Mb_Ma_mean = array_format(post_summ(post, sub_index("^phi_Mb_Ma[.+,LH_type,origin,pop]", pop = j, LH_type = 1, origin = 1))["mean",])[-1,1,1,j]
  
  # create vectors to predict survival at: for credible regions and mean curve
  L_Mb_seq = seq(min(L_Mb_mean), max(L_Mb_mean), length = 30)
  L_Mb_star_seq = seq(min(L_Mb_star_mean), max(L_Mb_star_mean), length = 30)
  
  # extract posteriors of coefficients
  tau0 = post_subset(post, sub_index("tau0[pop]", pop = j), matrix = TRUE)
  tau1 = post_subset(post, sub_index("tau1[pop]", pop = j), matrix = TRUE)
  
  # function to create predicted survival curve for one posterior sample
  pred_fn = function(i) {
    pred_phi_Mb_Ma = plogis(tau0[i] + tau1[i] * L_Mb_star_seq)
    names(pred_phi_Mb_Ma) = paste0("pred_phi_Mb_Ma[", 1:30, "]")
    pred_phi_Mb_Ma
  }
  
  # calculate predicted survivals for each posterior sample
  pred_phi_Mb_Ma = t(sapply(1:post_dim(post, "saved"), pred_fn))
  post_pred_phi_Mb_Ma = post_convert(cbind(postpack:::id_mat(post), pred_phi_Mb_Ma))
  
  # summarize predicted survival
  pred_phi_Mb_Ma_mean = post_summ(post_pred_phi_Mb_Ma, ".")["mean",]
  pred_phi_Mb_Ma_lwr = post_summ(post_pred_phi_Mb_Ma, ".")["2.5%",]
  pred_phi_Mb_Ma_upr = post_summ(post_pred_phi_Mb_Ma, ".")["97.5%",]
  
  # set xaxis limits
  if (is.null(xlim)) xlim = range(L_Mb_seq)
  
  # set yaxis limits
  if (is.null(ylim)) ylim = range(pred_phi_Mb_Ma_lwr, pred_phi_Mb_Ma_upr, phi_Mb_Ma_mean)
  
  # empty plot with correct dimensions
  plot(1,1, type = "n", 
       ylim = ylim,
       xlim = xlim,
       xlab = "", ylab = ""
  )
  
  # panel label
  panel_label(paste0("(", paste0("f", j), ")"), "topright", y_inp = 0.075)
  
  # add credible region for predicted curve
  polygon(c(L_Mb_seq, rev(L_Mb_seq)), c(pred_phi_Mb_Ma_lwr, rev(pred_phi_Mb_Ma_upr)), border = NA, col = tranp_col)
  
  # add posterior mean of predicted curve
  lines(pred_phi_Mb_Ma_mean ~ L_Mb_seq, col = solid_col2, lwd = 2)
  
  # add posterior mean of realized pairs
  points(phi_Mb_Ma_mean ~ L_Mb_mean, pch = 21, col = solid_col, bg = tranp_col, cex = pt_cex)
}

##### LARGE MULTI-PANEL PLOT #####

make_layout = function(nrow, ncol, show = FALSE) {
  m = NULL
  for (r in 1:nrow) {
    if (r == 1) {
      f = 1
    } else {
      f = max(m) + 1
    }
    l = f+ncol-1
    tmp = rbind(f:l, rep(l+1, ncol))
    m = rbind(m, tmp)
  }
  layout(m, heights = rep(c(1, 0.25), nrow))
  if (show) layout.show(max(m))
}

g = function(ylim = c(-3,3), xlim = c(-3,3)) {
  x = rnorm(30)
  y = rnorm(30)
  plot(y ~ x, xlab = "", ylab = "", xlim = xlim, ylim = ylim)
}

yaxis_label = function(label_main, label_sub) {
  mtext(side = 2, outer = FALSE, label_main, line = 2.25, cex = 0.9)
  mtext(side = 2, outer = FALSE, label_sub, line = 1.25, cex = 0.75, font = 3)
}

xaxis_label = function(label) {
  plot(1, 1, type = "n", xlim = c(0,1), ylim = c(0,1), axes = FALSE)
  text(0.5, 0.3, label = label, cex = 1.25, xpd = TRUE)
}

# open a graphics device
dev.on(file.path(this_dir, "relationships"), width = 7.2, height = 9, format = fig_type)

# construct the layout
make_layout(nrow = 6, ncol = 4, show = FALSE)

# graphical parameters
par(mar = c(0.5,1.5,0,0), oma = c(0,2,1.5,1), tcl = -0.1, mgp = c(200,0.05,0),
    lend = "square", ljoin = "mitre")

# egg to parr survival vs. egg density
# ylim = make_lim(post_summ(post, "^phi_E_Pb[")["mean",])
ylim = c(0,0.4)
xlim = make_lim(post_summ(post, "^E[")["mean",]); xlim[1] = 0
junk = sapply(1:4, function(j) {
  phi_E_Pa_VS_E_per_WUL(j, xlim = xlim, ylim = ylim)
  if (j == 1) yaxis_label("Egg \u2192 Parr Survival", "")
  mtext(side = 3, outer = FALSE, line = 0, pops[j], font = 2)
}); xaxis_label("Total Egg Production (Millions)")

# parr recruitment vs. total egg production
ylim = make_lim(post_summ(post, "^Pb[")["mean",]); ylim[1] = 0
xlim = make_lim(post_summ(post, "^E[")["mean",]); xlim[1] = 0
junk = sapply(1:4, function(j) {
  Pb_VS_E(j, xlim = xlim, ylim = ylim)
  if (j == 1) yaxis_label("Parr Recruitment", "Thousands")
}); xaxis_label("Total Egg Production (Millions)")

# parr size vs. egg production
ylim = make_lim(post_summ(post, "^L_Pb[")["mean",])
xlim = make_lim(post_summ(post, "^E[")["mean",]/1e6); xlim[1] = 0
junk = sapply(1:4, function(j) {
  L_Pb_VS_E_per_WUL(j, xlim = xlim, ylim = ylim)
  if (j == 1) yaxis_label("Parr Mean Length", "mm")
}); xaxis_label("Total Egg Production (Millions)")

# overwinter survival vs. parr size
ylim = make_lim(post_summ(post, "^phi_Pa_Mb[")["mean",])
xlim = make_lim(post_summ(post, "^L_Pb[")["mean",])
junk = sapply(1:4, function(j) {
  phi_Pa_Mb_VS_L_Pb(j, 1, xlim = xlim, ylim = ylim)
  phi_Pa_Mb_VS_L_Pb(j, 2, xlim = xlim, ylim = ylim)
  if (j == 1) yaxis_label("Parr \u2192 Smolt Survival", "")
}); xaxis_label("Parr Mean Length (mm)")

# smolt size vs. parr size
xlim = make_lim(post_summ(post, "^L_Pb[")["mean",])
ylim = make_lim(post_summ(post, "^L_Mb[")["mean",])
junk = sapply(1:4, function(j) {
  L_Mb_VS_L_Pb(j, xlim = xlim, ylim = ylim)
  if (j == 1) yaxis_label("Smolt Mean Length", "mm")
}); xaxis_label("Parr Mean Length (mm)")

# migration to LGR vs. smolt size
xlim = make_lim(post_summ(post, "^L_Mb[")["mean",])
ylim = make_lim(post_summ(post, sub_index("^phi_Mb_Ma[year,LH_type,origin,pop]", year = ".+", LH_type = "1", origin = 1, pop = ".")))

junk = sapply(1:4, function(j) {
  phi_Mb_Ma_VS_L_Mb(j, xlim = xlim, ylim = ylim)
  if (j == 1) yaxis_label("Smolt \u2192 LGR Survival", "")
}); xaxis_label("Smolt Mean Length (mm)")

# close the device
dev.off()
