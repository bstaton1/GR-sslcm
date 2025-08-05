
# load model info if not already in workspace
if (!exists("post")) source(file.path(this_dir, "load-model-info.R"))

# function to extract the posterior samples for all parameters needed for the calculation
# for one population
prep_samples = function(j) {
  
  # the parameters needed to calculate process model expected values
  the_params = c("alpha[pop]", "beta[pop]", "^lambda$",
                 "omega0[pop]", "omega1[pop]", "mu_pi[1,pop]",
                 "gamma0[1,pop]", "gamma1[1,pop]",
                 "gamma0[2,pop]", "gamma1[2,pop]",
                 "theta0[pop]", "theta1[pop]",
                 "tau0[pop]", "tau1[pop]"
  )
  
  # extract the posterior samples for all of these parameters for this population
  the_samps = post_subset(post, sub_index(the_params, pop = j), matrix = TRUE)
  
  # handle the names
  the_names = postpack:::drop_index(names(the_samps[1,]))
  the_names[which(the_names == "gamma0")] = paste0("gamma0", c("_fall", "_spring"))
  the_names[which(the_names == "gamma1")] = paste0("gamma1", c("_fall", "_spring"))
  colnames(the_samps) = the_names
  
  return(the_samps)
}

# prepare one posterior sample for one population
prep_sample = function(the_samps, i, j) {
  c(the_samps[i,],
    WUL = unname(jags_data$wul[j]),
    E_scale = unname(jags_data$E_scale),
    L_Pb_scale = unname(jags_data$L_Pb_scale[j]),
    L_Pb_center = unname(jags_data$L_Pb_center[j]),
    L_Mb_scale = unname(jags_data$L_Mb_scale[j]),
    L_Mb_center = unname(jags_data$L_Mb_center[j])
  )
}

# get deterministic population quantity at a given number of eggs, one population, one posterior sample
# x is the parameter vector
# E is the number of eggs
# WUL change is the % change in WUL
# start is the starting term in calculating survival
# out_type is type of output: phi is survival from start to Ma; Ma is Ma abundance
get_phi_hat = function(x, E, WUL_change = 0, start = "E", out_type = "phi") {
  with(as.list(c(x, E = E, WUL_change = WUL_change)), {
    
    # original deterministic beta estimate at current WUL
    beta_hat = lambda * WUL
    
    # updated WUL
    WUL_new = WUL * (1 + WUL_change)
    
    # updated deterministic beta estimate at updated WUL
    beta_hat_new = lambda * WUL_new
    
    # updated realized beta estimate -- uses the noise term from the original realized beta
    beta_new = beta_hat_new * (beta/beta_hat)
    
    # egg to parr survival
    phi_E_Pb = 1/(1/alpha + E/beta_new)
    
    # parr abundance
    Pb = E * phi_E_Pb
    
    # parr mean length
    L_Pb = exp(omega0 + omega1 * log((E/E_scale)/WUL_new))  # in mm
    L_Pb_star = (L_Pb - L_Pb_center)/L_Pb_scale             # scaled and centered version
    
    # apportion to migratory type
    Pa_fall = Pb * mu_pi
    Pa_spring = Pb * (1 - mu_pi)
    
    # overwinter survival
    phi_Pa_Mb_fall = plogis(gamma0_fall + gamma1_fall * L_Pb_star)
    phi_Pa_Mb_spring = plogis(gamma0_spring + gamma1_spring * L_Pb_star)
    
    # smolt abundance, in tributary
    Mb_fall = Pa_fall * phi_Pa_Mb_fall
    Mb_spring = Pa_spring * phi_Pa_Mb_spring
    
    # growth factor
    Delta = exp(theta0 + theta1 * L_Pb_star)
    
    # smolt mean length
    L_Mb = L_Pb * Delta
    L_Mb_star = (L_Mb - L_Mb_center)/L_Mb_scale
    
    # migration survival to LGR
    phi_Mb_Ma = plogis(tau0 + tau1 * L_Mb_star)
    
    # move smolt to LGR
    Ma_fall = Mb_fall * phi_Mb_Ma
    Ma_spring = Mb_spring * phi_Mb_Ma
    Ma = Ma_fall + Ma_spring
    
    # calculate output to return:
    # survival from a starting point to Ma if out_type == "phi"
    # Ma if out_type == "Ma"
    if (out_type == "phi") {
      if (start == "E") out = Ma/E
      if (start == "Pb") out = Ma/Pb
    } 
    
    if (out_type == "Ma") {
      out = Ma
    }
    out
  })
}

# wrapper around get_phi_hat for all posterior samples and a range of E values
get_phi_hat_post = function(j, E_pred, WUL_change = 0, start = "E", out_type = "phi") {
  
  # prepare all posterior samples for this population
  the_samps = prep_samples(j = j)
  
  # loop over posterior samples
  out = sapply(1:nrow(the_samps), function(i) {
    
    # loop over E values
    x = prep_sample(the_samps, i = i, j = j)
    sapply(E_pred, function(e) get_phi_hat(x, E = e, WUL_change = WUL_change, start = start, out_type = out_type))
  })
  
  # make an mcmc.list object
  out = t(out)
  colnames(out) = paste0("var", 1:ncol(out))
  postpack::post_convert(cbind(postpack:::id_mat(post), out))
}

# get the realized population quanitites
get_phi_real_post = function(j, start = "E", out_type) {
  
  if (start == "E") start_param = "^E[.+,pop]"
  if (start == "Pb") start_param = "^Pb[.+,pop]"
  
  Ma_fall = post_subset(post, sub_index("^Ma[.+,1,1,pop]", pop = j), matrix = TRUE)
  Ma_spring = post_subset(post, sub_index("^Ma[.+,2,1,pop]", pop = j), matrix = TRUE)
  start_out = post_subset(post, sub_index(start_param, pop = j), matrix = TRUE)
  
  if (out_type == "phi") {
    out = (Ma_fall + Ma_spring)/start_out
  } else {
    out = Ma_fall + Ma_spring
  }
  
  colnames(out) = paste0("var", 1:ncol(out))
  postpack::post_convert(cbind(postpack:::id_mat(post), out))
}

plot_f = function(j, start, out_type) {
  
  WUL_changes = c(0.5, 1, 2)
  lwd = 1.5
  
  # get eggs per spawner for this population
  eps = mean(post_summ(post, sub_index("E_per_Sa[year,pop]", year = ".+", pop = j))["mean",])
  
  # get the realized population value each year
  phi_real_post = get_phi_real_post(j = j, start = start, out_type = out_type)
  phi_real = post_summ(phi_real_post, ".")
  
  # get realized egg abundances each year
  E_real = post_summ(post, sub_index("^E[.+,pop]", pop = j))
  
  # create a vector of E values to calculate deterministic population value at
  E_pred = seq(eps * 100, max(E_real["mean",]), length = 30)
  
  # the years to keep
  keep_yrs = stringr::str_extract(colnames(E_real), "[:digit:]+,") |> 
    stringr::str_remove(",")
  
  # number of spawners for the predicted and actual values
  S_pred = E_pred/eps
  S_real = post_summ(post, sub_index("^Sa_tot[year,pop]", year = keep_yrs, pop = j))
  
  # obtain deterministic population values: current WUL
  phi_hat_post = get_phi_hat_post(j = j, E_pred, start = start, out_type = out_type)
  phi_hat = post_summ(phi_hat_post, ".")
  
  # obtain deterministic population values: hypothetical WUL
  phi_hat_new1_post = get_phi_hat_post(j = j, E_pred, WUL_change = WUL_changes[1], start = start, out_type = out_type)
  phi_hat_new1 = post_summ(phi_hat_new1_post, ".")
  phi_hat_new2_post = get_phi_hat_post(j = j, E_pred, WUL_change = WUL_changes[2], start = start, out_type = out_type)
  phi_hat_new2 = post_summ(phi_hat_new2_post, ".")
  phi_hat_new3_post = get_phi_hat_post(j = j, E_pred, WUL_change = WUL_changes[3], start = start, out_type = out_type)
  phi_hat_new3 = post_summ(phi_hat_new3_post, ".")
  
  # empty plot
  plot(phi_real["mean",] ~ S_real["mean",], type = "n",
       ylim = make_lim(0, phi_hat["97.5%",], phi_hat_new3["mean",]),
       # phi_real["mean",]),
       # ylim = make_lim(0, ifelse(start == "E", 0.075, 0.3)),
       xlim = c(0, max(S_real["mean",], S_pred)),
       xaxt = "n", yaxt = "n")
  
  # draw on realized points
  # points(phi_real["mean",] ~ S_real["mean",], pch = 21,
  #        col = solid_col2, bg = tranp_col, cex = 0.8)
  
  # draw uncertainty intervals
  polygon(x = c(S_pred, rev(S_pred)), 
          y = c(phi_hat["2.5%",], rev(phi_hat["97.5%",])),
          col = tranp_col, border = NA)
  
  # draw posterior mean curves
  lines(phi_hat["mean",] ~ S_pred, lwd = lwd)
  lines(phi_hat_new1["mean",] ~ S_pred, lty = 3, lwd = lwd)
  lines(phi_hat_new2["mean",] ~ S_pred, lty = 2, lwd = lwd)
  lines(phi_hat_new3["mean",] ~ S_pred, lty = 5, lwd = lwd)
  
  # draw x-axis
  at_x = axisTicks(par("usr")[1:2], log = FALSE)
  axis(side = 1, at = at_x, labels = at_x/1000)
  
  # draw y-axis
  at_y = axisTicks(par("usr")[3:4], log = FALSE)
  axis(side = 2, at = at_y, labels = at_y/ifelse(out_type == "Ma", 1000, 1))
  
  # draw panel label
  base = ifelse(start == "E", 0, 4)
  base = ifelse(out_type == "Ma", 8, base)
  panel_label(paste0("(", letters[base + j], ") "))
  
  # draw axes labels
  if (j == 1) {
    if (start == "E" & out_type == "phi") {
      axis_labels(ylab = "Egg \u2192 LGR Survival", outer = FALSE, yline = 2)
      legend("topright", title = "Weighted Habitat", legend = c("Current", paste0("\u2191", WUL_changes * 100, "%")), lty = c(1,3,2,5), seg.len = 2.7, bty = "n", lwd = lwd, cex = 0.9)
    } 
    if (start == "Pb" & out_type == "phi") axis_labels(ylab = "Parr \u2192 LGR Survival", outer = FALSE, yline = 2)
    if (out_type == "Ma") {
      axis_labels(ylab = "LGR Smolt", outer = FALSE, yline = 2)
      axis_labels(ylab = "Thousands", outer = FALSE, yline = 1, cex = 0.8, font = 3)
    } 
  }
  
  # add a population label 
  if (start == "E" & out_type == "phi") {
    mtext(side = 3, line = 0.15, font = 2, text = pops[j], cex = 1)
  }
}

# open a graphics device
dev.on(file.path(this_dir, "wul-change"), width = 7.2, height = 6, format = fig_type)
mypar(mfrow = c(3,4), oma = c(2.5,2.5,1,0), col.axis = "black", tcl = -0.1)
junk = sapply(1:4, function(j) plot_f(j, "E", "phi"))
junk = sapply(1:4, function(j) plot_f(j, "Pb", "phi"))
junk = sapply(1:4, function(j) plot_f(j, "Pb", "Ma"))
axis_labels(xlab = "Spawner Abundance")
axis_labels(xlab = "Thousands", outer = TRUE, xline = 1.5, cex = 0.8, font = 3)
dev.off()
