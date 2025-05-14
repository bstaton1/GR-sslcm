
# load model info if not already in workspace
if (!exists("post")) source(file.path(this_dir, "load-model-info.R"))

# initialize an empty list
plot_dat = list()

# prep values to plot: fall trap passage
plot_dat$Pa_fall = list(
  obs = jags_data$Pa_obs[,i_fall,],
  fitted = {
    tmp = post_summ(post, sub_index("^Pa[.+,LH_type,.]", LH_type = i_fall))
    list(
      mn = array_format(tmp["mean",])[,i_fall,],
      lwr = array_format(tmp["2.5%",])[,i_fall,],
      upr = array_format(tmp["97.5%",])[,i_fall,]
    )
  }
)

# prep values to plot: spring trap passage
plot_dat$Mb_spring = list(
  obs = jags_data$Mb_obs[,i_spring,o_nor,],
  fitted = {
    tmp = post_summ(post, sub_index("^Mb[.+,LH_type,origin,.]", LH_type = i_spring, origin = o_nor))
    list(
      mn = array_format(tmp["mean",])[,i_spring,o_nor,],
      lwr = array_format(tmp["2.5%",])[,i_spring,o_nor,],
      upr = array_format(tmp["97.5%",])[,i_spring,o_nor,]
    )
  }
)

# prep values to plot: adult return to river
plot_dat$Ra_tot = list(
  obs = jags_data$Ra_obs,
  fitted = {
    tmp = post_summ(post, "^Ra_tot[")
    list(
      mn = array_format(tmp["mean",]),
      lwr = array_format(tmp["2.5%",]),
      upr = array_format(tmp["97.5%",])
    )
  }
)

# prep values to plot: summer to LGR survival
plot_dat$phi_Pb_Ma = list(
  obs = plogis(jags_data$Lphi_obs_Pb_Ma),
  fitted = {
    tmp = post_summ(post, "^phi_Pb_Ma[")
    list(
      mn = array_format(tmp["mean",]),
      lwr = array_format(tmp["2.5%",]),
      upr = array_format(tmp["97.5%",])
    )
  }
)

# prep values to plot: fall to LGR survival
plot_dat$phi_Pa_Ma_fall = list(
  obs = plogis(jags_data$Lphi_obs_Pa_Ma[,i_fall,]),
  fitted = {
    tmp = post_summ(post, sub_index("^phi_Pa_Ma[year,LH_type,pop]", year = ".+", LH_type = i_fall, pop = ".+"))
    list(
      mn = array_format(tmp["mean",])[,i_fall,],
      lwr = array_format(tmp["2.5%",])[,i_fall,],
      upr = array_format(tmp["97.5%",])[,i_fall,]
    )
  }
)

# prep values to plot: winter to LGR survival
plot_dat$phi_Pa_Ma_winter = list(
  obs = plogis(jags_data$Lphi_obs_Pa_Ma[,i_spring,]),
  fitted = {
    tmp = post_summ(post, sub_index("^phi_Pa_Ma[year,LH_type,pop]", year = ".+", LH_type = i_spring, pop = ".+"))
    list(
      mn = array_format(tmp["mean",])[,i_spring,],
      lwr = array_format(tmp["2.5%",])[,i_spring,],
      upr = array_format(tmp["97.5%",])[,i_spring,]
    )
  }
)

# prep values to plot: spring to LGR survival (NOR)
plot_dat$phi_Mb_Ma_nor = list(
  obs = plogis(jags_data$Lphi_obs_Mb_Ma[,i_spring,o_nor,]),
  fitted = {
    tmp = post_summ(post, sub_index("^phi_Mb_Ma[year,LH_type,origin,pop]", year = ".+", LH_type = i_spring, origin = o_nor, pop = ".+"))
    list(
      mn = array_format(tmp["mean",])[,i_spring,o_nor,],
      lwr = array_format(tmp["2.5%",])[,i_spring,o_nor,],
      upr = array_format(tmp["97.5%",])[,i_spring,o_nor,]
    )
  }
)

# prep values to plot: spring to LGR survival (HOR)
plot_dat$phi_Mb_Ma_hor = list(
  obs = plogis(jags_data$Lphi_obs_Mb_Ma[,i_spring,o_hor,]),
  fitted = {
    tmp = post_summ(post, sub_index("^phi_Mb_Ma[year,LH_type,origin,pop]", year = ".+", LH_type = i_spring, origin = o_hor, pop = ".+"))
    list(
      mn = array_format(tmp["mean",])[,i_spring,o_hor,],
      lwr = array_format(tmp["2.5%",])[,i_spring,o_hor,],
      upr = array_format(tmp["97.5%",])[,i_spring,o_hor,]
    )
  }
)

# prep values to plot: spring to LGR survival (NOR & HOR both)
plot_dat$phi_Mb_Ma = list(
  obs = plogis(jags_data$Lphi_obs_Mb_Ma[,i_spring,,]),
  fitted = {
    tmp = post_summ(post, sub_index("^phi_Mb_Ma[year,LH_type,.,pop]", year = ".+", LH_type = i_spring, origin = o_hor, pop = ".+"))
    list(
      mn = array_format(tmp["mean",])[,i_spring,,],
      lwr = array_format(tmp["2.5%",])[,i_spring,,],
      upr = array_format(tmp["97.5%",])[,i_spring,,]
    )
  }
)

# prep values to plot: LGR to BON (NOR & HOR both)
plot_dat$phi_Ma_O0 = list(
  obs = plogis(jags_data$Lphi_obs_Ma_O0),
  fitted = {
    tmp = post_summ(post, "^phi_Ma_O0[")
    list(
      mn = array_format(tmp["mean",]),
      lwr = array_format(tmp["2.5%",]),
      upr = array_format(tmp["97.5%",])
    )
  }
)

# prep values to plot: BON to LGR (NOR & HOR both)
plot_dat$phi_Rb_Ra = list(
  obs = jags_data$x_LGR/jags_data$x_BON,
  fitted = {
    tmp = post_summ(post, "^phi_Rb_Ra[")
    list(
      mn = array_format(tmp["mean",]),
      lwr = array_format(tmp["2.5%",]),
      upr = array_format(tmp["97.5%",])
    )
  }
)

#### CALCULATE AGGREGATED AGE/ORIGIN COMPS: MODEL
# extract the full posteriors of the two composition sets
p_Ra_full = post_subset(post, "p_Ra", matrix = TRUE)
p_Sa_prime_full = post_subset(post, "p_Sa", matrix = TRUE)

# primary containers: will store final recalculated output
p_Ra_age = p_Sa_prime_age = p_Ra_origin = p_Sa_prime_origin = NULL

# number of posterior samples
n = post_dim(post, "saved")

# loop through posterior samples, calculate various aggregates of composition by type
for (i in 1:n) {
  
  # format the posterior draw as arrays -- easier to subset
  p_Ra = array_format(p_Ra_full[i,])
  p_Sa_prime = array_format(p_Sa_prime_full[i,])
  
  ## BY AGE, AGGREGATED ACROSS ORIGINS
  
  # containers
  p_Ra_age_tmp = p_Sa_prime_age_tmp = array(NA, dim = c(jags_data$ny, jags_data$nk, jags_data$nj))
  inds_1_age = inds_2_age = inds_3_age = p_Ra_age_tmp
  for (j in 1:jags_data$nj) {
    for (k in 1:jags_data$nk) {
      # calculate aggregate proportions for each age
      p_Ra_age_tmp[,k,j] = rowSums(p_Ra[,ko_age[[k]],j])
      p_Sa_prime_age_tmp[,k,j] = rowSums(p_Sa_prime[,ko_age[[k]],j])
      
      # build the element identifiers - for building names later
      inds_1_age[,k,j] = 1:jags_data$ny
      inds_2_age[,k,j] = k
      inds_3_age[,k,j] = j
    }
  }
  
  # vectorize, add names, and remove NAs
  p_Ra_age_v = as.numeric(p_Ra_age_tmp)
  p_Sa_prime_age_v = as.numeric(p_Sa_prime_age_tmp)
  names(p_Ra_age_v) = paste0("p_Ra_age[", as.numeric(inds_1_age), ",", as.numeric(inds_2_age), ",", as.numeric(inds_3_age), "]")
  names(p_Sa_prime_age_v) = paste0("p_Sa_prime_age[", as.numeric(inds_1_age), ",", as.numeric(inds_2_age), ",", as.numeric(inds_3_age), "]")
  p_Ra_age_v = p_Ra_age_v[!is.na(p_Ra_age_v)]
  p_Sa_prime_age_v = p_Sa_prime_age_v[!is.na(p_Sa_prime_age_v)]
  
  ## BY ORIGIN, AGGREGATED ACROSS AGES
  
  # containers
  p_Ra_origin_tmp = p_Sa_prime_origin_tmp = array(NA, dim = c(jags_data$ny, jags_data$no, jags_data$nj))
  inds_1_origin = inds_2_origin = inds_3_origin = p_Ra_origin_tmp
  for (j in 1:jags_data$nj) {
    for (o in 1:jags_data$no) {
      # calculate aggregate proportions for each origin
      p_Ra_origin_tmp[,o,j] = rowSums(p_Ra[,ko_origin[[o]],j])
      p_Sa_prime_origin_tmp[,o,j] = rowSums(p_Sa_prime[,ko_origin[[o]],j])
      
      # build the element identifiers - for building names later
      inds_1_origin[,o,j] = 1:jags_data$ny
      inds_2_origin[,o,j] = o
      inds_3_origin[,o,j] = j
    }
  }
  
  # vectorize, add names, and remove NAs
  p_Ra_origin_v = as.numeric(p_Ra_origin_tmp)
  p_Sa_prime_origin_v = as.numeric(p_Sa_prime_origin_tmp)
  names(p_Ra_origin_v) = paste0("p_Ra_origin[", as.numeric(inds_1_origin), ",", as.numeric(inds_2_origin), ",", as.numeric(inds_3_origin), "]")
  names(p_Sa_prime_origin_v) = paste0("p_Sa_prime_origin[", as.numeric(inds_1_origin), ",", as.numeric(inds_2_origin), ",", as.numeric(inds_3_origin), "]")
  p_Ra_origin_v = p_Ra_origin_v[!is.na(p_Ra_origin_v)]
  p_Sa_prime_origin_v = p_Sa_prime_origin_v[!is.na(p_Sa_prime_origin_v)]
  
  ## COMBINE NEW CALCULATED QUANTITIES ACROSS POSTERIOR SAMPLES
  
  p_Ra_age = rbind(p_Ra_age, p_Ra_age_v)
  p_Sa_prime_age = rbind(p_Sa_prime_age, p_Sa_prime_age_v)
  p_Ra_origin = rbind(p_Ra_origin, p_Ra_origin_v)
  p_Sa_prime_origin = rbind(p_Sa_prime_origin, p_Sa_prime_origin_v)
  
}

# combine new calculated quantities into a big matrix
q_new = cbind(p_Ra_age, p_Sa_prime_age, p_Ra_origin, p_Sa_prime_origin)

# combine with the rest of the posterior samples
post = post_bind(post, q_new)

#### CALCULATE AGGREGATED AGE/ORIGIN COMPS: DATA

# containers
x_Ra_age_obs = x_Sa_prime_age_obs = array(NA, dim = c(jags_data$ny_obs, jags_data$nk, jags_data$nj))
x_Ra_origin_obs = x_Sa_prime_origin_obs = array(NA, dim = c(jags_data$ny_obs, jags_data$no, jags_data$nj))

for (j in 1:jags_data$nj) {
  # calculate aggregate proportions for each age
  for (k in 1:jags_data$nk) {
    x_Ra_age_obs[,k,j] = rowSums(jags_data$x_Ra[,ko_age[[k]],j])
    x_Sa_prime_age_obs[,k,j] = rowSums(jags_data$x_Sa_prime[,ko_age[[k]],j])
  }
  
  # calculate aggregate proportions for each origin
  for (o in 1:jags_data$no) {
    x_Ra_origin_obs[,o,j] = rowSums(jags_data$x_Ra[,ko_origin[[o]],j])
    x_Sa_prime_origin_obs[,o,j] = rowSums(jags_data$x_Sa_prime[,ko_origin[[o]],j])
  }
}

dimnames(x_Ra_age_obs)[[1]] = dimnames(x_Ra_origin_obs)[[1]] = dimnames(x_Sa_prime_age_obs)[[1]] = dimnames(x_Sa_prime_origin_obs)[[1]] = dimnames(jags_data$x_Ra)[[1]]

x_Ra_new_obs = list(
  x_Ra_age_obs = x_Ra_age_obs,
  x_Ra_origin_obs = x_Ra_origin_obs
)

x_Sa_prime_new_obs = list(
  x_Sa_prime_age_obs = x_Sa_prime_age_obs,
  x_Sa_prime_origin_obs = x_Sa_prime_origin_obs
)

# prep values to plot: weir compositon age-3
plot_dat$p_Ra_age3 = list(
  obs = x_Ra_new_obs$x_Ra_age_obs[,k_3,]/apply(x_Ra_new_obs$x_Ra_age_obs, 3, rowSums),
  fitted = {
    tmp = post_summ(post, sub_index("^p_Ra_age[year,1,pop]", year = ".+", pop = ".+"))
    list(
      mn = array_format(tmp["mean",])[,k_3,],
      lwr = array_format(tmp["2.5%",])[,k_3,],
      upr = array_format(tmp["97.5%",])[,k_3,]
    )
  }
)

# prep values to plot: weir compositon age-4
plot_dat$p_Ra_age4 = list(
  obs = x_Ra_new_obs$x_Ra_age_obs[,k_4,]/apply(x_Ra_new_obs$x_Ra_age_obs, 3, rowSums),
  fitted = {
    tmp = post_summ(post, sub_index("^p_Ra_age[year,2,pop]", year = ".+", pop = ".+"))
    list(
      mn = array_format(tmp["mean",])[,k_4,],
      lwr = array_format(tmp["2.5%",])[,k_4,],
      upr = array_format(tmp["97.5%",])[,k_4,]
    )
  }
)

# prep values to plot: weir compositon age-5
plot_dat$p_Ra_age5 = list(
  obs = x_Ra_new_obs$x_Ra_age_obs[,k_5,]/apply(x_Ra_new_obs$x_Ra_age_obs, 3, rowSums),
  fitted = {
    tmp = post_summ(post, sub_index("^p_Ra_age[year,3,pop]", year = ".+", pop = ".+"))
    list(
      mn = array_format(tmp["mean",])[,k_5,],
      lwr = array_format(tmp["2.5%",])[,k_5,],
      upr = array_format(tmp["97.5%",])[,k_5,]
    )
  }
)

# prep values to plot: weir compositon origin
plot_dat$p_Ra_nor = list(
  obs = x_Ra_new_obs$x_Ra_origin_obs[,o_nor,]/apply(x_Ra_new_obs$x_Ra_origin_obs, 3, rowSums),
  fitted = {
    tmp = post_summ(post, sub_index("^p_Ra_origin[year,1,pop]", year = ".+", pop = ".+"))
    list(
      mn = array_format(tmp["mean",])[,o_nor,],
      lwr = array_format(tmp["2.5%",])[,o_nor,],
      upr = array_format(tmp["97.5%",])[,o_nor,]
    )
  }
)

# prep values to plot: carcass compositon age-3
plot_dat$p_Sa_prime_age3 = list(
  obs = x_Sa_prime_new_obs$x_Sa_prime_age_obs[,k_3,]/apply(x_Sa_prime_new_obs$x_Sa_prime_age_obs, 3, rowSums),
  fitted = {
    tmp = post_summ(post, sub_index("^p_Sa_prime_age[year,1,pop]", year = ".+", pop = ".+"))
    list(
      mn = array_format(tmp["mean",])[,k_3,],
      lwr = array_format(tmp["2.5%",])[,k_3,],
      upr = array_format(tmp["97.5%",])[,k_3,]
    )
  }
)

# prep values to plot: carcass compositon age-4
plot_dat$p_Sa_prime_age4 = list(
  obs = x_Sa_prime_new_obs$x_Sa_prime_age_obs[,k_4,]/apply(x_Sa_prime_new_obs$x_Sa_prime_age_obs, 3, rowSums),
  fitted = {
    tmp = post_summ(post, sub_index("^p_Sa_prime_age[year,2,pop]", year = ".+", pop = ".+"))
    list(
      mn = array_format(tmp["mean",])[,k_4,],
      lwr = array_format(tmp["2.5%",])[,k_4,],
      upr = array_format(tmp["97.5%",])[,k_4,]
    )
  }
)

# prep values to plot: weir compositon age-5
plot_dat$p_Sa_prime_age5 = list(
  obs = x_Sa_prime_new_obs$x_Sa_prime_age_obs[,k_5,]/apply(x_Sa_prime_new_obs$x_Sa_prime_age_obs, 3, rowSums),
  fitted = {
    tmp = post_summ(post, sub_index("^p_Sa_prime_age[year,3,pop]", year = ".+", pop = ".+"))
    list(
      mn = array_format(tmp["mean",])[,k_5,],
      lwr = array_format(tmp["2.5%",])[,k_5,],
      upr = array_format(tmp["97.5%",])[,k_5,]
    )
  }
)

# prep values to plot: weir compositon origin
plot_dat$p_Sa_prime_nor = list(
  obs = x_Sa_prime_new_obs$x_Sa_prime_origin_obs[,o_nor,]/apply(x_Sa_prime_new_obs$x_Sa_prime_origin_obs, 3, rowSums),
  fitted = {
    tmp = post_summ(post, sub_index("^p_Sa_prime_origin[year,1,pop]", year = ".+", pop = ".+"))
    list(
      mn = array_format(tmp["mean",])[,o_nor,],
      lwr = array_format(tmp["2.5%",])[,o_nor,],
      upr = array_format(tmp["97.5%",])[,o_nor,]
    )
  }
)

# function to create plot for one data type
plot_f = function(dat_list, f, label = "", adj_label = FALSE) {
  
  with(dat_list, {
    
    # handle outlying values
    # obs[obs == 0] = NA
    # obs[obs == 1] = NA
    obs[obs == "NaN"] = NA
    
    # apply transformations
    obs = f(obs)
    mn = f(fitted$mn)
    lwr = f(fitted$lwr)
    upr = f(fitted$upr)
    
    # set fitted value to NA if no data value
    mn[is.na(obs)] = NA
    lwr[is.na(obs)] = NA
    upr[is.na(obs)] = NA
    
    # obtain axis limits
    lim = range(mn, lwr, upr, obs, na.rm = TRUE)
    
    # create the plot
    plot(mn ~ obs, xlim = lim, ylim = lim, pch = 21, col = solid_col, bg = tranp_col, cex = 1.5)
    
    # draw fitted value CRIs
    # segments(obs, lwr, obs, upr)
    
    # draw 1:1 line
    abline(0,1, lty = 2)
    
    # add text
    panel_label(label, x_inp = ifelse(adj_label, 0.05, 0.015), cex = 0.95)
  })
}

# identity function
ifun = function(x) x

# open the device
dev.on(file.path(this_dir, "obs-v-fit"), height = 7, width = 7.2, format = fig_type)

# graphical parameters
mypar(mfrow = c(1,1), oma = c(1.5,1.5,0,1.5), col.axis = "black", mar = c(1,0.75,0.5,0.25))

# set up the layout
r1 = rep(1:3, each = 4); r2 = rep(4:9, each = 2); r3 = rep(10:13, each = 3); r4 = rep(14:17, each = 3)
layout(rbind(r1, r2, r3, r4))

# abundance plots
plot_f(plot_dat[["Pa_fall"]], ifun, "(a) Fall Trap Passage")
plot_f(plot_dat[["Mb_spring"]], ifun, "(b) Spring Trap Passage")
plot_f(plot_dat[["Ra_tot"]], ifun, "(c) Adult Return to River"); mtext(side = 4, text = "Abundance", outer = FALSE, line = 0.55)

# survival plots
plot_f(plot_dat[["phi_Pb_Ma"]], ifun, "(d) Summer \u2192 LGR", TRUE)
plot_f(plot_dat[["phi_Pa_Ma_fall"]], ifun, "(e) Fall \u2192 LGR", TRUE)
plot_f(plot_dat[["phi_Pa_Ma_winter"]], ifun, "(f) Winter \u2192 LGR", TRUE)
plot_f(plot_dat[["phi_Mb_Ma"]], ifun, "(g) Spring \u2192 LGR", TRUE)#; mtext(side = 4, text = "Survival", outer = FALSE, line = 0.75)
plot_f(plot_dat[["phi_Ma_O0"]], ifun, "(h) LGR \u2192 BON", TRUE)
plot_f(plot_dat[["phi_Rb_Ra"]], ifun, "(i) BON \u2192 LGR", TRUE); mtext(side = 4, text = "Survival", outer = FALSE, line = 0.55)

# composition plots
plot_f(plot_dat[["p_Ra_age3"]], ifun, "(j) Proportion Age-3")
plot_f(plot_dat[["p_Ra_age4"]], ifun, "(k) Proportion Age-4")
plot_f(plot_dat[["p_Ra_age5"]], ifun, "(l) Proportion Age-5")
plot_f(plot_dat[["p_Ra_nor"]], ifun, "(m) Proportion NOR"); mtext(side = 4, text = "Weir Composition", outer = FALSE, line = 0.55)

# composition plots
plot_f(plot_dat[["p_Sa_prime_age3"]], ifun, "(n) Proportion Age-3")
plot_f(plot_dat[["p_Sa_prime_age4"]], ifun, "(o) Proportion Age-4")
plot_f(plot_dat[["p_Sa_prime_age5"]], ifun, "(p) Proportion Age-5")
plot_f(plot_dat[["p_Sa_prime_nor"]], ifun, "(q) Proportion NOR"); mtext(side = 4, text = "Carcass Composition", outer = FALSE, line = 0.55)

# axis labels
axis_labels("Observed Value", "Fitted Value", yline = 0.25)

# close the device
dev.off()
