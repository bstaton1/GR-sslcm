
# load model info if not already in workspace
if (!exists("post")) source(file.path(this_dir, "load-model-info.R"))

#### SET INFO FOR COMPARE PLOTS ####

pch = c(21, 23)
cex_main = 2.5
cex_sub = 1.1

f_horiz = function(ests1, ests2, lim, label) {
  # prepare the first group of estimates
  ests1 = ests1[-1,]; x1 = colMeans(ests1)
  
  # figure out the correct axis labels
  if (ncol(ests1) == 4) tick_labels = pops else tick_labels = c("NOR", "HOR")
  jitter = runif(nrow(ests1), -0.25, 0.25)
  
  if (is.null(ests2)) {
    # empty barplot with correct dimensions
    mp = barplot(x1, col = "white", border = "white", xlim = lim, ylim = c(ncol(ests1) + 0.5, 0.05), horiz = TRUE)
    
    # figure out location of jittered points
    mp_all = t(sapply(1:nrow(ests1), function(i) mp))
    mp_all = apply(mp_all, 2, function(x) x + jitter)
    
    # draw year-specific estimates
    points(ests1, mp_all, col = solid_col, bg = tranp_col, pch = pch[1], cex = cex_sub)
    
    # draw year averages
    points(x1, mp, pch = pch[1], col = "white", bg = "black", cex = cex_main, lwd = 1.5)
    
    # draw axis
    omgp = par("mgp"); nmgp = omgp; nmgp[2] = 0.15; otcl = par("tcl")
    par(mgp = nmgp, tcl = 0)
    usr = par("usr"); xdiff = diff(usr[1:2]); ydiff = diff(usr[3:4])
    segments(usr[1], usr[3], usr[1], usr[4], xpd = TRUE, col = par("col.axis"))
    axis(side = 2, at = mp, labels = tick_labels, las = 1)
    par(mgp = omgp, tcl = otcl)
    
  } else {
    # prepare the second group of estimates if supplied
    ests2 = ests2[-1,]; x2 = colMeans(ests2)
    
    # empty barplot with correct dimensions
    mp = barplot(rbind(x1, x2), beside = TRUE, space = c(0.25,1), col = "white", border = "white", xlim = lim, ylim = c(ncol(ests1) * 3 + 1, 0.15), horiz = TRUE)
    mp1 = mp[1,]; mp2 = mp[2,]
    
    # figure out location of jittered points
    mp1_all = t(sapply(1:nrow(ests1), function(i) mp1))
    mp1_all = apply(mp1_all, 2, function(x) x + jitter)
    mp2_all = t(sapply(1:nrow(ests2), function(i) mp2))
    mp2_all = apply(mp2_all, 2, function(x) x + jitter)
    
    # draw year-specific estimates
    points(ests1, mp1_all, col = solid_col, bg = tranp_col, pch = pch[1], cex = cex_sub)
    points(ests2, mp2_all, col = solid_col, bg = tranp_col, pch = pch[2], cex = cex_sub)
    
    # draw year averages
    points(x1, mp1, pch = pch[1], cex = cex_main, bg = "black", col = "white", lwd = 1.5)
    points(x2, mp2, pch = pch[2], cex = cex_main, bg = "black", col = "white", lwd = 1.5)
    
    # draw axis
    omgp = par("mgp"); nmgp = omgp; nmgp[2] = 0.15; otcl = par("tcl")
    par(mgp = nmgp, tcl = 0)
    usr = par("usr"); xdiff = diff(usr[1:2]); ydiff = diff(usr[3:4])
    segments(usr[1], usr[3], usr[1], usr[4], xpd = TRUE, col = par("col.axis"))
    axis(side = 2, at = (mp1 + mp2)/2, labels = tick_labels, las = 1)
    par(mgp = omgp, tcl = otcl)
  }
  
  # draw panel label
  panel_label(label, y_inp = ifelse(ncol(ests1) == 4, 0.05, 0.1), cex = 1.2)
}

legend_f = function(title, grp1, grp2, loc = "topright") {
  legend(loc, title = title, legend = c(grp1, grp2), pch = pch, pt.cex = 1.3, cex = 0.9, bty = "n",
         text.col = par("col.axis"), col = "black", pt.bg = "black")
}

#### CREATE SURVIVAL PARAMETER COMPARE PLOT ####

# extract the posterior means to plot
phi_E_Pb = array_format(post_summ(post, sub_index("^phi_E_Pb[year,pop]", year = observable, pop = ".+"))["mean",])
phi_Pa_Mb_fall = array_format(post_summ(post, sub_index("^phi_Pa_Mb[year,LH_type,pop]", year = observable, pop = ".+", LH_type = i_fall))["mean",])[,i_fall,]
phi_Pa_Mb_spring = array_format(post_summ(post, sub_index("^phi_Pa_Mb[year,LH_type,pop]", year = observable, pop = ".+", LH_type = i_spring))["mean",])[,i_spring,]
phi_Mb_Ma_nor = array_format(post_summ(post, sub_index("^phi_Mb_Ma[year,LH_type,origin,pop]", year = observable, pop = ".+", LH_type = i_spring, origin = o_nor))["mean",])[,i_spring,o_nor,]
phi_Mb_Ma_hor = array_format(post_summ(post, sub_index("^phi_Mb_Ma[year,LH_type,origin,pop]", year = observable, pop = ".+", LH_type = i_spring, origin = o_hor))["mean",])[,i_spring,o_hor,]
phi_Ma_O0 = array_format(post_summ(post, sub_index("^phi_Ma_O0[year,origin]", year = observable, origin = ".+"))["mean",])
phi_O0_O1_nor = array_format(post_summ(post, sub_index("^phi_O0_O1[year,origin,pop]", year = observable, pop = ".+", origin = o_nor))["mean",])[,o_nor,]
phi_O0_O1_hor = array_format(post_summ(post, sub_index("^phi_O0_O1[year,origin,pop]", year = observable, pop = ".+", origin = o_hor))["mean",])[,o_hor,]
phi_Rb_Ra = array_format(post_summ(post, sub_index("^phi_Rb_Ra[year,origin]", year = observable, origin = ".+"))["mean",])

# set MIN to NA for HOR values
phi_Mb_Ma_hor[,j_min] = NA
phi_O0_O1_hor[,j_min] = NA

# open the graphics device
dev.on(file.path(this_dir, "surv-compare"), width = 3.4, height = 8, format = fig_type)

# graphical parameters/layout
mypar(mfrow = c(1,1), oma = c(1.5,1,0,1), col.axis = "black", tcl = -0.1)
layout(matrix(1:5, ncol = 1), height = c(1,1,1,0.5,1))

# egg survival panel
f_horiz(ests1 = phi_E_Pb, ests2 = NULL, lim = c(0, 0.5), label = "(a) Egg \u2192 Parr")

# overwinter survival panel
f_horiz(ests1 = phi_Pa_Mb_fall, ests2 = phi_Pa_Mb_spring, lim = c(0,1), label = "(b) Parr \u2192 Smolt")
legend_f("Mig. Type", "Fall", "Spring")

# outmigration survival panel
f_horiz(ests1 = phi_Mb_Ma_nor, ests2 = phi_Mb_Ma_hor, lim = c(0,1), label = "(c) Smolt \u2192 LGR")
legend_f("Origin", "NOR", "HOR")

# mainstem survival panel
f_horiz(ests1 = phi_Ma_O0, ests2 = phi_Rb_Ra, lim = c(0,1), label = "(d) Mainstem")
legend_f("Direction", "LGR \u2192 BON", "BON \u2192 LGR", "bottomleft")

# ocean survival panel
f_horiz(ests1 = phi_O0_O1_nor, ests2 = phi_O0_O1_hor, lim = c(0,0.4), label = "(e) Ocean-0 \u2192 Ocean-1")
legend_f("Origin", "NOR", "HOR", "bottomright")
axis_labels("Survival Probability")

# close the graphics device
dev.off()

#### CREATE APPORTION PARAMETER COMPARE PLOT ####

# extract the posterior means to plot
pi_fall = array_format(post_summ(post, sub_index("^pi[year,LH_type,pop]", year = observable, pop = ".+", LH_type = i_fall))["mean",])[,i_fall,]
psi_O1_nor = array_format(post_summ(post, sub_index("^psi_O1[year,origin,pop]", year = observable, pop = ".+", origin = o_nor))["mean",])[,o_nor,]
psi_O1_hor = array_format(post_summ(post, sub_index("^psi_O1[year,origin,pop]", year = observable, pop = ".+", origin = o_hor))["mean",])[,o_hor,]
psi_O2_nor = array_format(post_summ(post, sub_index("^psi_O2[year,origin,pop]", year = observable, pop = ".+", origin = o_nor))["mean",])[,o_nor,]
psi_O2_hor = array_format(post_summ(post, sub_index("^psi_O2[year,origin,pop]", year = observable, pop = ".+", origin = o_hor))["mean",])[,o_hor,]

# set MIN to NA for HOR values
psi_O1_hor[,j_min] = NA
psi_O2_hor[,j_min] = NA

# open the graphics device
dev.on(file.path(this_dir, "apportion-compare"), width = 3.4, height = 8 * (3/4.5), format = fig_type)

# set graphical parameters/layout
mypar(mfrow = c(1,1), oma = c(1.5,1,0,1), col.axis = "black", tcl = -0.1)
layout(matrix(1:3, ncol = 1), height = c(1,1,1))

# migratory strategy apportionment
f_horiz(ests1 = pi_fall, ests2 = NULL, lim = c(0, 0.6), label = "(a) Pr(Fall Migrant)")

# Pr(Mature Age-3)
f_horiz(ests1 = psi_O1_nor, ests2 = psi_O1_hor, lim = c(0,0.4), label = "(b) Pr(Mature Age-3)")
legend_f("Origin", "NOR", "HOR")

# Pr(Mature Age-4)
f_horiz(ests1 = psi_O2_nor, ests2 = psi_O2_hor, lim = c(0,1), label = "(c) Pr(Mature Age-4)")
legend_f("Origin", "NOR", "HOR", "bottomleft")

axis_labels("Transition Probability (Non-Survival)")

# close the graphics device
dev.off()
