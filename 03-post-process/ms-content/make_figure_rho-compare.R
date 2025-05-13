
# load model info if not already in workspace
if (!exists("post")) source(file.path(this_dir, "load-model-info.R"))

# get posterior summaries of the mean correlation across all pop pairs
# for each process separately
rho_params1 = match_params(post, "rho", type = "base_only")
rho_params2 = match_params(post, "rho_.+_pr", type = "base_only")
rho_params = rho_params1[!(rho_params1 %in% rho_params2)]

# function to extract posterior summaries of the rho parameters for a given process
f = function(p) {
  x = rowMeans(post_subset(post, paste0("^", p, "["), matrix = TRUE))
  c(mean = mean(x), sd = sd(x), quantile(x, c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)))
}

# apply it and order
rho_mean_out = t(sapply(rho_params, f))
rho_mean_out = rho_mean_out[order(rho_mean_out[,"mean"]),]

# assign names to the different rho parameters
rho_labels = c(
  "rho_Lphi_E_Pb" = "Egg \u2192 Parr",
  "rho_Lphi_O0_O1" = "Ocean-0 \u2192 Ocean-1",
  "rho_Lpsi_O1" = "Pr(Mature Age-3)",
  "rho_Lpsi_O2" = "Pr(Mature Age-4)",
  "rho_lDelta_L_Pb_Mb" = "\u0394 Mean Length",
  "rho_Lphi_Mb_Ma" = "Smolt \u2192 LGR",
  "rho_lL_Pb" = "Mean Parr Length",
  "rho_Lphi_Pa_Mb" = "Parr \u2192 Smolt",
  "rho_Lphi_Rb_Ra" = "BON \u2192 LGR",
  "rho_Lpi" = "Pr(Fall Migrant)",
  "rho_Lphi_Ma_O0" = "LGR \u2192 BON"
)

# open a graphics device
dev.on(file.path(this_dir, "rho-compare"), width = 3.4, height = 3.4, format = fig_type)

# set up empty plotting region
mypar(mar = c(2,6,0.75,0.75), yaxs = "i", mfrow = c(1,1), oma = c(0,0,0,0), mgp = c(2,0.1,0), col.axis = "black", cex.axis = 0.75)
mp = barplot(rho_mean_out[,"mean"], 
             horiz = TRUE, xlim = c(-1,1),
             names.arg = rho_labels[rownames(rho_mean_out)],
             las = 1, col = "white", border = "white", xaxt = "n")
usr = par("usr")

box_upr = mp - rbind(diff(mp)/2, NA)
box_lwr = mp + rbind(diff(mp)/2, NA)
box_upr[nrow(box_upr),] = box_upr[nrow(box_upr)-1,] + diff(box_upr[1:2])
box_lwr[nrow(box_lwr),] = box_lwr[nrow(box_upr)-1,] + diff(box_lwr[1:2])

pop_rhos = function(param, at_mp) {
  set.seed(1)
  # if the process is not one of the hydropower migration survival ones, do this
  # those processes are correlated across origins, because all populations share the same values each year
  if (!(param %in% c("rho_Lphi_Ma_O0", "rho_Lphi_Rb_Ra"))) {
    x = post_summ(post, paste0(param, vcov_inds), probs = c(0.025, 0.1, 0.25, 0.75, 0.9, 0.975))
    y = at_mp + runif(6, -0.35, 0.35)
    # segments(x["2.5%",], y, x["97.5%",], y, col = solid_col, bg = tranp_col)
    # points(y ~ x["mean",], pch = 21, col = "white", bg = "white")
    points(y ~ x["mean",], pch = 21, col = solid_col, bg = tranp_col)
  } 
}
abline(v = 0, lty = 1, lwd = 2, col = "grey")

junk = sapply(1:nrow(mp), function(i) pop_rhos(rownames(rho_mean_out)[i], mp[i,]))
segments(rho_mean_out[,"25%"], mp, rho_mean_out[,"75%"], mp, lwd = 5, col = "white")
segments(rho_mean_out[,"2.5%"], mp, rho_mean_out[,"97.5%"], mp, col = "white", lwd = 2)
points(mp ~ rho_mean_out[,"mean"], pch = 3, cex = 1, col = "white", lwd = 3)

segments(rho_mean_out[,"25%"], mp, rho_mean_out[,"75%"], mp, lwd = 4, col = "black")
segments(rho_mean_out[,"2.5%"], mp, rho_mean_out[,"97.5%"], mp, col = "black", lwd = 1)
points(mp ~ rho_mean_out[,"mean"], pch = 3, cex = 1, col = "black", lwd = 2)

abline(h = mp + rbind(diff(mp)/2, NA), xpd = FALSE, col = "grey", lty = 2)
box(col = par("col.axis"))
axis_labels(xlab = "Process Noise Correlation", xline = 1, outer = FALSE)

par(mgp = c(2,0,0))
axis(side = 1)

# close the device
dev.off()
