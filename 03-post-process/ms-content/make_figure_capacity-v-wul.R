
# load model info if not already in workspace
if (!exists("post")) source(file.path(this_dir, "load-model-info.R"))

# extract posterior summaries of the capacity parameters
beta_ests = post_summ(post, "^beta[", probs = c(0.025, 0.1, 0.5, 0.9, 0.975))

# extract posterior samples of the slope parameter
beta_per = post_subset(post, "lambda", TRUE)

# create a sequence of WUL values to predict capacity at
wul_seq = seq(0, max(jags_data$wul) * 1.1, length = 30)

# obtain the regression line for each posterior sample
pred_beta = t(sapply(beta_per, function(slope) wul_seq * slope))

# summarize posterior of regression line
pred_beta_summ = apply(pred_beta, 2, function(x) c(mean = mean(x), quantile(x, c(0.025, 0.1, 0.5, 0.9, 0.975))))

# open a graphics device
dev.on(file.path(this_dir, "capacity-v-wul"), width = 3.4, height = 3.4, format = fig_type)

# empty plot with correct labeling/dimensions
mypar(mfrow = c(1,1), col.axis = "black", tcl = -0.1)
plot(1,1, type = "n", xlim = range(wul_seq),
     # ylim = c(0, max(beta_ests["97.5%",], pred_beta_summ["97.5%",])),
     ylim = c(0, 3000) * 1000,
     xlab = "",
     ylab = "", yaxt = "n")

# draw prettier y-axis
at_y = c(0, 1000, 2000, 3000) * 1000
axis(side = 2, at = at_y, labels = at_y/1000)
usr = par("usr"); xdiff = diff(usr[1:2]); ydiff = diff(usr[3:4])

# draw the regression uncertainty
polygon(x = c(wul_seq, rev(wul_seq)), y = c(pred_beta_summ["2.5%",], rev(pred_beta_summ["97.5%",])), border = NA, col = tranp_col)
polygon(x = c(wul_seq, rev(wul_seq)), y = c(pred_beta_summ["10%",], rev(pred_beta_summ["90%",])), border = NA, col = tranp_col)

# draw the posterior mean regression line
lines(pred_beta_summ["mean",] ~ wul_seq, lwd = 2, lty = 2, col = solid_col2)

# draw the population-specific capacity estimates w/error bars and labels
segments(jags_data$wul, beta_ests["2.5%",], jags_data$wul, beta_ests["97.5%",], col = solid_col2)
segments(jags_data$wul, beta_ests["10%",], jags_data$wul, beta_ests["90%",], col = solid_col2, lwd = 3)
points(beta_ests["mean",] ~ jags_data$wul, cex = 1.5, pch = 21, col = solid_col2, bg = solid_col2)

# draw the population labels
xoff = 0.015
yoff = 0.03
text(x = jags_data$wul      + xdiff * c(1.5,-1,1,1) * xoff,
     y = beta_ests["mean",] + ydiff * c(1,-1,1,1) * yoff,
     labels = c("CAT", "LOS", "MIN", "UGR"), pos = c(2,4,2,2), font = 1, cex = 0.8)

# add axis labels
axis_labels("Weighted Usable Habitat Length (km)", "Parr Capacity (Thousands)")

# close the device
dev.off()
