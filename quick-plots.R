# load packages
source("00-packages.R")

# load all necessary functions
invisible(sapply(list.files(path = "01-functions", pattern = "\\.R$", full.names = T), source))

params = list(
  scenario = "fecund_long",
  keep_percent = 1
)

# set the input directory
in_dir = "02-model/model-output"

# read information from this model
cat("  Preparing Model Output for Analysis\n")
model_info = readRDS(file.path(in_dir, paste0("output_", params$scenario, ".rds")))

# extract the posterior samples
# post-MCMC thinning allowed to reduce run time of this document
if (params$keep_percent < 1) {
  post = post_thin(model_info$post, keep_percent = params$keep_percent)
} else {
  post = model_info$post
}

# extract the model data 
jags_data = model_info$jags_data

# the years included by the model
all_yrs = as.numeric(dimnames(jags_data$phi_SL)[[1]])

# which years are "observable" for at least some quantities?
observable = 2:jags_data$ny_obs

# population names
pops = colnames(jags_data$Ra_obs)

plot_f = function(y,j) {
  x = post_subset(post, sub_index(c("^phi_E_Pb[year,pop]", "^phi_Sb_Sa[year,pop]"), year = y, pop = j), matrix = TRUE)
  plot(x[,1] ~ x[,2], col = scales::alpha("grey20", 0.05), pch = 16,
       ylim = make_lim(x[,1]), cex = 0.8)
  panel_label(pops[j])
}

pdf("quick-plots.pdf", width = 4.5, height = 4.5)
for (y in observable) {
  cat("\r    Plotting Brood Year:", all_yrs[y])
  mypar(oma = c(1.5,1.5,1,0))
  sapply(1:4, plot_f, y = y)
  axis_labels(xlab = "Pre-Spawn Survival", ylab = "Egg-to-Parr Survival")
  mtext(side = 3, outer = TRUE, line = -0.25, paste("Brood Year:", all_yrs[y]), font = 2, adj = 0)
}
junk = dev.off()
