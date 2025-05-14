# load packages
source(file.path(proj_dir, "00-packages.R"))

# load all necessary functions
invisible(sapply(list.files(path = file.path(proj_dir, "01-functions"), pattern = "\\.R$", full.names = TRUE), source))

# read information from this model
model_info = readRDS(file.path(proj_dir, "02-model/model-output", paste0("output_", scenario, ".rds")))

# extract the posterior samples
# post-MCMC thinning allowed to reduce run time of this document
if (keep_percent < 1) {
  post = post_thin(model_info$post, keep_percent = keep_percent)
} else {
  post = model_info$post
}

# extract the model data 
jags_data = model_info$jags_data

# remove the model_info object: saves memory
rm(model_info)

# the years included in the model
all_yrs = as.numeric(dimnames(jags_data$Mb_obs)[[1]])

# which years are "observable" for at least some quantities?
observable = ts_yrs = 2:jags_data$ny_obs

# which years can be plotted for SARs? Drop off last kmax years
# since brood year adult returns not complete
sar_yrs = ts_yrs[-((length(ts_yrs) - jags_data$kmax + 1):length(ts_yrs))]

# population names
pops = colnames(jags_data$Ra_obs)

# origin names
origins = c("NOR", "HOR")

# dimension IDs
i_fall   = 1  # fall migrants are i = 1
i_spring = 2  # spring migrants are i = 2,
o_nor    = 1  # natural origin are o = 1,
o_hor    = 2  # hatchery origin are o = 2,
j_cat    = 1  # Catherine Creek is j = 1
j_los    = 2  # Lostine River is j = 2
j_min    = 3  # Minam River is j = 3
j_ugr    = 4  # Upper Grande Ronde River is j = 4
k_3      = 1  # age 3 is k = 1
k_4      = 2  # age 4 is k = 2
k_5      = 3  # age 5 is k = 3

# ko represents "age/origin" combo
# e.g., ko = 1 is age 3 NORs
# e.g., ko = 5 is age 4 HORs
# these objects specify which elements of ko are for different aggregations of these
ko_age = list(
  ko_3 = c(1, 4),
  ko_4 = c(2, 5),
  ko_5 = c(3, 6)
)

ko_origin = list(
  ko_nor = 1:3,
  ko_hor = 4:6
)

# make the labels for which elements of the rho matrices contain unique elements
dummy_cols = matrix(rep(1:jags_data$nj, each = jags_data$nj), jags_data$nj, jags_data$nj)
dummy_rows = matrix(rep(1:jags_data$nj, jags_data$nj), jags_data$nj, jags_data$nj)
vcov_cols = dummy_cols[lower.tri(dummy_cols)]
vcov_rows = dummy_rows[lower.tri(dummy_rows)]
vcov_inds = paste0("[", vcov_rows, ",", vcov_cols, "]")
vcov_labels = paste0(pops[vcov_rows], "-", pops[vcov_cols])

# plot settings
solid_col = "grey30"
solid_col2 = "black"
tranp_col = alpha(solid_col, 0.25)
pt_cex = 1.3
