
# means to compare 
# * p_HOS
# * Ra
# * psi_O2_Rb[,"NOR",LOS"] (to check scale of maturity problem)
# * BON to BON SAR (to check consequence of maturity problem)
# * BON to LGR SAR (to check consequence of maturity problem)
# * q_Ra
# * smolt per spawner (Mb_per_Sa_tot)

# SDs to compare
# * BON to BON SAR
# * BON to LGR SAR
# * smolt per spawner (Mb_per_Sa_tot)

# load all necessary functions
invisible(sapply(list.files(path = "01-functions", pattern = "\\.R$", full.names = T), source))

# read in the model input/output
model_info = readRDS("02-model/model-output/output-base-with-sim-medium.rds")

# extract the posterior samples
# post = model_info$post
post = post_thin(model_info$post, keep_percent = 0.1)

# extract the data object passed to JAGS
jags_data = model_info$jags_data

obs_y = with(jags_data, (kmax+1):ny_obs)
sim_y = with(jags_data, (ny_obs+1):max(ny))
obs_sar_y = with(jags_data, (kmax+1):(ny_obs-kmax))
sim_sar_y = with(jags_data, (ny_obs+1):max(ny-kmax))

# function to add p_hor to posterior samples
add_p_hor = function() {
  Ra = post_subset(post, "^Ra[", matrix = TRUE)
  
  out = NULL
  for (i in 1:nrow(Ra)) {
    Ra_i = array_format(Ra[i,])
    nor = Ra_i[,1,1,] + Ra_i[,2,1,] + Ra_i[,3,1,]
    hor = Ra_i[,1,2,] + Ra_i[,2,2,] + Ra_i[,3,2,]
    out = rbind(out, as.numeric(hor/(nor + hor)))
  }
  
  r_ind = rep(1:nrow(nor), ncol(nor))
  c_ind = rep(1:ncol(nor), each = nrow(nor))
  colnames(out) = paste0("p_hor[", r_ind, ",", c_ind, "]")
  na_ind = which(!is.na(out[1,]))
  out = out[,na_ind]
  
  post = post_bind(post, out)
  post
}

dim_names = function(LH_type = NULL, age = NULL, origin = NULL, pop = NULL) {
  
  out = list()
  
  if (!is.null(LH_type)) out = append(out, list(LH_type = c("fall_mig", "spring_mig")[LH_type]))
  if (!is.null(age)) out = append(out, list(age = c("3", "4", "5")[age]))
  if (!is.null(origin)) out = append(out, list(origin = c("NOR", "HOR")[origin]))
  if (!is.null(pop)) out = append(out, list(pop = c("CAT", "LOS", "MIN", "UGR")[pop]))
  
  return(out)
}

get_sim_to_obs_ratio = function(param, obs_y, sim_y, ...) {
  obs_param = sub_index(param, year = obs_y, ...)
  sim_param = sub_index(param, year = sim_y, ...)
  obs_x = post_subset(post, obs_param, matrix = TRUE)
  sim_x = post_subset(post, sim_param, matrix = TRUE)
  # cbind(
  #   cv_ratio = sapply(1:nrow(obs_x), function(i) (sd(sim_x[i,])/mean(sim_x[i,]))/(sd(obs_x[i,])/mean(obs_x[i,]))),
  #   mean_ratio = sapply(1:nrow(obs_x), function(i) mean(sim_x[i,])/mean(obs_x[i,]))
  # )
  
  named_param = do.call(sub_index, append(dim_names(...), list(x = param)))
  named_param = stringr::str_remove(named_param, "year,")
  named_param = stringr::str_extract(named_param, "\\[.+\\]")
  
  out = array(NA, dim = c(nrow(obs_x), 2, 1))
  dimnames(out) = list(iter = 1:nrow(obs_x), ratio_type = c("cv", "mean"), dims = named_param)
  out[,"cv",] = sapply(1:nrow(obs_x), function(i) (sd(sim_x[i,])/mean(sim_x[i,]))/(sd(obs_x[i,])/mean(obs_x[i,])))
  out[,"mean",] = sapply(1:nrow(obs_x), function(i) mean(sim_x[i,])/mean(obs_x[i,]))
  return(out)
}

my_plot = function(samps, ylab, ...) {
  
  summs = apply(samps, 2, function(x) quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE))
  
  if (ylab == "cv") {
    ylab = "CV(sim)/CV(obs)"
  } else{
    ylab = "mean(sim)/mean(obs)"
  }
  
  par(lend = "square", ljoin = "mitre")
  par(mar = c(5,3,1,1), mgp = c(2,0.35,0), tcl = -0.15)
  mp = barplot(summs["50%",], col = "white", ylab = ylab, border = "white", las = 2,
               ylim = c(0, max(summs, na.rm = TRUE)) * 1.1)
  usr = par("usr")
  axis(side = 1, at = mp, labels = F)
  rect(usr[1], 0.75, usr[2], 1.25, col = "grey90", border = NA)
  abline(h = 1, lty = 2, col = "grey")
  segments(mp, summs["2.5%",], mp, summs["97.5%",])
  segments(mp, summs["25%",], mp, summs["75%",], lwd = 4)
  # points(summs["50%",] ~ mp, pch = 23, col = "white", bg = "black", cex = 1.5, lwd = 2)
  points(summs["50%",] ~ mp, pch = 3, cex = 1.5)
  box()
}

post = add_p_hor()

##### Yr1 Ocean Survival #####
x1 = lapply(1:4, function(j) get_sim_to_obs_ratio("phi_O0_O1[year,origin,pop]", obs_y, sim_y, origin = 1, pop = j))
x2 = lapply(1:4, function(j) get_sim_to_obs_ratio("phi_O0_O1[year,origin,pop]", obs_y, sim_y, origin = 2, pop = j))
all_x = do.call(abind, append(x1, append(x2, list(along = 3))))

par(mfrow = c(1,2), oma = c(0,0,1,0))
my_plot(all_x[,"cv",], ylab = "cv")
my_plot(all_x[,"mean",], ylab = "mean")
mtext(side = 3, outer = TRUE, line = -0.5, "First Year Ocean Survival", cex = 1.25, font = 2)

##### SAR: BON to BON #####
x1 = lapply(1:4, function(j) get_sim_to_obs_ratio("phi_O0_Rb_BON[year,origin,pop]", obs_sar_y, sim_sar_y, pop = j, origin = 1))
x2 = lapply(1:4, function(j) get_sim_to_obs_ratio("phi_O0_Rb_BON[year,origin,pop]", obs_sar_y, sim_sar_y, pop = j, origin = 2))
all_x = do.call(abind, append(x1, append(x2, list(along = 3))))

par(mfrow = c(1,2), oma = c(0,0,1,0))
my_plot(all_x[,"cv",], ylab = "cv")
my_plot(all_x[,"mean",], ylab = "mean")
mtext(side = 3, outer = TRUE, line = -0.5, "BON to BON Survival Rate", cex = 1.25, font = 2)

##### Proportion Hatchery Returns #####
x = lapply(1:4, function(j) get_sim_to_obs_ratio("p_hor[year,pop]", obs_y, sim_y, pop = j))
all_x = do.call(abind, append(x, list(along = 3)))

par(mfrow = c(1,2), oma = c(0,0,1,0))
my_plot(all_x[,"cv",], ylab = "cv")
my_plot(all_x[,"mean",], ylab = "mean")
mtext(side = 3, outer = TRUE, line = -0.5, "Proportion Hatchery in Return", cex = 1.25, font = 2)

##### TOTAL ADULT RETURNS #####
x = lapply(1:4, function(j) get_sim_to_obs_ratio("^Ra_tot[year,pop]", obs_y, sim_y, pop = j))
all_x = do.call(abind, append(x, list(along = 3)))

par(mfrow = c(1,2), oma = c(0,0,1,0))
my_plot(all_x[,"cv",], ylab = "cv")
my_plot(all_x[,"mean",], ylab = "mean")
mtext(side = 3, outer = TRUE, line = -0.5, "Total Adult Return", cex = 1.25, font = 2)

##### TOTAL PARR RECRUITMENT #####
x = lapply(1:4, function(j) get_sim_to_obs_ratio("^Pb[year,pop]", obs_y, sim_y, pop = j))
all_x = do.call(abind, append(x, list(along = 3)))

par(mfrow = c(1,2), oma = c(0,0,1,0))
my_plot(all_x[,"cv",], ylab = "cv")
my_plot(all_x[,"mean",], ylab = "mean")
mtext(side = 3, outer = TRUE, line = -0.5, "Total Parr Recruitment", cex = 1.25, font = 2)

##### AGE-3 MATURITY #####
x1 = lapply(1:4, function(j) get_sim_to_obs_ratio("^psi_O1_Rb[year,origin,pop]", obs_y, sim_y, origin = 1, pop = j))
x2 = lapply(1:4, function(j) get_sim_to_obs_ratio("^psi_O1_Rb[year,origin,pop]", obs_y, sim_y, origin = 2, pop = j))
all_x = do.call(abind, append(x1, append(x2, list(along = 3))))

par(mfrow = c(1,2), oma = c(0,0,1,0))
my_plot(all_x[,"cv",], ylab = "cv")
my_plot(all_x[,"mean",], ylab = "mean")
mtext(side = 3, outer = TRUE, line = -0.5, "Age-3 Maturity", cex = 1.25, font = 2)

##### AGE-4 MATURITY #####
x1 = lapply(1:4, function(j) get_sim_to_obs_ratio("^psi_O2_Rb[year,origin,pop]", obs_y, sim_y, origin = 1, pop = j))
x2 = lapply(1:4, function(j) get_sim_to_obs_ratio("^psi_O2_Rb[year,origin,pop]", obs_y, sim_y, origin = 2, pop = j))
all_x = do.call(abind, append(x1, append(x2, list(along = 3))))

par(mfrow = c(1,2), oma = c(0,0,1,0))
my_plot(all_x[,"cv",], ylab = "cv")
my_plot(all_x[,"mean",], ylab = "mean")
mtext(side = 3, outer = TRUE, line = -0.5, "Age-4 Maturity", cex = 1.25, font = 2)


# x = lapply(1:4, function(j) get_sim_to_obs_ratio("Mb_per_Sa_tot[year,pop]", obs_y, sim_y, pop = j))
# x = lapply(1:4, function(j) get_sim_to_obs_ratio("Ra_tot[year,pop]", obs_y, sim_y, pop = j))
# x = lapply(1:4, function(j) get_sim_to_obs_ratio("Sa_tot[year,pop]", obs_y, sim_y, pop = j))
all_x = do.call(abind, append(x1, append(x2, list(along = 3))))

samps = all_x[,"mean",]

range(all_x)




boxplot(all_x[,"cv",], outline = F, las = 2); abline(h = 1)
boxplot(all_x[,"mean",], outline = F, las = 2); abline(h = 1)

mean(var_ratio)

hist(var_ratio)

boxplot(var_ratio); abline(h = 1, lty = 2)


# set up simulations
sims <- 1000000
c1 <- -0.9
c2 <- 0.9
x1 <- rnorm(sims,0,1)

# create negative correlation with x1
x2 <- c1*x1 + sqrt(1-c1^2)*rnorm(sims,0,1)

# create searies uncorrelated with x1
x3 <- rnorm(sims,0,1)

# create positive correlation with x1
x4 <- c2*x1 + sqrt(1-c2^2)*rnorm(sims,0,1)

s1 <- plogis(x1)
s2 <- plogis(x2)
s3 <- plogis(x3)
s4 <- plogis(x4)

cor(s1,s2)
cor(s1,s3)
cor(s1,s4)

sd(s1 * s2) # negative correlation
mean(s1 * s2) # negative correlation

sd(s1 * s3) # not correlated
mean(s1 * s3) # not correlated

sd(s1 * s4) # positive correlation
mean(s1 * s4) # positive correlation
