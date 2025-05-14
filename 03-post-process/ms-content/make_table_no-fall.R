
# load model info if not already in workspace
if (!exists("post")) source(file.path(this_dir, "load-model-info.R"))

# build a function to calculate coefficient of variation
cv = function(x) {sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE)}

##### STEP 1: EXTRACT POSTERIOR SAMPLES OF QUANTITIES NEEDED #####

# total parr recruitment (regardless of LH type)
Pb = post_subset(post, "^Pb[", matrix = TRUE)

# total smolt at LGR (separated by LH type)
Ma = post_subset(post, "^Ma[", matrix = TRUE)

# overwinter survival (separated by LH type)
phi_Pa_Mb = post_subset(post, "^phi_Pa_Mb[", matrix = TRUE)

# migration survival (separated by LH type, but assumed constant)
phi_Mb_Ma = post_subset(post, "^phi_Mb_Ma[", matrix = TRUE)

##### STEP 2: DEFINE A FUNCTION TO CALCULATE THE DESIRED STATS FOR A GIVEN POSTERIOR SAMPLE #####
f = function(i) {
  # put the values from this draw into the format used in the model
  Pb_i = array_format(Pb[i,])
  Ma_i = array_format(Ma[i,])
  phi_Pa_Mb_i = array_format(phi_Pa_Mb[i,])[,i_spring,]
  phi_Mb_Ma_i = array_format(phi_Mb_Ma[i,])[,i_spring,o_nor,]
  
  # calculate hypothetical NOR smolt reaching LGR if all parr were spring migrants
  Ma_tot_hyp_i = Pb_i * phi_Pa_Mb_i * phi_Mb_Ma_i
  
  # calculate actual NOR smolt reaching LGR
  Ma_tot_real_i = Ma_i[,i_fall,o_nor,] + Ma_i[,i_spring,o_nor,]
  
  # add columns storing the aggregate across populations
  Ma_tot_hyp_i = cbind(Ma_tot_hyp_i, rowSums(Ma_tot_hyp_i))
  Ma_tot_real_i = cbind(Ma_tot_real_i, rowSums(Ma_tot_real_i))
  
  # calculate the ratios between the actual and hypothetical
  mean_p_ratio_i = colMeans(Ma_tot_hyp_i/Ma_tot_real_i, na.rm = TRUE)
  cv_ratio_i = apply(Ma_tot_hyp_i, 2, cv)/apply(Ma_tot_real_i, 2, cv)
  
  # calculate interannual cv by population and scenario
  Ma_cv_hyp_i = apply(Ma_tot_hyp_i, 2, cv)
  Ma_cv_real_i = apply(Ma_tot_real_i, 2, cv)
  
  # calculate interannual mean by population and scenario
  Ma_mean_hyp_i = colMeans(Ma_tot_hyp_i, na.rm = TRUE)
  Ma_mean_real_i = colMeans(Ma_tot_real_i, na.rm = TRUE)
  
  # calculate percent changes
  Ma_mean_change_i = (Ma_mean_hyp_i - Ma_mean_real_i)/Ma_mean_real_i
  Ma_cv_change_i = (Ma_cv_hyp_i - Ma_cv_real_i)/Ma_cv_real_i
  
  # add names
  names(Ma_mean_hyp_i) = paste0("Ma_mean_hyp[", 1:5, "]")
  names(Ma_mean_real_i) = paste0("Ma_mean_real[", 1:5, "]")
  names(Ma_cv_hyp_i) = paste0("Ma_cv_hyp[", 1:5, "]")
  names(Ma_cv_real_i) = paste0("Ma_cv_real[", 1:5, "]")
  names(Ma_mean_change_i) = paste0("Ma_mean_change[", 1:5, "]")
  names(Ma_cv_change_i) = paste0("Ma_cv_change[", 1:5, "]")
  
  names(mean_p_ratio_i) = paste0("mean_p_ratio[", 1:5, "]")
  names(cv_ratio_i) = paste0("cv_ratio[", 1:5, "]")
  c(mean_p_ratio_i, cv_ratio_i, Ma_mean_hyp_i, Ma_mean_real_i, Ma_cv_hyp_i, Ma_cv_real_i, Ma_mean_change_i, Ma_cv_change_i)
}

# apply to each posterior sample
out = t(sapply(1:nrow(Pb), f))
out = post_convert(cbind(postpack:::id_mat(post), out))

# function to format the output
f = function(x, is_percent) {
  g = function(x, is_percent) {
    if (!is_percent) {
      round(x, 1)
    } else {
      paste0(round(x, 2) * 100, "%")
    }
  }
  paste0(g(x["mean",], is_percent), " (", g(x["sd",], is_percent), ")")
}

# build the table
tab = cbind(
  pop = c("CAT", "LOS", "MIN", "UGR", "Total"),
  Ma_with_fall = f(post_summ(out, "Ma_mean_real")/1000, FALSE),
  Ma_no_fall = f(post_summ(out, "Ma_mean_hyp")/1000, FALSE),
  Ma_pchange = f(post_summ(out, "Ma_mean_change"), TRUE),
  Ma_cv_with_fall = f(post_summ(out, "Ma_cv_real"), TRUE),
  Ma_cv_no_fall = f(post_summ(out, "Ma_cv_hyp"), TRUE)
)

# export the table to be read-in/kableExtra-ed in the msdown source
write.csv(tab, file.path(this_dir, "no-fall-tab.csv"), row.names = FALSE)
