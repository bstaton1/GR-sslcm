
##### SESSION SET UP #####

# load packages
source("00-packages.R")

# load all necessary functions
invisible(sapply(list.files(path = "01-functions", pattern = "\\.R$", full.names = T), source))
source("C:/Users/bstaton/Desktop/source_rmd_chunk.R")

# define scenario names
scenarios = c("trib-harv_long", "fecund_long", "fecund-no-PSvar_long")
scenario_names = c("const-fecund", "vary-fecund", "vary-fecund+noPSvar")
in_dir = "02-model/model-output"

# read in the model output
# model_info_list = lapply(scenarios, function(m) {
#   readRDS(file.path(in_dir, paste0("output_", m, ".rds")))
# })
# post_list_all = lapply(model_info_list, function(m) post_thin(m$post, 0.25))
# jags_data_list_all = lapply(model_info_list, function(m) m$jags_data)

keep_mods = c(2:3)

# extract the elements
post_list = post_list_all[keep_mods]
jags_data_list = jags_data_list_all[keep_mods]

names(post_list) = names(jags_data_list) = scenario_names[keep_mods]

# select main colors
main_cols = c(
  "model" = "salmon",
  "model2" = "salmon",
  "data" = "royalblue"
)[keep_mods]

# set the transparent versions based on main color
tran_cols = sapply(main_cols, scales::alpha, alpha = 0.25)

all_years = 1990:2019
ts_yrs = which(all_years > 1990)
pops = c("CAT", "LOS", "MIN", "UGR")

##### PLOTTING FUNCTIONS #####

compare_tseries = function(post_list, param, ylim = NULL, yrs = NULL, label_text = NULL, label_loc = "topleft", y_scale = 1, legend_loc, ...) {
  # ests = lapply(post_list, function(m) post_summ(m, sub_index(param, year = ts_yrs, pop = j)))
  ests = lapply(post_list, function(m) post_summ(m, sub_index(param, ...)))
  
  if (is.null(ylim)) {
    ylim = make_lim(unlist(lapply(ests, function(x) x[c("2.5%", "97.5%"),]/y_scale)))
  }
  
  plot_tseries(est = ests[[1]], ylim = ylim, yrs = all_years[ts_yrs], label_text = label_text, label_loc = label_loc, blank = TRUE, y_scale = y_scale)
  lapply(1:length(post_list), function(m) add_tseries(ests[[m]], yrs = all_years[ts_yrs], main_col = main_cols[m], tran_col = tran_cols[m], y_scale = y_scale))
  
  if (!is.na(legend_loc)) {
    legend(legend_loc, legend = names(post_list), pch = 22,
           col = main_cols, pt.bg = tran_cols, pt.cex = 2, bty = "n", cex = 0.8, text.col = par("col.axis"))
  }
}

compare_BH = function(post_list, pop) {
  # function to get a sequence of egg abundances to predict at
  # calculates a 30 element vector counting up from zero to the maximum 97.5% quantile of eggs by a given population across models
  get_pred_eggs = function(j) {
    max_eggs = max(unlist(lapply(post_list, function(post) max(post_summ(post, sub_index("^E[.+,pop]", pop = j))["97.5%",]))))
    pred_eggs = seq(0, max_eggs, length = 30)
    return(pred_eggs)
  }
  
  # function to posterior summaries of predicted parr at each egg abundance
  get_pred_parr = function(post, j) {
    alpha_post = post_subset(post, sub_index("alpha[pop]", pop = j), TRUE)
    beta_post = post_subset(post, sub_index("beta[pop]", pop = j), TRUE)
    pred_eggs = get_pred_eggs(j)
    pred_parr = sapply(1:nrow(alpha_post), function(i) BH(pred_eggs, alpha_post[i], beta_post[i]))
    apply(pred_parr, 1, function(x) c(mean = mean(x), quantile(x, c(0.025, 0.975))))
  }
  
  # get the sequence of egg abundances to predict parr abundances at
  pred_eggs = get_pred_eggs(pop)
  
  # get posterior of predicted parr abundances
  pred_parr = lapply(post_list, get_pred_parr, j = pop)
  
  plot(1,1, type = "n", xlim = range(pred_eggs),
       ylim = make_lim(0, max(unlist(lapply(pred_parr, max)))),
       main = "", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  usr = par("usr")
  panel_label(pops[pop], "topleft")
  junk = sapply(1:length(post_list), function(i) {
    polygon(x = c(pred_eggs, rev(pred_eggs)), y = c(pred_parr[[i]]["2.5%",], rev(pred_parr[[i]]["97.5%",])), col = tran_cols[i], border = NA)
    lines(pred_parr[[i]]["mean",] ~ pred_eggs, col = main_cols[i], lwd = 2)
  })
  
  # draw axes
  at_x = axisTicks(usr[1:2], log = FALSE)
  at_y = axisTicks(usr[3:4], log = FALSE)
  axis(side = 1, at = at_x, labels = at_x/1e6)
  axis(side = 2, at = at_y, labels = at_y/1000)
  
  if (pop == 1) {
    legend("bottomright", legend = names(post_list), pch = 22,
           col = main_cols, pt.bg = tran_cols, pt.cex = 2, bty = "n", cex = 0.8, text.col = par("col.axis"))  }
}

compare_param = function(post_list, param, ylim = NULL, pt_est = "mean", label_loc = NULL, label_text = NULL, legend_loc = "topright", y_scale = 1, ...) {
  # get identifier for each element that will be summarized
  param_ID = dim_IDs(param, ...)
  
  # get posterior summaries from all models
  ests = lapply(post_list, function(post) post_summ(post, sub_index(param, ...), probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))
  
  # convert summaries into an array format rather than a list format (models are third dimension)
  ests = do.call(abind, append(ests, list(along = 3)))
  dimnames(ests)[[2]] = param_ID
  
  # calculate the y-axis limits if not supplied
  if (is.null(ylim)) {
    y_range = make_lim(ests[stringr::str_detect(dimnames(ests)[[1]], "%"),,])
    y_diff = diff(y_range)
    ylim = y_range + c(-1,1) * 0.05 * y_diff
  }
  
  # empty plot
  mp = barplot(t(ests[pt_est,,]), beside = TRUE, col = "white", border = NA, ylim = ylim, yaxt = "n", main = "")
  usr = par("usr")
  segments(usr[1], usr[3], usr[2], usr[3], xpd = TRUE)
  segments(usr[1], usr[3], usr[1], usr[4], xpd = TRUE)
  if (!is.null(label_loc)) panel_label(label_text, label_loc)
  
  # draw uncertainty intervals
  segments(mp, t(ests["2.5%",,]), mp, t(ests["97.5%",,]), col = main_cols)
  segments(mp, t(ests["25%",,]), mp, t(ests["75%",,]), lwd = 6, col = main_cols)
  
  # draw point estimates
  points(x = mp, y = t(ests[pt_est,,]), cex = 1.5, pch = 21, col = "white", bg = "white")
  points(x = mp, y = t(ests[pt_est,,]), cex = 1.5, pch = 21,
         col = main_cols, bg = tran_cols)
  
  # draw axes
  at_y = axisTicks(usr[3:4], log = FALSE)
  axis(side = 1, at = colSums(mp)/nrow(mp), labels = FALSE)
  axis(side = 2, at = at_y, labels = at_y/y_scale)
  
  if (!is.na(legend_loc)) {
    legend(legend_loc, legend = names(post_list), pch = 22,
           col = main_cols, pt.bg = tran_cols, pt.cex = 2, bty = "n", cex = 0.8, text.col = par("col.axis"))
  }
}

pdf("plots.pdf", w = 6.5, h = 4.5)

##### BH parameters #####
mypar(oma = c(0,0,0,0))
compare_param(post_list, param = "^alpha[pop]", pop = 1:4, label_text = "BH Prod.", label_loc = "topleft", legend_loc = "topright")
compare_param(post_list, param = "^beta[pop]", pop = 1:4, label_text = "BH Cap.", label_loc = "topleft", legend_loc = NA, y_scale = 1000)
compare_param(post_list, param = "^sig_Lphi_E_Pb[pop]", pop = 1:4, label_text = "BH SD", label_loc = "topleft", legend_loc = NA)
compare_param(post_list, param = "^kappa_phi_E_Pb[pop]", pop = 1:4, label_text = "BH AR(1) Coef.", label_loc = "topleft", legend_loc = NA)

##### BH RELATIONSHIPS #####
mypar()
sapply(1:4, function(j) compare_BH(post_list, pop = j))
axis_labels("Total Egg Production (Millions)", "Total Parr Recruitment (Thousands)")

##### TOTAL SPAWNER ABUNDANCE #####
mypar()
sapply(1:4, function(j) {
  compare_tseries(post_list, param = "^Sa_tot[year,pop]", pop = j, year = ts_yrs,
                  label_text = pops[j], legend_loc = ifelse(j == 2, "topright", NA),
                  y_scale = 1e3)
})
axis_labels("Return Year", "Total Spawner Abundance (Thousands)")

##### TOTAL EGG PRODUCTION #####

mypar()
sapply(1:4, function(j) {
  compare_tseries(post_list, param = "^E[year,pop]", pop = j, year = ts_yrs, 
                  label_text = pops[j], legend_loc = ifelse(j == 2, "topright", NA),
                  y_scale = 1e6)
})
axis_labels("Return Year", "Total Egg Production (Millions)")

##### TOTAL PARR RECRUITMENT #####

mypar()
sapply(1:4, function(j) {
  compare_tseries(post_list, param = "^Pb[year,pop]", pop = j, year = ts_yrs, 
                  label_text = pops[j], legend_loc = ifelse(j == 2, "topright", NA),
                  y_scale = 1e3)
})
axis_labels("Brood Year", "Total Parr Recruitment (Thousands)")

##### EGG TO PARR SURVIVAL #####

mypar()
sapply(1:4, function(j) {
  compare_tseries(post_list, param = "^phi_E_Pb[year,pop]", pop = j, year = ts_yrs, 
                  label_text = pops[j], legend_loc = ifelse(j == 2, "topright", NA),
                  y_scale = 1)
})
axis_labels("Brood Year", "Egg to Parr Survival")

##### EGGS PER SPAWNER #####

mypar()
sapply(1:4, function(j) {
  compare_tseries(post_list, param = "^E_per_Sa[year,pop]", pop = j, year = ts_yrs, 
                  label_text = pops[j], legend_loc = ifelse(j == 2, "topright", NA),
                  y_scale = 1)
})
axis_labels("Brood Year", "Eggs Per Spawner")

##### EGG TO PARR SURVIVAL PROCESS NOISE TERM #####

# proc_noise_file = "03-post-process/output-plots-children/11-noise-term-diags.Rmd"
# source_rmd_chunk(proc_noise_file, "misc-resid-fns")
# source_rmd_chunk(proc_noise_file, "process-noise-term-fns")
# 
# Lphi_E_Pb_qresid = lapply(1:length(keep_mods), function(m) {
#   cat("\nModel:", names(post_list)[m], "\n")
#   postpack:::id_mat(post_list[[m]]) |>
#     cbind(make_qresid_all(real_param = "^phi_E_Pb[year,pop]",
#                           mean_param = "phi_E_Pb_dot2[year,pop]",
#                           sig_param = "sig_Lphi_E_Pb[pop]",
#                           type = "pop", post_use = post_list[[m]], jags_data_use = jags_data_list[[m]])) |>
#     postpack::post_convert()
# }); names(Lphi_E_Pb_qresid) = scenario_names[keep_mods]


# mypar()
# sapply(1:4, function(j) {
#   compare_tseries(Lphi_E_Pb_qresid, param = "^Lphi_E_Pb_qresid[year,pop]", pop = j, year = ts_yrs, 
#                   label_text = pops[j], legend_loc = ifelse(j == 3, "topright", NA),
#                   y_scale = 1)
# })
# axis_labels("Brood Year", "Egg to Parr Survival QResidual")

##### FECUNDITY DATA (NO MODEL) PLOTS #####

# plot_f = function(j) {
#   plot(1, 1, type = "n", xlim = range(years), ylim = make_lim(ylim))
#   ylim = range(jags_data_list[[2]]$f[,2:3,], na.rm = TRUE)
#   
#   years = as.numeric(dimnames(jags_data_list[[2]]$f)[[1]])
#   
#   a4_lty = 1
#   a5_lty = 2
#   
#   col_old = "salmon"
#   col_new = "black"
#   
#   abline(h = jags_data_list[[1]]$f[2:3], col = col_old, lty = c(a4_lty, a5_lty))
#   lines(jags_data_list[[2]]$f[,2,j] ~ years, col = col_new, lty = a4_lty)
#   lines(jags_data_list[[2]]$f[,3,j] ~ years, col = col_new, lty = a5_lty)
#   panel_label(pops[j])
#   
#   if (j == 1) legend("topright", legend = c("Age-4", "Age-5"), lty = c(a4_lty, a5_lty), col = "black", cex = 0.8, bty = "n")
# }
# 
# mypar()
# sapply(1:4, plot_f)
# axis_labels("Return Year", "Fecundity")

# mypar(mfrow = c(1,1))
# x = apply(jags_data_list[[2]]$f[,2:3,], 3, colMeans, na.rm = TRUE)
# barplot(x, beside = TRUE, border ="white", col = c("grey80", "grey50"))
# abline(h = c(3971, 4846), col = "salmon", lty = c(1,2))
# axis_labels(ylab = "Fecundity")

##### END #####

dev.off(); file.show("plots.pdf")

