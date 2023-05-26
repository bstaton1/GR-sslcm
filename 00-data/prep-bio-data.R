# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT TO PROCESS RAW BIOLOGICAL DATA FILES FOR CONSTRUCTION OF A MAIN DATAFRAME #
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #

# clear the workspace
rm(list = ls(all = T))

# weir assignments to origin are known to be problematic at UGR
# ~10% of fish assigned NOR are actually HOR (only 50% HOR are ad-clipped)
# but only applies to 2006 and later; see gibsonpp/GR-sslcm-data#30
# can be turned off by setting p = 0
p_missassign_nor_UGR = 0.097
p_missassign_nor_UGR_fyr = 2006

# load all necessary packages
source("00-packages.R")

# load all necessary functions
invisible(sapply(list.files(path = "01-functions", pattern = "\\.R$", full.names = T), source))

# specify where data files are found
data_dir = "../GR-sslcm-data/bio"

##### ADULT ABUNDANCE DATA #####

# read the data
tmp = read.csv(file.path(data_dir, "adult-abundance.csv"), stringsAsFactors = F)

# keep only relevant columns
tmp = tmp[,c("population", "return_year", "n_returned", "n_above_weir", "abund_cv")]

# convert cv to lognormal sd
tmp$abund_cv = cv2sig(tmp$abund_cv)

# update column names
# note: return_year renamed to "brood_year" to allow merging with other data sets
# and for consistent indexing in model.
# just note that for adults, brood_year is the year is the year adults SPAWNED, NOT the year they WERE SPAWNED
colnames(tmp) = c("population", "brood_year", "adults_at_weir", "adults_above_weir", "adult_log_se")

# rename the data frame, and remove "tmp" object
adult_abundance = tmp; rm(tmp)

##### ADULT COMPOSITION: CARCASS DATA #####

# read the data
tmp = read.csv(file.path(data_dir, "adult-indiv-carcass.csv"), stringsAsFactors = F)

# new age_best variable
# some of the age_best are NA, but there are age records in the other columns
# this will correct that issue and use an age if present, while prioritizing age methods
# but use the "age_best" by default
tmp$age_best2 = tmp$age_best
tmp$age_best2 = ifelse(is.na(tmp$age_best2) & !is.na(tmp$age_cwt), tmp$age_cwt, tmp$age_best2)
tmp$age_best2 = ifelse(is.na(tmp$age_best2) & !is.na(tmp$age_scale), tmp$age_scale, tmp$age_best2)
tmp$age_best2 = ifelse(is.na(tmp$age_best2) & !is.na(tmp$age_key), tmp$age_key, tmp$age_best2)
tmp$age_best2 = ifelse(is.na(tmp$age_best2) & !is.na(tmp$age_length), tmp$age_length, tmp$age_best2)

## THE FOLLOWING DISCARD LINES
## throw out about 10% of the records
## the big one is missing origin.
## would be nice to find a way to use these records

# discard records with unknown age
tmp = tmp[!is.na(tmp$age_best2),]

# discard records with age assigned as 2
tmp = tmp[tmp$age_best2 > 2,]

# discard records with unknown origin
tmp = tmp[tmp$origin != "Unk",]

# keep only relevant columns
tmp = tmp[,c("population", "year", "origin", "age_best2")]
tmp$count = 1

# calculate the sum of the counts by population, year, age, and origin
tmp = aggregate(count ~ population + year + origin + age_best2, data = tmp, FUN = sum)

# reformat: to wide
tmp = dcast(tmp, population + year ~ origin + age_best2, value.var = "count")

# update column names
# note: year renamed to "brood_year" to allow merging with other data sets
# and for consistent indexing in model.
# just note that for adults, brood_year is the year is the year adults SPAWNED, NOT the year they WERE SPAWNED
colnames(tmp)[2] = "brood_year"
colnames(tmp)[3:ncol(tmp)] = paste("carc", colnames(tmp)[3:ncol(tmp)], sep = "_")

# rename the data frame, and remove "tmp" objects
adult_carc_composition = tmp; rm(tmp)

##### ADULT COMPOSITION: WEIR DATA #####

# read the data
tmp = read.csv(file.path(data_dir, "adult-indiv-weir.csv"), stringsAsFactors = F)

# new age_best variable
# some of the age_best are NA, but there are age records in the other columns
# this will correct that issue and use an age if present, while prioritizing age methods
# but use the "age_best" by default
tmp$age_best2 = tmp$age_best
tmp$age_best2 = ifelse(is.na(tmp$age_best2) & !is.na(tmp$age_cwt), tmp$age_cwt, tmp$age_best2)
tmp$age_best2 = ifelse(is.na(tmp$age_best2) & !is.na(tmp$age_pit), tmp$age_pit, tmp$age_best2)
tmp$age_best2 = ifelse(is.na(tmp$age_best2) & !is.na(tmp$age_scale), tmp$age_scale, tmp$age_best2)
tmp$age_best2 = ifelse(is.na(tmp$age_best2) & !is.na(tmp$age_length), tmp$age_length, tmp$age_best2)

## THE FOLLOWING DISCARD LINES
## throw out about 7.5% of the records
## the big one is missing age.
## would be nice to find a way to use these records

# discard records with unknown age
tmp = tmp[!is.na(tmp$age_best2),]

# discard records with age assigned as 2
tmp = tmp[tmp$age_best2 > 2,]

# discard records with unknown origin
tmp = tmp[tmp$origin != "Unk",]

# discard records for recaptured fish
tmp = tmp[-which(tmp$recapture),]

# discard records for LOS in 2001 - 2008
# known sampling issues in these years
# however, records needed in data set
# to be used in weir removal calculations.
tmp = tmp[-which(tmp$population == "LOS" & tmp$trapyear %in% 2001:2008),]

# rename trap_year to year
colnames(tmp)[colnames(tmp) == "trapyear"] = "year"

# keep only relevant columns
tmp = tmp[,c("population", "year", "origin", "age_best2", "count")]

# calculate the sum of the counts by population, year, age, and origin
tmp = aggregate(count ~ population + year + origin + age_best2, data = tmp, FUN = sum)

# reformat: to wide
tmp = dcast(tmp, population + year ~ origin + age_best2, value.var = "count")

# handle misassignments to origin in UGR
# see gibsonpp/GR-sslcm-data#30
tmp_UGR = tmp[tmp$population == "UGR" & tmp$year >= p_missassign_nor_UGR_fyr,]
tmp_UGR[is.na(tmp_UGR)] = 0 # NA means zero fish of that age/origin were counted
tmp_UGR_fixed = data.frame(
  do.call(rbind, lapply(1:nrow(tmp_UGR), correct_origins, x = tmp_UGR, k = 1)),
  do.call(rbind, lapply(1:nrow(tmp_UGR), correct_origins, x = tmp_UGR, k = 2)),
  do.call(rbind, lapply(1:nrow(tmp_UGR), correct_origins, x = tmp_UGR, k = 3))
)
tmp_UGR_fixed = cbind(tmp_UGR[,c("population", "year")], tmp_UGR_fixed)
tmp_UGR_fixed = tmp_UGR_fixed[,colnames(tmp_UGR)]
tmp[tmp$population == "UGR" & tmp$year >= p_missassign_nor_UGR_fyr,] = tmp_UGR_fixed

# update column names
# note: year renamed to "brood_year" to allow merging with other data sets
# and for consistent indexing in model.
# just note that for adults, brood_year is the year is the year adults SPAWNED, NOT the year they WERE SPAWNED
colnames(tmp)[2] = "brood_year"
colnames(tmp)[3:ncol(tmp)] = paste("weir", colnames(tmp)[3:ncol(tmp)], sep = "_")

# rename the data frame, and remove "tmp" objects
adult_weir_composition = tmp; rm(tmp)

##### ADULT COMPOSITION: WEIR REMOVAL DATA #####

# read the data
tmp = read.csv(file.path(data_dir, "adult-indiv-weir.csv"), stringsAsFactors = F)

# new age_best variable
# some of the age_best are NA, but there are age records in the other columns
# this will correct that issue and use an age if present, while prioritizing age methods
# but use the "age_best" by default
tmp$age_best2 = tmp$age_best
tmp$age_best2 = ifelse(is.na(tmp$age_best2) & !is.na(tmp$age_cwt), tmp$age_cwt, tmp$age_best2)
tmp$age_best2 = ifelse(is.na(tmp$age_best2) & !is.na(tmp$age_pit), tmp$age_pit, tmp$age_best2)
tmp$age_best2 = ifelse(is.na(tmp$age_best2) & !is.na(tmp$age_scale), tmp$age_scale, tmp$age_best2)
tmp$age_best2 = ifelse(is.na(tmp$age_best2) & !is.na(tmp$age_length), tmp$age_length, tmp$age_best2)

## THE FOLLOWING DISCARD LINES
## throw out about 7.5% of the records
## the big one is missing age.
## would be nice to find a way to use these records

# discard records with unknown age
tmp = tmp[!is.na(tmp$age_best2),]

# discard records with age assigned as 2
tmp = tmp[tmp$age_best2 > 2,]

# discard records with unknown origin
tmp = tmp[tmp$origin != "Unk",]

# rename trap_year to year
colnames(tmp)[colnames(tmp) == "trapyear"] = "year"

# keep only fish that were removed
tmp = tmp[tmp$disposition == "removed",]

# keep only relevant columns
tmp = tmp[,c("population", "year", "origin", "age_best2", "count")]

# calculate the sum of the counts by population, year, age, and origin
tmp = aggregate(count ~ population + year + origin + age_best2, data = tmp, FUN = sum)

# reformat: to wide
tmp = dcast(tmp, population + year ~ origin + age_best2, value.var = "count")

# handle misassignments to origin in UGR
# see gibsonpp/GR-sslcm-data#30
tmp_UGR = tmp[tmp$population == "UGR" & tmp$year >= p_missassign_nor_UGR_fyr,]
tmp_UGR[is.na(tmp_UGR)] = 0 # NA means zero fish of that age/origin were counted
tmp_UGR_fixed = data.frame(
  do.call(rbind, lapply(1:nrow(tmp_UGR), correct_origins, x = tmp_UGR, k = 1)),
  do.call(rbind, lapply(1:nrow(tmp_UGR), correct_origins, x = tmp_UGR, k = 2)),
  do.call(rbind, lapply(1:nrow(tmp_UGR), correct_origins, x = tmp_UGR, k = 3))
)
tmp_UGR_fixed = cbind(tmp_UGR[,c("population", "year")], tmp_UGR_fixed)
tmp_UGR_fixed = tmp_UGR_fixed[,colnames(tmp_UGR)]
tmp[tmp$population == "UGR" & tmp$year >= p_missassign_nor_UGR_fyr,] = tmp_UGR_fixed

# update column names
# note: year renamed to "brood_year" to allow merging with other data sets
# and for consistent indexing in model.
# just note that for adults, brood_year is the year is the year adults SPAWNED, NOT the year they WERE SPAWNED
colnames(tmp)[2] = "brood_year"
colnames(tmp)[3:ncol(tmp)] = paste("rm", colnames(tmp)[3:ncol(tmp)], sep = "_")

# rename the data frame, and remove "tmp" objects
adult_rm_composition = tmp; rm(tmp)

##### ADULT PRE-SPAWN SURVIVAL: CARCASS DATA #####

# years with fewer carcasses sampled than this will be assigned
# NA values for numbers of carcasses sampled and found with spawn status
# intended to remove the effects of weakly informative data
status_count_threshold = 10

# read the data
tmp = read.csv(file.path(data_dir, "adult-indiv-carcass.csv"), stringsAsFactors = F)

# drop out non-known females
tmp = tmp[tmp$sex == "F",]

# drop out non-known spawn status
tmp = tmp[tmp$prespawn != "Unk" & !is.na(tmp$prespawn),]

# add a count variable: for use in summing records
tmp$count = 1

# aggregate data: total carcasses sampled by population and year
sampled = aggregate(count ~ population + year, tmp, sum)
colnames(sampled)[3] = "carcs_samp_for_status"

# aggregate data: total carcasses sampled that were spawned out by population and year
success = aggregate(count ~ population + year, subset(tmp, prespawn == "Spawned"), sum)
colnames(success)[3] = "carcs_status_spawned"

# merge these two data sets
tmp = merge(sampled, success, all = T)

# if carcasses were sampled in a year, but none were found successfully spawned,
# there will be an NA in that cell. These should be zeros
tmp$carcs_status_spawned[is.na(tmp$carcs_status_spawned)] = 0

# update column names
# note: year renamed to "brood_year" to allow merging with other data sets
# and for consistent indexing in model.
# just note that for adults, brood_year is the year is the year adults SPAWNED, NOT the year they WERE SPAWNED
colnames(tmp)[2] = "brood_year"

# discard years with too few fish sampled
tmp$carcs_status_spawned[tmp$carcs_samp_for_status <= status_count_threshold] = NA
tmp$carcs_samp_for_status[tmp$carcs_samp_for_status <= status_count_threshold] = NA

# rename the data frame, and remove "tmp" objects
adult_prespawn = tmp; rm(tmp)

##### ADULT SURVIVAL PAST SEA LIONS #####

# read the data
tmp = read.csv(file.path(data_dir, "07-sea-lion-survival.csv"))

# reformat to long format
tmp = melt(tmp, id.vars = "year", value.name = "surv_est_sea_lions", variable.name = "population")

# update column names
# note: return_year renamed to "brood_year" to allow merging with other data sets
# and for consistent indexing in model.
# just note that for adults, brood_year is the year is the year adults SPAWNED, NOT the year they WERE SPAWNED
colnames(tmp)[1] = "brood_year"

# rename object and remove "tmp" object
sea_lion_survival = tmp; rm(tmp)

##### ADULT HARVEST RATES BELOW BON #####

# read the data
tmp = read.csv(file.path(data_dir, "harvest-rates-below-BON.csv"))

# add a population variable
tmp$population = "ALL"

# update column names
colnames(tmp) = c("brood_year", "HR_NOR", "HR_HOR", "population")

# rename object and remove "tmp" object
harvest_rates = tmp; rm(tmp) 

##### ADULT SURVIVAL BON -> LGR #####

# read the data
tmp = read.csv(file.path(data_dir, "BON-LGR-adult-PIT-detections.csv"))

# reformat BON detection counts
tmp_BON = dcast(tmp[,c("year", "origin", "BON_adult_detections")], year ~ origin, value.var = "BON_adult_detections")
colnames(tmp_BON) = c("year", "HOR_BON_adults", "NOR_BON_adults")
tmp_LGR = dcast(tmp[,c("year", "origin", "LGR_adult_detections")], year ~ origin, value.var = "LGR_adult_detections")
colnames(tmp_LGR) = c("year", "HOR_LGR_adults", "NOR_LGR_adults")

# merge into one data set
tmp = merge(tmp_BON, tmp_LGR, by = "year")

# add a population column
tmp = cbind(population = "ALL", tmp)

# update column names
# note: year renamed to "brood_year" to allow merging with other data sets
# and for consistent indexing in model.
# just note that for adults, brood_year is the year is the year adults SPAWNED, NOT the year they WERE SPAWNED
colnames(tmp)[2] = "brood_year"

# rename the object and remove "tmp" object
dam_adult_counts = tmp; rm(tmp)

##### ADULT HARVEST IN TRIBUTARIES #####

# read in the raw data file
tmp = read.csv(file.path(data_dir, "tributary-harvest.csv"))

# sum over sexes
tmp = aggregate(count ~ population + year + origin + age, data = tmp, FUN = sum)

# relabel origin type
tmp$origin = ifelse(tmp$origin == "Hat", "HOR", ifelse(tmp$origin == "Nat", "NOR", NA))

# create "age/origin" combo variable
tmp$origin_age = paste0(tmp$origin, "_", tmp$age)

# reformat
tmp = dcast(tmp[,c("population", "year", "origin_age", "count")], population + year ~ origin_age, value.var = "count")

# update column names
# note: year renamed to "brood_year" to allow merging with other data sets
# and for consistent indexing in model.
# just note that for adults, brood_year is the year is the year adults SPAWNED, NOT the year they WERE SPAWNED
colnames(tmp)[2] = "brood_year"
colnames(tmp)[3:ncol(tmp)] = paste0("harv_", colnames(tmp)[3:ncol(tmp)])

# rename the data frame, and remove "tmp" objects
harv_composition = tmp; rm(tmp)

##### ADULT FECUNDITY #####

# read the data
tmp = read.csv(file.path(data_dir, "fecundity.csv"))

# reformat the data
tmp = dcast(tmp, pop + return_year ~ age, value.var = "fecundity")

# update column names
# note: year renamed to "brood_year" to allow merging with other data sets
# and for consistent indexing in model.
# just note that for adults, brood_year is the year is the year adults SPAWNED, NOT the year they WERE SPAWNED
colnames(tmp) = c("population", "brood_year", "fecund_4", "fecund_5")

# add in age-3 fecundity: assumed to be zero
tmp = cbind(tmp[,c("population", "brood_year")], "fecund_3" = 0, tmp[,c("fecund_4", "fecund_5")])

# rename the data frame, and remove "tmp" objects
fecundity = tmp; rm(tmp)

##### JUVENILE ABUNDANCE #####

# read the data
tmp = read.csv(file.path(data_dir, "juv-abundance.csv"), stringsAsFactors = F)

# remove any records that don't have a point estimate
tmp = tmp[!is.na(tmp$abund_est),]

# convert the mean estimate and standard error into a "CV"
tmp$abund_cv = tmp$abund_se/tmp$abund_est

# convert the "CV" into a lognormal standard deviation
tmp$abund_sigma = cv2sig(tmp$abund_cv)

# drop irrelevant columns
tmp = tmp[,c("population", "season", "brood_year", "abund_est", "abund_sigma")]

# give season levels: for ordering purposes only
tmp$season = factor(tmp$season, levels = c("fall", "spring"))

# create two data frames: one for the point ests and one for sigmas
tmp_est = tmp[,c("population", "season", "brood_year", "abund_est")]
tmp_se = tmp[,c("population", "season", "brood_year", "abund_sigma")]

# reshape these
tmp_est = dcast(tmp_est, population + brood_year ~ season, value.var = "abund_est")
tmp_se = dcast(tmp_se, population + brood_year ~ season, value.var = "abund_sigma")

# improve column names
colnames(tmp_est)[3:4] = paste0(colnames(tmp_est)[3:4], "_passage_est")
colnames(tmp_se)[3:4] = paste0(colnames(tmp_se)[3:4], "_passage_log_se")

# combine back into one data set
tmp = merge(tmp_est, tmp_se, by = c("population", "brood_year"))

# rename the data frame, and remove "tmp" objects
juvenile_abundance = tmp; rm(list = c("tmp", "tmp_se", "tmp_est"))

##### JUVENILE SURVIVAL DATA #####

# read the data
tmp = read.csv(file.path(data_dir, "juv-survival.csv"), stringsAsFactors = F)

# remove any records that don't have a point estimate
tmp = tmp[!is.na(tmp$surv_est),]

# calculate the logit-scale standard error
tmp$logit_surv_se = with(tmp, get_logit_se(surv_est, surv_se, surv_ci_low, surv_ci_high, alpha = 0.05))

# exclude some survival estimates: only keep the four main populations
tmp = tmp[tmp$population %in% c("CAT", "LOS", "MIN", "UGR"),]

# exclude some survival estimates: only keep summer, fall, and spring tagging to LGR
tmp = tmp[tmp$season %in% c("summer", "fall", "spring", "winter"),]

# add a brood_year column: two years prior to migration year
tmp$brood_year = tmp$mig_year - 2

# drop irrelevant columns
tmp = tmp[,c("population", "season", "brood_year", "surv_est", "logit_surv_se")]

# give season levels: for ordering purposes only
tmp$season = factor(tmp$season, levels = c("summer", "fall", "winter", "spring"))

# create two data frames: one for the point ests and one for SEs
tmp_est = tmp[,c("population", "season", "brood_year", "surv_est")]
tmp_se = tmp[,c("population", "season", "brood_year", "logit_surv_se")]

# reshape these
tmp_est = dcast(tmp_est, population + brood_year ~ season, value.var = "surv_est")
tmp_se = dcast(tmp_se, population + brood_year ~ season, value.var = "logit_surv_se")

# improve column names
colnames(tmp_est)[3:6] = paste0(colnames(tmp_est)[3:6], "_surv_est")
colnames(tmp_se)[3:6] = paste0(colnames(tmp_se)[3:6], "_surv_logit_se")

# combine back into one data set
tmp = merge(tmp_est, tmp_se, by = c("population", "brood_year"))

# rename the data frame and remove "tmp" objects
juvenile_survival = tmp; rm(list = c("tmp", "tmp_se", "tmp_est"))

##### JUVENILE LENGTH DATA #####

# read the data
tmp = read.csv(file.path(data_dir, "juv-mean-length.csv"))

# retain only rows corresponding to all fish measured (tagged and untagged)
tmp = tmp[tmp$fish_subset == "all",]

# set any mean lengths with sample size less than 30 to NA
tmp[tmp$n_length < 30,c("len_mean", "len_sd")] = NA

# keep only summer and spring length measurements
tmp = tmp[tmp$season %in% c("summer", "spring"),]
tmp$season = factor(tmp$season, levels = c("summer", "spring"))

# calculate standard error in mean length
tmp$len_se = tmp$len_sd/sqrt(tmp$n_length)

# ensure the data are ordered by population and year
tmp = tmp[order(tmp$population, tmp$season, tmp$brood_year),]

# calculate adjusted mean length for summer
# corrects for capture date variability among years
length_mean_adj = lapply(unique(tmp$population), function(pop) {
  tmp_sub = subset(tmp, population == pop & season == "summer")
  standardize_mean_length(tmp_sub$len_mean, tmp_sub$jday_med, resid_type = "mult")
})
tmp$len_mean_adj = NA
tmp$len_mean_adj[tmp$season == "summer"] = unlist(length_mean_adj)

# calculate standard error in adjusted mean length
# (assume the CV is the same)
tmp$len_se_adj = tmp$len_se/tmp$len_mean * tmp$len_mean_adj

# set the adjusted mean lengths & SEs for spring equal to the observed mean lengths & SEs
# allows using the same column name for extracting below; not adjusting the spring length data b/c passive capture
tmp[tmp$season == "spring",c("len_mean_adj", "len_se_adj")] = tmp[tmp$season == "spring",c("len_mean", "len_se")]

# retain only needed columns
tmp = tmp[,c("population", "season", "brood_year", "n_length", "len_mean", "len_se", "len_mean_adj", "len_se_adj")]

# extract and format mean length
length_mean = dcast(tmp, brood_year + population ~ season, value.var = "len_mean_adj")
colnames(length_mean)[3:ncol(length_mean)] = paste0("length_mean_", colnames(length_mean)[3:ncol(length_mean)])

# extract and format se mean length
length_se = dcast(tmp, brood_year + population ~ season, value.var = "len_se_adj")
colnames(length_se)[3:ncol(length_se)] = paste0("length_se_", colnames(length_se)[3:ncol(length_se)])

# combine into one data set
tmp = merge(length_mean, length_se, by = c("population", "brood_year"), all = TRUE)

# rename the data frame and remove "tmp" objects
juvenile_length = tmp; rm("tmp", "length_mean", "length_se", "length_mean_adj")

##### HATCHERY RELEASES OF SMOLTS AND SURVIVAL TO LGD #####

# read the data
tmp = read.csv(file.path(data_dir, "hatchery-juv-releases.csv"), stringsAsFactors = F)

# discard any records prior to 1995
# these releases were very early and did not follow the same
# procedures as in the rest of the "official" hatchery program years
tmp = tmp[tmp$brood_year > 1995,]

# note this record: some fish released in fall as parr
tmp[tmp$population == "UGR" & tmp$brood_year == 2000,"comments"]

# approach for now: ignore these fish, and just count smolt released as usual
tmp[tmp$population == "UGR" & tmp$brood_year == 2000,"n_smolt_released"] = 151443

# convert survival se into logit survival se
tmp$logit_surv_se = with(tmp, get_logit_se(surv_est, surv_se, NA, NA, alpha = 0.05))

# keep only relevant columns
tmp = tmp[,c("population", "brood_year", "n_smolt_released", "surv_est", "logit_surv_se")]

# rename columns
colnames(tmp)[colnames(tmp) == "n_smolt_released"] = "hatchery_smolt"
colnames(tmp)[colnames(tmp) == "surv_est"] = "hatchery_spring_surv_est"
colnames(tmp)[colnames(tmp) == "logit_surv_se"] = "hatchery_spring_surv_logit_se"

# rename the data frame, and remove "tmp" object
hatchery_release_survival = tmp; rm(tmp)

##### JUVENILE SURVIVAL: HYDROSYSTEM #####

# read the data: found in scratch folder for now
# until we decide on the best source for these data
tmp = read.csv(file.path(data_dir, "juv-survival-hydro.csv"), stringsAsFactors = F)

# convert migration year to brood year
tmp$brood_year = tmp$mig_year - 2

# calculate logit-scale SE of survival estimate: hatchery origin
tmp$logit_hat_surv_se = sapply(1:nrow(tmp), function(i) {
  get_logit_se(p_mean = tmp$hat_est[i],
               p_se = NA, 
               p_lwr = tmp$hat_lwr90[i], 
               p_upr = tmp$hat_upr90[i], 
               alpha = 0.1)
})

# calculate logit-scale SE of survival estimate: natural origin
tmp$logit_nat_surv_se = sapply(1:nrow(tmp), function(i) {
  get_logit_se(p_mean = tmp$nat_est[i],
               p_se = NA, 
               p_lwr = tmp$nat_lwr90[i], 
               p_upr = tmp$nat_upr90[i], 
               alpha = 0.1)
})

# retain only necessary columns
tmp = tmp[,c("brood_year", "nat_est", "logit_nat_surv_se", "hat_est", "logit_hat_surv_se")]

# add a "pop" column
tmp = cbind(data.frame(pop = "ALL"), tmp)

# change column names
colnames(tmp) = c("population", "brood_year", "nat_hydro_est", "nat_hydro_logit_se", "hat_hydro_est", "hat_hydro_logit_se")

# rename object and delete tmp
hydro_surv = tmp; rm(tmp)

##### COMBINE THESE DATA SOURCES INTO ONE DATA FRAME #####

# merge together the various data sets
bio_dat = merge(juvenile_abundance, juvenile_survival, by = c("population", "brood_year"), all = T)
bio_dat = merge(bio_dat, juvenile_length, by = c("population", "brood_year"), all = T)
bio_dat = merge(bio_dat, adult_abundance, by = c("population", "brood_year"), all = T)
bio_dat = merge(bio_dat, adult_carc_composition, by = c("population", "brood_year"), all = T)
bio_dat = merge(bio_dat, adult_weir_composition, by = c("population", "brood_year"), all = T)
bio_dat = merge(bio_dat, adult_rm_composition, by = c("population", "brood_year"), all = T)
bio_dat = merge(bio_dat, harv_composition, by = c("population", "brood_year"), all = T)
bio_dat = merge(bio_dat, fecundity, by = c("population", "brood_year"), all = T)
bio_dat = merge(bio_dat, adult_prespawn, by = c("population", "brood_year"), all = T)
bio_dat = merge(bio_dat, hatchery_release_survival, by = c("population", "brood_year"), all = T)
bio_dat = merge(bio_dat, hydro_surv, by = c("population", "brood_year"), all = T)
bio_dat = merge(bio_dat, dam_adult_counts, by = c("population", "brood_year"), all = T)
bio_dat = merge(bio_dat, harvest_rates, by = c("population", "brood_year"), all = T)

# create an empty data frame for merging
# this ensures all populations have rows for every year
empty_df = with(bio_dat, expand.grid(population = unique(population), brood_year = seq(min(brood_year), max(brood_year))))
bio_dat = merge(bio_dat, empty_df, by = c("population", "brood_year"), all = T)

# make hatchery releases be zero if NA
bio_dat$hatchery_smolt[is.na(bio_dat$hatchery_smolt) & bio_dat$population != "ALL"] = 0

colnames(bio_dat) = stringr::str_replace(colnames(bio_dat), "Nat", "NOR")
colnames(bio_dat) = stringr::str_replace(colnames(bio_dat), "Hat", "HOR")

# remove unnecessary objects from workspace, retain only bio_dat
rm(list = setdiff(ls(), "bio_dat"))
