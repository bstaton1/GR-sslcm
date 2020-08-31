# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT TO PROCESS RAW BIOLOGICAL DATA FILES FOR CONSTRUCTION OF A MAIN DATAFRAME #
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #

# clear the workspace
rm(list = ls(all = T))

# load all necessary packages
source("00-packages.R")

# load all necessary functions
invisible(sapply(list.files(path = "01-functions", pattern = "\\.R$", full.names = T), source))

# specify where data files are found
data_dir = "../GR-sslcm-data/bio"

##### ADULT ABUNDANCE DATA #####

# read the data
tmp = read.csv(file.path(data_dir, "01-adult-abundance.csv"), stringsAsFactors = F)

# keep only relevant columns
tmp = tmp[,c("population", "return_year", "n_returned", "n_above_weir")]

# update column names
# note: return_year renamed to "brood_year" to allow merging with other data sets
# and for consistent indexing in model.
# just note that for adults, brood_year is the year is the year adults SPAWNED, NOT the year they WERE SPAWNED
colnames(tmp) = c("population", "brood_year", "adults_at_weir", "adults_above_weir")

# rename the data frame, and remove "tmp" object
adult_abundance = tmp; rm(tmp)

##### ADULT COMPOSITION: CARCASS DATA #####

# read the data
tmp = read.csv(file.path(data_dir, "02a-adult-indiv-carcass.csv"), stringsAsFactors = F)

# new age_best variable
# some of the age_best are NA, but there are age records in the other columns
# this will correct that issue and use an age if present, while prioritizing age methods
tmp$age_best2 = NA
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

# discard records with unknown sex
tmp = tmp[tmp$sex != "Unk",]

# discard records with unknown origin
tmp = tmp[tmp$origin != "Unk",]

# keep only relevant columns
tmp = tmp[,c("population", "year", "sex", "origin", "age_best2")]
tmp$count = 1

# calculate the sum of the counts by population, year, sex, age, and origin
tmp = aggregate(count ~ population + year + sex + origin + age_best2, data = tmp, FUN = sum)

# reformat: to wide
tmp = dcast(tmp, population + year ~ origin + sex + age_best2, value.var = "count")

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
tmp = read.csv(file.path(data_dir, "02b-adult-indiv-weir.csv"), stringsAsFactors = F)

# new age_best variable
# some of the age_best are NA, but there are age records in the other columns
# this will correct that issue and use an age if present, while prioritizing age methods
tmp$age_best2 = NA
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

# discard records with unknown sex
tmp = tmp[tmp$sex != "Unk",]

# discard records with unknown origin
tmp = tmp[tmp$origin != "Unk",]

# rename trap_year to year
colnames(tmp)[colnames(tmp) == "trap_year"] = "year"

# keep only relevant columns
tmp = tmp[,c("population", "year", "sex", "origin", "age_best2", "count")]

# calculate the sum of the counts by population, year, sex, age, and origin
tmp = aggregate(count ~ population + year + sex + origin + age_best2, data = tmp, FUN = sum)

# reformat: to wide
tmp = dcast(tmp, population + year ~ origin + sex + age_best2, value.var = "count")

# update column names
# note: year renamed to "brood_year" to allow merging with other data sets
# and for consistent indexing in model.
# just note that for adults, brood_year is the year is the year adults SPAWNED, NOT the year they WERE SPAWNED
colnames(tmp)[2] = "brood_year"
colnames(tmp)[3:ncol(tmp)] = paste("weir", colnames(tmp)[3:ncol(tmp)], sep = "_")

# rename the data frame, and remove "tmp" objects
adult_weir_composition = tmp; rm(tmp)

##### JUVENILE ABUNDANCE #####

# read the data
tmp = read.csv(file.path(data_dir, "03-juv-abundance.csv"), stringsAsFactors = F)

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
tmp = read.csv(file.path(data_dir, "04-juv-survival.csv"), stringsAsFactors = F)

# remove any records that don't have a point estimate
tmp = tmp[!is.na(tmp$surv_est),]

# calculate the logit-scale standard error
tmp$logit_surv_se = with(tmp, get_logit_se(n_tagged, surv_est, surv_se, surv_ci_low, surv_ci_high))

# exclude some survival estimates: only keep the four main populations
tmp = tmp[tmp$population %in% c("CAT", "LOS", "MIN", "UGR"),]

# exclude some survival estimates: only keep summer, fall, and spring tagging to LGR
tmp = tmp[tmp$season %in% c("summer", "fall", "spring"),]

# add a brood_year column: two years prior to migration year
tmp$brood_year = tmp$mig_year - 2

# drop irrelevant columns
tmp = tmp[,c("population", "season", "brood_year", "surv_est", "logit_surv_se")]

# give season levels: for ordering purposes only
tmp$season = factor(tmp$season, levels = c("summer", "fall", "spring"))

# create two data frames: one for the point ests and one for SEs
tmp_est = tmp[,c("population", "season", "brood_year", "surv_est")]
tmp_se = tmp[,c("population", "season", "brood_year", "logit_surv_se")]

# reshape these
tmp_est = dcast(tmp_est, population + brood_year ~ season, value.var = "surv_est")
tmp_se = dcast(tmp_se, population + brood_year ~ season, value.var = "logit_surv_se")

# improve column names
colnames(tmp_est)[3:5] = paste0(colnames(tmp_est)[3:5], "_surv_est")
colnames(tmp_se)[3:5] = paste0(colnames(tmp_se)[3:5], "_surv_logit_se")

# combine back into one data set
tmp = merge(tmp_est, tmp_se, by = c("population", "brood_year"))

# rename the data frame and remove "tmp" objects
juvenile_survival = tmp; rm(list = c("tmp", "tmp_se", "tmp_est"))

