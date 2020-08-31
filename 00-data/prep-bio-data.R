# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT TO PROCESS RAW BIOLOGICAL DATA FILES FOR CONSTRUCTION OF A MAIN DATAFRAME #
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #

# clear the workspace
rm(list = ls(all = T))

# load all necessary packages
source("00-packages.R")

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

