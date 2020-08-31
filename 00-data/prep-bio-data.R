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

