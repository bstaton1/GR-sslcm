# ::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT TO HOUSE FUNCTIONS FOR UTILITY TASKS #
# ::::::::::::::::::::::::::::::::::::::::::: #

##### CONVERT NATURAL SCALE CV TO LOGNORMAL SD #####
cv2sig = function(cv) {
  sqrt(log((cv^2) + 1))
}

##### CONVERT LOGNORMAL SD TO NATURAL SCALE CV #####
sig2cv = function(sig) {
  sqrt(exp(sig^2) - 1)
}

##### LOGIT TRANFORMATION #####
# p should be a vector on the interval (0,1)

logit = function(p) {
  log(p/(1 - p))
}

##### INVERSE LOGIT TRANSFORMATION #####
# lp should be a vector on the interval (-Inf,Inf)
# a value in logit scale

expit = function(lp) {
  exp(lp)/(1 + exp(lp))
}

##### OPPOSITE OF %in% #####
# it is annoying to have to write !(a %in% b)
# to find cases where elements of a are NOT found in b

`%!in%` = function(x, y) {
  !(x %in% y)
}

##### FLEXIBLY OBTAIN OBJECT DIMENSIONS #####
# dim() only works for 2+ dimensional objects
# this function can be used on 1+ dimensional objects

# x must be a vector, matrix, or array
# returns the number of elements found in each dimension
xdim = function(x) {
  # error handler for correct object type
  if (!is.vector(x) & !is.matrix(x) & !is.array(x)) {
    stop ("'x' must be a vector, matrix, or array")
  }
  
  # use the right dimensions function depending on structure
  if (is.vector(x)) {
    out = length(x)
  } else {
    out = dim(x)
  }
  
  # return the output
  return(out)
}

##### FIND ELEMENT INDICES THAT ARE NOT NA VALUES #####

# creates an object that specifies the element indices in each dimension
# of an objec that are not missing.
# this will be useful for fitting a JAGS model to square objects
# that have missing values. We can instruct JAGS to only touch the elements in this object
# and skip over the rest for likelihood nodes

# EXAMPLE: 3D ARRAY
# x = array(rnorm(100), dim = c(25, 2, 2))
# x[,1,1] = NA # insert some NAs
# x[,2,1] = NA # insert some more NAs
# find_no_na_indices(x)

find_no_na_indices = function(x) {
  
  # obtain object dimensions
  dims = xdim(x)
  
  # obtain number of dimensions
  n_dims = length(dims)
  
  # error handler for if object is too large
  if (n_dims > 6) {
    stop ("'x' has more than 6 dimensions")
  }
  
  # obtain total elements in object (regardless of NA or not)
  n_ind = prod(dims)
  
  # depending on the dimensions, create a data frame storing the index of each element
  # as well as a logical vector indicating whether each element is NA
  if (n_dims == 1) {
    ind = 1:dims[1]
    not_na = sapply(1:n_ind, function(i) !is.na(x[ind[i]]))
  }
  if (n_dims == 2) {
    ind = expand.grid(d1 = 1:dims[1], d2 = 1:dims[2])
    not_na = sapply(1:n_ind, function(i) !is.na(x[ind$d1[i],ind$d2[i]]))
  }
  if (n_dims == 3) {
    ind = expand.grid(d1 = 1:dims[1], d2 = 1:dims[2], d3 = 1:dims[3])
    not_na = sapply(1:n_ind, function(i) !is.na(x[ind$d1[i],ind$d2[i],ind$d3[i]]))
  }
  if (n_dims == 4) {
    ind = expand.grid(d1 = 1:dims[1], d2 = 1:dims[2], d3 = 1:dims[3], d4 = 1:dims[4])
    not_na = sapply(1:n_ind, function(i) !is.na(x[ind$d1[i],ind$d2[i],ind$d3[i],ind$d4[i]]))
  }
  if (n_dims == 5) {
    ind = expand.grid(d1 = 1:dims[1], d2 = 1:dims[2], d3 = 1:dims[3], d4 = 1:dims[4], d5 = 1:dims[5])
    not_na = sapply(1:n_ind, function(i) !is.na(x[ind$d1[i],ind$d2[i],ind$d3[i],ind$d4[i],ind$d5[i]]))
  }
  if (n_dims == 6) {
    ind = expand.grid(d1 = 1:dims[1], d2 = 1:dims[2], d3 = 1:dims[3], d4 = 1:dims[4], d5 = 1:dims[5], d6 = 1:dims[6])
    not_na = sapply(1:n_ind, function(i) !is.na(x[ind$d1[i],ind$d2[i],ind$d3[i],ind$d4[i],ind$d5[i],ind$d6]))
  }
  
  # extract only the element indices that are not NA
  # dimensions will differ for 1D vs. 2+D objects
  if (n_dims == 1) {
    out = ind[not_na]
  } else {
    out = as.matrix(ind[not_na,])
    rownames(out) = NULL
  }
  
  # return this object
  return(out)
}

##### PRODUCE PREDICTION OUT OF A BEVERTON-HOLT MODEL #####
# x: the x-value (e.g., spawners)
# alpha: max recruits/spawner
# capacity: max recruits

BH = function(x, alpha, beta) {
  x/((1/alpha) + (x/beta))
}

##### REPLACE PLACEHOLDER INDEX #####
# x: a string like "psi_O1[year,sex,origin,pop]
# other arguments: a number to replace
# e.g., sub_index("psi_O1[year,sex,origin,pop]", year = ".+", sex = 1, origin = 2, pop = ".")
# returns "psi_O1[.+,1,2,.]"
# this reduces how much hardcoding is required in generating output summaries

sub_index = function(x, year = NULL, LH_type = NULL, age = NULL, origin = NULL, pop = NULL) {
  # create a duplicate object
  newx = x
  
  # if one of the dimensions was specified, perform replacement
  if (!is.null(year)) newx = stringr::str_replace(newx, "year", as.character(year))
  if (!is.null(age)) newx = stringr::str_replace(newx, "age", as.character(age))
  if (!is.null(LH_type)) newx = stringr::str_replace(newx, "LH_type", as.character(LH_type))
  if (!is.null(origin)) newx = stringr::str_replace(newx, "origin", as.character(origin))
  if (!is.null(pop)) newx = stringr::str_replace(newx, "pop", as.character(pop))
  
  # return the version with placeholders replaced by index values
  newx
}

##### OBTAIN PARAMETER DIMENSION ID #####
# E.G., dim_IDs("alpha[pop]", pop = 1) gives "CAT"
# used for labeling in figures

dim_IDs = function(param, ...) {
  
  # function to convert sub_index() output into named dimensions
  dim_names = function(LH_type = NULL, age = NULL, origin = NULL, pop = NULL) {
    
    # empty list
    out = list()
    
    # combine the dimention names for each type supplied
    if (!is.null(LH_type)) out = append(out, list(LH_type = c("fall", "spring")[LH_type]))
    if (!is.null(age)) out = append(out, list(age = c("age3", "age4", "age5")[age]))
    if (!is.null(origin)) out = append(out, list(origin = c("NOR", "HOR")[origin]))
    if (!is.null(pop)) out = append(out, list(pop = c("CAT", "LOS", "MIN", "UGR")[pop]))
    
    return(out)
  }
  
  named_param = do.call(sub_index, append(list(x = param), dim_names(...)))
  named_param = stringr::str_remove(named_param, "year,")
  named_param = stringr::str_extract(named_param, "\\[.+\\]")
  named_param = stringr::str_remove(named_param, "age[:digit:],")
  named_param = stringr::str_remove(named_param, "\\[")
  named_param = stringr::str_remove(named_param, "\\]")
  return(named_param)
}


##### COMBINE A LIST OF DATA FRAMES #####
# list: a list with data frames to be rbind-ed stored as elements

unlist_dfs = function(list) {
  # empty object
  output = NULL
  
  # loop through list elements, combining the data frame in each with all previous
  for (i in 1:length(list)) output = rbind(output, list[[i]])
  
  # return the output
  return(output)
}

##### ADD A NEW INDEX #####

# changes the name of quantities (param) stored in mcmc.list object (post)
# adds a new dimension and index value
# e.g., add_index(post, "alpha[1]", 2)
# turns "alpha[1]" node name to "alpha[1,2]"

add_index = function(post, param, index_value) {
  # extract the names of all quantities that match param
  matches = match_params(post, param)
  
  # convert samples to matrix format while retaining the chain and iter ID
  post_m = as.matrix(post, chains = TRUE, iters = TRUE)
  
  # determine which elements that match param
  which_matches = which(colnames(post_m) %in% matches)
  
  # extract them in the order they are found in the object
  name_matches = colnames(post_m)[which_matches]
  
  # append a new index dimension and value on the back of the quantity name
  new_names = stringr::str_replace(name_matches, "\\]$", paste0(",", index_value, "]"))
  
  # replace the old names with new names
  colnames(post_m)[which_matches] = new_names
  
  # convert back to mcmc.list format
  post_convert(post_m)
}

##### DROP AN INDEX #####

# changes the name of quantities (param) stored in mcmc.list object (post)
# removes the last dimension and index value
# e.g., rm_index(post, "Pb[1,1]")
# turns "Pb[1,1]" node name to "Pb[1]", and returns only posterior samples that match "Pb[1,1]"
# allows using vcov_decomp() on a covariance matrix stored as a >2d array

rm_index = function(post, param) {
  post_sub = post_subset(post, param, matrix = TRUE, chains = TRUE, iters = TRUE)
  colnames(post_sub) = stringr::str_replace(colnames(post_sub), ",[:digit:]+\\]$", "]")
  post_convert(post_sub)
}

##### CREATE A JAGS MODEL FILE FROM AN R FUNCTION #####

# this function does the same thing as postpack::write_model or R2OpenBUGS::write.model
# except that it retains the comments contained in the source function

# function to write function to file
write_model_code = function(fun_file, out_file) {
  
  # extract the function body, including comments
  # code = attr(fun, "srcref")
  # code = as.character(code)
  code = readLines(fun_file)
  
  # replace the first line
  code[1] = "model {"
  
  # replace the last line
  code[length(code)] = "}  # END OF MODEL"
  
  # remove any instances of "%_%"
  code = stringr::str_remove(code, "%_%\\s?")
  
  # write the code to a file
  writeLines(code, out_file)
}

##### TOGGLE THE CALCULATION OF DATA CHECKS #####

# the default JAGS model code includes data simulation for the observed period (for posterior predictive checks)
# and calculates log posterior predictive density (for WAIC)
# these are calculated separately for all observed stochastic nodes
# given these calculations may add some computation time, use this function to turn them off when fitting the model

toggle_data_diagnostics = function(do_lppd = FALSE, do_ppd_check = FALSE, jags_file = "02-model/model.txt") {
  
  # read in the existing jags model code (after running write_model_code())
  # has all data diagnostics toggled on
  model_lines = readLines(jags_file)
  
  # toggle off WAIC calculations if requested
  if (!do_lppd) {
    which_matches = stringr::str_which(model_lines, "_lppd")
    spaces = stringr::str_remove(stringr::str_extract(model_lines[which_matches], "^[:space:]+[:alpha:]"), "[:alpha:]")
    model_lines[which_matches] = paste0(spaces, "# WAIC Calculations Toggled Off")
  }
  
  # toggle off posterior predictive check calculations
  if (!do_ppd_check) {
    which_matches = stringr::str_which(model_lines, "expected_|_new|_dev")
    spaces = stringr::str_remove(stringr::str_extract(model_lines[which_matches], "^[:space:]+[:alpha:]"), "[:alpha:]")
    model_lines[which_matches] = paste0(spaces, "# Posterior Predictive Check Calculations Toggled Off")
  }
  
  # write over the old jags model code
  if (!do_lppd | !do_ppd_check) {
    writeLines(model_lines, jags_file)
  }
  
}

toggle_HOR_Rb_init = function(jags_file = "02-model/model.txt") {
  # read in the existing jags model code
  # has Rb[y,k,o_hor,j] <- 0 for the init years
  model_lines = readLines(jags_file)
  
  # find the line numbers where Rb in init years is defined
  which_matches = stringr::str_which(model_lines, "^[:space:]*Rb\\[y,k,o_hor,j\\] <- 0")
  
  # replace the "<- 0" with a prior to turn on the estimation of the Rb params in these years
  model_lines[which_matches] = stringr::str_replace(model_lines[which_matches], "<- 0", "~ dunif(0, max_Rb_init[k])")
  
  # write over the old jags model code
  writeLines(model_lines, jags_file)
}

##### PRINT A NICE MESSAGE TO CONSOLE #####

my_cat = function(label, value, total_width = 60, indent = 2, first = FALSE) {
  label = paste(c(rep(" ", indent), "| ", label), collapse = "")
  value = paste(c(rep(" ", indent), value, " |"), collapse = "")
  label_width = nchar(label)
  value_width = nchar(value)
  blank_width = total_width - label_width - value_width
  if (blank_width < 0) {stop ("total_width too narrow")}
  if (first) cat(paste(c(rep(" ", indent), "|-", rep("-", total_width - 3 - indent), "|"), collapse = ""), "\n", sep = "")
  cat(label, paste(rep(" ", blank_width), collapse = ""), value, "\n", sep = "")
  cat(paste(c(rep(" ", indent), "|-", rep("-", total_width - 3 - indent), "|"), collapse = ""), "\n", sep = "")
  Sys.sleep(0.25)
}
