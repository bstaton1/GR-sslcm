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
