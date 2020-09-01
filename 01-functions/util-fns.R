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
