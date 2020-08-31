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

