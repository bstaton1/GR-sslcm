# :::::::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT TO HOUSE FUNCTIONS FOR PREPARING RAW DATA #
# :::::::::::::::::::::::::::::::::::::::::::::::: #

# OBTAIN LOGIT-SCALE STANDARD ERROR OF A PROPORTION
# some estimates have p_mean and CIs
# some estimates have p_mean and p_se
# either way, obtain CI's of proportion, then convert to a standard error on the logit scale
get_logit_se = function(N, p_mean, p_se, p_lwr, p_upr) {
  # turn se's into ci if that is what is available
  p_lwr = ifelse(is.na(p_se), p_lwr, p_mean + qnorm(0.025) * sqrt((p_mean * (1 - p_mean))/N) - (0.5/N))
  p_upr = ifelse(is.na(p_se), p_upr, p_mean + qnorm(0.975) * sqrt((p_mean * (1 - p_mean))/N) + (0.5/N))
  
  # cap the CI
  p_lwr = ifelse(p_lwr <= 0, 0.001, p_lwr)
  p_upr = ifelse(p_upr >= 1, 0.999, p_upr)
  
  # convert CI into logit-normal standard error: M. Liermann's approximation
  (logit(p_upr) - logit(p_lwr))/(2 * qnorm(0.975))
}
