#' Likelihood functions for data.frame
#' @title Normal-distribution-based likelihood function for dichotomous data, used by xregControl and xreg
#' @param d_df data.frame provided by xreg.
#' @param theta_est Subsequent parameters define required formulas to generate necessary variables, such as standard deviations. If these are not manually provided (which could be useful to model heteroscedasticity, etc.), xreg will get them from the functions.
#' @description Likelihood function to use with xregControl
#' @export
dich_normal <- function(d_df, theta_est = 0) {
  if(theta_est == 1) return(list(formula = formula(theta_est ~ THETA, env = globalenv()),
                                 start_values = c(THETA = 1)))
  return(d_df$internal_lb*pnorm(d_df$Xb/d_df$theta_est, 0, 1, log.p = TRUE) + (1-d_df$internal_lb)*pnorm(-d_df$Xb/d_df$theta_est, 0, 1, log.p = TRUE))
}
