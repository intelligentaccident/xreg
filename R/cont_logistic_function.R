#' Likelihood functions for data.frame
#' @title Logistic-distribution-based likelihood function for xregControl and xreg
#' @param d_df data.frame provided by xreg.
#' @param sigma_est Subsequent parameters define required formulas to genrete necessary variables, such as standard deviations. If these are not manually provided (which could be useful to model heteroscedasticity, etc.), xreg will get them from the functions.
#' @description Likelihood function to use with xregControl
cont_logistic <- function(d_df, sigma_est = 0) {
  if(sigma_est == 1) return(list(formula = formula(sigma_est ~ SIGMA, env = globalenv()),
                                 start_values = c(SIGMA = 1)))  #lb == ub
  sel <- d_df$internal_ub == d_df$internal_lb
  d_df$p[sel] <- log(exp(-(d_df$internal_lb[sel]- d_df$Xb[sel])/(d_df$sigma_est[sel])) / ((d_df$sigma_est[sel])*(1+exp(-(d_df$internal_lb[sel]- d_df$Xb[sel])/d_df$sigma_est[sel]))^2))

  # ub != lb
  sel <- d_df$internal_ub != d_df$internal_lb
  d_df$p[sel] <- log((1/(1+exp(-(d_df$internal_ub[sel]- d_df$Xb[sel])/d_df$sigma_est[sel])))- (1/(1+exp(-(d_df$internal_lb[sel]- d_df$Xb[sel])/d_df$sigma_est[sel]))))

  return(d_df$p)
}
