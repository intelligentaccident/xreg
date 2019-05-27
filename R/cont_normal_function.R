#' Likelihood functions for data.frame
#' @title Normal-distribution-based likelihood function for xregControl and xreg
#' @param d_df data.frame provided by xreg.
#' @param sigma_est Subsequent parameters define required formulas to generate necessary variables, such as standard deviations. If these are not manually provided (which could be useful to model heteroscedasticity, etc.), xreg will get them from the functions.
#' @description Likelihood function to use with xregControl
#' @export
cont_normal <- function(d_df, sigma_est = 0) {
  if (sigma_est == 1) 
    return(list(formula = formula(sigma_est ~ exp(LN_SIGMA), env = globalenv()), 
                start_values = c(LN_SIGMA = 0)))
  if (min(d_df$sigma_est) <= 0) 
    return(NA)
  sel <- d_df$internal_ub == d_df$internal_lb
  d_df$p[sel] <- dnorm((d_df$Xb[sel] - d_df$internal_lb[sel]), 
                       0, d_df$sigma_est[sel], log = TRUE)
  sel <- d_df$internal_ub != d_df$internal_lb
  d_df$p[sel] <- log(pnorm((d_df$Xb[sel] - d_df$internal_lb[sel])/d_df$sigma_est[sel], 
                           0, 1) - pnorm((d_df$Xb[sel] - d_df$internal_ub[sel])/d_df$sigma_est[sel], 
                                         0, 1))
  d_df$p[d_df$p == -Inf] <- log(.Machine$double.xmin)
  return(d_df$p)
}


