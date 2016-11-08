#' Likelihood functions for data.frame
#' @title Likelihood functions for xregControl
#' @param d_df data.frame provided by xreg.
#' @param sigma_est Subsequent parameters define required formulas to genrete necessary variables, such as standard deviations. If these are not manually provided (which could be useful to model heteroscedasticity, etc.), xreg will get them from the functions.
#' @return \code{cont_normal} normal-distribution-based likelihood for continuous variables
#' @return \code{cont_logistic} logistic-distribution-based likelihood for continuous variables
#' @return \code{dich_normal} normal-distribution-based likelihood for dichotomous variables (probit)
#' @return \code{dich_logistic} logistic-distribution-based likelihood for dichotomous variables (logit)
#' @return \code{berkson_w} chi_square distance
#' @return \code{berkson_w_aggr} weighted aggregation for chi-square
#' @return \code{mse_p} means square error
#' @aliases cont_logistic, dich_normal, dich_logistic, berkson_w, berkson_w_aggr, mse_p
#' @description Likelihood functions to use with xregControl
cont_normal <- function(d_df, sigma_est = 0) {
  if(sigma_est == 1) return(list(formula = formula(sigma_est ~ SIGMA, env = globalenv()),
                                 start_values = c(SIGMA = 1)))

  if(min(d_df$sigma_est) <= 0) return(NA)
  #lb == ub
  sel <- d_df$internal_ub == d_df$internal_lb
  d_df$p[sel] <- dnorm((d_df$Xb[sel] - d_df$internal_lb[sel]), 0, d_df$sigma_est[sel], log = TRUE)

  # ub != lb
  sel <- d_df$internal_ub != d_df$internal_lb
  d_df$p[sel] <- log(pnorm((d_df$Xb[sel]-d_df$internal_lb[sel]) / d_df$sigma_est[sel], 0, 1)-pnorm((d_df$Xb[sel]-d_df$internal_ub[sel])/d_df$sigma_est[sel], 0, 1))

  return(d_df$p)
}

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

dich_normal <- function(d_df, theta_est = 0) {
  if(theta_est == 1) return(list(formula = formula(theta_est ~ THETA, env = globalenv()),
                                 start_values = c(THETA = 1)))
  return(d_df$internal_lb*pnorm(d_df$Xb/d_df$theta_est, 0, 1, log.p = TRUE) + (1-d_df$internal_lb)*pnorm(-d_df$Xb/d_df$theta_est, 0, 1, log.p = TRUE))
}

dich_logistic <- function(d_df, theta_est = 0) {
  if(theta_est == 1) return(list(formula = formula(theta_est ~ THETA, env = globalenv()),
                                 start_values = c(THETA = 1)))

  return(d_df$internal_lb*(-log(1+exp(-d_df$Xb*d_df$theta_est)))   + (1-d_df$internal_lb)*(-d_df$Xb*d_df$theta_est   - log(1+exp(-d_df$Xb*d_df$theta_est))))
}

berkson_w <- function(d_df) {
  with(d_df, {
    sel <- (choice>=1 | choice<=0)
    w <- 0
    w[!sel] <- (1/(choice[!sel]*(1-choice[!sel])))
    w[sel] <- (4*internal_count[sel]^2)/(2*internal_count[sel]-1)
    return(w)
  })
}


berkson_w_aggr <- function(d_df) {
  return(sum(d_df$internal_count * (d_df$internal_lb - d_df$Xb)^2*d_df$p))
}

mse_p <- function(d_df) {
  return(-(d_df$internal_lb - d_df$Xb)^2)
}
