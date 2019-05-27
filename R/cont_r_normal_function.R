#' Likelihood functions for data.frame
#' @title Normal-distribution-based random intercept likelihood function for xregControl and xreg
#' @param d_df data.frame provided by xreg.
#' @param sigma_est Within-unit variance parameter, may be specified as a function to vary between observations. Defaults to exp(LN_SIGMA)
#' @param omega_est Between-unit variance parameter, must be a fixed result for all observations, as only the first instance will be used. Defaults to exp(LN_OMEGA)
#' @description Likelihood function to use with xregControl
#' @export
cont_r_normal <- function (d_df, sigma_est = 0, omega_est = 0, internal_id = 0) {
if (sigma_est == 1)
  return(list(formula = formula(sigma_est ~ exp(LN_SIGMA), env = globalenv()),
              start_values = c(LN_SIGMA = 0)))
if (internal_id == 1)
  return(list(formula = formula(internal_id ~ id, env = globalenv())))
if (omega_est == 1)
  return(list(formula = formula(omega_est ~ exp(LN_OMEGA), env = globalenv()),
              start_values = c(LN_OMEGA = 0)))
if (min(d_df$omega_est) <= 0)
  return(NA)
if (min(d_df$sigma_est) <= 0)
  return(NA)

betwSD <- d_df$omega_est[1]

uid <- sort(unique(d_df$internal_id))
q <- NROW(uid)

d_df$p <- 0

order <- 20L

ws <- as.matrix(xreg_GHS[[order]]$w)
xs <- xreg_GHS[[order]]$x
mus <- sqrt(2) * betwSD * xs

tmpm <- matrix(NA, nrow = NROW(d_df), ncol = order)


sel <- (d_df$internal_ub == d_df$internal_lb)
sel2 <- (d_df$internal_ub != d_df$internal_lb)

tmpm[sel,] <- dnorm(d_df$internal_lb[sel], mean = outer(d_df$Xb[sel], mus, FUN = "+"), sd = d_df$sigma_est[sel], log = F)
tmpm[sel2,] <- pnorm((d_df$internal_lb[sel2]),outer(d_df$Xb[sel2], mus, FUN = "+"), d_df$sigma_est[sel2], lower.tail = F) -
  pnorm((d_df$internal_ub[sel2]), outer(d_df$Xb[sel2], mus, FUN = "+"), d_df$sigma_est[sel2], lower.tail = F)



d_df$p[match(uid, d_df$internal_id)] <- log(1/sqrt(pi)*as.matrix(aggregate(tmpm, by = list(d_df$internal_id), FUN = prod)[,-1]) %*% ws)
d_df$p[d_df$p == -Inf] <- log(.Machine$double.xmin)


return(d_df$p)
}
