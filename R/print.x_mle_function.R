#' Print x_mle-object
#'
#' @param x_mle x_mle object
#' @description Prints important results from an x_mle-object
#' @author Kim Rand-Hendriksen


print.x_mle <- function(x_mle) {

  # if(NROW(x_mle$fullcoef)) {
  #   cat(paste0("\nFitted coefficients:\n"))
  #   print(x_mle$fullcoef)
  # }
  # if(NROW(x_mle$fixed_values)) {
  #   cat(paste0("\nFixed values:\n"))
  #   #print(ndf(Fixed = x_mle$fixed_values))
  #   print(x_mle$fixed_values)
  #
  # }
  if(NROW(x_mle$pars)) {
    cat(paste0("\nParameters:\n"))
    #print(ndf(Fixed = x_mle$fixed_values))
    print(x_mle$pars)

  }
  cat("\nMinima:\n")
  print(x_mle$minima)
}
