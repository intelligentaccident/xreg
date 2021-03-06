% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loglik_functions.R
\name{cont_normal}
\alias{cont_normal}
\title{Normal-distribution-based likelihood function for xregControl and xreg}
\usage{
cont_normal(d_df, sigma_est = 0)
}
\arguments{
\item{d_df}{data.frame provided by xreg.}

\item{sigma_est}{Subsequent parameters define required formulas to generate necessary variables, such as standard deviations. If these are not manually provided (which could be useful to model heteroscedasticity, etc.), xreg will get them from the functions.}
}
\description{
Likelihood function to use with xregControl
}
\details{
Likelihood functions for data.frame
}
