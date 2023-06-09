% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc_helpers.R
\name{generate_random_start}
\alias{generate_random_start}
\title{Generate random starting values from a parameter table}
\usage{
generate_random_start(parTab, random_seed, f)
}
\arguments{
\item{parTab}{a parameter table data frame of the form needed by the \code{run_mcmc} function}

\item{random_seed}{random seed for generation
The info used here are the columns
parTab$fixed: 1 indicates the parameter is fixed, 0 indicates the parameter is to be fitted
parTab$lower_bound: lower bounds for parameters allowed for proposal steps
parTab$upper_bound: lower bounds for parameters allowed for proposal steps}
}
\value{
a data frame of the same for as parTab, but the values of the fitted
parameters are chosen from a uniform distribution with bounds \link{lower_bound, upper_bound}
}
\description{
\code{generate_random_start} generates random starting values for all parameters to be fitted
}
