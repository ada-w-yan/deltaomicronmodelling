% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_sol.R
\name{calc_sol_clare_model}
\alias{calc_sol_clare_model}
\title{calculate solution to ODEs}
\usage{
calc_sol_clare_model(pars, solving_time, model = "model_teiv_dual")
}
\arguments{
\item{pars}{named character vector of parameter values}

\item{solving_time}{vector of times for which to solve the ODEs}

\item{model}{name of ODE model file}
}
\value{
data frame with solution to ODE
}
\description{
calculate solution to ODEs
}
