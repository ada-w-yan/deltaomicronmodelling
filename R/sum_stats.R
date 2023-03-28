#' calculate viral infectivity for a given set of active pathways
#'
#' @param pars named list of parameter values
#' @param pathway "both", "endosomal", "tmprss2", "both_no_IFITM", "endosomal_no_IFITM", "tmprss2_no_IFITM"
#' @return numeric vector with values of beta_E and beta_T
#'
calc_beta_pathway <- function(pars, pathway) {
  if(!grepl("no_inhibition", pathway)) {
    pars$beta_E <- pars$beta_E * (1 - pars$f_E)
    pars$beta_T <- pars$beta_T * (1 - pars$f_T)
  }
  if(grepl("_dependent", pathway)) {
    beta1 <- c(pars$beta_T, 0)
  } else if(grepl("independent", pathway)) {
    beta1 <- rep(pars$beta_E, 2)
  } else {
    beta1 <- c(pars$beta_T, 0) + rep(pars$beta_E, 2)
  }
  beta1
}

#' calculate R_0 for all combinations of active pathways for model 2
#'
#' @param pars named list of parameter values
#' @return named numeric vector with values of R_0; names are the pathways
#'
calc_R_0_camostat <- function(pars) {
  calc_R_0_pathway <- function(pathway) {

    pars$beta <- calc_beta_pathway(pars, pathway)
    R_0 <- sum(pars$beta * pars$omega_Inf * pars$T_0 / pars$delta) /
      (pars$kappa_Inf + sum(pars$beta * pars$T_0))

    R_0
  }

  pathways <- get_pathways(pars)
  R_0 <- vnapply(pathways, calc_R_0_pathway)

  names(R_0) <- paste0("R_0_", names(R_0))
  R_0
}

#' calculate r for all combinations of active pathways for model 2
#'
#' @param pars named list of parameter values
#' @return named numeric vector with values of R_0; names are the pathways
#'
calc_r_camostat <- function(pars, warn_if_neg = FALSE) {
  calc_r_pathway <- function(pathway) {
    if(grepl("independent", pathway)) {
      pars$beta_T <- 0
    } else if(grepl("_dependent", pathway)) {
      pars$beta_E <- 0
    }

    if(!grepl("no_inhibition", pathway)) {
      pars$beta_E <- pars$beta_E * (1 - pars$f_E)
      pars$beta_T <- pars$beta_T * (1 - pars$f_T)
    }

    eigenvalue_mat <- matrix(c(-pars$tau_E, 0, 0, pars$beta_E * sum(pars$T_0),
                               0, -pars$tau_T, 0, pars$beta_T * pars$T_0[1],
                               pars$tau_E, pars$tau_T, -pars$delta, 0,
                               0, 0, pars$omega_Inf, -(pars$kappa_Inf + (pars$beta_E + pars$beta_T) * pars$T_0[1] +
                                                     pars$beta_E * pars$T_0[2])), ncol = 4)
    ev <- eigen(eigenvalue_mat)
    ev <- ev$values
    ev <- ev[Im(ev) == 0]
    stopifnot(length(ev) > 0)
    r <- max(as.double(ev))
    if(warn_if_neg && r <= 0) {
      warning("negative value of r")
    }
    r
  }

  pathways <- get_pathways(pars)
  r <- vnapply(pathways, calc_r_pathway)

  names(r) <- paste0("r_", names(r))
  r
}

#' create closure to calculate summary statistics for model 2 and provide bounds for them
#'
#' @param parTab data frame
#' @return a list with elements
#' "calc_summary": closure to calculate
gen_summary_statistics_fn_camostat <-
  function(parTab) {

    dummy_pars <- parTab$names
    names(dummy_pars) <- parTab$names
    pathways <- get_pathways(dummy_pars)

    sum_stat_names <- c("R_0", "r", "doubling_time")

    sum_stat_names_pathway <- outer(sum_stat_names, pathways, paste0) %>% t %>%
      as.character
    # sum_stat_names_pathway <- c(sum_stat_names_pathway, "R_0_tmprss2_on_R_0_both", "R_0_endosomal_on_R_0_both")

    # function to transform parameters to calculate summary statistics

    transform_pars <- transform_pars_wrapper_clare(parTab)

    # functions to calculate summary statistics given parameter values

    # closure to generate function to calculate summary statistics given parameter values
    calc_summary <- function(values) {

      values <- transform_pars(values)

      calc_funcs <- list(calc_R_0_camostat, calc_r_camostat)

      sum_stats <-
        vapply(calc_funcs, function(f)
          f(values), double(length(pathways))) %>%
        as.numeric

      names(sum_stats) <- sum_stat_names_pathway[seq_along(sum_stats)]

      doubling_times <- log(2) / sum_stats[paste0("r", pathways)]
      names(doubling_times) <- paste0("doubling_time", pathways)
      sum_stats <- c(sum_stats, doubling_times)

      sum_stats
    }

    dummy_values <- calc_summary(parTab$values)
    par_names_plot <- names(dummy_values)
    lower_bound <- rep(0, length(par_names_plot))
    upper_bound <- rep(c(100, 1, 100), length(pathways))

    names(par_names_plot) <- names(lower_bound) <- names(upper_bound) <- par_names_plot

    # nice parameter names for plotting

    list(
      calc_summary = calc_summary,
      lower_bound = lower_bound,
      upper_bound = upper_bound,
      par_names_plot = par_names_plot
    )
  }

#' create closure to calculate summary statistics for model 1 and provide bounds for them
#'
#' @param parTab data frame
#' @return a list with elements
#' "calc_summary": closure to calculate
gen_summary_statistics_fn_two_strain <-
  function(parTab, strains) {

    dummy_pars <- parTab$names
    names(dummy_pars) <- parTab$names
    pathways <- ""
    sum_stat_names <- c("R_0", "r", "doubling_time")

    sum_stat_names_pathway <- outer(sum_stat_names, pathways, paste0) %>% t %>%
      as.character

    # function to transform parameters to calculate summary statistics

    transform_pars <- transform_pars_wrapper_two_strain(parTab, strains)

    # functions to calculate summary statistics given parameter values

    # closure to generate function to calculate summary statistics given parameter values
    calc_summary <- function(values) {

      values <- transform_pars(values)

      calc_funcs <- list(calc_R_0_clare, calc_r_clare)

      calc_sum_stats_single_strain <- function(values, strain) {
        sum_stats <- vapply(calc_funcs, function(f)
          f(values), double(length(pathways))) %>%
          as.numeric
        names(sum_stats) <- sum_stat_names_pathway[seq_along(sum_stats)]

        doubling_times <- log(2) / sum_stats[paste0("r", pathways)]
        names(doubling_times) <- paste0("doubling_time", pathways)
        sum_stats <- c(sum_stats, doubling_times)
        # names(sum_stats) <- paste0(names(sum_stats), ".", strain)

        sum_stats
      }

      sum_stats <- Map(calc_sum_stats_single_strain, values, strains) %>%
        do.call(c, .)
      sum_stats
    }

    dummy_values <- calc_summary(parTab$values)
    par_names_plot <- names(dummy_values)
    lower_bound <- upper_bound <- rep(0, length(par_names_plot)) # not maintaining
    # the plotting functionality any more
    # upper_bound <- rep(c(100, 1, 100, 100), length(pathways))#, 1, 1)

    names(par_names_plot) <- names(lower_bound) <- names(upper_bound) <- par_names_plot

    # nice parameter names for plotting

    list(
      calc_summary = calc_summary,
      lower_bound = lower_bound,
      upper_bound = upper_bound,
      par_names_plot = par_names_plot
    )
  }

#' calculate R_0 for model 1
#'
#' @param pars named list of parameter values
#' @return value of R_0
#'
calc_R_0_clare <- function(pars) {
  pars$beta * pars$omega_Inf * pars$T_0 / pars$delta /
    (pars$kappa_Inf + pars$beta * pars$T_0)
}

#' calculate r for model 1
#'
#' @param pars named list of parameter values
#' @return value of r
#'
calc_r_clare <- function(pars) {
  eigenvalue_mat <- matrix(c(-pars$tau, pars$tau, 0,
                             0, -pars$delta, pars$omega_Inf,
                             pars$beta * pars$T_0, 0,
                             -(pars$beta * pars$T_0 + pars$kappa_Inf)),
                           ncol = 3)
  ev <- eigen(eigenvalue_mat)
  ev <- ev$values
  ev <- ev[Im(ev) == 0]
  stopifnot(length(ev) > 0)
  r <- max(as.double(ev))
  r
}

#' list pathways
#'
#' @return character vector of pathway names
#'
get_pathways <- function(pars) {
  c("tmprss2_dependent", "tmprss2_independent", "both",
    "tmprss2_dependent_no_inhibition", "tmprss2_independent_no_inhibition", "both_no_inhibition")
}
