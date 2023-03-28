#' fit model 1 to data
#'
#' @param save_dir string ending in "/": directory to save results in
#' @param Calu3 logical. if TRUE, fit model to Calu3 data, if FALSE, use hNEC data
#' @param length_run integer.  If length == 1, short run (for testing purposes); if length == 2, normal run; if length = 3, extra long run.
#' otherwise, run until convergence
#' @param run_flag logical.  if TRUE, run and postprocess; if FALSE, postprocess only
#' @return NULL (results saved to file)
#' @export
run_exp_two_strain <- function(save_dir,
                               Calu3,
                               length_run = 2,
                               run_flag = TRUE) {

  # number of parallel chains to run
  n_replicates <- 3
  # number of temperatures for parallel tempering
  n_temperatures <- 5
  # set random number generator seed
  seeds <- set_default_seeds(n_replicates, n_temperatures)

  filenames = paste0(save_dir, vapply(seeds, function(x) x[1], double(1)))


  n_replicates <- 3

  CREATE_PRIOR_FUNC <- function(x) {
    function(y) 0
  }

  # define parameters for MCMC
  if(length_run == 1) {
    mcmcPars <- gen_mcmcPars_timing()
  } else if(length_run == 2){
    mcmcPars <- gen_mcmcPars(temperature = seq_len(n_temperatures), parallel_tempering_iter = 10)
  } else if (length_run == 3) {
    mcmcPars <- gen_mcmcPars_long(temperature = seq_len(n_temperatures), parallel_tempering_iter = 10)
  } else {
    stop("unknown length run")
  }

  # define function to generate random starting values
  generate_start_values <- function(x, y, z, f) generate_random_start(x, y, f)

  strains <- c("omicron", "delta")
  get_data <- function() read_two_strain_data(strains, Calu3)

  specify_parameters <- function() specify_two_strain_fn(strains,
                                                         Calu3 = Calu3)

  gen_summary_statistics_fn <- \(x) gen_summary_statistics_fn_two_strain(x, strains)
  CREATE_POSTERIOR_FUNC_fn <- CREATE_POSTERIOR_FUNC_fn_two_strain

  # run MCMC
  run_all(starting_point_seed = seeds,
          filenames = filenames,
          run_parallel = TRUE,
          mvr = FALSE,
          specify_parameters = specify_parameters,
          get_data = get_data,
          CREATE_POSTERIOR_FUNC = CREATE_POSTERIOR_FUNC_fn,
          CREATE_PRIOR_FUNC = CREATE_PRIOR_FUNC,
          mcmcPars = mcmcPars,
          generate_start_values = generate_start_values,
          gen_summary_statistics = gen_summary_statistics_fn,
          calc_residuals = NULL,
          plot_flags = c(trace = FALSE,
                         diagnostics= FALSE,
                         priors = FALSE,
                         posteriors = FALSE,
                         bivariate = FALSE),
          plot_dynamics = NULL,
          plot_residuals = NULL,
          run = c("run" = run_flag, "process" = TRUE))
}

#' fit model 2 to data
#'
#' @param save_dir string ending in "/": directory to save results in
#' @param strain character vector.  name of strain.
#' @param Calu3 logical. if TRUE, fit model to Calu3 data, if FALSE, use hNEC data
#' @param tau_E fixed value of tau_E
#' @param tau_T fixed value of tau_T
#' @param length_run integer.  If length == 1, short run (for testing purposes); if length == 2, normal run; if length = 3, extra long run.
#' otherwise, run until convergence
#' @param run_flag logial.  if TRUE, run and postprocess; if FALSE, postprocess only
#' @return NULL (results saved to file)
#' @export
run_exp_camostat_amphoB <- function(save_dir,
                                    strain,
                                    Calu3,
                                    tau_E = 4,
                                    tau_T = 4,
                                    length_run = 2,
                                    run_flag = TRUE) {

  different_eclipse <- tau_E != tau_T

  # number of parallel chains to run
  n_replicates <- 3
  # number of temperatures for parallel tempering
  n_temperatures <- 5
  # set random number generator seed
  seeds <- set_default_seeds(n_replicates, n_temperatures)

  filenames = paste0(save_dir, vapply(seeds, function(x) x[1], double(1)))

  CREATE_PRIOR_FUNC <- function(x) {
    function(y) 0
  }

  # define parameters for MCMC
  if(length_run == 1) {
    mcmcPars <- gen_mcmcPars_timing()
  } else if(length_run == 2){
    mcmcPars <- gen_mcmcPars(temperature = seq_len(n_temperatures), parallel_tempering_iter = 10)
  } else if (length_run == 3) {
    mcmcPars <- gen_mcmcPars_long(temperature = seq_len(n_temperatures), parallel_tempering_iter = 10)
  } else {
    stop("unknown length run")
  }

  # define function to generate random starting values
  generate_start_values <- function(x, y, z, f) generate_random_start(x, y, f)

  get_data <- function() read_data(strain, Calu3)

  specify_parameters <- function() specify_clare_IFITM_fn(Calu3 = Calu3,
                                                          tau_E = tau_E,
                                                          tau_T = tau_T)


  gen_summary_statistics_fn <- gen_summary_statistics_fn_camostat
  CREATE_POSTERIOR_FUNC_fn <- CREATE_POSTERIOR_FUNC_fn_IFITM

  # run MCMC
  run_all(starting_point_seed = seeds,
          filenames = filenames,
          run_parallel = TRUE,
          mvr = FALSE,
          specify_parameters = specify_parameters,
          get_data = get_data,
          CREATE_POSTERIOR_FUNC = CREATE_POSTERIOR_FUNC_fn,
          CREATE_PRIOR_FUNC = CREATE_PRIOR_FUNC,
          mcmcPars = mcmcPars,
          generate_start_values = generate_start_values,
          gen_summary_statistics = gen_summary_statistics_fn,
          calc_residuals = NULL,
          plot_flags = c(trace = FALSE,
                         diagnostics= FALSE,
                         priors = FALSE,
                         posteriors = FALSE,
                         bivariate = FALSE),
          plot_dynamics = NULL,
          plot_residuals = NULL,
          run = c("run" = run_flag, "process" = TRUE))
}




