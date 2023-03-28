#' creates a closure to calculate the likelihood and model predictions for given
#' parameter values for model 1
#' @param parTab data frame of parameter values created by specify_two_strain_fn
#' @param data_df data frame of data created by read_two_strain_data
#' @param PRIOR_FUNC prior function f <- function(x) 0
#' @return closure to calculate likelihood and model predictions for given
#' parameter values
CREATE_POSTERIOR_FUNC_fn_two_strain <-
  function(parTab, data_df, PRIOR_FUNC) {

    strains <- unique(data_df$strain)

    transform_pars <- transform_pars_wrapper_two_strain(parTab, strains)

    prediction_compartments <- c("V", "V_tot")

    # function to calculate log likelihood
    f <- function(pars) {
      if(PRIOR_FUNC(pars) == -Inf) {
        return(list(lik = -Inf, misc = NA))
      }

      transformed_pars <- transform_pars(pars)

      calc_likelihood_strain <- function(transformed_pars, strain1) {

        transformed_pars <- transformed_pars[[strain1]]
        data_df <- data_df %>% filter(strain == strain1)

        sigma_Inf <- transformed_pars[["sigma_Inf"]]
        obs_threshold <- transformed_pars[["obs_threshold"]]
        sigma_RNA <- transformed_pars[["sigma_RNA"]]

        sigma_vec <- rep(c(sigma_Inf, sigma_RNA), each = length(prediction_compartments)/2)

        obs_threshold_vec <- rep(obs_threshold, length(sigma_vec)) # no PCR data points below threshold so can set to whatever

        solving_time <- make_solving_time(data_df$t)

        sol <- calc_sol_clare_model(transformed_pars, solving_time, "model_teiv_dual")

        sol <- sol[,c("t", prediction_compartments)]

        individual_lik <- calc_LL_log10normal_threshold_wrapper(
          data_df,
          prediction_compartments,
          sol,
          prediction_compartments,
          sigma_vec,
          obs_threshold_vec)

        misc <- sol %>%
          as.list %>%
          lapply(name_predictions) %>%
          unlist
        names(misc) <- sub(".", paste0(".", strain1, "."), names(misc), fixed = TRUE)

        lik <- sum(individual_lik)

        list(lik = lik, misc = misc)
      }

      lik_list <- lapply(strains, calc_likelihood_strain, transformed_pars = transformed_pars)
      lik <- vnapply(lik_list, \(x) x$lik) %>% sum
      misc <- lapply(lik_list, \(x) x$misc) %>% do.call(c, .)

      list(lik = lik, misc = misc)
    }
    f
  }

#' creates a closure to calculate the likelihood and model predictions for given
#' parameter values for model 2
#' @param parTab data frame of parameter values created by specify_clare_IFITM_fn
#' @param data_df data frame of data created by read_data
#' @param PRIOR_FUNC prior function f <- function(x) 0
#' @return closure to calculate likelihood and model predictions for given
#' parameter values
CREATE_POSTERIOR_FUNC_fn_IFITM <-
  function(parTab, data_df, PRIOR_FUNC) {

    transform_pars <- transform_pars_wrapper_clare(parTab)

    prediction_compartments <- c("V", "V_Camostat", "V_AmphoB", "V_Camostat_AmphoB")
    prediction_compartments <- intersect(prediction_compartments, colnames(data_df))
    prediction_compartments <- c(prediction_compartments,
                                 gsub("V", "V_tot", prediction_compartments))


    # function to calculate log likelihood
    f <- function(pars) {
      if(PRIOR_FUNC(pars) == -Inf) {
        return(list(lik = -Inf, misc = NA))
      }

      transformed_pars <- transform_pars(pars)

      sigma_Inf <- transformed_pars[["sigma_Inf"]]
      obs_threshold <- transformed_pars[["obs_threshold"]]
      sigma_RNA <- transformed_pars[["sigma_RNA"]]

      sigma_vec <- rep(c(sigma_Inf, sigma_RNA), each = length(prediction_compartments)/2)

      obs_threshold_vec <- rep(obs_threshold, length(sigma_vec)) # no PCR data points below threshold so can set to whatever

      solving_time <- make_solving_time(data_df$t)

      model_name <- "model_different_eclipse_dual"
      # calculate solution for no drugs
      transformed_pars_no_drug <- transformed_pars
      transformed_pars_no_drug$beta_T <- transformed_pars_no_drug$beta_T * (1 - transformed_pars$f_T)

      transformed_pars_no_drug$beta_E <- transformed_pars_no_drug$beta_E * (1 - transformed_pars$f_E)

      sol <- calc_sol_clare_model(transformed_pars_no_drug, solving_time, model_name)

      # calculate solution for Camostat
      transformed_pars_camostat <- transformed_pars_no_drug
      transformed_pars_camostat$beta_T <- 0

      transformed_pars_amphoB <- transformed_pars

      if("V_Camostat_AmphoB" %in% prediction_compartments) {
        transformed_pars_camostat_amphoB <- transformed_pars
        transformed_pars_camostat_amphoB$beta_T <- 0
      }


      sol_Camostat <- calc_sol_clare_model(transformed_pars_camostat, solving_time, model_name)

      # calculate solution for AmphoB

      sol_AmphoB <- calc_sol_clare_model(transformed_pars_amphoB, solving_time, model_name)

      if("V_Camostat_AmphoB" %in% prediction_compartments) {
        sol_Camostat_AmphoB <- calc_sol_clare_model(transformed_pars_camostat_amphoB, solving_time, model_name)
      }

      sol <- data.frame(t = sol$t,
                        V = sol$V,
                        V_Camostat = sol_Camostat$V,
                        V_AmphoB = sol_AmphoB$V,
                        V_tot = sol$V_tot,
                        V_tot_Camostat = sol_Camostat$V_tot,
                        V_tot_AmphoB = sol_AmphoB$V_tot)

      if("V_Camostat_AmphoB" %in% prediction_compartments) {
        sol <- cbind(sol,
                     data.frame(V_Camostat_AmphoB = sol_Camostat_AmphoB$V,
                                V_tot_Camostat_AmphoB = sol_Camostat_AmphoB$V_tot))
      }

      individual_lik <- calc_LL_log10normal_threshold_wrapper(
        data_df,
        prediction_compartments,
        sol,
        prediction_compartments,
        sigma_vec,
        obs_threshold_vec)

      misc <- sol %>%
        as.list %>%
        lapply(name_predictions) %>%
        unlist

      lik <- sum(individual_lik)

      list(lik = lik, misc = misc)

    }
    f
  }

# function to name a vector of predictions
name_predictions <- function(prediction) {
  names(prediction) <- as.character(seq_along(prediction))
  prediction
}
