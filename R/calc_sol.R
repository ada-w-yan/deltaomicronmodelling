#' calculate solution to ODEs
#'
#' @param pars named character vector of parameter values
#' @param solving_time vector of times for which to solve the ODEs
#' @param model name of ODE model file
#' @return data frame with solution to ODE
#'
calc_sol_clare_model <- function(pars, solving_time, model = "model_teiv_dual") {

  model_filename <- paste0(model, ".R")

  gen <- compile_model(model_filename)
  par_names_mod <- get_mod_info(model_filename)$pars

  pars <- pars[par_names_mod]
  mod <- gen$new(user = pars)

  sol <- solve_ODE_error_handling(mod, solving_time)

  sol <- as.data.frame(sol)
  sol <- sol[,c("t", "V", "V_tot")]

  sol
}

#' calculate solution to ODEs for each pathway for model 2, and extract median, 95\% CI and max LL
#'
#' @param filename .RData file produced by lazymcmc
#' @param model string: "both", "endosomal", "tmprss2", "both_no_IFITM", "endosomal_no_IFITM", "tmprss2_no_IFITM"
#' @param ODE_model_filename name of ODE model file
#' @param solving_time_max time to solve until.   Default = 3 dpi.
#' @return tibble with columns model, t, V, T1... and probs, which can take the values 0.025, 0.5, 0.975 or "max_LL"
#'
#' @import dplyr
#' @import odin
calc_sol_chain <- function(filename, model, ODE_model_filename, solving_time_max = 3) {

  # load posterior samples
  chain <- get_MCMC_chain(filename)
  n_samples <- 1000 # number of samples to use to calculate 95% CI
  max_LL_params <- get_max_LL_params(filename) # get maximum likelihood parameters
  # sample uniformly from rest of chain
  chain <- chain[seq(1, nrow(chain), length.out = n_samples),] %>%
    as_tibble %>%
    select(names(max_LL_params))

  # times for which to solve ODE
  solving_time <- seq(0, solving_time_max, by = 1/24)

  gen <- compile_model(ODE_model_filename)
  par_names_mod <- get_mod_info(ODE_model_filename)$pars

  parTab <- load_into_list(filename, "parTab")$parTab
  transform_pars <- transform_pars_wrapper_clare(parTab)

  calc_sol <- function(pars) {

    pars <- transform_pars(pars)

    if(grepl("independent", model)) {
      pars$beta_T <- 0
    } else if(grepl("_dependent", model)) {
      pars$beta_E <- 0
    }

    if(!grepl("no_inhibition", model)) {
      pars$beta_E <- pars$beta_E * (1 - pars$f_E)
      pars$beta_T <- pars$beta_T * (1 - pars$f_T)
    }

    pars <- pars[par_names_mod]

    mod <- gen$new(user = pars)

    sol <- solve_ODE_error_handling(mod, solving_time)
    sol <- as.data.frame(sol)
    sol <- sol[,c("t", "V", "V_tot")]

    sol
  }

  probs <- c(0.025, .5, .975)

  prctiles <- apply(chain, 1, calc_sol) %>%
    bind_rows %>%
    group_by(t) %>%
    summarise(across(everything(), quantile, probs = probs)) %>%
    mutate(probs = c("lower", "median", "upper"))

  max_LL_trajectory <- calc_sol(max_LL_params) %>%
    mutate(probs = "max_LL")

  prctiles <- prctiles %>%
    bind_rows(max_LL_trajectory) %>%
    mutate(model = model)
  ungroup

  prctiles
}

#' calculate solution to ODEs for model 1 and extract median, 95\% CI and max LL
#'
#' @param filename .RData file produced by lazymcmc
#' @param ODE_model_filename name of ODE model file
#' @param solving_time_max time to solve until.   Default = 3 dpi
#' @return tibble with columns model, t, V, T1... and probs, which can take the values 0.025, 0.5, 0.975 or "max_LL"
#'
#' @import dplyr
#' @import odin
calc_sol_two_strain <- function(filename, ODE_model_filename, solving_time_max = 3) {

  # load posterior samples
  chain <- get_MCMC_chain(filename)
  n_samples <- 1000 # number of samples to use to calculate 95% CI

  max_LL_params <- get_max_LL_params(filename) # get maximum likelihood parameters
  # sample uniformly from rest of chain
  chain <- chain[seq(1, nrow(chain), length.out = n_samples),] %>%
    as_tibble %>%
    select(names(max_LL_params))# %>%
  # rbind(max_LL_params)

  # times for which to solve ODE
  solving_time <- seq(0, solving_time_max, by = 1/24)

  gen <- compile_model(ODE_model_filename)
  par_names_mod <- get_mod_info(ODE_model_filename)$pars

  parTab <- load_into_list(filename, "parTab")$parTab
  data_df <- load_into_list(filename, "data_df")$data_df
  strains <- unique(data_df$strain)
  transform_pars <- transform_pars_wrapper_two_strain(parTab, strains)

  calc_sol_strain <- function(pars, strain1) {

    pars <- transform_pars(pars)
    pars <- pars[[strain1]]

    pars <- pars[par_names_mod]

    mod <- gen$new(user = pars)

    sol <- solve_ODE_error_handling(mod, solving_time)
    sol <- as.data.frame(sol)
    sol <- sol[,c("t", "V", "V_tot")]
    sol
  }

  calc_prctiles_strain <- function(strain1) {
    probs <- c(0.025, .5, .975)

    prctiles <- apply(chain, 1, \(x) calc_sol_strain(x, strain1)) %>%
      bind_rows %>%
      group_by(t) %>%
      summarise(across(everything(), quantile, probs = probs)) %>%
      mutate(probs = c("lower", "median", "upper"))

    max_LL_trajectory <- calc_sol_strain(max_LL_params, strain1) %>%
      mutate(probs = "max_LL")

    prctiles <- prctiles %>%
      bind_rows(max_LL_trajectory) %>%
      mutate(strain = strain1) %>%
      ungroup

    prctiles
  }

  sol <- lapply(strains, calc_prctiles_strain) %>% bind_rows()
  sol
}

#' wrapper for calculating median and 95% CI trajectories from model fits
#' @param filename .RData file produced by lazymcmc
#' @param ODE_model_filename name of ODE model file
#' @param two_strain TRUE for model 1, FALSE for model 2
#'
#' @return tibble of median and 95% CI trajectories from model fits (for each pathway for model 2)
#'
calc_trajectories <- function(filename, ODE_model_filename, two_strain = FALSE) {
  parTab <- load_into_list(filename)$parTab
  pars <- parTab$values
  names(pars) <- parTab$names

  if(two_strain) {
    trajectories <- calc_sol_two_strain(filename, ODE_model_filename)
  } else {
    pathways <- get_pathways(pars)
    trajectories <- lapply(pathways, calc_sol_chain, ODE_model_filename = ODE_model_filename,
                           filename = filename) %>%
      bind_rows()
  }

  trajectories
}

#' plot model 1 fits
#'
#' @param strain1 "omicron" or "Delta"
#' @param Calu3 if TRUE plot Calu3 trajectories; if FALSE plot hNEC
#' @param trajectory_filename filename which trajectories were saved in
#'
#' @return filenames which plots are saved to
#'
plot_trajectories <- function(strain1, Calu3, trajectory_filename) {

  tmp <- strsplit(trajectory_filename, "/", fixed = TRUE)[[1]]
  dir_name <- paste0(tmp[-length(tmp)], "/", collapse = "")

  trajectories <- readRDS(trajectory_filename)

  data_df <- read_data(strain1, Calu3)

  trajectories <- readRDS(trajectory_filename)

  # if trajectories were calculated for more than one strain,
  # filter on strain we are plotting
  if("strain" %in% colnames(trajectories)) {
    trajectories <- trajectories %>%
      filter(strain == strain1)
  }

  data_df <- data_df %>%
    filter(t >= 0)
  plot_inner <- function(data_df, trajectories, PCR) {

    if(PCR) {
      trajectories <- trajectories %>%
        select(t, V_tot, probs) %>%
        pivot_wider(names_from = probs, values_from = V_tot) %>%
        filter(t >= 0)
      y_label = "RNA copy number/mL"
      y_max <- 1e12
    } else {
      trajectories <- trajectories %>%
        select(t, V, probs) %>%
        pivot_wider(names_from = probs, values_from = V) %>%
        filter(t >= 0)
      y_label = "Viral load (pfu/mL)"
      y_max <- 1e10
    }

    g <- ggplot(data_df) +
      geom_ribbon(data = trajectories, aes(x = t, ymin = lower, ymax = upper),
                  fill = "#D8D8D8") +
      geom_point(aes_string(x = "t", y = ifelse(PCR, "V_tot", "V"))) +
      geom_line(data = trajectories, aes(x = t, y = max_LL)) +
      geom_hline(yintercept = 10, linetype = "dashed") +
      scale_y_log10(y_label) +
      coord_cartesian(ylim = c(1, y_max)) +
      theme_bw() +
      xlab("Time (days)")
    g
  }

  plot_grid <- tibble(PCR = c(FALSE, TRUE)) %>%
    mutate(filename = paste0(dir_name, "model_predictions_", strain1,
                             ifelse(PCR, "_PCR", ""), ".png"))
  g <- lapply(plot_grid$PCR, function(x) plot_inner(data_df, trajectories, x))
  Map(function(x, y) ggsave(x, y, width = 3, height = 3),
      plot_grid$filename, g)
}

#' plot model 2 fits
#'
#' @param strain1 "omicron" or "Delta"
#' @param Calu3 if TRUE plot Calu3 trajectories; if FALSE plot hNEC
#' @param trajectory_filename filename which trajectories were saved in
#'
#' @return filenames which plots are saved to
#'
plot_trajectories_camostat_amphoB <- function(strain1, Calu3 = FALSE, trajectory_filename) {

  lod <- 10

  tmp <- strsplit(trajectory_filename, "/", fixed = TRUE)[[1]]
  dir_name <- paste0(tmp[-length(tmp)], "/", collapse = "")

  trajectories <- readRDS(trajectory_filename)

  data_df <- read_data(strain = strain1, Calu3 = Calu3)

  if(Calu3) {
    drugs <- c("No drug", "Camostat", "Amphotericin B", "Camostat + AmphoB")
  } else {
    drugs <- c("No drug", "Camostat", "Amphotericin B")
  }

  data_df <- data_df %>%
    pivot_longer(starts_with("V"), names_to = "model", values_to = "V") %>%
    mutate(PCR = grepl("tot", model),
           model = gsub("_tot", "", model),
           model = recode(model, V = "No drug",
                          V_Camostat = "Camostat",
                          V_AmphoB = "Amphotericin B",
                          V_Camostat_AmphoB = "Camostat + AmphoB")) %>%
    filter(model %in% drugs, t >= 0) %>%
    mutate(model = factor(model, levels = drugs),
           V = pmax(V, lod))

  plot_inner <- function(data_df, trajectories, PCR, facet1 = FALSE) {

    if(PCR) {
      trajectories <- trajectories %>%
        select(t, V_tot, probs, model) %>%
        pivot_wider(names_from = probs, values_from = V_tot) %>%
        filter(t >= 0)
      y_label = "RNA copy number/mL"
      y_max <- 1e12
      data_df <- data_df %>%
        filter(PCR)
    } else {
      trajectories <- trajectories %>%
        select(t, V, probs, model) %>%
        pivot_wider(names_from = probs, values_from = V) %>%
        filter(t >= 0)
      if("PCR" %in% colnames(data_df)) {
        data_df <- data_df %>% filter(PCR == FALSE)
      }
      y_label = "Viral load (pfu/mL)"
      y_max <- 1e10
    }

    trajectories <- trajectories %>%
      filter(model %in% c("endosomal", "both", "both_no_IFITM", "endosomal_no_IFITM")) %>%
      mutate(model = recode(model, endosomal = "Camostat",
                            both = "No drug",
                            both_no_IFITM = "Amphotericin B",
                            endosomal_no_IFITM = "Camostat + AmphoB")) %>%
      filter(model %in% drugs) %>%
      mutate(model = factor(model, levels = drugs))
    g <- ggplot(data_df) +
      geom_ribbon(data = trajectories, aes(x = t, ymin = lower, ymax = upper, fill = model, group = model),
                  alpha = .3) +
      geom_point(aes(x = t, y = V,
                     color = model, group = model)) +
      geom_line(data = trajectories, aes(x = t, y = max_LL, color = model, group = model)) +
      scale_y_log10(y_label) +
      coord_cartesian(ylim = c(1, y_max)) +
      theme_bw() +
      xlab("Time (days)") +
      theme(legend.position = "none") +
      geom_hline(yintercept = 10, linetype = "dashed")

    if(facet1) {
      g <- g +
        facet_wrap(~model, ncol = 1)
    }
    g
  }

  plot_grid <- expand_grid(PCR = c(FALSE, TRUE), facet1 = c(FALSE, TRUE)) %>%
    mutate(filename = paste0(dir_name, "model_predictions_", strain1,
                             ifelse(PCR, "_PCR", ""),
                             ifelse(facet1, "_facet", ""), ".png"))
  g <- Map(function(x, y) plot_inner(data_df, trajectories, x, y), plot_grid$PCR, plot_grid$facet1)

  Map(function(x, y, z) ggsave(x, y, width = ifelse(z, 2, 3), height = ifelse(z, 4, 3)),
      plot_grid$filename, g, plot_grid$facet1)


}
