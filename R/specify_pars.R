#' make parameter table for model 1
#' @param strains vector of strain names
#' @param Calu3 if TRUE, parameter values for Calu-3 cells; if FALSE, values for hNECs
#' @return data_df frame with columns
#' values: default parameter value
#' names: internal name for parameter
#' fixed:  if 0, fit parameter, otherwise fix parameter
#' lower_bound: lower bound of parameter for fitting
#' upper_bound: upper bound of parameter for fitting
#' steps: width of proposal distribution in parameter space
#' names_plot: character vector which which to label parameter distribution plot
specify_two_strain_fn <- function(strains, Calu3) {

  shared_pars <- c("log10_s_inf", "log10_s_rna")

  # prop ace2 pos cells
  T_0 <- if(Calu3) 1.2e6 else (5e5 * (0.0376 + 0.0158 + 0.0043 + 0.0471) / 4)
  parTab <- specify_clare_parameters_fn(Calu3 = FALSE)
  parTab[parTab$names == "T_0", "values"] <- T_0

  parTab_fixed_or_shared <- parTab %>% filter(fixed == 1 | names %in% shared_pars)
  parTab_fitted_separate <- parTab %>% filter(!(fixed == 1 | names %in% shared_pars))
  append_strain_name <- function(parTab, strain) {
    parTab %>%
      mutate(names = paste0(names, ".", strain))
  }

  parTab_fitted_separate_two_strain <- lapply(strains, append_strain_name, parTab = parTab_fitted_separate) %>%
    bind_rows

  parTab <- bind_rows(parTab_fixed_or_shared, parTab_fitted_separate_two_strain) %>%
    as.data.frame
  parTab
}

#' make parameter table for model 2
#' @param Calu3 if TRUE, parameter values for Calu-3 cells; if FALSE, values for hNECs
#' @param tau_E value of tau_E
#' @param tau_T value of tau_T
#' @return data_df frame with columns
#' values: default parameter value
#' names: internal name for parameter
#' fixed:  if 0, fit parameter, otherwise fix parameter
#' lower_bound: lower bound of parameter for fitting
#' upper_bound: upper bound of parameter for fitting
#' steps: width of proposal distribution in parameter space
#' names_plot: character vector which which to label parameter distribution plot
specify_clare_IFITM_fn <- function(Calu3, tau_E, tau_T) {

  parTab <- specify_clare_parameters_fn(Calu3)

  default_step <- 0.1

  parTab <-
    rbind(
      parTab,
      data.frame(
        # Table 1, weight all studies equally
        values = ifelse(Calu3, 1, (0.0376 + 0.0158 + 0.0043 + 0.0471) / 4),
        names = "prop_ace2_pos",
        fixed = 1,
        lower_bound = -Inf,
        upper_bound = Inf,
        steps = default_step,
        names_plot = "$ACE2^+$"
      )
    )

  parTab <-
    rbind(
      parTab,
      data.frame(
        # Table 1, weight all studies equally
        values = ifelse(Calu3, 1, (1.09 / 3.76 + 0.73 / 1.58 + 0.06 / 0.43 + 2.52 / 4.71) / 4),
        # proportion tmprss2+ out of ace2+
        names = "prop_tmprss2_pos",
        fixed = 1,
        lower_bound = -Inf,
        upper_bound = Inf,
        steps = default_step,
        names_plot = "$TMPRSS2^+$"
      )
    )

  parTab <-
    rbind(
      parTab,
      data.frame(
        values = -5,
        names = "log10_beta_E",
        fixed = 0,
        lower_bound = -12,
        upper_bound = -2,
        steps = default_step,
        names_plot = "$\\log_{10}\\beta_E$"
      )
    )

  parTab <-
    rbind(
      parTab,
      data.frame(
        values = -5,
        names = "log10_beta_T",
        fixed = 0,
        lower_bound = -12,
        upper_bound = -2,
        steps = default_step,
        names_plot = "$\\log_{10}\\beta_T$"
      )
    )

  different_eclipse <- tau_E != tau_T

  if(different_eclipse) {
    parTab <-
      rbind(
        parTab,
        data.frame(
          values = tau_E,
          names = "tau_E",
          fixed = 1,
          lower_bound = -Inf,
          upper_bound = Inf,
          steps = default_step,
          names_plot = "$\\tau_E$"
        )
      )

    parTab <-
      rbind(
        parTab,
        data.frame(
          values = tau_T,
          names = "tau_T",
          fixed = 1,
          lower_bound = -Inf,
          upper_bound = Inf,
          steps = default_step,
          names_plot = "$\\tau_T$"
        )
      )

    parTab <- parTab[parTab$names != "tau",]
  }

  parTab <-
    rbind(
      parTab,
      data.frame(
        values = .5,
        names = "f_E",
        fixed = 0,
        lower_bound = 0,
        upper_bound = 1,
        steps = default_step,
        names_plot = "f_E"
      )
    )

    parTab <-
      rbind(
        parTab,
        data.frame(
          values = .5,
          names = "f_T",
          fixed = 0,
          lower_bound = 0,
          upper_bound = 1,
          steps = default_step,
          names_plot = "f_T"
        )
      )
  parTab <- parTab[parTab$names != "log10_beta",]
  parTab
}

#' make parameter table
#'
#' @param Calu3 if TRUE, parameter values for Calu-3 cells
#' @return data_df frame with columns
#' values: default parameter value
#' names: internal name for parameter
#' fixed:  if 0, fit parameter, otherwise fix parameter
#' lower_bound: lower bound of parameter for fitting
#' upper_bound: upper bound of parameter for fitting
#' steps: width of proposal distribution in parameter space
#' names_plot: character vector which which to label parameter distribution plot
specify_clare_parameters_fn <- function(Calu3) {

  default_step <- 0.1

  # inefficient but easier to read
  # check lower and upper bounds for everything

  # initial number of target cells
  T_0 <- ifelse(Calu3, 1.2e6, 5e5)

  parTab <- data.frame(
    values = log10(10),
    names = "log10_kappa_Inf",
    fixed = 1,
    lower_bound = -Inf,
    upper_bound = Inf,
    steps = default_step,
    names_plot = "$\\log_{10}\\kappa_{inf}$",
    stringsAsFactors = FALSE
  )

  parTab <- rbind(parTab, data.frame(
    values = log10(0.05 *24),
    names = "log10_kappa_RNA",
    fixed = 1,
    lower_bound = -Inf,
    upper_bound = Inf,
    steps = default_step,
    names_plot = "$\\log_{10}\\kappa_{RNA}$",
    stringsAsFactors = FALSE
  ))

  parTab <-
    rbind(
      parTab,
      data.frame(
        values = -5,
        names = "log10_beta",
        fixed = 0,
        lower_bound = -12,
        upper_bound = -2,
        steps = default_step,
        names_plot = "$\\log_{10}\\beta$"
      )
    )

  parTab <-
    rbind(
      parTab,
      data.frame(
        values = 4,
        names = "tau",
        fixed = 1,
        lower_bound = -Inf,
        upper_bound = Inf,
        steps = default_step,
        names_plot = "$\\tau$"
      )
    )

  parTab <-
    rbind(
      parTab,
      data.frame(
        values = log10(1.7),
        names = "log10_delta",
        fixed = 1,
        lower_bound = -Inf,
        upper_bound = Inf,
        steps = default_step,
        names_plot = "$\\log_{10}\\delta$"
      )
    )

  parTab <- rbind(parTab,
                  data.frame(values = 0,
                             names = "log10_omega_Inf",
                             fixed = 0,
                             lower_bound = 0,
                             upper_bound = 10,
                             steps = default_step,
                             names_plot = "$\\log_{10}\\omega_{Inf}$"))

  ## initial conditions

  parTab <- rbind(
    parTab,
    data.frame(
      values = ifelse(Calu3, 1.2e3, 2.5e4),
      names = "V_0",
      fixed = 1,
      lower_bound = -Inf,
      upper_bound = Inf,
      steps = default_step,
      names_plot = "V0"
    )
  )

  if(Calu3) {
    s_inf <- 4.95e-3
    s_rna <- 17.85
  } else {
    s_inf <- 3.86e-2
    s_rna <- 334.64
  }
  parTab <- rbind(
    parTab,
    data.frame(
      values = log10(s_inf),
      names = "log10_s_inf",
      fixed = 0,
      lower_bound = -20,
      upper_bound = 1,
      steps = default_step,
      names_plot = "log10_s_inf"
    )
  )

  parTab <- rbind(
    parTab,
    data.frame(
      values = log10(s_rna),
      names = "log10_s_rna",
      fixed = 0,
      lower_bound = -20,
      upper_bound = 4,
      steps = default_step,
      names_plot = "log10_s_rna"
    )
  )

  parTab <- rbind(
    parTab,
    data.frame(
      values = 0.5,
      names = "log10_omega_RNA",
      fixed = 0,
      lower_bound = 0,
      upper_bound = 10,
      steps = default_step,
      names_plot = "log10_omega_RNA"
    )
  )

  parTab <-
    rbind(
      parTab,
      data.frame(
        values = T_0,
        names = "T_0",
        fixed = 1,
        lower_bound = -Inf,
        upper_bound = Inf,
        steps = default_step,
        names_plot = "$T_0$"
      )
    )

  parTab <-
    rbind(
      parTab,
      data.frame(
        values = .24,
        names = "sigma_Inf",
        fixed = 1,
        lower_bound = -Inf,
        upper_bound = Inf,
        steps = default_step,
        names_plot = "$\\sigma_Inf$"
      )
    )

  parTab <-
    rbind(
      parTab,
      data.frame(
        values = .22,
        names = "sigma_RNA",
        fixed = 1,
        lower_bound = -Inf,
        upper_bound = Inf,
        steps = default_step,
        names_plot = "$\\sigma_{RNA}$"
      )
    )

  parTab <-
    rbind(
      parTab,
      data.frame(
        values = 10,
        names = "obs_threshold",
        fixed = 1,
        lower_bound = -Inf,
        upper_bound = Inf,
        steps = default_step,
        names_plot = "obs threshold"
      )
    )


  # observation threshold of 10 pfu.  treat values below threshold as censored

  # convert logicals to numerics
  parTab$fixed <- as.numeric(parTab$fixed)
  # set bounds for fixed parameters to c(-Inf Inf)
  parTab[parTab$fixed == 1, "lower_bound"] <- -Inf
  parTab[parTab$fixed == 1, "upper_bound"] <- Inf

  parTab
}
