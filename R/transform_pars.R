#' create a closure for model 2 which transforms an unnamed numeric vector into
#' lists of named parameters to calculate the likelihood
#'
#' @param parTab data frame containing priors
#' @return a function which transforms an unnamed numeric vector into
#' lists of named parameters to calculate the likelihood
transform_pars_wrapper_clare <- function(parTab) {

  par_names <- parTab$names
  log10_ind <- grep("log10", par_names)
  par_names <- sub("log10_", "", par_names)

  transform_pars <- function(pars) {
    pars[log10_ind] <- 10 ^(pars[log10_ind])
    names(pars) <- par_names
    pars <- as.list(pars)
    pars$V_tot0 <- pars$V_0 * pars$s_rna
    pars$V_0 <- pars$V_0 * pars$s_inf
    pars$L_0 <- pars$I_0 <- 0

    # for model 2 only
    if("prop_ace2_pos" %in% names(pars)) {
      if(!("tau_E" %in% names(pars))) {
        pars$tau_E <- pars$tau_T <- pars$tau
      }
      pars$L_0 <- double(3)
      pars$I_0 <- double(2)

      pars$T_0 <- pars$T_0 * pars$prop_ace2_pos * c(pars$prop_tmprss2_pos, 1 - pars$prop_tmprss2_pos)
      # cell type 1 is tmprss2+
    }

    pars
  }
  transform_pars
}

#' create a closure for model 1 which transforms an unnamed numeric vector into
#' lists of named parameters to calculate the likelihood
#'
#' @param parTab data frame containing priors
#' @return a function which transforms an unnamed numeric vector into
#' lists of named parameters to calculate the likelihood
transform_pars_wrapper_two_strain <- function(parTab, strains) {
  par_names <- parTab$names
  log10_ind <- grep("log10", par_names)
  par_names <- sub("log10_", "", par_names)

  transform_pars_one_strain <- transform_pars_wrapper_clare(parTab)

  transform_pars <- function(pars) {
    pars <- transform_pars_one_strain(pars)
    get_pars_strain <- function(pars, strain) {
      exc_strain <- setdiff(strains, strain)
      # the dot is needed so as to not exclude log10_delta when strain == "delta"
      pars <- pars[!grepl(paste0(".", exc_strain), names(pars), fixed = TRUE)]
      names(pars) <- sub(paste0(".", strain), "", names(pars))
      pars
    }
    pars_two_strain <- lapply(strains, get_pars_strain, pars = pars)
    names(pars_two_strain) <- strains
    pars_two_strain
  }
  transform_pars
}
