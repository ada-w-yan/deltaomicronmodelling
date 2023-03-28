#' read data
#'
#' @param strain1 "omicron" or "delta"
#' @param Calu3 logical. if TRUE, read Calu3 data, if FALSE, read hNEC data
#' @import tidyr
#' @import dplyr
#' @return a data frame
read_data <- function(strain1, Calu3) {
  filename <- paste0("data/PCR_", strain1, ifelse(Calu3, "_Calu3", ""), ".csv")

  data_df <- readr::read_csv(filename) %>%
    rename_with(~ gsub("_no_drug", "", .x, fixed = TRUE)) %>%
    mutate(t = t/24) %>%
    filter(t >= 0) %>%
    as.data.frame()
  return(data_df)

}

#' read data for two strains
#'
#' @param strains vector of strain names
#' @param Calu3 logical. if TRUE, read Calu3 data, if FALSE, read hNEC data
#' @import tidyr
#' @import dplyr
#' @return a data frame
read_two_strain_data <- function(strains, Calu3) {

  read_wrapper <- function(strain1) {
    read_data(strain1, Calu3) %>% mutate(strain = strain1)
  }

  lapply(strains, read_wrapper) %>%
    bind_rows()
}
