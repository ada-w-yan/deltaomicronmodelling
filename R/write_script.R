#' write script to run job on cluster
#'
#' @param R_filename string ending in .R with R script to be run
#' @return file with same name as R_filename but without extension
write_script <- function(R_filename) {
  script_filename <- strsplit(R_filename, ".", fixed = TRUE)[[1]][[1]]
  script_lines <- c(
    "#PBS -lselect=1:ncpus=3:mem=32gb",
    "#PBS -lwalltime=72:00:0",
    "",
    "# Load modules for any applications",
    "",
    "module load anaconda3/personal",
  "source activate r410",
  "",
  "# Change to the submission directory",
  "",
  "cd $PBS_O_WORKDIR",
  "",
  "# Run program",
  "",
  "export MC_CORES=3",
  "",
  paste0("R CMD BATCH --slave $HOME/git_repos/deltaomicron/",
         R_filename, " $HOME/git_repos/deltaomicron/", script_filename, ".out"),
  "",
  "mkdir $WORK/$PBS_JOBID",
  "cp * $WORK/$PBS_JOBID"
  )
  fileConn<-file(script_filename)
  writeLines(script_lines, fileConn)
  close(fileConn)
}

#' write scripts for omicron and delta on hNECs and Calu3
#'
#' @param base_name string ending in .R which forms basis of filename
#' @return scripts
write_script_wrapper <- function(base_name) {
  par_grid <- expand_grid(strain = c("omicron", "delta"),
                          cells = c("hNEC", "Calu3")) %>%
    mutate(R_filename = paste0("scripts/", strain, "_", cells, "_", base_name))
  lapply(par_grid$R_filename, write_script)
}
