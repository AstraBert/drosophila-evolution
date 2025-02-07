# Generate simulate new sumstats from a matrix of parameters
library(here)
#######################################################################################################################
## Functions
#######################################################################################################################

# DIYABC_CMD <- "/save_home/pbastide/diyabc/build/src-JMC-C++/general -p ./"
# DIYABC_CMD <- "diyabc -p ./"
DIYABC_CMD <- "./diyabc-RF-linux-v1.1.54 -p ./"

# function to match parameter orders (si pas le bon ordre de parametres dans headerRF.txt)
get_parameter_order <- function(scenario = paste0(1:6)) {
  #scenario <- paste0(scen.ref)
  param_order <- switch(scenario,
                        "1" = c("scenario", "N1", "N2", "N3", "N4", "raa", "ta", "t1", "N5", "t2", "Nanc"),
                        "2" = c("scenario", "N1", "N2", "N3", "N4", "raa", "ta", "t1", "N5", "t2", "Nanc"),
                        "3" = c("scenario", "N1", "N2", "N3", "N4", "raa", "ta", "t1", "N5", "t2", "Nanc"),
                        "4" = c("scenario", "N1", "N2", "N3", "N4", "raa", "ta", "t1", "N5", "t2", "Nanc"),
                        "5" = c("scenario", "N1", "N2", "N3", "N4", "raa", "ta", "t1", "N5", "t2", "Nanc"),
                        "6" = c("scenario", "N1", "N2", "N3", "N4", "raa", "ta", "t1", "N5", "t2", "Nanc"))
  return(param_order)
}

# # function to reorder parameters
# reorder_param_table <- function(params_table, scenario) {
#   new_params <- params_table[, match(get_parameter_order(paste0(scenario)), colnames(params_table))]
#   return(new_params)
# }

# function to reorder parameters (si pas le bon ordre de parametres dans headerRF.txt)
reorder_param_table <- function(params_table, scenario) {
  #new_params <- params_table[, match(get_parameter_order(paste0(scenario)), colnames(params_table))]
  #return(new_params)
  return(params_table)
}

setup_sim <- function(path_to_headers, datestamp, seed, ncores) {
  cwd <- getwd()
  on.exit(setwd(cwd))
  path_to_sim <- file.path(cwd,"sim")
  if (!dir.exists(path_to_sim)) {
    dir.create(path_to_sim)
    setwd(path_to_sim)
    # copy headers and maf
    files_to_copy <- list.files(path_to_headers, full.names = TRUE)
    file.copy(files_to_copy, path_to_sim)
    # Generate RNG
    system(paste0(DIYABC_CMD, " -n \"t:", ncores, ";c:1;s:", seed, ";f:f\""))

  }
  return(path_to_sim)
}

simulate_swd <- function(params_table, scenario, path_to_sim, ncores) {
  # re-order table
  new_params_table <- reorder_param_table(params_table, scenario)
  # tmpext <- sub("/", "", tempfile("", "", ""))
  # newparamname <- paste0("reordered_param_human_", tmpext)
  newparamname <- "reordered_param_swd.txt"
  write.table(file = here(path_to_sim, newparamname), new_params_table, quote = F, col.names = T, row.names = F)
  #write.table(file = paste0(path_to_sim,"/", newparamname), new_params_table, quote = F, col.names = T, row.names = F)

  # do sim
  cwd <- getwd()
  on.exit(setwd(cwd))
  setwd(path_to_sim)
  # outputsumstats <- paste0("generated_sumstats_human.", tmpext)
  outputsumstats <- "generated_sumstats_swd.txt"
  system(paste0(DIYABC_CMD, " -o ", newparamname, " -i ", outputsumstats, " -g 100 -m -t ", ncores))
  # read result
  new_sumstats <- read.table(outputsumstats, header = FALSE)
  # delete temp files
  #unlink(c(newparamname, outputsumstats))
  return(new_sumstats)
}







# #######################################################################################################################
# ## Test
# #######################################################################################################################
# library(here)
#
# # Files names
# inputfile <- here("data", "dummy", "test_human_params_scenario_2.txt")
# path_to_headers <- here("data", "2024-03-18_human_resim")
# seed <- 1289
# ncores <- 4
#
# # setup simulations
# path_to_sim <- setup_sim(path_to_headers, seed, ncores)
#
# # parameters
# scenario <- 2
# params_table <- read.table(inputfile, header = TRUE)
#
# # sim
# new_sumstats <- simulate_human(params_table, scenario, path_to_sim, ncores)
#
# new_sumstats
