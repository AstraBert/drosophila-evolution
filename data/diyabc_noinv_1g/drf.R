# AE: DRF and RF  nalysis - DESIGN FOR FOR GHOST POP ANALYSIS 01-2025 from version 05-02-2023 - with RMSE + HDP mono et biplot + energy score
# Key paper for DRF: Distributional Random Forest; D Cevid, L Michel, J Näf, P Bühlmann, N Meinshausen 2022 Journal of Machine Learning Research 23 (333), 1-79


#### Initial cleaning and general options #################
# Clean all objects in the global environment
rm(list = ls())
# Free memory used by deleted objects
gc()
# Clean all open graphics (if using RStudio)
if (!is.null(dev.list())) dev.off()
# Reset options
options(default = TRUE)
#options(prompt = FALSE)
options(ask = FALSE)
# Clean history (optional)
cat("", file = ".Rhistory")

######## Libraries ##########################
library(doParallel)
library(drf)
library(abcrf)
help(package = "drf")
help(package = "abcrf")

#### Various options
options(max.print=10000) # To print tables until 10000 lines
# Number of cores available for your computer
n_cores <- detectCores()
# Nbre of CPU cores I want to use for parallel computation (if not define then all available CPU cores will be used)
nbre.threads.used.for.computation = n_cores-2

# Saving results
save <- TRUE  # Set save to TRUE or FALSE
if (save == TRUE) {
  output_file_name = "DRF_vs_RF_astra_EstimParam_ta_raa_t1000_sim20000_bis.tx"
  #output_file_name = "test1.txt"
}

# Key file names
# File names
file.reference.sim = "reftableRF.bin"
file.header = "headerRF.txt"
name.reference.table <- "reftableRF.bin"
name.header.file <- "headerRF.txt"
#name.observed.dataset <- "statobsRF_all_statobs_vectors.txt" # Note that this observed dataset may includes several lines (i.e. several vectors of observed data/sumstats) 
#name.observed.dataset <- "statobsRF_target_dataset.mss.dataNANNE.dat_US-Wat_US-Sok_US-Haw_JP-Sap_CN-Lan_BR-Poa.txt"
name.observed.dataset <- "statobsRF.txt"

# Selected model, dimension of n.pods & DRF - RF running parameters
selected.model=1
n.pods=1 ## If statobs contain a single row of statobs then n.pods=1 
# n.tree and n.train configurations to be explored
# vector.n.train=c(500,2000,4000,10000,20000,50000)
# vector.n.tree <- list(1000, 2000)
# vector.n.train <- list(10000)
# vector.n.tree <- list(500, 1000, 1500, 2000, 3000, 5000)
vector.n.train <- list(11000)
vector.n.tree <- list(2000)
# N.load.reftable = total number of simuations one wants to load from the reference table
N.load.reftable = 12000
# N.train = 2000 # A affiner

# Computation options
DRF=TRUE
RF=FALSE
param_original = TRUE
param_compound = FALSE
CORRELATIONS_COMPILATION_GRAPH=TRUE
addoc_graphs = FALSE
PRECISION.METRICS.FOR.ind.PODS.COMPUTATIONS=FALSE
PRECISION.METRICS.FOR.MEAN.COMPUTATIONS=FALSE
ANNEX.COMPUTATION.FOR.DRF.RF.EVALUATION.ala.NAF=FALSE
POST_TREATMENTS_HPD_MONO_plus_BI_PLOT_and_more =FALSE
PLS=FALSE

#### Reading pods for fixed-values parameter sets
# # pods.fixed.parameter.values <- rep(c(7000, 2000, 1000, 1000, 3000, 4000, 0.3), time = 10000) # set1
# pods.fixed.parameter.values <- rep(c(1000, 500, 2000, 300, 3000, 5000, 0.2), time = 10000)  # set2
# #file.name.fixed.pods = "reftableRF - pods_set1.bin"
# file.name.fixed.pods = "reftableRF - pods_set2.bin"
# param.pods <- matrix(pods.fixed.parameter.values, nrow = 10000, ncol = 7, byrow = TRUE)
# param.pods <- as.data.frame(param.pods)
# colnames(param.pods) <- c("N12", "N3", "N4", "t124", "t23", "t12", "ra")
# param.pods = param.pods[1:n.pods,]
# #dim(param.pods)
# #head(param.pods)
# data.pods <- readRefTable(filename = file.name.fixed.pods, header=file.header, N=n.pods)
# stats.pods = as.data.frame(data.pods$stats)
# # dim(stats.pods)
# # head(stats.pods,n=5)
# #stats.pods = read.table("pods_param_set1_first_records_of_the_reference_table_0.txt", header=TRUE, skip=7)
# #stats.pods <- stat.pods[1:n.pods,-1]
# #dim(stats.pods)
# #head(stats.pods,n=5)

# ########## Reading pods data from the reftable (initial code)###############
# # data.pods <- readRefTable(filename = file.reference.sim, header=file.header, N=n.pods)
# # Loading the statobs file
# # Check the presence of the element "fin" in the file name.observed.dataset
# if (file.exists(name.observed.dataset)) {
#   # Read all lines from the file
#   lines <- readLines(name.observed.dataset)
#   # Check if "fin" is present in the file
#   fin_index <- which(grepl("fin", lines))
#   if (length(fin_index) > 0) {
#     # Skip all lines before and including the line containing "fin"
#     data_start <- fin_index[length(fin_index)] + 1
#     # Read the following lines with header = TRUE
#     stat.obs <- read.table(text = paste(lines[data_start:length(lines)], collapse = "\n"), header = TRUE)
#     #stat.obs <- read.table(text = paste(lines[data_start:length(lines)], collapse = "\n"), header = TRUE,  nrows = 1)
#   } else {
#     # If "fin" is not found, read the entire file with header = TRUE
#     stat.obs <- read.table(name.observed.dataset, header = TRUE)
#     #stat.obs <- read.table(name.observed.dataset, header = TRUE, nrows = 1)
#   }
# } else {
#   stop(paste("The file", name.observed.dataset, "does not exist."))
# }

# Loading statObs and reference.table files
# N=N.load.reftable = total number of simuations one wants to load from the reference table
reference.table <- readRefTable(filename = name.reference.table, header=name.header.file, N=N.load.reftable)
indexesModel <- which(reference.table$scenarios == selected.model) #### Choix du scenario pour lequel on estime le parametre
# Loading the statobs file
# Check the presence of the element "fin" in the file name.observed.dataset
if (file.exists(name.observed.dataset)) {
  # Read all lines from the file
  lines <- readLines(name.observed.dataset)
  # Check if "fin" is present in the file
  fin_index <- which(grepl("fin", lines))
  if (length(fin_index) > 0) {
    # Skip all lines before and including the line containing "fin"
    data_start <- fin_index[length(fin_index)] + 1
    # Read the following lines with header = TRUE
    stat.obs <- read.table(text = paste(lines[data_start:length(lines)], collapse = "\n"), header = TRUE)
    #stat.obs <- read.table(text = paste(lines[data_start:length(lines)], collapse = "\n"), header = TRUE,  nrows = 1)
  } else {
    # If "fin" is not found, read the entire file with header = TRUE
    stat.obs <- read.table(name.observed.dataset, header = TRUE)
    #stat.obs <- read.table(name.observed.dataset, header = TRUE, nrows = 1)
  }
} else {
  stop(paste("The file", name.observed.dataset, "does not exist."))
}

### Pre-processing and sectioning of data = indexing the addoc number of data, focussing on the selected model, and cleaning NA's parameter colunms
reference.table$scenarios <- reference.table$scenarios[indexesModel]
reference.table$stats <- reference.table$stats[indexesModel,]
reference.table$params <- reference.table$params[indexesModel,]
# WARNING = il faut virer les param avec des NA (param des autres scenarios que celui selected)
reference.table$params <- reference.table$params[, colSums(is.na(reference.table$params)) == 0] # Virer les param avec des NA (autre scenario que celui selected)
colnames(reference.table$params)

# # Complements (optional)
# # Combien on a de données au final
# N <- length(reference.table$scenarios) 
# # Combien de stats
# N.stats <- ncol(reference.table$stats)
# # Combien de param
# N.param <- ncol(reference.table$params)
# # Nombre de données potentielles pour faire des tests
# N.test_possible <- N - N.train # datatests = potentiellement toutes les donnees n'ayant pas servis pour entrainement
# 
# # Randomization of data after selecting those of a given scenario - set.seed(1) # If one want some controle
# indicesTrain <- sample(1:N, N.train, replace=FALSE) # indices d'entrainement
# indicesTest <- c(1:N)[-indicesTrain] # indices de test

########## Connecting the INITIAL SCRIPT AND NEW SCRIPT training set data and StatObs data ###############
PARAMS.FULL <- reference.table$params 
#dim(PARAMS.FULL)
STATS.FULL <- reference.table$stats
#dim(STATS.FULL)
N.rec <- nrow(PARAMS.FULL) # the number simulations
N.rec
colnames(PARAMS.FULL)
#head(PARAMS.FULL, n = 10)
#head(STATS.FULL, n = 10)
stats.pods = as.data.frame(stat.obs)
###########################################

# ########## Reading pods data from the reftable ###############
data.pods <- readRefTable(filename = file.reference.sim, header=file.header, N=n.pods)
# dim(data.pods$stats)
# dim(data.pods$params)
# head(data.pods$params, n = 10)
# head(data.pods$stats, n = 10)
param.pods = as.data.frame(data.pods$params)
# # File names
# file.reference.sim = "reftableRF.bin"
# file.header = "headerRF.txt"
# stats.pods = as.data.frame(data.pods$stats)
# # File names
# file.reference.sim = "reftableRF.bin"
# file.header = "headerRF.txt"

# ########## Reading training set data ###############
# reftable <- readRefTable(filename = file.reference.sim, header=file.header)
# PARAMS.FULL <- as.data.frame(reftable$params[, ]) 
# #dim(PARAMS.FULL)
# STATS.FULL <- as.data.frame(reftable$stats[, ])
# #dim(STATS.FULL)
# N.rec <- reftable$nrec # the number simulations
# #N.rec
# #colnames(PARAMS.FULL)
# #head(PARAMS.FULL, n = 10)
# #head(STATS.FULL, n = 10)

################################################################################
### Compound parameters  #######################################################
################################################################################

# list_compound_parameters <- "None" # Compound parameters for admixture rates
# if (param_compound == TRUE) {
#   param1 <- PARAMS.FULL[, "raaas1"]
#   param2 <- PARAMS.FULL[, "raaas2"]
#   # Définir les objets composés
#   ra.Wat.SA <- as.data.frame(1 - param1)
#   ra.Gan2.SA <- as.data.frame(param2 * param1)
#   ra.nc.SA <- as.data.frame((1 - param2) * param1)
#   ra.watGan2.SA <- as.data.frame(ra.Wat.SA + ra.Gan2.SA)
#   
#   # Liste des paramètres composés
#   colnames(ra.Wat.SA)="ra.Wat.SA"
#   colnames(ra.Gan2.SA)="ra.Gan2.SA"
#   colnames(ra.nc.SA)="ra.nc.SA"
#   colnames(ra.watGan2.SA)="ra.watGan2.SA"
#   list_compound_parameters <- c("ra.Wat.SA", "ra.Gan2.SA", "ra.nc.SA", "ra.watGan2.SA")
#   
#   PARAMS.FULL = cbind(PARAMS.FULL, ra.Wat.SA, ra.Gan2.SA, ra.nc.SA, ra.watGan2.SA)
#   param.list = colnames(PARAMS.FULL)
#   
#   # dim(PARAMS.FULL)
#   # head(PARAMS.FULL, n=3)
# }

########## List of original Parameters to estimate
list_original_parameters = "NONE"
if (param_original==TRUE) {
  
  # Code addoc pour liste automatique a partir du reftable
  # Liste des colonnes initiales
  #colnames(reference.table$params)
  # Créer la liste à partir des noms de colonnes de reference.table$params
  list_original_parameters <- as.list(colnames(reference.table$params))
  # Combiner le contenu de la liste en un vecteur
  list_original_parameters <- unlist(list_original_parameters, use.names = FALSE)
  #print(list_original_parameters)
  
  #list_original_parameters = c("raan1","raan2", "raawat", "raasd", "raaas1", "raaas2")
  # list_original_parameters = c("DBc1", "NBc1", "DBc2", "NBc2", "DBc3", "NBc3", "DBarg", "NBarg",
  #                              "DBbra", "NBbra", "DBas1", "NBas1", "DBnc", "NBnc", "DBwis", "NBwis",
  #                              "DBgen", "NBgen", "DBcol", "NBcol", "DBan3", "NBan3", "DBsd", "NBsd",
  #                              "DBsok", "NBsok", "DBan2", "NBan2","DBan1", "NBan1", "DBh", "NBh", "DBGhw", "NBGhw")
  # list_original_parameters <- c(
  #   "Nwat", "Nsok", "Nhw", "Nsap", "Nlia", "Nsd", "Nnc", "Nwis", "Ngen",
  #   "Ncol", "Nbra", "Narg", "Ncl1", "Ncl2", "Ncl3", "NGan2", "NGan1", "NGhw", "NGan3", "NGas1", "NGas2",
  #   "NAC", "tc1", "tc2", "tc3", "DBc1", "NBc1", "DBc2", "NBc2", "DBc3", "NBc3", "targ", "DBarg", 
  #   "NBarg", "tbra", "DBbra", "NBbra", "tas1", "DBas1", "NBas1", "raaas1", "raaas2", "tnc", "DBnc",
  #   "NBnc", "twis", "DBwis", "NBwis", "tgen", "DBgen", "NBgen", "tcol", "DBcol", "NBcol", "tan3",
  #   "DBan3", "NBan3", "tsd", "DBsd", "NBsd", "raasd", "twat", "DBwat", "NBwat", "raawat", "tsok",
  #   "DBsok", "NBsok", "tan2", "DBan2", "NBan2", "raan2", "tan1", "DBan1", "NBan1", "raan1", "th",
  #   "DBh", "NBh", "tGhw", "DBGhw", "NBGhw", "tj", "tc", "µmic_1", "pmic_1", "snimic_1"
  # )
}

########## Compound parameters for Bottleneck Intensity #################
list_compound_parameters <- "None"
if (param_compound == TRUE) {
  list_selected_original_params_for_compound_params <- 
    c("DBc1", "NBc1", "DBc2", "NBc2", "DBc3", "NBc3", "DBarg", "NBarg",
    "DBbra", "NBbra", "DBas1", "NBas1", "DBnc", "NBnc", "DBwis", "NBwis",
    "DBgen", "NBgen", "DBcol", "NBcol", "DBan3", "NBan3", "DBsd", "NBsd",
    "DBsok", "NBsok", "DBan2", "NBan2","DBan1", "NBan1", "DBh", "NBh", "DBGhw", "NBGhw")
  
  # Nombre total de binômes (chaque binôme correspond à deux paramètres consécutifs)
  num_BI_binomes <- length(list_selected_original_params_for_compound_params) / 2
  # Initialiser un dataframe vide avec le bon nombre de colonnes et lignes
  BI_binomes_dataframe <- data.frame(matrix(nrow = nrow(PARAMS.FULL), ncol = num_BI_binomes))
  colnames(BI_binomes_dataframe) <- rep("", num_BI_binomes)  # Initialement vide pour les noms de colonnes
  # Générer dynamiquement les binômes et les remplir dans le dataframe
  col_index <- 1
  for (i in seq(1, length(list_selected_original_params_for_compound_params), by = 2)) {
    db_param <- list_selected_original_params_for_compound_params[i]
    nb_param <- list_selected_original_params_for_compound_params[i + 1]
    bi_name <- paste0("BI", sub("DB", "", db_param))
    
    # Vérifier si les colonnes existent dans PARAMS.FULL avant de procéder
    if (exists("PARAMS.FULL") && all(c(db_param, nb_param) %in% colnames(PARAMS.FULL))) {
      # Calculer le binôme
      binome_col <- PARAMS.FULL[, db_param] / PARAMS.FULL[, nb_param]
      # Ajouter au dataframe à l'index approprié
      BI_binomes_dataframe[, col_index] <- binome_col
      colnames(BI_binomes_dataframe)[col_index] <- bi_name  # Nommer la colonne
      col_index <- col_index + 1  # Passer à la colonne suivante
    } else {
      cat("Les colonnes", db_param, "et/ou", nb_param, "n'existent pas dans PARAMS.FULL\n")
    }
  }
  
  # Supprimer les colonnes vides, le cas échéant (si des colonnes manquent)
  BI_binomes_dataframe <- BI_binomes_dataframe[, colnames(BI_binomes_dataframe) != ""]
  # Afficher le dataframe final contenant tous les binômes
  # dim(BI_binomes_dataframe)
  # head(BI_binomes_dataframe, n=3)
  # colnames(BI_binomes_dataframe)
  
  list_compound_parameters <- colnames(BI_binomes_dataframe)
  PARAMS.FULL = cbind(PARAMS.FULL,BI_binomes_dataframe)
  # dim(PARAMS.FULL)
  # head(PARAMS.FULL, n=3)
  # colnames(PARAMS.FULL)
}

# Combiner les deux listes original.param et compound.param en une seule si necessaire
if ((param_original == TRUE)&(param_compound == TRUE)) param.list <- c(list_original_parameters, list_compound_parameters)
if ((param_original == TRUE)&(param_compound == FALSE)) param.list <- c(list_original_parameters)
if ((param_original == FALSE)&(param_compound == TRUE)) param.list <- c(list_compound_parameters)

#################################################################
# NEW CODE FOR A SINGLE STATOBS: Compile RF results
# Buiding result table for all parameters
# Liste combinée des paramètres (remplacez avec vos données réelles)
# Dimensions du dataframe
n_rows <- length(param.list)
col_names <- c("Parameter", "Mean", "Median", "Q5", "Q95", 
               "post.NMAE.median", "post.NMAE.mean", 
               "prior.NMAE.median", "prior.NMAE.mean", "Coverage_OOB")

# Initialisation du dataframe avec NA
RESULT_ESTIM_PARAM_RF <- data.frame(matrix(NA, nrow = n_rows, ncol = length(col_names)))
colnames(RESULT_ESTIM_PARAM_RF) <- col_names
# Remplir la colonne "Parameter" avec les noms des paramètres
RESULT_ESTIM_PARAM_RF$Parameter <- param.list
##################################################

# INITIAL CODE: MetaDataframe for each RF parameter
dim.param <- length(param.list)
param.pred.rf = vector("list", dim.param)

# INITIAL CODE: Creating RESULTS dataframe (for mean computation if large number of pods)
dim.n.train <- as.numeric(length(vector.n.train))
dim.n.tree <- as.numeric(length(vector.n.tree))
dim.n.train.n.tree <- as.numeric(dim.n.train * dim.n.tree)
RESULTS<- as.data.frame(matrix(0, nrow = dim.n.train.n.tree, ncol = 2+length(param.list)*6+4))
colnames(RESULTS)[1] = "n.tree"
colnames(RESULTS)[2] = "n.train"

# INITIAL CODE: Créer une liste RESULTS.list.n.pod composée de n.pod dataframes (to store all n.pods computation)
# Définir les dimensions
dim.n.train <- as.numeric(length(vector.n.train))
dim.n.tree <- as.numeric(length(vector.n.tree))
dim.n.train.n.tree <- as.numeric(dim.n.train * dim.n.tree)
# Créer un dataframe modèle
RESULTS.single.pod <- as.data.frame(matrix(0, nrow = dim.n.train.n.tree, ncol = 2 + length(param.list) * 6 + 4))
colnames(RESULTS.single.pod)[1] = "n.tree"
colnames(RESULTS.single.pod)[2] = "n.train"
# Créer une liste de dataframes
RESULTS.list.n.pods <- vector("list", n.pods)
# Remplir chaque élément de la liste avec une copie du dataframe modèle
for (i in 1:n.pods) RESULTS.list.n.pods[[i]] = RESULTS.single.pod

### Here is a bench of useful functions
calculate_NMAE_range <- function(true_values, point_estimates) {
  absolute_errors <- abs(true_values - point_estimates)
  normalized_errors <- absolute_errors / (max(true_values) - min(true_values))
  nmae.range <- mean(normalized_errors)
  return(nmae.range)
}
calculate_NMAE_mean <- function(true_values, point_estimates) {
  absolute_errors <- abs(true_values - point_estimates)
  normalized_errors <- absolute_errors / mean(true_values)
  nmae.mean <- mean(normalized_errors)
  return(nmae.mean)
}
calculate_dif.drf.rf.mean <- function(true_values, point_estimates.drf, point_estimates.rf) {
  dif.drf = abs(true_values - point_estimates.drf)*1/true_values
  dif.rf = abs(true_values - point_estimates.rf)*1/true_values
  dif.drf.rf = dif.drf - dif.rf
  return(dif.drf.rf)
}
calculate_dif.mean <- function(true_values, point_estimates) {
  dif = abs(true_values - point_estimates)*1/true_values
  return(dif)
}
calculate_dif.drf.rf.sd <- function(true_values, point_estimates.drf, point_estimates.rf) {
  dif.drf.rf = (point_estimates.drf - point_estimates.rf)*1/true_values
  return(dif.drf.rf)
}
calculate_90CI <- function(true_values, Q5_point_estimates, Q95_point_estimates) {
  # Calculer la proportion de valeurs dans l'intervalle [Q5, Q95]
  CI90 <- sum(true_values >= Q5_point_estimates & true_values <= Q95_point_estimates) / length(true_values)
  return(CI90)
}
calculate_dif.drf.rf.lengthCI90 <- function(true_values, Q5_point_estimates.drf, Q95_point_estimates.drf, Q5_point_estimates.rf, Q95_point_estimates.rf) {
  lengthCI90.drf = (Q95_point_estimates.drf - Q5_point_estimates.drf)/true_values
  lengthCI90.rf = (Q95_point_estimates.rf - Q5_point_estimates.rf)/true_values 
  dif.drf.rf.lengthCI90 = lengthCI90.drf - lengthCI90.rf
  return(dif.drf.rf.lengthCI90)
}

# Calculer le CI90 pour intersection des deux variables V1 et V2 
# c'est-à-dire la proportion de valeurs où les deux variables se trouvent simultanement dans leurs intervalles respectifs [Q5, Q95]
# On attends 0.9*0.9 = 0.81 ????????????
# calculate_90CI_intersection <- function(V1_true_values, V1_Q5, V1_Q95, V2_true_values, V2_Q5, V2_Q95) {
#   # Calculer la proportion de valeurs où les deux variables sont dans leurs intervalles respectifs
#   CI90_intersection <- sum(V1_true_values >= V1_Q5 & V1_true_values <= V1_Q95 & 
#                              V2_true_values >= V2_Q5 & V2_true_values <= V2_Q95) / length(V1_true_values)
#   return(CI90_intersection)
#}

sink(file = output_file_name, split = TRUE)
cat(" ############## PARAMETER ESTIMATION using DRF & RF ############### ","\n")
cat("\n")
cat("Name of the statobs file:", name.observed.dataset,"\n")
cat("Number of statobs in the statobs file =", nrow(stat.obs),"\n")
cat("Name of the reference table file:", name.reference.table,"\n")
cat("Name of the headerfile file:", name.header.file,"\n")
cat("Number of simulations loaded from the reference table =",reference.table$nrec, "\n")
cat("Number of scenarios (i.e. models) in the reference table =",reference.table$nscen, "\n")
cat("Number of simulations available for each scenario from the loaded reference table =",reference.table$nrecscen,"\n")
cat("Number of parameters in the reference table =",reference.table$nparam, "\n")
cat("Number of summary statistics in the reference table (without LDA axes) =",ncol(reference.table$stats), "\n")
cat("Selected model for parameter estimation = ",selected.model, "\n")
cat("Number of cores available = ",n_cores, "\n")
cat("Number of cores used for computation =",nbre.threads.used.for.computation, "\n")
if (param_original == TRUE) cat("Original parameters estimated = ",list_original_parameters, "\n")
if (param_compound == TRUE) cat("Compound parameters estimated = ",list_compound_parameters, "\n")
cat("Parameters in fine estimated = ",param.list, "\n")
cat("DRF analysis = ",DRF, "\n")
cat("RF analysis = ",RF, "\n")
cat("\n")

###########################################################################################################
####################### STARTING LOOP OVER n.tree and n.train   ###########################################
###########################################################################################################
k=0
for (i in 1:dim.n.tree) {
  for (j in 1:dim.n.train) {
    n.tree = as.numeric(vector.n.tree[i])
    n.train = as.numeric(vector.n.train[j])
    k = k+1
    cat("\n")
    cat("ANALYSIS with n.tree =",n.tree , "n.train =",n.train ,"\n")
    cat(" ","\n")
    PARAMS <- PARAMS.FULL[(n.pods+1):(n.pods+n.train), ]
    STATS <- STATS.FULL[(n.pods+1):(n.pods+n.train), ]
    # dim(STATS)
    # dim(PARAMS)
    # dim(stats.pods)
    # str(STATS)
    # str(PARAMS)
    # str(stats.pods)
    # head(STATS, n = 2)
    # head(PARAMS, n = 2)
    # head(stats.pods, n=2)

    ####### RF and DRF computation ###################################################################
    if (RF==TRUE) {
    ####### RF computation ###################################################################
    param.nbr=0
    for (i.param in 1:length(param.list)) {
      param.nbr = param.nbr+1
      param.name <- param.list[i.param]
      cat("\n") 
      cat("RF analysis - parameter ",param.name, "\n")
      cat("Parameter #",param.nbr," over a total of",length(param.list), "parameters to estimate", "\n")
      target.param = as.data.frame(PARAMS[,param.name])
      colnames(target.param) = "target"
      #dim(target.param)
      #head(target.param, n=5)
      data.RF <- data.frame(target.param, STATS)
      stats.pods = as.data.frame(stats.pods)
      #dim(data.RF)
      #head(data.RF, N=5)
      model.rf <- regAbcrf(target~., data.RF, ntree = n.tree, min.node.size = 5,
                    paral=TRUE, ncores=nbre.threads.used.for.computation)
      pred.rf <- predict(model.rf, stats.pods, data.RF, quantiles = c(0.05,0.95), 
                         prior.err.med = TRUE, post.err.med = TRUE, paral=TRUE, 
                         ncores=nbre.threads.used.for.computation)
      param.pred.rf[[i.param]]=pred.rf
      
      # Extraire les informations de pref.rf for storage
      output <- capture.output(pred.rf) # Capture les sorties affichées en console
      # Identifier et extraire les lignes contenant des informations particulières
      NMAE_prior_mean <- as.numeric(sub(".*mean: ", "", grep("Prior out-of-bag normalized mean absolute error computed with mean", output, value = TRUE)))
      NMAE_prior_median <- as.numeric(sub(".*median: ", "", grep("Prior out-of-bag normalized mean absolute error computed with median", output, value = TRUE)))
      prior_OOB_coverage <- as.numeric(sub(".*coverage: ", "", grep("Prior out-of-bag credible interval coverage", output, value = TRUE)))
      # Verify
      # cat("Prior out-of-bag normalized mean absolute error (mean):", NMAE_prior_mean, "\n")
      # cat("Prior out-of-bag normalized mean absolute error (median):", NMAE_prior_median, "\n")
      # cat("Prior out-of-bag credible interval coverage:", prior_OOB_coverage, "\n")
      
      pred.rf.key.results <- c(pred.rf$expectation,
                             pred.rf$med,
                             pred.rf$quantiles[1],
                             pred.rf$quantiles[2],
                             pred.rf$post.NMAE.med,
                             pred.rf$post.NMAE.mean,
                             NMAE_prior_median,
                             NMAE_prior_mean, 
                             prior_OOB_coverage)
      ### Dynamic Building of result Table with the RF key results for each parameter
      RESULT_ESTIM_PARAM_RF[param.nbr,2:ncol(RESULT_ESTIM_PARAM_RF)] = pred.rf.key.results
      print(RESULT_ESTIM_PARAM_RF[param.nbr,])
      cat("\n") 
     }
    }
    ##########################################################################
    
    if (DRF==TRUE) {
    ####### DRF computation ###################################################################
    # WARNING:  Pour DRF on travaille avec des matrices (et pas des dzataframe avec RF) !!!
    #  La structure de stats.pods (dataframe) ne correspond pas à celle de STATS (matrice)...et donc
    cat(" ", "\n")
    cat("DRF analysis", "\n")
    stats.pods <- as.matrix(stats.pods)
    stats.pods <- stats.pods[, colnames(STATS), drop = FALSE] ### Verif que colonnes stats.pods et STATS sont dans le meme ordre
    
    # model.drf <- drf(X = STATS, Y = PARAMS, num.trees = n.tree, honesty = TRUE, min.node.size = 5, 
    #                   sample.fraction = 0.5, ci.group.size = 2, num.threads =nbre.threads.used.for.computation, splitting.rule = "FourierMMD")
    model.drf <- drf(X = STATS, Y = PARAMS, num.trees = n.tree, honesty = FALSE, min.node.size = 1,
                     sample.fraction = 2/3, ci.group.size = 1, num.threads =nbre.threads.used.for.computation, splitting.rule = "FourierMMD")
    
    # weights
    W.drf <- predict(model.drf, newdata = stats.pods, ntree = n.tree, num.threads = nbre.threads.used.for.computation)
    # mean
    pred.drf.mean <- predict(model.drf, newdata = stats.pods, ntree = n.tree, num.threads = nbre.threads.used.for.computation, functional = "mean")
    # median
    pred.drf.median <- predict(model.drf, newdata = stats.pods, num.threads =nbre.threads.used.for.computation, functional = "quantile", quantiles = c(0.5))
    # quantiles
    pred.drf.quant <- predict(model.drf, newdata = stats.pods,ntree = n.tree, num.threads =nbre.threads.used.for.computation, functional = "quantile", quantiles = c(0.05, 0.95))
    # sd
    pred.drf.sd <- predict(model.drf, newdata = stats.pods, ntree = n.tree,num.threads =nbre.threads.used.for.computation, functional = "sd")
    # cor
    pred.drf.cor <- predict(model.drf, newdata = stats.pods,ntree = n.tree, num.threads =nbre.threads.used.for.computation, functional = "cor")
    
    pred.drf.mean = as.data.frame(pred.drf.mean)
    pred.drf.median = as.data.frame(pred.drf.median)
    pred.drf.quant = as.data.frame(pred.drf.quant)
    pred.drf.sd = as.data.frame(pred.drf.sd)
    pred.drf.cor = as.data.frame(pred.drf.cor)
    # dim(pred.drf.mean)
    # head(pred.drf.mean, n=10)
    # dim(pred.drf.quant)
    # head(pred.drf.quant, n=10)
    # dim(pred.drf.sd)
    # head(pred.drf.sd, n=10)
    # dim(pred.drf.cor)
    # head(pred.drf.cor, n=10)
    
    ######## DRF RESULTS POUR REAL STATOBS #############
    # Convertir explicitement les objets en data frames...si nécessaire
    results.drf.mean <- as.data.frame(pred.drf.mean)
    results.drf.median <- as.data.frame(pred.drf.median)
    results.drf.quant <- as.data.frame(pred.drf.quant)
    results.drf.Q5 <- as.data.frame(pred.drf.quant[, 1:ncol(PARAMS)])
    results.drf.Q95 <- as.data.frame(pred.drf.quant[, (ncol(PARAMS) + 1):(2 * ncol(PARAMS))])
    # Appliquer les noms de colonnes
    colnames(results.drf.mean) <- colnames(PARAMS)
    colnames(results.drf.median) <- colnames(PARAMS)
    colnames(results.drf.Q5) <- colnames(PARAMS)
    colnames(results.drf.Q95) <- colnames(PARAMS)
    results.drf.global = rbind(results.drf.mean, results.drf.median, results.drf.Q5, results.drf.Q95)
    row.names(results.drf.global) = c("Mean", "Median", "Q5", "Q95"  )
    #dim(results.drf.global)
    results.drf.global.inverted = t(results.drf.global)
    target.list = c(list_original_parameters,list_compound_parameters)
    # Filtrer les lignes dont les noms sont présents dans target.list
    results.drf.global.inverted.target.parameters <- results.drf.global.inverted[rownames(results.drf.global.inverted) %in% target.list, ]
    ##########################################################################
    }
    
    ###### AFFICHAGE / ECRITURE RESULTS
    cat("\n")
    cat("######################################################################################################################","\n")
    cat("################################### ALL PARAMETER ESTIMATION RESULTS (NO NOISE - NO PLS)         ####################","\n")
    cat("######################################################################################################################","\n")
    if (DRF==TRUE) {
      cat("\n")
      cat("DRF Estimation","\n")
      cat("\n")
      print(results.drf.global.inverted.target.parameters, row.names=FALSE)
      cat("\n")
    }
    if (RF==TRUE) {
    cat("\n")
    cat("RF Estimation","\n")
    cat("\n")
    print(RESULT_ESTIM_PARAM_RF, row.names=FALSE)
    cat("\n")
    }
    cat("######################################################################################################################","\n")
    cat("######################################################################################################################","\n")
    
    

    if (CORRELATIONS_COMPILATION_GRAPH==TRUE){
    
    # CORRELATION GRAPHS
    # Convertir toutes les colonnes en un vecteur
    correlation_values <- unlist(pred.drf.cor, use.names = FALSE)
    
    # Convertir pred.drf.cor en une matrice 91x91
    correlation_matrix <- matrix(unlist(pred.drf.cor, use.names = FALSE), nrow = length(param.list), ncol = length(param.list), byrow = TRUE)
    dim(correlation_matrix)
    # Supprimer les valeurs de la diagonale
    correlation_values <- correlation_matrix[lower.tri(correlation_matrix) | upper.tri(correlation_matrix)]
    length(correlation_values)  # Vérifier le nombre de valeurs restantes
    # Supprimer les valeurs proches de -1.0000
    tolerance <- .Machine$double.eps^0.5  # Tolérance pour la précision numérique
    correlation_values <- correlation_values[abs(correlation_values - (-1)) > tolerance]
    # Sauvegarder le graphique dans un fichier PNG
    png(filename = paste0(output_file_name, "_correlation.graph.png"), width = 800, height = 600)
    # Créer le graphique des correlations
    plot(correlation_values, 
         type = "p", 
         main = "Correlation coefficient for all pairs of parameters", 
         xlab = "Index of the pairs of parameters", 
         ylab = "Correlation", 
         col = "blue", 
         pch = 1,  # Rond vide
         ylim = c(-1, 1),
         cex.axis = 0.6,  # Réduire la taille des caractères des axes
         las = 1)  # Orientation des étiquettes des axes à l'horizontale
    # Ajouter une échelle des y avec des pas de 0.1
    axis(2, at = seq(-1, 1, by = 0.1), cex.axis = 0.6, las = 1)  # Taille réduite pour l'axe Y
    # Ajouter une ligne horizontale en pointillé à y = 0
    abline(h = 0, lty = 2, col = "red")  # lty = 2 pour pointillé, col = "red" pour visibilité
    # Fermer le fichier PNG
    dev.off()
    
    # Extraire les paires de parametres avec de fortes correlations
    
    # Assigner les noms des lignes et des colonnes à correlation_matrix
    rownames(correlation_matrix) <- colnames(PARAMS.FULL)
    colnames(correlation_matrix) <- colnames(PARAMS.FULL)
    # Trouver les indices des valeurs > 0.5 ou < -0.5
    selected_indices <- which(correlation_matrix > 0.5 | correlation_matrix < -0.5, arr.ind = TRUE)
    # Extraire les labels des lignes et colonnes correspondants
    selected_values <- data.frame(
      Row = rownames(correlation_matrix)[selected_indices[, 1]],
      Column = colnames(correlation_matrix)[selected_indices[, 2]],
      Value = correlation_matrix[selected_indices]
    )
    # Éliminer les lignes où Row et Column sont identiques
    selected_values_filtered <- selected_values[selected_values$Row != selected_values$Column, ]
    # Utiliser un seuil de tolérance pour vérifier les valeurs égales à -1
    selected_values_filtered <- selected_values_filtered[abs(selected_values_filtered$Value - (-1)) > .Machine$double.eps, ]
    # Trier par ordre décroissant de Value
    selected_sorted_corr_values_filtered <- selected_values_filtered[order(selected_values_filtered$Value, decreasing = TRUE), ]
    cat("\n")
    cat(" ############## DRF CORRELATIONS PAIRS OF PARAMETERS OF INTEREST ############### ","\n")
    print(selected_sorted_corr_values_filtered, row.names=FALSE)
  }
    
    sink()
    
  }  
}  
  ###########################################################################################################
  ####################### ENDING LOOP OVER n.tree AND n.train ###############################################
  ###########################################################################################################   
    
if (addoc_graphs == TRUE) {
# ####### Graphiques addoc
# 
# Bottleneck Intensity (sorted on Mean values)

# Charger la bibliothèque ggplot2
library(ggplot2)

# Données au format texte
# data_text <- "BIcl1     0.2974968 1.442886e-01 1.304348e-02    0.7507375
# BIcl2     0.4147124 1.838073e-01 2.073227e-02    1.2741935
# BIcl3     0.2994772 1.616022e-01 1.294964e-02    0.8446602
# BIarg    0.1429912 8.097051e-02 9.300095e-03    0.4281833
# BIbra    0.3681970 1.057674e-01 1.181102e-02    0.9886755
# BIGas1    0.2565786 1.033464e-01 1.078932e-02    0.8689143
# BInca     0.3041404 9.307212e-02 1.021863e-02    0.7260080
# BIwis    0.2692379 9.555189e-02 1.019620e-02    0.8378378
# BIgen    0.1552850 9.356822e-02 1.025163e-02    0.4848485
# BIcol    0.1370445 7.746895e-02 7.874546e-03    0.2930789
# BIGan3    0.2711061 9.865208e-02 1.011810e-02    0.6283962
# BIsdi     0.2538763 8.995233e-02 9.324009e-03    0.5986175
# BIsok    0.5693245 1.931705e-01 2.288335e-02    1.7101449
# BIGan1    0.6696274 2.071130e-01 2.484472e-02    2.3085271
# BIhaw     0.9289952 3.342697e-01 5.225134e-02    2.9750000
#BIGh    0.6758741 2.051357e-01 2.240732e-02    2.2435897"


data_text <- "BIcl1     0.2855232 1.681416e-01 1.554404e-02    0.8768473
BIcl2     0.3392199 2.038297e-01 3.701800e-02    0.9816268
BIcl3     0.1981089 1.347926e-01 1.095290e-02    0.5601841
BIarg    0.1773621 8.399444e-02 8.056395e-03    0.5642857
BIbra    0.2990730 9.870142e-02 9.461426e-03    0.8291019
BIGas1   0.2242034 9.443867e-02 1.139321e-02    0.6682238
BInca    0.2169323 8.247423e-02 6.642255e-03    0.5011251
BIwis    0.3482098 1.028440e-01 1.239669e-02    0.8152174
BIgen    0.2377231 9.141561e-02 9.909363e-03    0.6511476
BIcol    0.2507043 8.295001e-02 7.042254e-03    0.4928367
BIGan3   0.2407510 9.746835e-02 9.922861e-03    0.6203008
BIsdi    0.1808891 9.176606e-02 1.073734e-02    0.5266693
BIsok    0.2982627 1.636712e-01 1.363229e-02    0.7863018
BIGan2   0.8775933 2.854080e-01 3.348730e-02    3.2413793
BIGan1    0.4140476 1.870187e-01 1.695323e-02    1.2488889
BIhaw    0.6585131 2.393986e-01 3.825137e-02    1.9147050
BIGh    0.4880664 2.003082e-01 2.443177e-02    1.3617986"


# Lire les données
data <- read.table(text = data_text, header = FALSE, col.names = c("Population", "Mean", "Median", "Q5", "Q95"), row.names = 1)
# Trier les données par la colonne Median en ordre croissant
data_sorted <- data[order(data$Median), ]
# Extraire les lettres après les deux premières lettres des noms de population
data_sorted$PopulationShort <- substr(row.names(data_sorted), 3, nchar(row.names(data_sorted)))

######### Créer le graphique avec ggplot2  - version 1 sans limites sur les axes
ggplot(data_sorted, aes(x = reorder(PopulationShort, Median))) +
  # Points pour Mean et Median
  geom_point(aes(y = Mean, color = "Mean"), size = 3, shape = 19) +
  geom_point(aes(y = Median, color = "Median"), size = 3, shape = 17) +
  # Bandes pour les intervalles de crédibilité
  geom_errorbar(aes(ymin = Q5, ymax = Q95), width = 0.2, color = "black") +
  # Personnalisation des axes et des légendes
  labs(
    x = "Population",
    y = "Bottleneck Intensity (BI)",
    color = "Legend"
  ) +
  # Personnalisation des couleurs
  scale_color_manual(values = c("Mean" = "blue", "Median" = "red")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

######### Créer le graphique avec ggplot2  - version 2 AVEC limites sur les axes
ggplot(data_sorted, aes(x = reorder(PopulationShort, Median))) +
  # Points pour Mean et Median
  geom_point(aes(y = Mean, color = "Mean"), size = 3, shape = 19) +
  geom_point(aes(y = Median, color = "Median"), size = 3, shape = 17) +
  # Bandes pour les intervalles de crédibilité
  geom_errorbar(aes(ymin = Q5, ymax = Q95), width = 0.2, color = "black") +
  # Personnalisation des axes et des légendes
  labs(
    x = "Population",
    y = "Bottleneck Intensity (BI)",
    color = "Legend"
  ) +
  # Personnalisation des couleurs
  scale_color_manual(values = c("Mean" = "blue", "Median" = "red")) +
  # Limiter l'axe des y
  coord_cartesian(ylim = c(0, 1)) +
  # Appliquer un thème minimal
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

######### Créer le graphique avec ggplot2  - version 3 AVEC limites sur les axes et ssi mediane
# Créer le graphique avec ggplot2
ggplot(data_sorted, aes(x = reorder(PopulationShort, Median))) +
  # Points pour Median
  geom_point(aes(y = Median, color = "Median"), size = 3, shape = 17) +
  # Personnalisation des axes et des légendes
  labs(
    x = "Population",
    y = "Bottleneck Intensity (BI)",
    color = "Legend"
  ) +
  # Personnalisation des couleurs
  scale_color_manual(values = c("Median" = "red")) +
  # Limiter l'axe des y avec des graduations spécifiques
  scale_y_continuous(limits = c(0, 0.3), breaks = seq(0, 0.3, by = 0.05)) +
  # Appliquer un thème minimal
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


}




    
  if (POST_TREATMENTS_HPD_MONO_plus_BI_PLOT_and_more ==TRUE){
    
    ####################################################################################################
    ############## START PRECISION METRICS COMPUTATION WITHOUT KEEPING iPODS VALUES  ##################
    ####################################################################################################
    RESULTS[k,1] <- n.tree
    RESULTS[k,2] <- n.train
    col.target = 2
    
    #############################################################
    ############ HPD mono-param computation for drf and rf #################
    #############################################################
    library(dplyr)
    library(coda)
    sample_size <- 100
    
    # Objets utiles pour HPD drf et rf
    pod.weighted.sample.rf = vector("list", n.pods)
    pod.weighted.sample.drf = vector("list", n.pods)
    pod.weighted.sample.rf<- lapply(pod.weighted.sample.rf, function(x) NULL) # Vider la liste
    pod.weighted.sample.drf<- lapply(pod.weighted.sample.drf, function(x) NULL) # Vider la liste
    # Création d'une liste avec replicate contenant hpd_results
    hpd_result = matrix(0, nrow = length(param.list)+1, ncol = 2)
    HPD90drf.lower.upper <- replicate(n.pods, hpd_result, simplify = FALSE)
    # HPD90drf.lower.upper <- lapply(1:n.pods, function(x) hpd_result) # Ca marche aussi pour creer la meme liste
    HPD90drf.lower.upper<- lapply(pod.weighted.sample.drf, function(x) NULL) # Vider la liste
    # rf = liste de liste: Initialiser la liste principale
    HPD90rf.lower.upper <- vector("list", length(param.list))
    # Remplir la liste principale avec les sous-listes
    for (i in 1:length(param.list)) {HPD90rf.lower.upper[[i]] <- vector("list", n.pods)}
    # On a alors une liste de liste que l'on manipule via HPD90rf.lower.upper[[i]][[j]]
    # str(HPD90rf.lower.upper)
    # Créer une matrice avec des valeurs initiales = zéro + Convertir la matrice en dataframe
    initial_matrix <- matrix(0, nrow = n.pods, ncol = length(param.list))
    HPD90drf.lower.for.each.param.all.pods <- as.data.frame(initial_matrix)
    HPD90drf.upper.for.each.param.all.pods <- as.data.frame(initial_matrix)
    colnames(HPD90drf.lower.for.each.param.all.pods) = colnames(PARAMS)
    colnames(HPD90drf.upper.for.each.param.all.pods) = colnames(PARAMS)
    HPD90rf.lower.for.each.param.all.pods <- as.data.frame(initial_matrix)
    HPD90rf.upper.for.each.param.all.pods <- as.data.frame(initial_matrix)
    colnames(HPD90rf.lower.for.each.param.all.pods) = colnames(PARAMS)
    colnames(HPD90rf.upper.for.each.param.all.pods) = colnames(PARAMS)
    
    # Weighted sampling for each pod and HPD computation
    for (i.pod in 1:n.pods) # for DRF
    {
      # drf: Extraction des weights correspondant a pod numéro i.pod
      weights.drf = as.vector(W.drf$weights[i.pod, ])
      #str(weights.drf)
      # sum(weights.drf)
      # mean(weights.drf)
      # drf: Création d'un data frame qui combine 'Ne', 'µ', 'Neµ' et les poids de i.pod
      param.W.drf = data.frame(PARAMS, weights = weights.drf)
      # head(param.W.drf)
      # dim(param.W.drf)
      pod.weighted.sample.drf[[i.pod]] <- param.W.drf %>%
        slice_sample(n = sample_size, replace = TRUE, weight_by = weights)
      # head(weighted.sample.drf, n=10)
      # dim(weighted.sample.drf)
      HPD90drf.lower.upper[[i.pod]] <- HPDinterval(as.mcmc(pod.weighted.sample.drf[[i.pod]]), prob = 0.9)
      # # Autre méthode donnant exactement les memes resultats que avec coda
      # #Librairie bayestestR : Fournit une fonction hdi() pour calculer l'intervalle de densité le plus élevé. 
      # # C'est une bonne alternative pour ceux qui travaillent dans un contexte bayésien et cherchent des fonctions 
      # # faciles à utiliser pour l'interprétation des résultats.
      # # Installer bayestestR si nécessaire
      # install.packages("bayestestR")
      # # Utiliser hdi() pour calculer l'intervalle HPD
      # library(bayestestR)
      # library(bayestestR)
      # result <- hdi(pod.weighted.sample.drf[[i.pod]], ci = 0.9)
      # # Affichage de CI pour les differents parametres
      # result$CI_low
      # result$CI_high
      # # Affichage de CI_low pour µ
      # result$CI_low[result$Parameter == "µ"]
      # result$CI_high[result$Parameter == "µ"]
      
    #   for (i.param in 1:length(param.list)) # For RF
    #   {
    #     selected.weights.rf <- as.data.frame(W.rf.param.pod[[i.param]][, i.pod])
    #     selected.trainingset.param.values = PARAMS[i.param]
    #     param.W.rf <- cbind(selected.trainingset.param.values, selected.weights.rf)
    #     colnames (param.W.rf) = c(colnames(PARAMS[i.param]),"weights")
    #     # head(param.W.rf)
    #     # dim(param.W.rf)
    #     pod.weighted.sample.rf[[i.param]][[i.pod]] <- param.W.rf %>%
    #       slice_sample(n = sample_size, replace = TRUE, weight_by = weights)
    #     # head(weighted.sample.rf, n=10)
    #     # dim(weighted.sample.rf)
    #     HPD90rf.lower.upper[[i.param]][[i.pod]] <- HPDinterval(as.mcmc(pod.weighted.sample.rf[[i.param]][[i.pod]]), prob = 0.9)
    #   }
    # }
    
    for (i.param in 1:length(param.list)) 
    {
      for (i.pod in 1:n.pods) 
      {
        HPD90drf.lower.for.each.param.all.pods[i.pod, i.param] = HPD90drf.lower.upper[[i.pod]][i.param, "lower"]
        HPD90drf.upper.for.each.param.all.pods[i.pod, i.param] = HPD90drf.lower.upper[[i.pod]][i.param, "upper"]
        # HPD90rf.lower.for.each.param.all.pods[i.pod, i.param] = HPD90rf.lower.upper[[i.param]][[i.pod]][1, "lower"]
        # HPD90rf.upper.for.each.param.all.pods[i.pod, i.param] = HPD90rf.lower.upper[[i.param]][[i.pod]][1, "upper"]
      }
    }
    HPD90drf.lower.upper.for.each.param.all.pods = rbind(HPD90drf.lower.for.each.param.all.pods,HPD90drf.upper.for.each.param.all.pods)
    t.HPD90drf.lower.upper.for.each.param.all.pods =t(HPD90drf.lower.upper.for.each.param.all.pods)
    # dim(HPD90drf.lower.for.each.param.all.pods)
    # head(HPD90drf.lower.for.each.param.all.pods)
    # dim(HPD90drf.upper.for.each.param.all.pods)
    # head(HPD90drf.upper.for.each.param.all.pods)
    # dim(t.HPD90drf.lower.upper.for.each.param.all.pods)
    # head(t.HPD90drf.lower.upper.for.each.param.all.pods)
    print(t.HPD90drf.lower.upper.for.each.param.all.pods)
    # dim(HPD90rf.lower.for.each.param.all.pods)
    # head(HPD90rf.lower.for.each.param.all.pods)
    # dim(HPD90rf.upper.for.each.param.all.pods)
    # head(HPD90rf.upper.for.each.param.all.pods)
    
################## FIN HPD90 mono param #################################
    
    ############################################################################### 
    ### GRAPH BIPLOT ASSOCIE ########## NEW 01-01-2025 ############################ Ajusted for raa and ta
    ############################################################################### 
    
    # Charger le package parallel
    library(parallel)
    # Définir le nombre de cœurs à utiliser
    n.cores <- detectCores() - 1  # Utiliser tous les cœurs sauf un
    # Créer le cluster
    cl <- makeCluster(n.cores)
    # Vérifier que le cluster est actif
    print(cl)
    
    
    
    # HPD levels correspond à votre conception (1-HPD) :
    # Ce que vous décrivez (HPD90 = 90% des valeurs estimées se trouvent à l'intérieur de la ligne HPD90) correspond 
    # effectivement à une probabilité cumulative. Dans hdr.2d, les niveaux HPD (falpha) fonctionnent différemment, 
    # indiquant les niveaux de densité associés. 
    # Il faut donc reformuler les contours pour refléter votre interprétation, en inversant l'ordre des niveaux.
    
    ############## NEW ############# DRF (adapted 03-01-2024 !!! for Ghost project !!!):  DBan2.NBan2 (Ghost USA2)
    colnames(pod.weighted.sample.drf[[i.pod]])
    
    # PARALLELISATION Methode drf-hdr-poly
    calcul_HPD_drf_hdr_poly <- function(i.pod) {
      hpd.levels <- hdrcde::hdr.2d(x = pod.weighted.sample.drf[[i.pod]]$ta, y = pod.weighted.sample.drf[[i.pod]]$raa, prob = c(0.1, 0.25, 0.5, 0.75, 0.9))
      drf.hdr <- contourLines(x = hpd.levels$den$x, y = hpd.levels$den$y, z = hpd.levels$den$z, levels = hpd.levels$falpha[1])
      answer = sp::point.in.polygon(param.pods[i.pod,"ta"], param.pods[i.pod,"raa"], drf.hdr[[1]]$x, drf.hdr[[1]]$y)
      if (answer == 1) {
        return(1/n.pods)
      } else {
        return(0)
      }
    }
    
    # Exporter les objets nécessaires au cluster
    clusterExport(cl, varlist = c("calcul_HPD_drf_hdr_poly", "pod.weighted.sample.drf", "param.pods", "n.pods"))
    # Exécuter en parallèle
    results <- parLapply(cl, 1:n.pods, calcul_HPD_drf_hdr_poly)
    # Calculer le total
    HPD90drf.DBan2.NBan2.hdr.poly <- sum(unlist(results))
    gc()
    
    # # Results: visualisation and saving: methode hdr
    # cat("HPD90rf.NBan2.NBan2.hdr.poly =",HPD90rf.NBan2.NBan2.hdr.poly,"\n")
    # col.target = col.target+1
    # RESULTS[k,col.target]=HPD90rf.NBan2.NBan2.hdr.poly
    # colnames(RESULTS)[col.target] <- "HPD90rf.NBan2.NBan2.hdr.poly"
    # cat("HPD90drf.NBan2.NBan2.hdr.poly =",HPD90drf.NBan2.NBan2.hdr.poly,"\n")
    # col.target = col.target+1
    # RESULTS[k,col.target]=HPD90drf.NBan2.NBan2.hdr.poly
    # colnames(RESULTS)[col.target] <- "HPD90drf.NBan2.NBan2.hdr.poly"
    
    
    ### GRAPH BIPLOT ASSOCIE ######################################
    # HPD levels correspond à votre conception (1-HPD) :
    # Ce que vous décrivez (HPD90 = 90% des valeurs estimées se trouvent à l'intérieur de la ligne HPD90) correspond 
    # effectivement à une probabilité cumulative. Dans hdr.2d, les niveaux HPD (falpha) fonctionnent différemment, 
    # indiquant les niveaux de densité associés. 
    # Il faut donc reformuler les contours pour refléter votre interprétation, en inversant l'ordre des niveaux.
    
    # Charger les bibliothèques nécessaires
    library(hdrcde)
    library(ggplot2)
    
    # Calcul des HPD bivariés via hdr.2d
    hpd.levels <- hdr.2d(
      x = pod.weighted.sample.drf[[i.pod]]$ta,
      y = pod.weighted.sample.drf[[i.pod]]$raa,
      prob = c(0.9, 0.75, 0.5, 0.25, 0.1)  # Notez l'inversion pour refléter 1-HPD
    )
    
    # Extraction des lignes de contours pour chaque niveau HPD
    drf.hdr <- contourLines(
      x = hpd.levels$den$x,
      y = hpd.levels$den$y,
      z = hpd.levels$den$z,
      levels = hpd.levels$falpha
    )
    
    # Convertir les contours en données exploitables pour ggplot2
    contour_data <- do.call(rbind, lapply(seq_along(drf.hdr), function(i) {
      data.frame(x = drf.hdr[[i]]$x, y = drf.hdr[[i]]$y, level = 1 - hpd.levels$alpha[i])  # Reflète 1-HPD
    }))
    
    # Supprimer les éventuelles lignes NA
    contour_data <- contour_data[!is.na(contour_data$level), ]
    
    # Charger les bibliothèques nécessaires
    # Charger les bibliothèques nécessaires
    library(hdrcde)
    library(ggplot2)
    
    # Calcul des HPD bivariés via hdr.2d
    hpd.levels <- hdr.2d(
      x = pod.weighted.sample.drf[[i.pod]]$DBan2,
      y = pod.weighted.sample.drf[[i.pod]]$NBan2,
      prob = c(0.9, 0.75, 0.5, 0.25, 0.1)  # Notez l'inversion pour refléter 1-HPD
    )
    
    # Extraction des lignes de contours pour chaque niveau HPD
    drf.hdr <- contourLines(
      x = hpd.levels$den$x,
      y = hpd.levels$den$y,
      z = hpd.levels$den$z,
      levels = hpd.levels$falpha
    )
    
    # Convertir les contours en données exploitables pour ggplot2
    contour_data <- do.call(rbind, lapply(seq_along(drf.hdr), function(i) {
      data.frame(x = drf.hdr[[i]]$x, y = drf.hdr[[i]]$y, level = 1 - hpd.levels$alpha[i])  # Reflète 1-HPD
    }))
    
    # Supprimer les éventuelles lignes NA
    contour_data <- contour_data[!is.na(contour_data$level), ]
    
    # Créer le graphique avec ggplot2
    ggplot() +
      # Ajouter les courbes de niveau pour chaque niveau HPD
      geom_path(data = contour_data, aes(x = x, y = y, group = level, color = as.factor(level)), size = 1) +
      # Personnalisation des axes et légendes
      labs(
        title = "DRF - HPD Bivariate Plot for ta and raa (introgression of East-Europe from Asia t=1000 sim=2000)" ,
        x = "ta",
        y = "raa",
        color = "HPD Level"
      ) +
      # Palette de couleurs inversée
      scale_color_manual(
        #values = c("blue", "green", "yellow", "orange", "red"), 
        values = c("red", "orange",  "yellow",  "green", "blue"),
        #labels = paste0(c("90%", "75%", "50%", "25%", "10%"), " HPD")
        labels = paste0(c("10%", "25%", "50%", "75%", "90%"), " HPD")
      ) +
      # Rendre les axes plus visibles
      theme_minimal() +
      theme(
        legend.position = "right",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        panel.grid = element_blank(),
        axis.line = element_line(size = 1.2, color = "black"),  # Épaisseur et couleur des axes
        axis.ticks = element_line(size = 1, color = "black"),   # Taille et couleur des ticks
        axis.text = element_text(size = 10, color = "black"),   # Couleur et taille des labels
        axis.title = element_text(size = 12, face = "bold")     # Style des titres
      )
    
    ###############################################################################################    
    
    ############## NEW ############# DRF (adapted 03-01-2024 !!! for Ghost project !!!):  DBnc.NBnc (Ghost USA2)
    colnames(pod.weighted.sample.drf[[i.pod]])
    
    # PARALLELISATION Methode drf-hdr-poly
    calcul_HPD_drf_hdr_poly <- function(i.pod) {
      hpd.levels <- hdrcde::hdr.2d(x = pod.weighted.sample.drf[[i.pod]]$DBnc, y = pod.weighted.sample.drf[[i.pod]]$NBnc, prob = c(0.1, 0.25, 0.5, 0.75, 0.9))
      drf.hdr <- contourLines(x = hpd.levels$den$x, y = hpd.levels$den$y, z = hpd.levels$den$z, levels = hpd.levels$falpha[1])
      answer = sp::point.in.polygon(param.pods[i.pod,"DBnc"], param.pods[i.pod,"DBnc"], drf.hdr[[1]]$x, drf.hdr[[1]]$y)
      if (answer == 1) {
        return(1/n.pods)
      } else {
        return(0)
      }
    }
    
    # Exporter les objets nécessaires au cluster
    clusterExport(cl, varlist = c("calcul_HPD_drf_hdr_poly", "pod.weighted.sample.drf", "param.pods", "n.pods"))
    # Exécuter en parallèle
    results <- parLapply(cl, 1:n.pods, calcul_HPD_drf_hdr_poly)
    # Calculer le total
    HPD90drf.DBnc.NBnc.hdr.poly <- sum(unlist(results))
    gc()
    
    # # Results: visualisation and saving: methode hdr
    # cat("HPD90rf.NBnc.NBnc.hdr.poly =",HPD90rf.NBnc.NBnc.hdr.poly,"\n")
    # col.target = col.target+1
    # RESULTS[k,col.target]=HPD90rf.NBnc.NBnc.hdr.poly
    # colnames(RESULTS)[col.target] <- "HPD90rf.NBnc.NBnc.hdr.poly"
    # cat("HPD90drf.NBnc.NBnc.hdr.poly =",HPD90drf.NBnc.NBnc.hdr.poly,"\n")
    # col.target = col.target+1
    # RESULTS[k,col.target]=HPD90drf.NBnc.NBnc.hdr.poly
    # colnames(RESULTS)[col.target] <- "HPD90drf.NBnc.NBnc.hdr.poly"
    
    
    ### GRAPH BIPLOT ASSOCIE ######################################
    # HPD levels correspond à votre conception (1-HPD) :
    # Ce que vous décrivez (HPD90 = 90% des valeurs estimées se trouvent à l'intérieur de la ligne HPD90) correspond 
    # effectivement à une probabilité cumulative. Dans hdr.2d, les niveaux HPD (falpha) fonctionnent différemment, 
    # indiquant les niveaux de densité associés. 
    # Il faut donc reformuler les contours pour refléter votre interprétation, en inversant l'ordre des niveaux.
    
    # Charger les bibliothèques nécessaires
    library(hdrcde)
    library(ggplot2)
    
    # Calcul des HPD bivariés via hdr.2d
    hpd.levels <- hdr.2d(
      x = pod.weighted.sample.drf[[i.pod]]$DBnc,
      y = pod.weighted.sample.drf[[i.pod]]$NBnc,
      prob = c(0.9, 0.75, 0.5, 0.25, 0.1)  # Notez l'inversion pour refléter 1-HPD
    )
    
    # Extraction des lignes de contours pour chaque niveau HPD
    drf.hdr <- contourLines(
      x = hpd.levels$den$x,
      y = hpd.levels$den$y,
      z = hpd.levels$den$z,
      levels = hpd.levels$falpha
    )
    
    # Convertir les contours en données exploitables pour ggplot2
    contour_data <- do.call(rbind, lapply(seq_along(drf.hdr), function(i) {
      data.frame(x = drf.hdr[[i]]$x, y = drf.hdr[[i]]$y, level = 1 - hpd.levels$alpha[i])  # Reflète 1-HPD
    }))
    
    # Supprimer les éventuelles lignes NA
    contour_data <- contour_data[!is.na(contour_data$level), ]
    
    # Charger les bibliothèques nécessaires
    # Charger les bibliothèques nécessaires
    library(hdrcde)
    library(ggplot2)
    
    # Calcul des HPD bivariés via hdr.2d
    hpd.levels <- hdr.2d(
      x = pod.weighted.sample.drf[[i.pod]]$DBnc,
      y = pod.weighted.sample.drf[[i.pod]]$NBnc,
      prob = c(0.9, 0.75, 0.5, 0.25, 0.1)  # Notez l'inversion pour refléter 1-HPD
    )
    
    # Extraction des lignes de contours pour chaque niveau HPD
    drf.hdr <- contourLines(
      x = hpd.levels$den$x,
      y = hpd.levels$den$y,
      z = hpd.levels$den$z,
      levels = hpd.levels$falpha
    )
    
    # Convertir les contours en données exploitables pour ggplot2
    contour_data <- do.call(rbind, lapply(seq_along(drf.hdr), function(i) {
      data.frame(x = drf.hdr[[i]]$x, y = drf.hdr[[i]]$y, level = 1 - hpd.levels$alpha[i])  # Reflète 1-HPD
    }))
    
    # Supprimer les éventuelles lignes NA
    contour_data <- contour_data[!is.na(contour_data$level), ]
    
    # Créer le graphique avec ggplot2
    ggplot() +
      # Ajouter les courbes de niveau pour chaque niveau HPD
      geom_path(data = contour_data, aes(x = x, y = y, group = level, color = as.factor(level)), size = 1) +
      # Personnalisation des axes et légendes
      labs(
        title = "HPD Bivariate Plot for DBnca and NBnca" ,
        x = "DBnca",
        y = "NBnca",
        color = "HPD Level"
      ) +
      # Palette de couleurs inversée
      scale_color_manual(
        #values = c("blue", "green", "yellow", "orange", "red"), 
        values = c("red", "orange",  "yellow",  "green", "blue"),
        #labels = paste0(c("90%", "75%", "50%", "25%", "10%"), " HPD")
        labels = paste0(c("10%", "25%", "50%", "75%", "90%"), " HPD")
      ) +
      # Rendre les axes plus visibles
      theme_minimal() +
      theme(
        legend.position = "right",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        panel.grid = element_blank(),
        axis.line = element_line(size = 1.2, color = "black"),  # Épaisseur et couleur des axes
        axis.ticks = element_line(size = 1, color = "black"),   # Taille et couleur des ticks
        axis.text = element_text(size = 10, color = "black"),   # Couleur et taille des labels
        axis.title = element_text(size = 12, face = "bold")     # Style des titres
      )
    
    ###############################################################################################    
    
    ############## NEW ############# DRF (adapted 03-01-2024 !!! for Ghost project !!!):  DBh.NBh (Ghost USA2)
    colnames(pod.weighted.sample.drf[[i.pod]])
    
    # PARALLELISATION Methode drf-hdr-poly
    calcul_HPD_drf_hdr_poly <- function(i.pod) {
      hpd.levels <- hdrcde::hdr.2d(x = pod.weighted.sample.drf[[i.pod]]$DBh, y = pod.weighted.sample.drf[[i.pod]]$NBh, prob = c(0.1, 0.25, 0.5, 0.75, 0.9))
      drf.hdr <- contourLines(x = hpd.levels$den$x, y = hpd.levels$den$y, z = hpd.levels$den$z, levels = hpd.levels$falpha[1])
      answer = sp::point.in.polygon(param.pods[i.pod,"DBh"], param.pods[i.pod,"DBh"], drf.hdr[[1]]$x, drf.hdr[[1]]$y)
      if (answer == 1) {
        return(1/n.pods)
      } else {
        return(0)
      }
    }
    
    # Exporter les objets nécessaires au cluster
    clusterExport(cl, varlist = c("calcul_HPD_drf_hdr_poly", "pod.weighted.sample.drf", "param.pods", "n.pods"))
    # Exécuter en parallèle
    results <- parLapply(cl, 1:n.pods, calcul_HPD_drf_hdr_poly)
    # Calculer le total
    HPD90drf.DBh.NBh.hdr.poly <- sum(unlist(results))
    gc()
    
    # # Results: visualisation and saving: methode hdr
    # cat("HPD90rf.NBh.NBh.hdr.poly =",HPD90rf.NBh.NBh.hdr.poly,"\n")
    # col.target = col.target+1
    # RESULTS[k,col.target]=HPD90rf.NBh.NBh.hdr.poly
    # colnames(RESULTS)[col.target] <- "HPD90rf.NBh.NBh.hdr.poly"
    # cat("HPD90drf.NBh.NBh.hdr.poly =",HPD90drf.NBh.NBh.hdr.poly,"\n")
    # col.target = col.target+1
    # RESULTS[k,col.target]=HPD90drf.NBh.NBh.hdr.poly
    # colnames(RESULTS)[col.target] <- "HPD90drf.NBh.NBh.hdr.poly"
    
    ### GRAPH BIPLOT ASSOCIE ######################################
    # HPD levels correspond à votre conception (1-HPD) :
    # Ce que vous décrivez (HPD90 = 90% des valeurs estimées se trouvent à l'intérieur de la ligne HPD90) correspond 
    # effectivement à une probabilité cumulative. Dans hdr.2d, les niveaux HPD (falpha) fonctionnent différemment, 
    # indiquant les niveaux de densité associés. 
    # Il faut donc reformuler les contours pour refléter votre interprétation, en inversant l'ordre des niveaux.
    
    # Charger les bibliothèques nécessaires
    library(hdrcde)
    library(ggplot2)
    
    # Calcul des HPD bivariés via hdr.2d
    hpd.levels <- hdr.2d(
      x = pod.weighted.sample.drf[[i.pod]]$DBh,
      y = pod.weighted.sample.drf[[i.pod]]$NBh,
      prob = c(0.9, 0.75, 0.5, 0.25, 0.1)  # Notez l'inversion pour refléter 1-HPD
    )
    
    # Extraction des lignes de contours pour chaque niveau HPD
    drf.hdr <- contourLines(
      x = hpd.levels$den$x,
      y = hpd.levels$den$y,
      z = hpd.levels$den$z,
      levels = hpd.levels$falpha
    )
    
    # Convertir les contours en données exploitables pour ggplot2
    contour_data <- do.call(rbind, lapply(seq_along(drf.hdr), function(i) {
      data.frame(x = drf.hdr[[i]]$x, y = drf.hdr[[i]]$y, level = 1 - hpd.levels$alpha[i])  # Reflète 1-HPD
    }))
    
    # Supprimer les éventuelles lignes NA
    contour_data <- contour_data[!is.na(contour_data$level), ]
    
    # Charger les bibliothèques nécessaires
    # Charger les bibliothèques nécessaires
    library(hdrcde)
    library(ggplot2)
    
    # Calcul des HPD bivariés via hdr.2d
    hpd.levels <- hdr.2d(
      x = pod.weighted.sample.drf[[i.pod]]$DBh,
      y = pod.weighted.sample.drf[[i.pod]]$NBh,
      prob = c(0.9, 0.75, 0.5, 0.25, 0.1)  # Notez l'inversion pour refléter 1-HPD
    )
    
    # Extraction des lignes de contours pour chaque niveau HPD
    drf.hdr <- contourLines(
      x = hpd.levels$den$x,
      y = hpd.levels$den$y,
      z = hpd.levels$den$z,
      levels = hpd.levels$falpha
    )
    
    # Convertir les contours en données exploitables pour ggplot2
    contour_data <- do.call(rbind, lapply(seq_along(drf.hdr), function(i) {
      data.frame(x = drf.hdr[[i]]$x, y = drf.hdr[[i]]$y, level = 1 - hpd.levels$alpha[i])  # Reflète 1-HPD
    }))
    
    # Supprimer les éventuelles lignes NA
    contour_data <- contour_data[!is.na(contour_data$level), ]
    
    # Créer le graphique avec ggplot2
    ggplot() +
      # Ajouter les courbes de niveau pour chaque niveau HPD
      geom_path(data = contour_data, aes(x = x, y = y, group = level, color = as.factor(level)), size = 1) +
      # Personnalisation des axes et légendes
      labs(
        title = "HPD Bivariate Plot for DBhaw and NBhaw" ,
        x = "DBhaw",
        y = "NBhaw",
        color = "HPD Level"
      ) +
      # Palette de couleurs inversée
      scale_color_manual(
        #values = c("blue", "green", "yellow", "orange", "red"), 
        values = c("red", "orange",  "yellow",  "green", "blue"),
        #labels = paste0(c("90%", "75%", "50%", "25%", "10%"), " HPD")
        labels = paste0(c("10%", "25%", "50%", "75%", "90%"), " HPD")
      ) +
      # Rendre les axes plus visibles
      theme_minimal() +
      theme(
        legend.position = "right",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        panel.grid = element_blank(),
        axis.line = element_line(size = 1.2, color = "black"),  # Épaisseur et couleur des axes
        axis.ticks = element_line(size = 1, color = "black"),   # Taille et couleur des ticks
        axis.text = element_text(size = 10, color = "black"),   # Couleur et taille des labels
        axis.title = element_text(size = 12, face = "bold")     # Style des titres
      )
    
    ###############################################################################################    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    ############################# CE QUI SUIT = fait avant le 01-01-2025 
    
    for (i.param in 1:length(param.list)) # Je calcule toutes les metrics d'intéret
    {
      single.param.pods = as.vector(param.pods[,i.param])
      
      # NMAE normalized by mean
      nom.result <- paste0(param.list[i.param], ".NMAE.mean.mean.rf")
      NMAE.mean.mean.rf = calculate_NMAE_mean(true_values = single.param.pods , point_estimates = param.pred.rf[[i.param]]$expectation)
      assign(nom.result, NMAE.mean.mean.rf)
      cat(nom.result, "=",NMAE.mean.mean.rf ,"\n")
      col.target = col.target+1
      RESULTS[k,col.target]=NMAE.mean.mean.rf
      colnames(RESULTS)[col.target] <- nom.result
      nom.result <- paste0(param.list[i.param], ".NMAE.mean.mean.drf")
      NMAE.mean.mean.drf = calculate_NMAE_mean(true_values = single.param.pods , point_estimates = pred.drf.mean[,i.param])
      assign(nom.result, NMAE.mean.mean.drf)
      cat(nom.result, "=",NMAE.mean.mean.drf ,"\n")
      col.target = col.target+1
      RESULTS[k,col.target]=NMAE.mean.mean.drf
      colnames(RESULTS)[col.target] <- nom.result
      
      # RMSE
      nom.result <- paste0(param.list[i.param], ".RMSE.rf")
      RMSE.rf = calculate_RMSE(true_values = single.param.pods , point_estimates = param.pred.rf[[i.param]]$expectation)
      assign(nom.result, RMSE.rf)
      cat(nom.result, "=",RMSE.rf ,"\n")
      col.target = col.target+1
      RESULTS[k,col.target]=RMSE.rf
      colnames(RESULTS)[col.target] <- nom.result
      nom.result <- paste0(param.list[i.param], ".RMSE.drf")
      RMSE.drf = calculate_RMSE(true_values = single.param.pods , point_estimates = pred.drf.mean[,i.param])
      assign(nom.result, RMSE.drf)
      cat(nom.result, "=",RMSE.drf ,"\n")
      col.target = col.target+1
      RESULTS[k,col.target]=RMSE.drf
      colnames(RESULTS)[col.target] <- nom.result
      
      # mean(SD)   
      nom.result <- paste0(param.list[i.param], ".sd.mean.rf")
      sd.mean.rf = mean(sqrt(param.pred.rf[[i.param]]$variance.cdf))
      assign(nom.result, sd.mean.rf)
      cat(nom.result, "=",sd.mean.rf,"\n")
      col.target = col.target+1
      RESULTS[k,col.target]=sd.mean.rf
      colnames(RESULTS)[col.target] <- nom.result
      nom.result <- paste0(param.list[i.param], ".sd.mean.drf")
      sd.mean.drf = mean(pred.drf.sd[,i.param])
      assign(nom.result, sd.mean.drf)
      cat(nom.result, "=",sd.mean.drf,"\n")
      col.target = col.target+1
      RESULTS[k,col.target]=sd.mean.drf
      colnames(RESULTS)[col.target] <- nom.result
      
      # 90% CI
      nom.result <- paste0(param.list[i.param], ".CI90.rf")
      CI90.rf <- calculate_90CI(true_values = single.param.pods, 
                                Q5_point_estimates = param.pred.rf[[i.param]][["quantiles"]][,1], Q95_point_estimates = param.pred.rf[[i.param]][["quantiles"]][,2])
      assign(nom.result,CI90.rf)
      cat(nom.result, "=",CI90.rf ,"\n")
      col.target = col.target+1
      RESULTS[k,col.target]=CI90.rf
      colnames(RESULTS)[col.target] <- nom.result
      nom.result <- paste0(param.list[i.param], ".CI90.drf")
      CI90.drf <- calculate_90CI(true_values = single.param.pods, 
                                 Q5_point_estimates = pred.drf.quant[, i.param], Q95_point_estimates = pred.drf.quant[, i.param+length(param.list)])
      assign(nom.result,CI90.drf)
      cat(nom.result, "=",CI90.drf ,"\n")
      col.target = col.target+1
      RESULTS[k,col.target]=CI90.drf
      colnames(RESULTS)[col.target] <- nom.result
      
      # Length 90% CI
      nom.result <- paste0(param.list[i.param], ".lengthCI90.rf")
      lengthCI90.rf <- mean(calculate_lengthCI90(true_values = single.param.pods, 
                                                 Q5_point_estimates = param.pred.rf[[i.param]][["quantiles"]][,1], Q95_point_estimates = param.pred.rf[[i.param]][["quantiles"]][,2]))
      assign(nom.result,lengthCI90.rf)
      cat(nom.result, "=",lengthCI90.rf ,"\n")
      col.target = col.target+1
      RESULTS[k,col.target]=lengthCI90.rf
      colnames(RESULTS)[col.target] <- nom.result
      nom.result <- paste0(param.list[i.param], ".lengthCI90.drf")
      lengthCI90.drf <- mean(calculate_lengthCI90(true_values = single.param.pods, 
                                                  Q5_point_estimates = pred.drf.quant[, i.param], Q95_point_estimates = pred.drf.quant[, i.param+length(param.list)]))
      assign(nom.result,lengthCI90.drf)
      cat(nom.result, "=",lengthCI90.drf ,"\n")
      col.target = col.target+1
      RESULTS[k,col.target]=lengthCI90.drf
      colnames(RESULTS)[col.target] <- nom.result
      
      # HPD90%.monoparam
      nom.result <- paste0(param.list[i.param], ".HPD90.monoparam.rf")
      HPD90.monoparam.rf <- calculate_HPD90(true_values = param.pods[,i.param], lower = HPD90rf.lower.for.each.param.all.pods[,i.param],
                                            upper = HPD90rf.upper.for.each.param.all.pods[,i.param])
      assign(nom.result,HPD90.monoparam.rf)
      cat(nom.result, "=",HPD90.monoparam.rf ,"\n")
      col.target = col.target+1
      RESULTS[k,col.target]=HPD90.monoparam.rf
      colnames(RESULTS)[col.target] <- nom.result
      
      nom.result <- paste0(param.list[i.param], ".HPD90.monoparam.drf")
      HPD90.monoparam.drf <- calculate_HPD90(true_values = param.pods[,i.param], lower = HPD90drf.lower.for.each.param.all.pods[,i.param],
                                             upper = HPD90drf.upper.for.each.param.all.pods[,i.param])
      assign(nom.result,HPD90.monoparam.drf)
      cat(nom.result, "=",HPD90.monoparam.drf ,"\n")
      col.target = col.target+1
      RESULTS[k,col.target]=HPD90.monoparam.drf
      colnames(RESULTS)[col.target] <- nom.result
      
      # Length HPD90%
      nom.result <- paste0(param.list[i.param], ".lengthHPD90.rf")
      lengthHPD90.rf <- mean(calculate_lengthHPD90(true_values = param.pods[,i.param], lower = HPD90rf.lower.for.each.param.all.pods[,i.param],
                                                   upper = HPD90rf.upper.for.each.param.all.pods[,i.param]))
      assign(nom.result,lengthHPD90.rf)
      cat(nom.result, "=",lengthHPD90.rf ,"\n")
      col.target = col.target+1
      RESULTS[k,col.target]=lengthHPD90.rf
      colnames(RESULTS)[col.target] <- nom.result
      
      nom.result <- paste0(param.list[i.param], ".lengthHPD90.drf")
      lengthHPD90.drf <- mean(calculate_lengthHPD90(true_values = param.pods[,i.param], lower = HPD90drf.lower.for.each.param.all.pods[,i.param],
                                                    upper = HPD90drf.upper.for.each.param.all.pods[,i.param]))    
      assign(nom.result,lengthHPD90.drf)
      cat(nom.result, "=",lengthHPD90.drf ,"\n")
      col.target = col.target+1
      RESULTS[k,col.target]=lengthHPD90.drf
      colnames(RESULTS)[col.target] <- nom.result
      
    }
    
    ##CORRELATIONS #######################################
    ##INTERSECTION.CI90.drf = calculate_90CI_intersection(V1_true_values = N11.TARGET.param.pods.vector, V1_Q5 = all.pred.drf.quant.target$quantile.8.q.0.05, V1_Q95 = all.pred.drf.quant.target$quantile.8.q.0.95,
    ##                                                                                  V2_true_values = t.TARGET.param.pods.vector, V2_Q5 = all.pred.drf.quant.target$quantile.5.q.0.05, V2_Q95 = all.pred.drf.quant.target$quantile.5.q.0.95)
    #SCENARIO 6 POPS
    #colnames(param.pods)
    # [1] "N123" "N4"   "N6"   "t124" "Dbn3" "Nbn3" "t12"  "t15"  "ra"  
    #cat("INTERSECTION.CI90.drf =",INTERSECTION.CI90.drf,"\n")
    # cor N123.t12
    mean.cor = mean(pred.drf.cor$cor.1.7)
    cat("t12.N123.mean.cor =",mean.cor,"\n")
    col.target = col.target+1
    RESULTS[k,col.target]=mean.cor
    colnames(RESULTS)[col.target] <- "t12.N123.mean.cor"
    sd.cor = sqrt(var(pred.drf.cor$cor.1.7))
    cat("t12.N123.sd.cor =",sd.cor,"\n")
    col.target = col.target+1
    RESULTS[k,col.target]=sd.cor
    colnames(RESULTS)[col.target] <- "t12.N123.sd.cor"
    
    # cor Dbn3.Nbn3
    mean.cor = mean(pred.drf.cor$cor.5.6)
    cat("Dbn3.Nbn3.mean.cor =",mean.cor,"\n")
    col.target = col.target+1
    RESULTS[k,col.target]=mean.cor
    colnames(RESULTS)[col.target] <- "Dbn3.Nbn3.mean.cor"
    sd.cor = sqrt(var(pred.drf.cor$cor.5.6))
    cat("Dbn3.Nbn3.sd.cor =",sd.cor,"\n")
    col.target = col.target+1
    RESULTS[k,col.target]=sd.cor
    colnames(RESULTS)[col.target] <- "Dbn3.Nbn3.sd.cor"
    
    # cor N6.ra
    mean.cor = mean(pred.drf.cor$cor.3.9)
    cat("N6.ra.mean.cor =",mean.cor,"\n")
    col.target = col.target+1
    RESULTS[k,col.target]=mean.cor
    colnames(RESULTS)[col.target] <- "N6.ra.mean.cor"
    sd.cor = sqrt(var(pred.drf.cor$cor.3.9))
    cat("N6.ra.sd.cor =",sd.cor,"\n")
    col.target = col.target+1
    RESULTS[k,col.target]=sd.cor
    colnames(RESULTS)[col.target] <- "N6.ra.sd.cor"
    
    # #SCENARIO Microsat 4 pops t.fixed
    # #colnames(param.pods) --> "Ne"  "µ"   "Neµ"
    # # cor Ne.µ
    # mean.cor = mean(pred.drf.cor$cor.1.2)
    # cat("Ne.µ.mean.cor =",mean.cor,"\n")
    # col.target = col.target+1
    # RESULTS[k,col.target]=mean.cor
    # colnames(RESULTS)[col.target] <- "Ne.µ.mean.cor"
    # sd.cor = sqrt(var(pred.drf.cor$cor.1.2))
    # cat("Ne.µ.sd.cor =",sd.cor,"\n")
    # col.target = col.target+1
    # RESULTS[k,col.target]=sd.cor
    # colnames(RESULTS)[col.target] <- "Ne.µ.sd.cor"
    
    ##############################################################################################
    ##############################################################################################
    ## HPD biparam Dbn3.Nbn3 (param 5 et 6) ----------- Methode "Paul" = hdr
    ##############################################################################################
    ##############################################################################################
    
    # Initialisation du cluster (valable pour tous les calculs HPD biparam --> a la fin des calculs paralleles faire stopCluster(cl) Fermer le cluster après utilisation  )
    #library(parallel)
    numCores <- detectCores() - 1
    cl <- makeCluster(numCores)
    # Charger les bibliothèques nécessaires sur chaque nœud du cluster
    clusterEvalQ(cl, {
      library(hdrcde)
      library(MASS)
      library(splancs)  # Pour areapl
      library(sp)  # Pour point.in.polygon
    })
    
    # ########################### RF
    # 
    # # PARALLELISATION: Methode rf-hdr-poly
    # # Définition de la fonction pour calculer HPD et vérifier les points
    # calcul_HPD_rf_hdr_poly <- function(i.pod) {
    #   hpd.levels <- hdrcde::hdr.2d(x = pod.weighted.sample.rf[[5]][[i.pod]]$Dbn3, y = pod.weighted.sample.rf[[6]][[i.pod]]$Nbn3, prob = c(0.1, 0.25, 0.5, 0.75, 0.9))
    #   rf.hdr <- contourLines(x = hpd.levels$den$x, y = hpd.levels$den$y,z = hpd.levels$den$z,levels = hpd.levels$falpha[1])
    #   answer = sp::point.in.polygon(param.pods[i.pod,"Dbn3"], param.pods[i.pod,"Nbn3"], rf.hdr[[1]]$x, rf.hdr[[1]]$y)
    #   if (answer == 1) {
    #     return(1)
    #   } else {
    #     return(0)
    #   }
    # }
    
    # Exporter les objets nécessaires au cluster
    clusterExport(cl, varlist = c("calcul_HPD_rf_hdr_poly", "pod.weighted.sample.rf", "param.pods"))
    # Exécuter en parallèle
    results <- parLapply(cl, 1:n.pods, calcul_HPD_rf_hdr_poly)
    # Calculer le total et nettoyage
    HPD90rf.Dbn3.Nbn3.hdr.poly <- sum(unlist(results)) / n.pods
    gc()
    
    ############## NEW ############# DRF (adapted 03-01-2024 !!! for Ghost project !!!):  DBan2.NBan2 (Ghost USA2)
    
    # PARALLELISATION Methode drf-hdr-poly
    calcul_HPD_drf_hdr_poly <- function(i.pod) {
      hpd.levels <- hdrcde::hdr.2d(x = pod.weighted.sample.drf[[i.pod]]$DBan2, y = pod.weighted.sample.drf[[i.pod]]$NBan2, prob = c(0.1, 0.25, 0.5, 0.75, 0.9))
      drf.hdr <- contourLines(x = hpd.levels$den$x, y = hpd.levels$den$y, z = hpd.levels$den$z, levels = hpd.levels$falpha[1])
      answer = sp::point.in.polygon(param.pods[i.pod,"DBan2"], param.pods[i.pod,"DBan2"], drf.hdr[[1]]$x, drf.hdr[[1]]$y)
      if (answer == 1) {
        return(1/n.pods)
      } else {
        return(0)
      }
    }
    
    # Exporter les objets nécessaires au cluster
    clusterExport(cl, varlist = c("calcul_HPD_drf_hdr_poly", "pod.weighted.sample.drf", "param.pods", "n.pods"))
    # Exécuter en parallèle
    results <- parLapply(cl, 1:n.pods, calcul_HPD_drf_hdr_poly)
    # Calculer le total
    HPD90drf.DBan2.NBan2.hdr.poly <- sum(unlist(results))
    gc()
    
    # Results: visualisation and saving: methode hdr
    cat("HPD90rf.NBan2.NBan2.hdr.poly =",HPD90rf.NBan2.NBan2.hdr.poly,"\n")
    col.target = col.target+1
    RESULTS[k,col.target]=HPD90rf.NBan2.NBan2.hdr.poly
    colnames(RESULTS)[col.target] <- "HPD90rf.NBan2.NBan2.hdr.poly"
    cat("HPD90drf.NBan2.NBan2.hdr.poly =",HPD90drf.NBan2.NBan2.hdr.poly,"\n")
    col.target = col.target+1
    RESULTS[k,col.target]=HPD90drf.NBan2.NBan2.hdr.poly
    colnames(RESULTS)[col.target] <- "HPD90drf.NBan2.NBan2.hdr.poly"
    
    
    
    
  ############################################################################### 
  ### GRAPH BIPLOT ASSOCIE ########## NEW 01-01-2025 ############################
  ############################################################################### 
    
    # HPD levels correspond à votre conception (1-HPD) :
    # Ce que vous décrivez (HPD90 = 90% des valeurs estimées se trouvent à l'intérieur de la ligne HPD90) correspond 
    # effectivement à une probabilité cumulative. Dans hdr.2d, les niveaux HPD (falpha) fonctionnent différemment, 
    # indiquant les niveaux de densité associés. 
    # Il faut donc reformuler les contours pour refléter votre interprétation, en inversant l'ordre des niveaux.
    
    # Charger les bibliothèques nécessaires
    library(hdrcde)
    library(ggplot2)
    
    # Calcul des HPD bivariés via hdr.2d
    hpd.levels <- hdr.2d(
      x = pod.weighted.sample.drf[[i.pod]]$DBan2,
      y = pod.weighted.sample.drf[[i.pod]]$NBan2,
      prob = c(0.9, 0.75, 0.5, 0.25, 0.1)  # Notez l'inversion pour refléter 1-HPD
    )
    
    # Extraction des lignes de contours pour chaque niveau HPD
    drf.hdr <- contourLines(
      x = hpd.levels$den$x,
      y = hpd.levels$den$y,
      z = hpd.levels$den$z,
      levels = hpd.levels$falpha
    )
    
    # Convertir les contours en données exploitables pour ggplot2
    contour_data <- do.call(rbind, lapply(seq_along(drf.hdr), function(i) {
      data.frame(x = drf.hdr[[i]]$x, y = drf.hdr[[i]]$y, level = 1 - hpd.levels$alpha[i])  # Reflète 1-HPD
    }))
    
    # Supprimer les éventuelles lignes NA
    contour_data <- contour_data[!is.na(contour_data$level), ]
    
    # Charger les bibliothèques nécessaires
    # Charger les bibliothèques nécessaires
    library(hdrcde)
    library(ggplot2)
    
    # Calcul des HPD bivariés via hdr.2d
    hpd.levels <- hdr.2d(
      x = pod.weighted.sample.drf[[i.pod]]$DBan2,
      y = pod.weighted.sample.drf[[i.pod]]$NBan2,
      prob = c(0.9, 0.75, 0.5, 0.25, 0.1)  # Notez l'inversion pour refléter 1-HPD
    )
    
    # Extraction des lignes de contours pour chaque niveau HPD
    drf.hdr <- contourLines(
      x = hpd.levels$den$x,
      y = hpd.levels$den$y,
      z = hpd.levels$den$z,
      levels = hpd.levels$falpha
    )
    
    # Convertir les contours en données exploitables pour ggplot2
    contour_data <- do.call(rbind, lapply(seq_along(drf.hdr), function(i) {
      data.frame(x = drf.hdr[[i]]$x, y = drf.hdr[[i]]$y, level = 1 - hpd.levels$alpha[i])  # Reflète 1-HPD
    }))
    
    # Supprimer les éventuelles lignes NA
    contour_data <- contour_data[!is.na(contour_data$level), ]
    
    # Créer le graphique avec ggplot2
    ggplot() +
      # Ajouter les courbes de niveau pour chaque niveau HPD
      geom_path(data = contour_data, aes(x = x, y = y, group = level, color = as.factor(level)), size = 1) +
      # Personnalisation des axes et légendes
      labs(
        title = "HPD Bivariate Plot for DBan2 and NBan2" ,
        x = "DBan2",
        y = "NBan2",
        color = "HPD Level"
      ) +
      # Palette de couleurs inversée
      scale_color_manual(
        #values = c("blue", "green", "yellow", "orange", "red"), 
        values = c("red", "orange",  "yellow",  "green", "blue"),
        #labels = paste0(c("90%", "75%", "50%", "25%", "10%"), " HPD")
        labels = paste0(c("10%", "25%", "50%", "75%", "90%"), " HPD")
      ) +
      # Rendre les axes plus visibles
      theme_minimal() +
      theme(
        legend.position = "right",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        panel.grid = element_blank(),
        axis.line = element_line(size = 1.2, color = "black"),  # Épaisseur et couleur des axes
        axis.ticks = element_line(size = 1, color = "black"),   # Taille et couleur des ticks
        axis.text = element_text(size = 10, color = "black"),   # Couleur et taille des labels
        axis.title = element_text(size = 12, face = "bold")     # Style des titres
      )
    
    #######################################################################################################
    
    ############## NEW ############# DRF (adapted 03-01-2024 !!! for Ghost project !!!):  DBan2.NBan2 (Ghost USA2)
    
    # PARALLELISATION Methode drf-hdr-poly
    calcul_HPD_drf_hdr_poly <- function(i.pod) {
      hpd.levels <- hdrcde::hdr.2d(x = pod.weighted.sample.drf[[i.pod]]$DBan2, y = pod.weighted.sample.drf[[i.pod]]$NBan2, prob = c(0.1, 0.25, 0.5, 0.75, 0.9))
      drf.hdr <- contourLines(x = hpd.levels$den$x, y = hpd.levels$den$y, z = hpd.levels$den$z, levels = hpd.levels$falpha[1])
      answer = sp::point.in.polygon(param.pods[i.pod,"DBan2"], param.pods[i.pod,"DBan2"], drf.hdr[[1]]$x, drf.hdr[[1]]$y)
      if (answer == 1) {
        return(1/n.pods)
      } else {
        return(0)
      }
    }
    
    # Exporter les objets nécessaires au cluster
    clusterExport(cl, varlist = c("calcul_HPD_drf_hdr_poly", "pod.weighted.sample.drf", "param.pods", "n.pods"))
    # Exécuter en parallèle
    results <- parLapply(cl, 1:n.pods, calcul_HPD_drf_hdr_poly)
    # Calculer le total
    HPD90drf.DBan2.NBan2.hdr.poly <- sum(unlist(results))
    gc()
    
    # Results: visualisation and saving: methode hdr
    cat("HPD90rf.NBan2.NBan2.hdr.poly =",HPD90rf.NBan2.NBan2.hdr.poly,"\n")
    col.target = col.target+1
    RESULTS[k,col.target]=HPD90rf.NBan2.NBan2.hdr.poly
    colnames(RESULTS)[col.target] <- "HPD90rf.NBan2.NBan2.hdr.poly"
    cat("HPD90drf.NBan2.NBan2.hdr.poly =",HPD90drf.NBan2.NBan2.hdr.poly,"\n")
    col.target = col.target+1
    RESULTS[k,col.target]=HPD90drf.NBan2.NBan2.hdr.poly
    colnames(RESULTS)[col.target] <- "HPD90drf.NBan2.NBan2.hdr.poly"
    
    
    ### GRAPH BIPLOT ASSOCIE ######################################
    # HPD levels correspond à votre conception (1-HPD) :
    # Ce que vous décrivez (HPD90 = 90% des valeurs estimées se trouvent à l'intérieur de la ligne HPD90) correspond 
    # effectivement à une probabilité cumulative. Dans hdr.2d, les niveaux HPD (falpha) fonctionnent différemment, 
    # indiquant les niveaux de densité associés. 
    # Il faut donc reformuler les contours pour refléter votre interprétation, en inversant l'ordre des niveaux.
    
    # Charger les bibliothèques nécessaires
    library(hdrcde)
    library(ggplot2)
    
    # Calcul des HPD bivariés via hdr.2d
    hpd.levels <- hdr.2d(
      x = pod.weighted.sample.drf[[i.pod]]$DBan2,
      y = pod.weighted.sample.drf[[i.pod]]$NBan2,
      prob = c(0.9, 0.75, 0.5, 0.25, 0.1)  # Notez l'inversion pour refléter 1-HPD
    )
    
    # Extraction des lignes de contours pour chaque niveau HPD
    drf.hdr <- contourLines(
      x = hpd.levels$den$x,
      y = hpd.levels$den$y,
      z = hpd.levels$den$z,
      levels = hpd.levels$falpha
    )
    
    # Convertir les contours en données exploitables pour ggplot2
    contour_data <- do.call(rbind, lapply(seq_along(drf.hdr), function(i) {
      data.frame(x = drf.hdr[[i]]$x, y = drf.hdr[[i]]$y, level = 1 - hpd.levels$alpha[i])  # Reflète 1-HPD
    }))
    
    # Supprimer les éventuelles lignes NA
    contour_data <- contour_data[!is.na(contour_data$level), ]
    
    # Charger les bibliothèques nécessaires
    # Charger les bibliothèques nécessaires
    library(hdrcde)
    library(ggplot2)
    
    # Calcul des HPD bivariés via hdr.2d
    hpd.levels <- hdr.2d(
      x = pod.weighted.sample.drf[[i.pod]]$DBan2,
      y = pod.weighted.sample.drf[[i.pod]]$NBan2,
      prob = c(0.9, 0.75, 0.5, 0.25, 0.1)  # Notez l'inversion pour refléter 1-HPD
    )
    
    # Extraction des lignes de contours pour chaque niveau HPD
    drf.hdr <- contourLines(
      x = hpd.levels$den$x,
      y = hpd.levels$den$y,
      z = hpd.levels$den$z,
      levels = hpd.levels$falpha
    )
    
    # Convertir les contours en données exploitables pour ggplot2
    contour_data <- do.call(rbind, lapply(seq_along(drf.hdr), function(i) {
      data.frame(x = drf.hdr[[i]]$x, y = drf.hdr[[i]]$y, level = 1 - hpd.levels$alpha[i])  # Reflète 1-HPD
    }))
    
    # Supprimer les éventuelles lignes NA
    contour_data <- contour_data[!is.na(contour_data$level), ]
    
    # Créer le graphique avec ggplot2
    ggplot() +
      # Ajouter les courbes de niveau pour chaque niveau HPD
      geom_path(data = contour_data, aes(x = x, y = y, group = level, color = as.factor(level)), size = 1) +
      # Personnalisation des axes et légendes
      labs(
        title = "HPD Bivariate Plot for DBan2 and NBan2" ,
        x = "DBan2",
        y = "NBan2",
        color = "HPD Level"
      ) +
      # Palette de couleurs inversée
      scale_color_manual(
        #values = c("blue", "green", "yellow", "orange", "red"), 
        values = c("red", "orange",  "yellow",  "green", "blue"),
        #labels = paste0(c("90%", "75%", "50%", "25%", "10%"), " HPD")
        labels = paste0(c("10%", "25%", "50%", "75%", "90%"), " HPD")
      ) +
      # Rendre les axes plus visibles
      theme_minimal() +
      theme(
        legend.position = "right",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        panel.grid = element_blank(),
        axis.line = element_line(size = 1.2, color = "black"),  # Épaisseur et couleur des axes
        axis.ticks = element_line(size = 1, color = "black"),   # Taille et couleur des ticks
        axis.text = element_text(size = 10, color = "black"),   # Couleur et taille des labels
        axis.title = element_text(size = 12, face = "bold")     # Style des titres
      )

###############################################################################################    
    
    ############## NEW ############# DRF (adapted 03-01-2024 !!! for Ghost project !!!):  DBan2.NBan2 (Ghost USA2)
    
    # PARALLELISATION Methode drf-hdr-poly
    calcul_HPD_drf_hdr_poly <- function(i.pod) {
      hpd.levels <- hdrcde::hdr.2d(x = pod.weighted.sample.drf[[i.pod]]$DBan2, y = pod.weighted.sample.drf[[i.pod]]$NBan2, prob = c(0.1, 0.25, 0.5, 0.75, 0.9))
      drf.hdr <- contourLines(x = hpd.levels$den$x, y = hpd.levels$den$y, z = hpd.levels$den$z, levels = hpd.levels$falpha[1])
      answer = sp::point.in.polygon(param.pods[i.pod,"DBan2"], param.pods[i.pod,"DBan2"], drf.hdr[[1]]$x, drf.hdr[[1]]$y)
      if (answer == 1) {
        return(1/n.pods)
      } else {
        return(0)
      }
    }
    
    # Exporter les objets nécessaires au cluster
    clusterExport(cl, varlist = c("calcul_HPD_drf_hdr_poly", "pod.weighted.sample.drf", "param.pods", "n.pods"))
    # Exécuter en parallèle
    results <- parLapply(cl, 1:n.pods, calcul_HPD_drf_hdr_poly)
    # Calculer le total
    HPD90drf.DBan2.NBan2.hdr.poly <- sum(unlist(results))
    gc()
    
    # Results: visualisation and saving: methode hdr
    cat("HPD90rf.NBan2.NBan2.hdr.poly =",HPD90rf.NBan2.NBan2.hdr.poly,"\n")
    col.target = col.target+1
    RESULTS[k,col.target]=HPD90rf.NBan2.NBan2.hdr.poly
    colnames(RESULTS)[col.target] <- "HPD90rf.NBan2.NBan2.hdr.poly"
    cat("HPD90drf.NBan2.NBan2.hdr.poly =",HPD90drf.NBan2.NBan2.hdr.poly,"\n")
    col.target = col.target+1
    RESULTS[k,col.target]=HPD90drf.NBan2.NBan2.hdr.poly
    colnames(RESULTS)[col.target] <- "HPD90drf.NBan2.NBan2.hdr.poly"
    
    
    ### GRAPH BIPLOT ASSOCIE ######################################
    # HPD levels correspond à votre conception (1-HPD) :
    # Ce que vous décrivez (HPD90 = 90% des valeurs estimées se trouvent à l'intérieur de la ligne HPD90) correspond 
    # effectivement à une probabilité cumulative. Dans hdr.2d, les niveaux HPD (falpha) fonctionnent différemment, 
    # indiquant les niveaux de densité associés. 
    # Il faut donc reformuler les contours pour refléter votre interprétation, en inversant l'ordre des niveaux.
    
    # Charger les bibliothèques nécessaires
    library(hdrcde)
    library(ggplot2)
    
    # Calcul des HPD bivariés via hdr.2d
    hpd.levels <- hdr.2d(
      x = pod.weighted.sample.drf[[i.pod]]$DBan2,
      y = pod.weighted.sample.drf[[i.pod]]$NBan2,
      prob = c(0.9, 0.75, 0.5, 0.25, 0.1)  # Notez l'inversion pour refléter 1-HPD
    )
    
    # Extraction des lignes de contours pour chaque niveau HPD
    drf.hdr <- contourLines(
      x = hpd.levels$den$x,
      y = hpd.levels$den$y,
      z = hpd.levels$den$z,
      levels = hpd.levels$falpha
    )
    
    # Convertir les contours en données exploitables pour ggplot2
    contour_data <- do.call(rbind, lapply(seq_along(drf.hdr), function(i) {
      data.frame(x = drf.hdr[[i]]$x, y = drf.hdr[[i]]$y, level = 1 - hpd.levels$alpha[i])  # Reflète 1-HPD
    }))
    
    # Supprimer les éventuelles lignes NA
    contour_data <- contour_data[!is.na(contour_data$level), ]
    
    # Charger les bibliothèques nécessaires
    # Charger les bibliothèques nécessaires
    library(hdrcde)
    library(ggplot2)
    
    # Calcul des HPD bivariés via hdr.2d
    hpd.levels <- hdr.2d(
      x = pod.weighted.sample.drf[[i.pod]]$DBan2,
      y = pod.weighted.sample.drf[[i.pod]]$NBan2,
      prob = c(0.9, 0.75, 0.5, 0.25, 0.1)  # Notez l'inversion pour refléter 1-HPD
    )
    
    # Extraction des lignes de contours pour chaque niveau HPD
    drf.hdr <- contourLines(
      x = hpd.levels$den$x,
      y = hpd.levels$den$y,
      z = hpd.levels$den$z,
      levels = hpd.levels$falpha
    )
    
    # Convertir les contours en données exploitables pour ggplot2
    contour_data <- do.call(rbind, lapply(seq_along(drf.hdr), function(i) {
      data.frame(x = drf.hdr[[i]]$x, y = drf.hdr[[i]]$y, level = 1 - hpd.levels$alpha[i])  # Reflète 1-HPD
    }))
    
    # Supprimer les éventuelles lignes NA
    contour_data <- contour_data[!is.na(contour_data$level), ]
    
    # Créer le graphique avec ggplot2
    ggplot() +
      # Ajouter les courbes de niveau pour chaque niveau HPD
      geom_path(data = contour_data, aes(x = x, y = y, group = level, color = as.factor(level)), size = 1) +
      # Personnalisation des axes et légendes
      labs(
        title = "HPD Bivariate Plot for DBan2 and NBan2" ,
        x = "DBan2",
        y = "NBan2",
        color = "HPD Level"
      ) +
      # Palette de couleurs inversée
      scale_color_manual(
        #values = c("blue", "green", "yellow", "orange", "red"), 
        values = c("red", "orange",  "yellow",  "green", "blue"),
        #labels = paste0(c("90%", "75%", "50%", "25%", "10%"), " HPD")
        labels = paste0(c("10%", "25%", "50%", "75%", "90%"), " HPD")
      ) +
      # Rendre les axes plus visibles
      theme_minimal() +
      theme(
        legend.position = "right",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        panel.grid = element_blank(),
        axis.line = element_line(size = 1.2, color = "black"),  # Épaisseur et couleur des axes
        axis.ticks = element_line(size = 1, color = "black"),   # Taille et couleur des ticks
        axis.text = element_text(size = 10, color = "black"),   # Couleur et taille des labels
        axis.title = element_text(size = 12, face = "bold")     # Style des titres
      )
    
    ###############################################################################################    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
        
    
    
    
    ##############################################################################################
    ##############################################################################################
    ## HPD biparam N123.t12 (param 1 et 7) ----------- Methode "Paul" = hdr
    ##############################################################################################
    ##############################################################################################
    
    ########################### RF
    
    # PARALLELISATION: Methode rf-hdr-poly
    # Définition de la fonction pour calculer HPD et vérifier les points
    calcul_HPD_rf_hdr_poly <- function(i.pod) {
      hpd.levels <- hdrcde::hdr.2d(x = pod.weighted.sample.rf[[1]][[i.pod]]$N123, y = pod.weighted.sample.rf[[7]][[i.pod]]$t12, prob = c(0.1, 0.25, 0.5, 0.75, 0.9))
      rf.hdr <- contourLines(x = hpd.levels$den$x, y = hpd.levels$den$y,z = hpd.levels$den$z,levels = hpd.levels$falpha[1])
      answer = sp::point.in.polygon(param.pods[i.pod,"N123"], param.pods[i.pod,"t12"], rf.hdr[[1]]$x, rf.hdr[[1]]$y)
      if (answer == 1) {
        return(1)
      } else {
        return(0)
      }
    }
    
    # Exporter les objets nécessaires au cluster
    clusterExport(cl, varlist = c("calcul_HPD_rf_hdr_poly", "pod.weighted.sample.rf", "param.pods"))
    # Exécuter en parallèle
    results <- parLapply(cl, 1:n.pods, calcul_HPD_rf_hdr_poly)
    # Calculer le total et nettoyage
    HPD90rf.N123.t12.hdr.poly <- sum(unlist(results)) / n.pods
    gc()
    
    ########################### DRF
    
    # PARALLELISATION Methode drf-hdr-poly
    calcul_HPD_drf_hdr_poly <- function(i.pod) {
      hpd.levels <- hdrcde::hdr.2d(x = pod.weighted.sample.drf[[i.pod]]$N123, y = pod.weighted.sample.drf[[i.pod]]$t12, prob = c(0.1, 0.25, 0.5, 0.75, 0.9))
      drf.hdr <- contourLines(x = hpd.levels$den$x, y = hpd.levels$den$y, z = hpd.levels$den$z, levels = hpd.levels$falpha[1])
      answer = sp::point.in.polygon(param.pods[i.pod,"N123"], param.pods[i.pod,"t12"], drf.hdr[[1]]$x, drf.hdr[[1]]$y)
      if (answer == 1) {
        return(1/n.pods)
      } else {
        return(0)
      }
    }
    
    # Exporter les objets nécessaires au cluster
    clusterExport(cl, varlist = c("calcul_HPD_drf_hdr_poly", "pod.weighted.sample.drf", "param.pods", "n.pods"))
    # Exécuter en parallèle
    results <- parLapply(cl, 1:n.pods, calcul_HPD_drf_hdr_poly)
    # Calculer le total
    HPD90drf.N123.t12.hdr.poly <- sum(unlist(results))
    gc()
    
    # Results: visualisation and saving: methode hdr
    cat("HPD90rf.N123.t12.hdr.poly =",HPD90rf.N123.t12.hdr.poly,"\n")
    col.target = col.target+1
    RESULTS[k,col.target]=HPD90rf.N123.t12.hdr.poly
    colnames(RESULTS)[col.target] <- "HPD90rf.N123.t12.hdr.poly"
    cat("HPD90drf.N123.t12.hdr.poly =",HPD90drf.N123.t12.hdr.poly,"\n")
    col.target = col.target+1
    RESULTS[k,col.target]=HPD90drf.N123.t12.hdr.poly
    colnames(RESULTS)[col.target] <- "HPD90drf.N123.t12.hdr.poly" 
    
    ##############################################################################################
    ## HPD biparam N6.ra (param 3 and 9)----------- Methode "Paul" = hdr
    ##############################################################################################
    ##############################################################################################
    
    ########################### RF
    
    # PARALLELISATION: Methode rf-hdr-poly
    # Définition de la fonction pour calculer HPD et vérifier les points
    calcul_HPD_rf_hdr_poly <- function(i.pod) {
      hpd.levels <- hdrcde::hdr.2d(x = pod.weighted.sample.rf[[3]][[i.pod]]$N6, y = pod.weighted.sample.rf[[9]][[i.pod]]$ra, prob = c(0.1, 0.25, 0.5, 0.75, 0.9))
      rf.hdr <- contourLines(x = hpd.levels$den$x, y = hpd.levels$den$y,z = hpd.levels$den$z,levels = hpd.levels$falpha[1])
      answer = sp::point.in.polygon(param.pods[i.pod,"N6"], param.pods[i.pod,"ra"], rf.hdr[[1]]$x, rf.hdr[[1]]$y)
      if (answer == 1) {
        return(1)
      } else {
        return(0)
      }
    }
    
    # Exporter les objets nécessaires au cluster
    clusterExport(cl, varlist = c("calcul_HPD_rf_hdr_poly", "pod.weighted.sample.rf", "param.pods"))
    # Exécuter en parallèle
    results <- parLapply(cl, 1:n.pods, calcul_HPD_rf_hdr_poly)
    # Calculer le total et nettoyage
    HPD90rf.N6.ra.hdr.poly <- sum(unlist(results)) / n.pods
    gc()
    
    ########################### DRF
    
    # PARALLELISATION Methode drf-hdr-poly
    calcul_HPD_drf_hdr_poly <- function(i.pod) {
      hpd.levels <- hdrcde::hdr.2d(x = pod.weighted.sample.drf[[i.pod]]$N6, y = pod.weighted.sample.drf[[i.pod]]$ra, prob = c(0.1, 0.25, 0.5, 0.75, 0.9))
      drf.hdr <- contourLines(x = hpd.levels$den$x, y = hpd.levels$den$y, z = hpd.levels$den$z, levels = hpd.levels$falpha[1])
      answer = sp::point.in.polygon(param.pods[i.pod,"N6"], param.pods[i.pod,"ra"], drf.hdr[[1]]$x, drf.hdr[[1]]$y)
      if (answer == 1) {
        return(1/n.pods)
      } else {
        return(0)
      }
    }
    
    # Exporter les objets nécessaires au cluster
    clusterExport(cl, varlist = c("calcul_HPD_drf_hdr_poly", "pod.weighted.sample.drf", "param.pods", "n.pods"))
    # Exécuter en parallèle
    results <- parLapply(cl, 1:n.pods, calcul_HPD_drf_hdr_poly)
    # Calculer le total
    HPD90drf.N6.ra.hdr.poly <- sum(unlist(results))
    stopCluster(cl)
    gc()
    
    # Results: visualisation and saving: methode hdr (poly & near)
    cat("HPD90rf.N6.ra.hdr.poly =",HPD90rf.N6.ra.hdr.poly,"\n")
    col.target = col.target+1
    RESULTS[k,col.target]=HPD90rf.N6.ra.hdr.poly
    colnames(RESULTS)[col.target] <- "HPD90rf.N6.ra.hdr.poly"
    cat("HPD90drf.N6.ra.hdr.poly =",HPD90drf.N6.ra.hdr.poly,"\n")
    col.target = col.target+1
    RESULTS[k,col.target]=HPD90drf.N6.ra.hdr.poly
    
    #################################################################
    ####################### START COMPUTATION OF "ENERGY" SCORE #####  
    #################################################################
    # JEFF: below = code for various energy scores for RF and for DRF = es_sample or vs_sample (with p=0.5 = default or p=1.0 ) or mmds_sample
    # JEFF: different ways of standardisations available cf. 2 standardisations possibles pour objet1 et objet2
    
    ####################### RF ###################################################################
    cat("Computing energy score RF","\n")
    
    ##################### Creer un objet2.rf correctement configuré ####################################
    # Objet 2 = posterior distributions
    objet2.rf = vector("list", n.pods)
    objet2.rf<- lapply(objet2.rf, function(x) NULL) # Vider la liste
    # single.pod.all.param.weighted.sample.rf = matrix(data = NA, nrow = length(param.list), ncol = sample_size)
    # rownames(single.pod.all.param.weighted.sample.rf) = colnames(PARAMS)
    # # Création d'un sous-objet pod.weighted.#note: sample.without.weights.col.rf sans le composant $weights
    #note: pod.weighted.sample.rf[[i.param]][[i.pod]]$nomparam
    #head(pod.weighted.sample.rf[[1]][[1000]])
    #str(pod.weighted.sample.rf)
    pod.weighted.sample.without.weights.col.rf <- lapply(pod.weighted.sample.rf, function(x) {
      x <- lapply(x, function(df) {
        df$weights <- NULL  # Suppression de la colonne weights
        return(df)
      })
      return(x)
    })
    #note: pod.weighted.sample.without.weights.col.rf[[i.param]][[i.pod]]$nomparam
    #head(pod.weighted.sample.without.weights.col.rf[[10]][[100]])
    
    # Initialisation de l'objet2 comme une liste vide
    objet2.rf <- vector("list", length = length(pod.weighted.sample.without.weights.col.rf[[1]]))
    # Boucle pour remplir objet2 avec les matrices nécessaires
    for (i in seq_along(objet2.rf)) {
      # Initialiser une liste pour stocker les vecteurs de chaque paramètre
      param_vectors <- vector("list", ncol(param.pods))
      # Boucle sur chaque paramètre
      for (j in seq_len(ncol(param.pods))) {
        # Extraire les données pour chaque simulation ipod et chaque paramètre
        param_vectors[[j]] <- t(pod.weighted.sample.without.weights.col.rf[[j]][[i]])
      }
      # Création de la matrice pour la simulation ipod
      # et ajout à objet2. Chaque ligne contient les valeurs pour un paramètre
      objet2.rf[[i]] <- do.call(rbind, param_vectors)
    }
    # str(objet2.rf[[100]])
    # head(objet2.rf[[100]])
    
    ##################################### START STANDARDISATIONS OBJET1 ET OBJET2 RF ####################################
    # Facteurs standardisation issus de objet2.rf et objet1.rf
    
    ######## Calcul des facteurs de standardise issus des i.pod posterior distrib de objets2.rf
    # Calcul du nombre de paramètres et initialisation du data.frame pour les facteurs de standardise
    nombre_parametres <- ncol(param.pods)
    noms_parametres <- rownames(objet2.rf[[1]])
    # Initialiser un data.frame pour stocker les facteurs de standardise
    standardise.factors.objet2.rf <- data.frame(matrix(nrow = n.pods, ncol = nombre_parametres * 2))
    #str(standardise.factors.objet2.rf)
    # Nommer les colonnes du data.frame
    colnames(standardise.factors.objet2.rf) <- c(rbind(paste0("mean_", noms_parametres), paste0("sd_", noms_parametres)))
    # Boucle pour remplir le data.frame avec les facteurs de standardise pour chaque simulation et chaque paramètre
    for (i.pod in 1:n.pods) {
      mat <- objet2.rf[[i.pod]]
      for (j in 1:nombre_parametres) {
        param_name <- noms_parametres[j]
        # Calculer et stocker la moyenne et l'écart type pour chaque paramètre
        standardise.factors.objet2.rf[i.pod, paste0("mean_", param_name)] <- mean(mat[j, ])
        standardise.factors.objet2.rf[i.pod, paste0("sd_", param_name)] <- sd(mat[j, ])
      }
    }
    #head(standardise.factors.objet2.rf)
    #dim(standardise.factors.objet2.rf)
    
    ###### Calcul des facteurs de standardisation issus des i.pod de param.pods = objet.1
    standardise.factors.param.pods = matrix(nrow=2, ncol=ncol(param.pods))
    colnames(standardise.factors.param.pods)=colnames(param.pods)
    rownames(standardise.factors.param.pods)=c("mean","sd")
    for (i.param in 1:ncol(param.pods)) {
      standardise.factors.param.pods[1,i.param] = mean(param.pods[,i.param])
      standardise.factors.param.pods[2,i.param] = sd(param.pods[,i.param])
    }
    
    ## 2 standardisations possibles pour objet1 et objet2
    
    ###### objet1.standardise.SFO1.rf = standardisation directe avec mean et sd de param.pods = objet1
    objet1.rf = param.pods
    objet1.standardise.SFO1.rf <- as.data.frame(scale(objet1.rf))
    # head(objet1.standardise.rf)
    # str(objet1.standardise.rf)
    #summary(objet1.standardise.rf)
    
    ##### objet1.standardise.SFO2.rf = standardise des colonnes de objet1 a partir des facteurs de standardises de objet2.rf
    # #head(objet1.rf)
    objet1.standardise.SFO2.rf <- matrix(0, nrow = n.pods, ncol = ncol(param.pods))
    colnames(objet1.standardise.SFO2.rf)=paste0(colnames(param.pods),".std")
    for (i.pod in 1:n.pods) {
      for (i.param in 1:ncol(objet1.rf)) {
        objet1.standardise.SFO2.rf[i.pod,i.param] = (objet1.rf[i.pod, i.param]-standardise.factors.objet2.rf[i.pod, 2*i.param-1])/standardise.factors.objet2.rf[i.pod, 2*i.param]
      }
    }
    
    ###### objet2.standardise.SFO2.rf = standardisation direct avec mean et sd de objet2
    objet2.standardise.SFO2.rf <- lapply(objet2.rf, function(mat) {
      # Application de la standardise à chaque ligne de la matrice
      apply(mat, 1, function(row) scale(row))
    })
    for (i.pod in 1:n.pods) {objet2.standardise.SFO2.rf[[i.pod]] = t(objet2.standardise.SFO2.rf[[i.pod]])}
    # head(objet2.rf[[1000]])
    # str(objet2.standardise.rf[[1000]])
    # head(objet2.standardise.rf[[1000]])
    # summary(objet2.standardise.rf[[1000]])
    
    ###### objet2.standardise.SFO1.rf = standardisation avec mean et sd de objet1 = param.pods
    objet2.standardise.SFO1.rf <- objet2.rf
    #head(objet2.rf[[1000]])
    for (i.pod in 1:n.pods) {
      for (i.param in 1:ncol(param.pods)) {
        objet2.standardise.SFO1.rf[[i.pod]][i.param,] = (objet2.rf[[i.pod]][i.param,]-standardise.factors.param.pods[1,i.param])/standardise.factors.param.pods[2,i.param]
      }
    }
    # head(objet2.standardise.rf[[1000]])
    # str(objet2.standardise.rf)
    # summary(objet2.standardise.rf[[1000]])
    
    ##################################### COMPUTE SCORE RF ####################################
    # JEFF: here is the code for computing RF scores for raw data (objects) = score.rf.ipod or standardize data (objects) = score.standardise.rf.ipod
    # JEFF: note that if you initially (before DRF and RF) have already standardise the data then score.rf.ipod = standardized score !
    # JEFF: and score.standardise.rf.ipod is then actually "doubled-standardised" which makes little sense so ignore score.standardise.rf.ipod in this case
    ######### CHOIX STANDARDISATIONS
    objet1.standardise.rf = objet1.standardise.SFO1.rf
    #objet1.standardise.rf = objet1.standardise.SFO2.rf
    #objet2.standardise.rf = objet2.standardise.SFO2.rf
    objet2.standardise.rf = objet2.standardise.SFO1.rf
    ##########################################################################################
    
    score.rf.ipod <- vector("numeric", length = n.pods)
    score.rf.ipod <- rep(0, n.pods)
    score.rf.ipod = sapply(1:n.pods, function(i.pod){es_sample(y = as.numeric(objet1.rf[i.pod,]), dat = objet2.rf[[i.pod]]  ) })
    #score.rf.ipod = sapply(1:n.pods, function(i.pod){vs_sample(y = as.numeric(objet1.rf[i.pod,3]), dat = objet2.uniparam.rf[[i.pod]]  ) })
    #score.rf.ipod = sapply(1:n.pods, function(i.pod){mmds_sample(y = as.numeric(objet1.rf[i.pod,3]), dat = objet2.uniparam.rf[[i.pod]]  ) })
    mean.score.rf<-mean(score.rf.ipod)
    sd.score.rf<-sd(score.rf.ipod)
    
    score.standardise.rf.ipod <- vector("numeric", length = n.pods)
    score.standardise.rf.ipod <- rep(0, n.pods)
    score.standardise.rf.ipod = sapply(1:n.pods, function(i.pod){es_sample(y = as.numeric(objet1.standardise.rf[i.pod,]), dat = objet2.standardise.rf[[i.pod]]  ) })
    #score.standardise.rf.ipod = sapply(1:n.pods, function(i.pod){vs_sample(y = as.numeric(objet1.standardise.rf[i.pod,]), dat = objet2.standardise.rf[[i.pod]]  ) })
    #score.standardise.rf.ipod = sapply(1:n.pods, function(i.pod){mmds_sample(y = as.numeric(objet1.standardise.rf[i.pod,]), dat = objet2.standardise.rf[[i.pod]]  ) })
    mean.score.standardise.rf<-mean(score.standardise.rf.ipod)
    sd.score.standardise.rf<-sd(score.standardise.rf.ipod)
    
    ####################### DRF ###################################################################
    cat("Computing energy score DRF","\n")
    
    ##################### Creer un objet2.drf correctement configuré ####################################
    # Objet 2: posterior samples des pods
    # Supprimer la colonne "weight" de pod.weighted.sample.drf
    pod.weighted.sample.drf.without.weights <- lapply(pod.weighted.sample.drf, function(df) {
      df[, !(names(df) %in% "weights")]
    })
    # Transformer en matrix et standardiser les valeurs
    objet2.drf=pod.weighted.sample.drf.without.weights
    # str(objet2.drf)
    # str(objet2.drf[[1000]])
    # head(objet2.drf[[1000]])
    
    ##################################### START STANDARDISATIONS OBJET1 ET OBJET2 DRF ####################################
    # Facteurs standardisation issus de objet2.drf et objet1.drf
    
    # Calcul des facteurs de standardise issus des i.pod posterior distrib de objets2.drf
    # Calcul du nombre de paramètres et initialisation du data.frame pour les facteurs de standardise
    nombre_parametres <- ncol(param.pods)
    noms_parametres <- colnames(param.pods)
    # Initialiser un data.frame pour stocker les facteurs de standardise
    standardise.factors.objet2.drf <- data.frame(matrix(nrow = n.pods, ncol = nombre_parametres * 2))
    # Nommer les colonnes du data.frame
    colnames(standardise.factors.objet2.drf) <- c(rbind(paste0("mean_", noms_parametres), paste0("sd_", noms_parametres)))
    # Boucle pour remplir le data.frame avec les facteurs de standardise pour chaque simulation et chaque paramètre
    for (i.pod in 1:n.pods) {
      mat <- objet2.drf[[i.pod]]
      for (j in 1:nombre_parametres) {
        param_name <- noms_parametres[j]
        # Calculer et stocker la moyenne et l'écart type pour chaque paramètre
        standardise.factors.objet2.drf[i.pod, paste0("mean_", param_name)] <- mean(mat[, j])
        standardise.factors.objet2.drf[i.pod, paste0("sd_", param_name)] <- sd(mat[, j])
      }
    }
    #head(standardise.factors.objet2.drf)
    
    ###### Calcul des facteurs de standardisation issus des i.pod de param.pods = objet.1
    standardise.factors.param.pods = matrix(nrow=2, ncol=ncol(param.pods))
    colnames(standardise.factors.param.pods)=colnames(param.pods)
    rownames(standardise.factors.param.pods)=c("mean","sd")
    for (i.param in 1:ncol(param.pods)) {
      standardise.factors.param.pods[1,i.param] = mean(param.pods[,i.param])
      standardise.factors.param.pods[2,i.param] = sd(param.pods[,i.param])
    }
    
    ## 2 standardisations possibles pour objet1.drf et objet2.drf
    
    ###### objet1.standardise.SFO1.drf = standardisation directe avec mean et sd de param.pods = objet1.drf
    objet1.drf = param.pods
    objet1.standardise.SFO1.drf <- as.data.frame(scale(objet1.drf))
    # head(objet1.standardise.SFO1.drf)
    # str(objet1.standardise.SFO1.drf)
    #summary(objet1.standardise.SFO1.drf)
    
    ###### objet1.standardise.SFO2.drf = standardisation avec mean et sd de objet2.drf
    # #head(objet1.drf)
    objet1.standardise.SFO2.drf <- matrix(0, nrow = n.pods, ncol = ncol(param.pods))
    colnames(objet1.standardise.SFO2.drf)=paste0(colnames(param.pods),".std")
    for (i.pod in 1:n.pods) {
      for (i.param in 1:ncol(objet1.drf)) {
        objet1.standardise.SFO2.drf[i.pod,i.param] = (objet1.drf[i.pod, i.param]-standardise.factors.objet2.drf[i.pod, 2*i.param-1])/standardise.factors.objet2.drf[i.pod, 2*i.param]
      }
    }
    
    ###### objet2.standardise.SFO2.drf = standardisation direct avec mean et sd de objet2.drf
    #head(objet2.drf[[1000]])
    # Copie de l'objet original pour conserver la structure
    objet2.standardise.SFO2.drf <- objet2.drf
    # Parcourir chaque data.frame dans la liste et standardiser chaque colonne
    for (i in seq_along(objet2.drf)) {
      # Application de la fonction scale à chaque colonne du data.frame
      objet2.standardise.SFO2.drf[[i]] <- as.data.frame(lapply(objet2.drf[[i]], scale))
    }
    
    ###### objet2.standardise.SFO1.drf = standardisation avec mean et sd de objet1.drf = param.pods
    objet2.standardise.SFO1.drf <- objet2.drf
    #head(objet2.drf[[1000]])
    for (i.pod in 1:n.pods) {
      for (i.param in 1:ncol(param.pods)) {
        objet2.standardise.SFO1.drf[[i.pod]][,i.param] = (objet2.drf[[i.pod]][,i.param]-standardise.factors.param.pods[1,i.param])/standardise.factors.param.pods[2,i.param]
      }
    }
    # head(objet2.standardise.SFO1.drf[[1000]])
    # str(objet2.standardise.SFO1.drf)
    # summary(objet2.standardise.SFO1.drf[[1000]])
    
    ##################################### COMPUTE SCORE DRF ####################################
    # JEFF: here is the code for computing DRF scores for raw data (objects) = score.drf.ipod or standardize data (objects) = score.standardise.drf.ipod
    # JEFF: note that if you initially (before DRF and RF) have already standardise the data then score.drf.ipod = standardized score !
    # JEFF: and score.standardise.rf.ipod is then actually "doubled-standardised" which makes little sense so ignore score.standardise.drf.ipod in this case
    ######### CHOIX STANDARDISATIONS
    objet1.standardise.drf = objet1.standardise.SFO1.drf
    #objet1.standardise.drf = objet1.standardise.SFO2.drf
    #objet2.standardise.drf = objet2.standardise.SFO2.drf
    objet2.standardise.drf = objet2.standardise.SFO1.drf
    ##########################################################################################
    
    score.drf.ipod <- vector("numeric", length = n.pods)
    score.drf.ipod <- rep(0, n.pods)
    score.drf.ipod = sapply(1:n.pods, function(i.pod){es_sample(y = as.numeric(objet1.drf[i.pod,]), dat = t(objet2.drf[[i.pod]] ) ) })
    #score.drf.ipod = sapply(1:n.pods, function(i.pod){vs_sample(y = as.numeric(objet1.drf[i.pod,]), dat = t(objet2.drf[[i.pod]] ) ) })
    #score.drf.ipod = sapply(1:n.pods, function(i.pod){mmds_sample(y = as.numeric(objet1.drf[i.pod,]), dat = t(objet2.drf[[i.pod]] ) ) })
    mean.score.drf<-mean(score.drf.ipod)
    sd.score.drf<-sd(score.drf.ipod)
    
    score.standardise.drf.ipod <- vector("numeric", length = n.pods)
    score.standardise.drf.ipod <- rep(0, n.pods)
    score.standardise.drf.ipod = sapply(1:n.pods, function(i.pod){es_sample(y = as.numeric(objet1.standardise.drf[i.pod,]), dat = t(objet2.standardise.drf[[i.pod]])  ) })
    #score.standardise.drf.ipod = sapply(1:n.pods, function(i.pod){vs_sample(y = as.numeric(objet1.standardise.drf[i.pod,]), dat = t(objet2.standardise.drf[[i.pod]])  ) })
    #score.standardise.drf.ipod = sapply(1:n.pods, function(i.pod){mmds_sample(y = as.numeric(objet1.standardise.drf[i.pod,]), dat = t(objet2.standardise.drf[[i.pod]])  ) })
    
    mean.score.standardise.drf<-mean(score.standardise.drf.ipod)
    sd.score.standardise.drf<-sd(score.standardise.drf.ipod)
    
    # ## Some stats
    summary(score.rf.ipod)
    summary(score.drf.ipod)
    summary(score.standardise.rf.ipod)
    summary(score.standardise.drf.ipod)
    # 
    # # Figures simples energy scores
    # # 1/ Création d'un data frame avec les data
    # data <- data.frame(score_standardise_drf_ipod = score.standardise.drf.ipod,
    #                    score_standardise_rf_ipod = score.standardise.rf.ipod)
    # ggplot(data, aes(x = score_standardise_drf_ipod, y = score_standardise_rf_ipod)) +
    #   geom_point(alpha = 0.5) +  # Ajouter les points de données avec une certaine transparence
    #   geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +  # Ajouter la ligne y=x
    #   labs(title = "        Complex population genetics model (9 original variable parameters)
    #        Standardise parameter values
    #        n.tree=1000, n.train=10000, npods=1000 test datasets",
    #        x = "Standardised Energy Scores DRF",
    #        y = "Standardised Energy Scores RF") +
    #   theme_minimal()  # Utiliser un thème minimal pour le graphique  
    # 
    # # 2/ Création d'un data frame avec les data
    # data <- data.frame(score_drf_ipod = score.drf.ipod,
    #                    score_rf_ipod = score.rf.ipod)
    # ggplot(data, aes(x = score_drf_ipod, y = score_rf_ipod)) +
    #   geom_point(alpha = 0.5) +  # Ajouter les points de données avec une certaine transparence
    #   geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +  # Ajouter la ligne y=x
    #   labs(title = "        Complex population genetics model (9 original variable parameters)
    #        Raw parameter values
    #        n.tree=1000, n.train=10000, npods=1000 test datasets",
    #        x = "Energy Scores DRF",
    #        y = "Energy Scores RF") +
    #   theme_minimal()  # Utiliser un thème minimal pour le graphique
    
    ## Affichage et sauvegarde energy score
    cat("mean.score.rf =",mean.score.rf,"\n")
    col.target = col.target+1
    RESULTS[k,col.target]=mean.score.rf
    colnames(RESULTS)[col.target] <- "mean.score.rf"
    cat("sd.score.rf =",sd.score.rf,"\n")
    col.target = col.target+1
    RESULTS[k,col.target]=sd.score.rf
    colnames(RESULTS)[col.target] <- "sd.score.rf"
    cat("mean.score.drf =",mean.score.drf,"\n")
    col.target = col.target+1
    RESULTS[k,col.target]=mean.score.drf
    colnames(RESULTS)[col.target] <- "mean.score.drf"
    cat("sd.score.drf =",sd.score.drf,"\n")
    col.target = col.target+1
    RESULTS[k,col.target]=sd.score.drf
    colnames(RESULTS)[col.target] <- "sd.score.drf"
    
    cat("mean.score.standardise.rf =",mean.score.standardise.rf,"\n")
    col.target = col.target+1
    RESULTS[k,col.target]=mean.score.standardise.rf
    colnames(RESULTS)[col.target] <- "mean.score.standardise.rf"
    cat("sd.score.standardise.rf =",sd.score.standardise.rf,"\n")
    col.target = col.target+1
    RESULTS[k,col.target]=sd.score.standardise.rf
    colnames(RESULTS)[col.target] <- "sd.score.standardise.rf"
    cat("mean.score.standardise.drf =",mean.score.standardise.drf,"\n")
    col.target = col.target+1
    RESULTS[k,col.target]=mean.score.standardise.drf
    colnames(RESULTS)[col.target] <- "mean.score.standardise.drf"
    cat("sd.score.standardise.drf =",sd.score.standardise.drf,"\n")
    col.target = col.target+1
    RESULTS[k,col.target]=sd.score.standardise.drf
    colnames(RESULTS)[col.target] <- "sd.score.standardise.drf"
    
    ### Code for ENERGY SCORES ala AE
    # JEFF: here is my second personnal code for computing scores for raw data (objects) 
    # JEFF: note that if you initially (before DRF and RF) have already standardise the data then the scores correspond to that of the corresponding standardized data
    
    # ######## RF 
    # Initialisation d'un vecteur pour stocker les scores d'énergie pour chaque pod
    n.pods <- nrow(param.pods) 
    energy.scores.AE.ipod.rf <- numeric(n.pods) # Vecteur pour stocker les scores d'énergie individuels
    # Boucle sur chaque pod pour calculer le score d'énergie
    for (i.pod in 1:n.pods) {
      observations_standardized <- as.numeric(objet1.rf[i.pod, ])  # Observations pour le pod actuel
      predictions_standardized <- as.matrix(t(objet2.rf[[i.pod]]))  # Prédictions pour le pod actuel
      # observations_standardized <- as.numeric(objet1.standardise.rf[i.pod, ])  # Observations pour le pod actuel
      # predictions_standardized <- as.matrix(t(objet2.standardise.rf[[i.pod]]))  # Prédictions pour le pod actuel
      # Calcul du score d'énergie pour le pod actuel
      es_value <- energy_score(observations_standardized, predictions_standardized)
      energy.scores.AE.ipod.rf[i.pod] <- es_value # Stocker le score d'énergie pour le pod actuel
    }
    # Calculer le score d'énergie moyen sur tous les pods
    mean.energy.score.AE.rf <- mean(energy.scores.AE.ipod.rf)
    
    # ######## DRF 
    # Initialisation d'un vecteur pour stocker les scores d'énergie pour chaque pod
    n.pods <- nrow(param.pods) 
    energy.scores.AE.ipod.drf <- numeric(n.pods) # Vecteur pour stocker les scores d'énergie individuels
    # Boucle sur chaque pod pour calculer le score d'énergie
    for (i.pod in 1:n.pods) {
      observations_standardized <- as.numeric(objet1.drf[i.pod, ])  # Observations pour le pod actuel
      predictions_standardized <- as.matrix(objet2.drf[[i.pod]])  # Prédictions pour le pod actuel
      # observations_standardized <- as.numeric(objet1.standardise.drf[i.pod, ])  # Observations pour le pod actuel
      # predictions_standardized <- as.matrix(objet2.standardise.drf[[i.pod]])  # Prédictions pour le pod actuel
      # Calcul du score d'énergie pour le pod actuel
      es_value <- energy_score(observations_standardized, predictions_standardized)
      energy.scores.AE.ipod.drf[i.pod] <- es_value # Stocker le score d'énergie pour le pod actuel
    }
    # Calculer le score d'énergie moyen sur tous les pods
    mean.energy.score.AE.drf <- mean(energy.scores.AE.ipod.drf)
    
    ## Some stats
    summary(energy.scores.AE.ipod.rf)
    summary(energy.scores.AE.ipod.drf)
    
    # Figure simples energy scores AE
    data <- data.frame(score_standardise_drf_ipod = energy.scores.AE.ipod.drf,
                       score_standardise_rf_ipod = energy.scores.AE.ipod.rf)
    ggplot(data, aes(x = score_standardise_drf_ipod, y = score_standardise_rf_ipod)) +
      geom_point(alpha = 0.5) +  # Ajouter les points de données avec une certaine transparence
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +  # Ajouter la ligne y=x
      labs(title = "        Complex population genetics model (9 original variable parameters)
           Standardise parameter values
           n.tree=1000, n.train=10000, npods=1000 test datasets",
           x = "AE Standardised Energy Scores DRF",
           y = "AE Standardised Energy Scores RF") +
      theme_minimal()  # Utiliser un thème minimal pour le graphique
    
    ## Affichage et sauvegarde energy score AE
    cat("mean.energy.score.AE.rf =",mean.energy.score.AE.rf,"\n")
    col.target = col.target+1
    RESULTS[k,col.target]=mean.score.rf
    colnames(RESULTS)[col.target] <- "mean.energy.score.AE.rf"
    cat("mean.energy.score.AE.drf =",mean.energy.score.AE.drf,"\n")
    col.target = col.target+1
    RESULTS[k,col.target]=mean.score.rf
    colnames(RESULTS)[col.target] <- "mean.energy.score.AE.drf"
    
    # ##### VARIANCE CO-VARIANCE ##########################
    # head(objet1.drf)
    # prior_cov = cov(objet1.drf)
    # for (i.pod in seq_along(objet2.rf)) {post_cov[[i.pod]] <- cov(t(objet2.rf[[i.pod]]))}
    # # Calculer la matrice de variance-covariance pour objet1.drf
    # var_cov_matrix <- cov(objet2.standardise.drf[[1000]])
    # var_cov_matrix <- cov(param.pods)
    # # Afficher la matrice de variance-covariance
    # print(var_cov_matrix)
    # ###############################
    
    ####### Adapted MAHALANOBIS distance pour comparer les distributions a posteriori par rapport aux valeurs de test ####
    # JEFF: here is the code for computing the Adapted MAHALANOBIS distance for raw data = mahalanobis.distances.raw.ipod or standardise data = mahalanobis.distances.standardise.ipod
    # JEFF: note that if you initially (before DRF and RF) have already standardise the data then mahalanobis.distances.raw.ipod = standardized distance !
    # JEFF: and mahalanobis.distances.standardise.ipod is then actually "doubled-standardised" which makes little sense so ignore mahalanobis.distances.standardise.ipod in this case
    # Supposons que 'objet' soit un dataframe ou une matrice où chaque colonne représente les échantillons pondérés pour un paramètre
    # et chaque ligne un échantillon de ces paramètres.
    
    ### RF - Raw
    # Construire param.mean.ipod.pred.raw.rf
    param.mean.ipod.pred.rf <- as.data.frame(matrix(NA, nrow = ncol(param.pods), ncol = n.pods))
    rownames(param.mean.ipod.pred.rf) <- colnames(param.pods)
    #head(param.mean.ipod.pred.rf)
    #str(param.mean.ipod.pred.rf)
    #dim(param.mean.ipod.pred.rf)
    for (i.param in 1:ncol(param.pods)) param.mean.ipod.pred.rf[i.param,] = param.pred.rf[[i.param]]$expectation
    param.mean.ipod.pred.raw.rf= t(param.mean.ipod.pred.rf)
    # head(param.mean.ipod.pred.raw.rf)
    # str(param.mean.ipod.pred.raw.rf)
    # dim(param.mean.ipod.pred.raw.rf)
    # summary(param.mean.ipod.pred.raw.rf)
    
    # Calcul de la matrice de covariance postérieure (i.e. des échantillons pondérés) a partir des post-samples obtenus par bootstrap
    post_cov <- vector("list", n.pods)
    # Remplir chaque élément de 'post_cov' avec NA
    for(i in 1:n.pods) {
      post_cov[[i]] <- NA  # Ceci place un simple NA dans chaque élément de la liste
    }
    # head(t(objet2.rf[[1]]))
    # dim(t(objet2.rf[[1]]))
    # summary(t(objet2.rf[[1]]))
    for (i.pod in seq_along(objet2.rf)) {post_cov[[i.pod]] <- cov(t(objet2.rf[[i.pod]]))}
    #post_cov[[1]]
    #param.pod = dataframe de valeurs vraies
    # mean.posterior.estimation = dataframe des moyennes a posteriori
    # post_cov_list = liste de i.pod matrices de covariance
    #library(MASS) # Pour la fonction mahalanobis
    # Initialisation d'un vecteur pour stocker les distances de Mahalanobis
    mahalanobis.distances.raw.ipod.rf <- numeric(length = n.pods) 
    mahalanobis.distances.raw.ipod.rf <- rep(0, n.pods)
    # Boucle sur chaque ensemble de paramètres
    for (i.pod in 1:n.pods) {   
      # Sélection des valeurs vraies et des moyennes a posteriori pour i.pod (sous forme de vecteurs !!!)
      true.param.values <- as.vector(t(objet1.drf[i.pod, ]))
      post.param.means <- as.vector(param.mean.ipod.pred.raw.rf[i.pod,])
      # Calcul de la distance de Mahalanobis pour i.pod
      mahalanobis.distances.raw.ipod.rf[i.pod] <- mahalanobis(true.param.values, center = post.param.means, cov = post_cov[[i.pod]])
    }
    # Calcul de la métrique globale comme la moyenne des distances de Mahalanobis RF
    mean.mahalanobis.distance.raw.rf <- mean(mahalanobis.distances.raw.ipod.rf) # On peut aussi calculer un sd
    sd.mahalanobis.distance.raw.rf <- sd(mahalanobis.distances.raw.ipod.rf)
    mean.mahalanobis.distance.raw.rf
    sd.mahalanobis.distance.raw.rf
    
    ### DRF - raw
    #head(objet2.drf[[1000]])
    # Calcul de la matrice de covariance postérieure (i.e. des échantillons pondérés) a partir des post-samples obtenus par bootstrap
    post_cov <- vector("list", n.pods)
    # Remplir chaque élément de 'post_cov' avec NA
    for(i in 1:n.pods) {
      post_cov[[i]] <- NA  # Ceci place un simple NA dans chaque élément de la liste
    }
    for (i.pod in seq_along(objet2.drf)) post_cov[[i.pod]] <- cov(objet2.drf[[i.pod]])
    #post_cov[[1]]
    #param.pod = dataframe de valeurs vraies
    # mean.posterior.estimation = dataframe des moyennes a posteriori
    # post_cov_list = liste de i.pod matrices de covariance
    #library(MASS) # Pour la fonction mahalanobis
    # Initialisation d'un vecteur pour stocker les distances de Mahalanobis
    mahalanobis.distances.raw.ipod.drf <- numeric(length = n.pods) 
    mahalanobis.distances.raw.ipod.drf <- rep(0, n.pods)
    # Boucle sur chaque vecteur de paramètres
    for (i.pod in 1:n.pods) {   
      # Sélection des valeurs vraies et des moyennes a posteriori pour i.pod (sous forme de vecteurs !!!)
      true.param.values <- as.vector(t(param.pods[i.pod, ]))
      post.param.means <- as.vector(t(pred.drf.mean[i.pod,]))
      # Calcul de la distance de Mahalanobis pour i.pod
      mahalanobis.distances.raw.ipod.drf[i.pod] <- mahalanobis(true.param.values, center = post.param.means, cov = post_cov[[i.pod]])
    }
    # Calcul de la métrique globale comme la moyenne des distances de Mahalanobis DRF
    mean.mahalanobis.distance.raw.drf <- mean(mahalanobis.distances.raw.ipod.drf)
    sd.mahalanobis.distance.raw.drf <- sd(mahalanobis.distances.raw.ipod.drf)
    mean.mahalanobis.distance.raw.drf
    sd.mahalanobis.distance.raw.drf
    
    ### RF - Standardise
    # Construire param.mean.ipod.pred.standardise.rf
    param.mean.ipod.pred.rf <- as.data.frame(matrix(NA, nrow = ncol(param.pods), ncol = n.pods))
    rownames(param.mean.ipod.pred.rf) <- colnames(param.pods)
    #head(param.mean.ipod.pred.rf)
    #str(param.mean.ipod.pred.rf)
    #dim(param.mean.ipod.pred.rf)
    for (i.param in 1:ncol(param.pods)) param.mean.ipod.pred.rf[i.param,] = param.pred.rf[[i.param]]$expectation
    param.mean.ipod.pred.rf= t(param.mean.ipod.pred.rf)
    #head(param.mean.ipod.pred.rf)
    param.mean.ipod.pred.standardise.rf = as.data.frame(scale(param.mean.ipod.pred.rf))
    # head(param.mean.ipod.pred.standardise.rf)
    # str(param.mean.ipod.pred.standardise.rf)
    # dim(param.mean.ipod.pred.standardise.rf)
    # summary(param.mean.ipod.pred.standardise.rf)
    
    # Calcul de la matrice de covariance postérieure (i.e. des échantillons pondérés) a partir des post-samples obtenus par bootstrap
    post_cov <- vector("list", n.pods)
    # Remplir chaque élément de 'post_cov' avec NA
    for(i in 1:n.pods) {
      post_cov[[i]] <- NA  # Ceci place un simple NA dans chaque élément de la liste
    }
    # head(t(objet2.rf[[1]]))
    # dim(t(objet2.rf[[1]]))
    # summary(t(objet2.rf[[1]]))
    for (i.pod in seq_along(objet2.standardise.rf)) {post_cov[[i.pod]] <- cov(t(objet2.standardise.rf[[i.pod]]))}
    #post_cov[[1]]
    #param.pod = dataframe de valeurs vraies
    # mean.posterior.estimation = dataframe des moyennes a posteriori
    # post_cov_list = liste de i.pod matrices de covariance
    #library(MASS) # Pour la fonction mahalanobis
    # Initialisation d'un vecteur pour stocker les distances de Mahalanobis
    mahalanobis.distances.standardise.ipod.rf <- numeric(length = n.pods) 
    mahalanobis.distances.standardise.ipod.rf <- rep(0, n.pods)
    # Boucle sur chaque ensemble de paramètres
    for (i.pod in 1:n.pods) {   
      # Sélection des valeurs vraies et des moyennes a posteriori pour i.pod (sous forme de vecteurs !!!)
      true.param.values <- as.vector(t(objet1.standardise.rf[i.pod, ]))
      post.param.means <- as.vector(t(param.mean.ipod.pred.standardise.rf[i.pod,]))
      # Calcul de la distance de Mahalanobis pour i.pod
      mahalanobis.distances.standardise.ipod.rf[i.pod] <- mahalanobis(true.param.values, center = post.param.means, cov = post_cov[[i.pod]])
    }
    # Calcul de la métrique globale comme la moyenne des distances de Mahalanobis RF
    mean.mahalanobis.distance.standardise.rf <- mean(mahalanobis.distances.standardise.ipod.rf) # On peut aussi calculer un sd
    sd.mahalanobis.distance.standardise.rf <- sd(mahalanobis.distances.standardise.ipod.rf)
    mean.mahalanobis.distance.standardise.rf
    sd.mahalanobis.distance.standardise.rf
    
    ### DRF - Standardise
    #head(pred.drf.mean)
    #summary(pred.drf.mean)
    pred.drf.mean.standardise = as.data.frame(scale(pred.drf.mean))
    colnames(pred.drf.mean.standardise) = colnames(param.pods)
    #head(pred.drf.mean.standardise)
    #summary(pred.drf.mean.standardise)
    # Calcul de la matrice de covariance postérieure (i.e. des échantillons pondérés) a partir des post-samples obtenus par bootstrap
    post_cov <- vector("list", n.pods)
    # Remplir chaque élément de 'post_cov' avec NA
    for(i in 1:n.pods) {
      post_cov[[i]] <- NA  # Ceci place un simple NA dans chaque élément de la liste
    }
    #head(objet2.standardise.drf[[1]])
    for (i.pod in seq_along(objet2.drf)) {post_cov[[i.pod]] <- cov(objet2.standardise.drf[[i.pod]])}
    #head(post_cov[[1000]])
    # Supposons que 'param.pod' soit votre dataframe de valeurs vraies,
    # 'mean.posterior.estimation' votre dataframe des moyennes a posteriori,
    # et 'post_cov_list' votre liste de 1000 matrices de covariance.
    #library(MASS) # Pour la fonction mahalanobis
    # Initialisation d'un vecteur pour stocker les distances de Mahalanobis
    mahalanobis.distances.standardise.ipod.drf <- numeric(length = n.pods) 
    mahalanobis.distances.standardise.ipod.drf <- rep(0, n.pods)
    # Boucle sur chaque vecteur de paramètres
    for (i.pod in 1:n.pods) {   
      # Sélection des valeurs vraies et des moyennes a posteriori pour i.pod (sous forme de vecteurs !!!)
      true.param.values <- as.vector(t(objet1.standardise.drf[i.pod, ]))
      post.param.means <- as.vector(t(pred.drf.mean.standardise[i.pod,]))
      # Calcul de la distance de Mahalanobis pour i.pod
      mahalanobis.distances.standardise.ipod.drf[i.pod] <- mahalanobis(true.param.values, center = post.param.means, cov = post_cov[[i.pod]])
    }
    # Calcul de la métrique globale comme la moyenne des distances de Mahalanobis DRF
    mean.mahalanobis.distance.standardise.drf <- mean(mahalanobis.distances.standardise.ipod.drf) # On peut aussi calculer un sd
    sd.mahalanobis.distance.standardise.drf <- sd(mahalanobis.distances.standardise.ipod.drf)
    mean.mahalanobis.distance.standardise.drf
    sd.mahalanobis.distance.standardise.drf
    
    ## Some stats and figures
    summary(mahalanobis.distances.raw.ipod.rf)
    summary(mahalanobis.distances.raw.ipod.drf)
    summary(mahalanobis.distances.standardise.ipod.rf)
    summary(mahalanobis.distances.standardise.ipod.drf)
    
    # Figures simples mahalanobis
    # 1/ Création d'un data frame avec les data raw
    data <- data.frame(mahalanobis_standardise_drf_ipod = mahalanobis.distances.raw.ipod.drf,
                       mahalanobis_standardise_rf_ipod = mahalanobis.distances.raw.ipod.rf)
    ggplot(data, aes(x = mahalanobis_standardise_drf_ipod, y = mahalanobis_standardise_rf_ipod)) +
      geom_point(alpha = 0.5) +  # Ajouter les points de données avec une certaine transparence
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +  # Ajouter la ligne y=x
      labs(title = "        COMPLEX MODEL
           Raw parameter values
           n.tree=1000, n.train=10000, npods=1000 test datasets",
           x = "raw mahalanobis DRF",
           y = "raw mahalanobis RF") +
      theme_minimal()  # Utiliser un thème minimal pour le graphique
    
    # 2/ Création d'un data frame avec les data standardise
    data <- data.frame(mahalanobis_drf_ipod = mahalanobis.distances.standardise.ipod.drf,
                       mahalanobis_rf_ipod = mahalanobis.distances.standardise.ipod.rf)
    ggplot(data, aes(x = mahalanobis_drf_ipod, y = mahalanobis_rf_ipod)) +
      geom_point(alpha = 0.5) +  # Ajouter les points de données avec une certaine transparence
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +  # Ajouter la ligne y=x
      labs(title = "        COMPLEX MODEL
           Standardise parameter values
           n.tree=1000, n.train=10000, npods=1000 test datasets",
           x = "standardise mahalanobis DRF",
           y = "standardise mahalanobis RF") +
      theme_minimal()  # Utiliser un thème minimal pour le graphique
    
    ## Affichage et sauvegarde mahalanobis
    cat("mean.mahalanobis.distance.raw.rf =",mean.mahalanobis.distance.raw.rf,"\n")
    col.target = col.target+1
    RESULTS[k,col.target]=mean.mahalanobis.distance.raw.rf
    colnames(RESULTS)[col.target] <- "mean.mahalanobis.distance.raw.rf"
    
    cat("sd.mahalanobis.distance.raw.rf =",sd.mahalanobis.distance.raw.rf,"\n")
    col.target = col.target+1
    RESULTS[k,col.target]=sd.mahalanobis.distance.raw.rf
    colnames(RESULTS)[col.target] <- "sd.mahalanobis.distance.raw.drf"
    
    cat("mean.mahalanobis.distance.raw.drf =",mean.mahalanobis.distance.raw.drf,"\n")
    col.target = col.target+1
    RESULTS[k,col.target]=mean.mahalanobis.distance.raw.drf
    colnames(RESULTS)[col.target] <- "mean.mahalanobis.distance.raw.drf"
    
    cat("sd.mahalanobis.distance.raw.drf =",sd.mahalanobis.distance.raw.drf,"\n")
    col.target = col.target+1
    RESULTS[k,col.target]=sd.mahalanobis.distance.raw.drf
    colnames(RESULTS)[col.target] <- "sd.mahalanobis.distance.raw.drf"
    
    cat("mean.mahalanobis.distance.standardise.rf =",mean.mahalanobis.distance.standardise.rf,"\n")
    col.target = col.target+1
    RESULTS[k,col.target]=mean.mahalanobis.distance.standardise.rf
    colnames(RESULTS)[col.target] <- "mean.mahalanobis.distance.standardise.rf"
    
    cat("sd.mahalanobis.distance.standardise.rf =",sd.mahalanobis.distance.standardise.rf,"\n")
    col.target = col.target+1
    RESULTS[k,col.target]=sd.mahalanobis.distance.standardise.rf
    colnames(RESULTS)[col.target] <- "sd.mahalanobis.distance.standardise.rf"
    
    cat("mean.mahalanobis.distance.standardise.drf =",mean.mahalanobis.distance.standardise.drf,"\n")
    col.target = col.target+1
    RESULTS[k,col.target]=mean.mahalanobis.distance.standardise.drf
    colnames(RESULTS)[col.target] <- "mean.mahalanobis.distance.standardise.drf"
    
    cat("sd.mahalanobis.distance.standardise.drf =",sd.mahalanobis.distance.standardise.drf,"\n")
    col.target = col.target+1
    RESULTS[k,col.target]=sd.mahalanobis.distance.standardise.drf
    colnames(RESULTS)[col.target] <- "sd.mahalanobis.distance.standardise.drf"
    
    ####################################################################################################
    ############## END COMPUTATION OF "ENERGY" SCORES                              #####################
    ####################################################################################################
    
    ####################################################################################################
    ############## END PRECISION METRICS COMPUTATION WITHOUT KEEPING iPODS VALUES  #####################
    ####################################################################################################
    
    ####################################################################################################
    ############## START PRECISION METRICS FOR recording each ind-PODS COMPUTATION #####################
    ####################################################################################################
    col.target = 2
    for (i.pods in 1:n.pods) {
      RESULTS.list.n.pods[[i.pods]][k,1] = n.tree
      RESULTS.list.n.pods[[i.pods]][k,2] = n.train
    }
    
    for (i.param in 1:length(param.list)) 
    {
      single.param.pods = as.vector(param.pods[,i.param])
      
      # dif.drf.rf for estimation = mean
      nom.result <- paste0(param.list[i.param], ".dif.drf.rf.mean")
      col.target = col.target+1
      dif.drf.rf.mean = calculate_dif.drf.rf.mean(true_values = single.param.pods , point_estimates.drf = pred.drf.mean[,i.param], point_estimates.rf = param.pred.rf[[i.param]]$expectation)
      for (i.pods in 1:n.pods) {
        colnames(RESULTS.list.n.pods[[i.pods]])[col.target]=nom.result
        RESULTS.list.n.pods[[i.pods]][k,col.target] = dif.drf.rf.mean[i.pods]
      }
      sum_dif.drf.rf.mean = summary(dif.drf.rf.mean)
      ### Affichage
      summary_str <- paste("Min:", format(sum_dif.drf.rf.mean["Min."], digits = 4, nsmall = 3),
                           "1st Qu.:", format(sum_dif.drf.rf.mean["1st Qu."], digits = 4, nsmall = 3),
                           "Median:", format(sum_dif.drf.rf.mean["Median"], digits = 4, nsmall = 3),
                           "Mean:", format(sum_dif.drf.rf.mean["Mean"], digits = 4, nsmall = 3),
                           "3rd Qu.:", format(sum_dif.drf.rf.mean["3rd Qu."], digits = 4, nsmall = 3),
                           "Max:", format(sum_dif.drf.rf.mean["Max."], digits = 4, nsmall = 3), sep = " ")
      cat("Computation done for n.pods of", nom.result, "-->", summary_str, "\n")
      ###
      
      # Relative dif drf for estimation = mean
      nom.result <- paste0(param.list[i.param], ".dif.drf.mean")
      col.target = col.target+1
      dif.drf.mean = calculate_dif.mean(true_values = single.param.pods , point_estimates = pred.drf.mean[,i.param])
      for (i.pods in 1:n.pods) {
        colnames(RESULTS.list.n.pods[[i.pods]])[col.target]=nom.result
        RESULTS.list.n.pods[[i.pods]][k,col.target] = dif.drf.mean[i.pods]
      }
      sum_dif.drf.mean = summary(dif.drf.mean)
      ### Affichage
      summary_str <- paste("Min:", format(sum_dif.drf.mean["Min."], digits = 4, nsmall = 3),
                           "1st Qu.:", format(sum_dif.drf.mean["1st Qu."], digits = 4, nsmall = 3),
                           "Median:", format(sum_dif.drf.mean["Median"], digits = 4, nsmall = 3),
                           "Mean:", format(sum_dif.drf.mean["Mean"], digits = 4, nsmall = 3),
                           "3rd Qu.:", format(sum_dif.drf.mean["3rd Qu."], digits = 4, nsmall = 3),
                           "Max:", format(sum_dif.drf.mean["Max."], digits = 4, nsmall = 3), sep = " ")
      cat("Computation done for n.pods of", nom.result, "-->", summary_str, "\n")
      ###
      
      # Relative dif rf for estimation = mean
      nom.result <- paste0(param.list[i.param], ".dif.rf.mean")
      col.target = col.target+1
      dif.rf.mean = calculate_dif.mean(true_values = single.param.pods , point_estimates = param.pred.rf[[i.param]]$expectation)
      for (i.pods in 1:n.pods) {
        colnames(RESULTS.list.n.pods[[i.pods]])[col.target]=nom.result
        RESULTS.list.n.pods[[i.pods]][k,col.target] = dif.rf.mean[i.pods]
      }
      sum_dif.rf.mean = summary(dif.rf.mean)
      ### Affichage
      summary_str <- paste("Min:", format(sum_dif.rf.mean["Min."], digits = 4, nsmall = 3),
                           "1st Qu.:", format(sum_dif.rf.mean["1st Qu."], digits = 4, nsmall = 3),
                           "Median:", format(sum_dif.rf.mean["Median"], digits = 4, nsmall = 3),
                           "Mean:", format(sum_dif.rf.mean["Mean"], digits = 4, nsmall = 3),
                           "3rd Qu.:", format(sum_dif.rf.mean["3rd Qu."], digits = 4, nsmall = 3),
                           "Max:", format(sum_dif.rf.mean["Max."], digits = 4, nsmall = 3), sep = " ")
      cat("Computation done for n.pods of", nom.result, "-->", summary_str, "\n")
      ###
      
      # dif.drf.rf for estimation = sd
      nom.result <- paste0(param.list[i.param], ".dif.drf.rf.sd")
      col.target = col.target+1
      dif.drf.rf.sd = calculate_dif.drf.rf.sd(true_values = single.param.pods , point_estimates.drf = pred.drf.sd[,i.param], point_estimates.rf = sqrt(param.pred.rf[[i.param]]$variance.cdf))
      for (i.pods in 1:n.pods) {
        colnames(RESULTS.list.n.pods[[i.pods]])[col.target]=nom.result
        RESULTS.list.n.pods[[i.pods]][k,col.target] = dif.drf.rf.sd[i.pods]
      }
      sum_dif.drf.rf.sd = summary(dif.drf.rf.sd)
      ### Affichage
      summary_str <- paste("Min:", format(sum_dif.drf.rf.sd["Min."], digits = 4, nsmall = 3),
                           "1st Qu.:", format(sum_dif.drf.rf.sd["1st Qu."], digits = 4, nsmall = 3),
                           "Median:", format(sum_dif.drf.rf.sd["Median"], digits = 4, nsmall = 3),
                           "Mean:", format(sum_dif.drf.rf.sd["Mean"], digits = 4, nsmall = 3),
                           "3rd Qu.:", format(sum_dif.drf.rf.sd["3rd Qu."], digits = 4, nsmall = 3),
                           "Max:", format(sum_dif.drf.rf.sd["Max."], digits = 4, nsmall = 3), sep = " ")
      cat("Computation done for n.pods of", nom.result, "-->", summary_str, "\n")
      
      # dif.drf.rf for estimation = lengthCI90
      nom.result <- paste0(param.list[i.param], ".dif.drf.rf.lengthCI90")
      col.target = col.target+1
      dif.drf.rf.lengthCI90 <- calculate_dif.drf.rf.lengthCI90(true_values = single.param.pods, Q5_point_estimates.drf = pred.drf.quant[, i.param],
                                                               Q95_point_estimates.drf = pred.drf.quant[, i.param+length(param.list)],
                                                               Q5_point_estimates.rf = param.pred.rf[[i.param]][["quantiles"]][,1], 
                                                               Q95_point_estimates.rf = param.pred.rf[[i.param]][["quantiles"]][,2])
      for (i.pods in 1:n.pods) {
        colnames(RESULTS.list.n.pods[[i.pods]])[col.target]=nom.result
        RESULTS.list.n.pods[[i.pods]][k,col.target] = dif.drf.rf.lengthCI90[i.pods]
      }
      sum_dif.drf.rf.lengthCI90 = summary(dif.drf.rf.lengthCI90)
      ### Affichage
      summary_str <- paste("Min:", format(sum_dif.drf.rf.lengthCI90["Min."], digits = 4, nsmall = 3),
                           "1st Qu.:", format(sum_dif.drf.rf.lengthCI90["1st Qu."], digits = 4, nsmall = 3),
                           "Median:", format(sum_dif.drf.rf.lengthCI90["Median"], digits = 4, nsmall = 3),
                           "Mean:", format(sum_dif.drf.rf.lengthCI90["Mean"], digits = 4, nsmall = 3),
                           "3rd Qu.:", format(sum_dif.drf.rf.lengthCI90["3rd Qu."], digits = 4, nsmall = 3),
                           "Max:", format(sum_dif.drf.rf.lengthCI90["Max."], digits = 4, nsmall = 3), sep = " ")
      cat("Computation done for n.pods of", nom.result, "-->", summary_str, "\n")
    }
    
    ###CORRELATIONS #####################################
    # SCENARIO 6 POPS
    #colnames(param.pods)
    # [1] "N123" "N4"   "N6"   "t124" "Dbn3" "Nbn3" "t12"  "t15"  "ra"  
    # i.pods values for estimation = cor t12.N123
    nom.result <- "cor.t12.N123"
    col.target = col.target+1
    cor.1.7 = pred.drf.cor$cor.1.7
    for (i.pods in 1:n.pods) {
      colnames(RESULTS.list.n.pods[[i.pods]])[col.target]=nom.result
      RESULTS.list.n.pods[[i.pods]][k,col.target] = cor.1.7[i.pods]
    }
    sum_cor.1.7 = summary(cor.1.7)
    ### Affichage
    summary_str <- paste("Min:", format(sum_cor.1.7["Min."], digits = 4, nsmall = 3),
                         "1st Qu.:", format(sum_cor.1.7["1st Qu."], digits = 4, nsmall = 3),
                         "Median:", format(sum_cor.1.7["Median"], digits = 4, nsmall = 3),
                         "Mean:", format(sum_cor.1.7["Mean"], digits = 4, nsmall = 3),
                         "3rd Qu.:", format(sum_cor.1.7["3rd Qu."], digits = 4, nsmall = 3),
                         "Max:", format(sum_cor.1.7["Max."], digits = 4, nsmall = 3), sep = " ")
    cat("Computation done for n.pods of", nom.result, "-->", summary_str, "\n")
    
    # i.pods values for estimation = cor Dbn3.Nbn3
    nom.result <- "cor.Dbn3.Nbn3 "
    col.target = col.target+1
    cor.5.6 = pred.drf.cor$cor.5.6
    for (i.pods in 1:n.pods) {
      colnames(RESULTS.list.n.pods[[i.pods]])[col.target]=nom.result
      RESULTS.list.n.pods[[i.pods]][k,col.target] = cor.5.6[i.pods]
    }
    sum_cor.5.6 = summary(cor.5.6)
    ### Affichage
    summary_str <- paste("Min:", format(sum_cor.5.6["Min."], digits = 4, nsmall = 3),
                         "1st Qu.:", format(sum_cor.5.6["1st Qu."], digits = 4, nsmall = 3),
                         "Median:", format(sum_cor.5.6["Median"], digits = 4, nsmall = 3),
                         "Mean:", format(sum_cor.5.6["Mean"], digits = 4, nsmall = 3),
                         "3rd Qu.:", format(sum_cor.5.6["3rd Qu."], digits = 4, nsmall = 3),
                         "Max:", format(sum_cor.5.6["Max."], digits = 4, nsmall = 3), sep = " ")
    cat("Computation done for n.pods of", nom.result, "-->", summary_str, "\n")
    
    # i.pods values for estimation = cor N6.ra
    nom.result <- "cor.N6.ra"
    col.target = col.target+1
    cor.3.9 = pred.drf.cor$cor.3.9
    for (i.pods in 1:n.pods) {
      colnames(RESULTS.list.n.pods[[i.pods]])[col.target]=nom.result
      RESULTS.list.n.pods[[i.pods]][k,col.target] = cor.3.9[i.pods]
    }
    sum_cor.3.9 = summary(cor.3.9)
    ### Affichage
    summary_str <- paste("Min:", format(sum_cor.3.9["Min."], digits = 4, nsmall = 3),
                         "1st Qu.:", format(sum_cor.3.9["1st Qu."], digits = 4, nsmall = 3),
                         "Median:", format(sum_cor.3.9["Median"], digits = 4, nsmall = 3),
                         "Mean:", format(sum_cor.3.9["Mean"], digits = 4, nsmall = 3),
                         "3rd Qu.:", format(sum_cor.3.9["3rd Qu."], digits = 4, nsmall = 3),
                         "Max:", format(sum_cor.3.9["Max."], digits = 4, nsmall = 3), sep = " ")
    cat("Computation done for n.pods of", nom.result, "-->", summary_str, "\n")
    
    # # SCENARIO MICROSAT 4 POPS t.fixed #######################################
    # # i.pods values for estimation = cor Neµ.µ 
    # nom.result <- "cor.Ne.µ"
    # col.target = col.target+1
    # cor.1.2 = pred.drf.cor$cor.1.2
    # for (i.pods in 1:n.pods) {
    #   colnames(RESULTS.list.n.pods[[i.pods]])[col.target]=nom.result
    #   RESULTS.list.n.pods[[i.pods]][k,col.target] = cor.1.2[i.pods]
    # }
    # sum_cor.1.2 = summary(cor.1.2)
    # ### Affichage
    # summary_str <- paste("Min:", format(sum_cor.1.2["Min."], digits = 4, nsmall = 3),
    #                      "1st Qu.:", format(sum_cor.1.2["1st Qu."], digits = 4, nsmall = 3),
    #                      "Median:", format(sum_cor.1.2["Median"], digits = 4, nsmall = 3),
    #                      "Mean:", format(sum_cor.1.2["Mean"], digits = 4, nsmall = 3),
    #                      "3rd Qu.:", format(sum_cor.1.2["3rd Qu."], digits = 4, nsmall = 3),
    #                      "Max:", format(sum_cor.1.2["Max."], digits = 4, nsmall = 3), sep = " ")
    # cat("Computation done for n.pods of", nom.result, "-->", summary_str, "\n")
    # 
    # Saving all i.pod ENERGY SCORE values
    col.target = col.target+1
    for (i.pod in 1:n.pods) {
      colnames(RESULTS.list.n.pods[[i.pod]])[col.target]="score.rf.ipod"
      RESULTS.list.n.pods[[i.pods]][k,col.target] = score.rf.ipod[i.pod]
    }
    col.target = col.target+1
    for (i.pod in 1:n.pods) {
      colnames(RESULTS.list.n.pods[[i.pod]])[col.target]="score.drf.ipod"
      RESULTS.list.n.pods[[i.pods]][k,col.target] = score.drf.ipod[i.pod]
    }
    col.target = col.target+1
    for (i.pod in 1:n.pods) {
      colnames(RESULTS.list.n.pods[[i.pod]])[col.target]="score.standardise.rf.ipod"
      RESULTS.list.n.pods[[i.pods]][k,col.target] = score.standardise.rf.ipod[i.pod]
    }
    col.target = col.target+1
    for (i.pod in 1:n.pods) {
      colnames(RESULTS.list.n.pods[[i.pod]])[col.target]="score.standardise.drf.ipod"
      RESULTS.list.n.pods[[i.pods]][k,col.target] = score.standardise.drf.ipod[i.pod]
    }
    
    
    ####################################################################################################
    ############## END PRECISION METRICS FOR recording each ind-PODS COMPUTATION   #####################
    ####################################################################################################
    
    # Recurrent saving as a RDS object (each n.tree*n.train pair)
    # saveRDS(RESULTS, file = "RESULTS.mean.SCENARIO_5POPS_indSNPs_ADMIXTURE_13parameters_SNP-Param.Pods.PRIOR_various.ntree.ntrain.rds")
    # saveRDS(RESULTS.list.n.pods, file = "RESULTS.list.n.pods.SCENARIO_5POPS_indSNPs_ADMIXTURE_13parameters_SNP-Param.Pods.PRIOR_various.ntree.ntrain.rds")
    
    # saveRDS(RESULTS, file = "RESULTS.mean.SCENARIO.ONEpop.microsat.Pods.PRIOR_various.ntree.ntrain.rds")
    # saveRDS(RESULTS.list.n.pods, file = "RESULTS.list.n.pods.SCENARIO.ONEpop.microsat.Pods.PRIOR_various.ntree.ntrain.rds")
    # 
    saveRDS(RESULTS, file = "RESULTS.6POP.RMSE.HDP.SCORE.rds")
    saveRDS(RESULTS.list.n.pods, file = "RESULTS.list.n.pods.6POP.RMSE.HDP.SCORE.rds")
    # JEFF: here are the two output files = one with punctual mean and sd values (RESULTS.6POP.RMSE.HDP.SCORE.rds)
    # JEFF: and one with all ipods values recorded (RESULTS.list.n.pods.6POP.RMSE.HDP.SCORE.rds)
    
  }
}
###########################################################################################################
####################### END LOOP OVER n.tree AND n.train ###############################################
###########################################################################################################


# RESULTS <- readRDS("RESULTS.6POP.RMSE.HDP.SCORE.rds")
# RESULTS.list.n.pods <- readRDS("RESULTS.list.n.pods.6POP.RMSE.HDP.SCORE.rds")

}



    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    if (PRECISION.METRICS.FOR.MEAN.COMPUTATIONS==TRUE){
    ####################################################################################################
    ############## START PRECISION METRICS OR MEAN COMPUTATIONS             ############################
    ####################################################################################################
    RESULTS[k,1] <- n.tree
    RESULTS[k,2] <- n.train
    col.target = 2
    
    for (i.param in 1:length(param.list)) 
  {
    single.param.pods = as.vector(param.pods[,i.param])

    # NMAE normalized by mean
    nom.result <- paste0(param.list[i.param], ".NMAE.mean.mean.rf")
    NMAE.mean.mean.rf = calculate_NMAE_mean(true_values = single.param.pods , point_estimates = param.pred.rf[[i.param]]$expectation)
    assign(nom.result, NMAE.mean.mean.rf)
    cat(nom.result, "=",NMAE.mean.mean.rf ,"\n")
    col.target = col.target+1
    RESULTS[k,col.target]=NMAE.mean.mean.rf
    colnames(RESULTS)[col.target] <- nom.result
    nom.result <- paste0(param.list[i.param], ".NMAE.mean.mean.drf")
    NMAE.mean.mean.drf = calculate_NMAE_mean(true_values = single.param.pods , point_estimates = pred.drf.mean[,i.param])
    assign(nom.result, NMAE.mean.mean.drf)
    cat(nom.result, "=",NMAE.mean.mean.drf ,"\n")
    col.target = col.target+1
    RESULTS[k,col.target]=NMAE.mean.mean.drf
    colnames(RESULTS)[col.target] <- nom.result

     # mean(SD)   
    nom.result <- paste0(param.list[i.param], ".sd.mean.rf")
    sd.mean.rf = mean(sqrt(param.pred.rf[[i.param]]$variance.cdf))
    assign(nom.result, sd.mean.rf)
    cat(nom.result, "=",sd.mean.rf,"\n")
    col.target = col.target+1
    RESULTS[k,col.target]=sd.mean.rf
    colnames(RESULTS)[col.target] <- nom.result
    nom.result <- paste0(param.list[i.param], ".sd.mean.drf")
    sd.mean.drf = mean(pred.drf.sd[,i.param])
    assign(nom.result, sd.mean.drf)
    cat(nom.result, "=",sd.mean.drf,"\n")
    col.target = col.target+1
    RESULTS[k,col.target]=sd.mean.drf
    colnames(RESULTS)[col.target] <- nom.result
   
    # 90% CI
    nom.result <- paste0(param.list[i.param], ".CI90.rf")
    CI90.rf <- calculate_90CI(true_values = single.param.pods, 
               Q5_point_estimates = param.pred.rf[[i.param]][["quantiles"]][,1], Q95_point_estimates = param.pred.rf[[i.param]][["quantiles"]][,2])
    assign(nom.result,CI90.rf)
    cat(nom.result, "=",CI90.rf ,"\n")
    col.target = col.target+1
    RESULTS[k,col.target]=CI90.rf
    colnames(RESULTS)[col.target] <- nom.result
    nom.result <- paste0(param.list[i.param], ".CI90.drf")
    CI90.drf <- calculate_90CI(true_values = single.param.pods, 
               Q5_point_estimates = pred.drf.quant[, i.param], Q95_point_estimates = pred.drf.quant[, i.param+length(param.list)])
    assign(nom.result,CI90.drf)
    cat(nom.result, "=",CI90.drf ,"\n")
    col.target = col.target+1
    RESULTS[k,col.target]=CI90.drf
    colnames(RESULTS)[col.target] <- nom.result
    
  }
    
    #INTERSECTION.CI90.drf = calculate_90CI_intersection(V1_true_values = N11.TARGET.param.pods.vector, V1_Q5 = all.pred.drf.quant.target$quantile.8.q.0.05, V1_Q95 = all.pred.drf.quant.target$quantile.8.q.0.95,
    #                                                                                   V2_true_values = t.TARGET.param.pods.vector, V2_Q5 = all.pred.drf.quant.target$quantile.5.q.0.05, V2_Q95 = all.pred.drf.quant.target$quantile.5.q.0.95)
    #colnames(param.pods)
    # "N123"     "N4"       "N5"       "t124"     "Dbn3"     "Nbn3"     "t13"      "t12"      "t15"      "ra"       "t12N123"  "Dbn3Nbn3" "t15N5"   
    #cat("INTERSECTION.CI90.drf =",INTERSECTION.CI90.drf,"\n")
    # cor t12.N123
    mean.cor = mean(pred.drf.cor$cor.1.8)
    cat("t12.N123.mean.cor =",mean.cor,"\n")
    col.target = col.target+1
    RESULTS[k,col.target]=mean.cor
    colnames(RESULTS)[col.target] <- "t12.N123.mean.cor"
    sd.cor = sqrt(var(pred.drf.cor$cor.1.8))
    cat("t12.N123.sd.cor =",sd.cor,"\n")
    col.target = col.target+1
    RESULTS[k,col.target]=sd.cor
    colnames(RESULTS)[col.target] <- "t12.N123.sd.cor"
    
    # cor Dbn3.Nbn3
    mean.cor = mean(pred.drf.cor$cor.5.6)
    cat("Dbn3.Nbn3.mean.cor =",mean.cor,"\n")
    col.target = col.target+1
    RESULTS[k,col.target]=mean.cor
    colnames(RESULTS)[col.target] <- "Dbn3.Nbn3.mean.cor"
    sd.cor = sqrt(var(pred.drf.cor$cor.5.6))
    cat("Dbn3.Nbn3.sd.cor =",sd.cor,"\n")
    col.target = col.target+1
    RESULTS[k,col.target]=sd.cor
    colnames(RESULTS)[col.target] <- "Dbn3.Nbn3.sd.cor"
    
    # cor t15.N5
    mean.core = mean(pred.drf.cor$cor.3.9)
    cat("t15.N5.mean.cor =",mean.cor,"\n")
    col.target = col.target+1
    RESULTS[k,col.target]=mean.core
    colnames(RESULTS)[col.target] <- "t15.N5.mean.cor"
    sd.cor = sqrt(var(pred.drf.cor$cor.3.9))
    cat("t15.N5.sd.cor =",sd.cor,"\n")
    col.target = col.target+1
    RESULTS[k,col.target]=sd.cor
    colnames(RESULTS)[col.target] <- "t15.N5.sd.cor"
   
    # cor N5.ra 
    mean.core = mean(pred.drf.cor$cor.3.10)
    cat("N5.ra.mean.cor =",mean.cor,"\n")
    col.target = col.target+1
    RESULTS[k,col.target]=mean.core
    colnames(RESULTS)[col.target] <- "N5.ra.mean.cor"
    sd.cor = sqrt(var(pred.drf.cor$cor.3.10))
    cat("N5.ra.sd.cor =",sd.cor,"\n")
    col.target = col.target+1
    RESULTS[k,col.target]=sd.cor
    colnames(RESULTS)[col.target] <- "N5.ra.sd.cor"
    
    ####################################################################################################
    ############## END PRECISION METRICS FOR MEAN COMPUTATIONS                   #######################
    ####################################################################################################
  }
    
    if (PRECISION.METRICS.FOR.ind.PODS.COMPUTATIONS==TRUE){
    ####################################################################################################
    ############## START PRECISION METRICS FOR ind-PODS COMPUTATIONS             #######################
    ####################################################################################################
    col.target = 2
    for (i.pods in 1:n.pods) {
      RESULTS.list.n.pods[[i.pods]][k,1] = n.tree
      RESULTS.list.n.pods[[i.pods]][k,2] = n.train
    }

    for (i.param in 1:length(param.list)) 
 {
      single.param.pods = as.vector(param.pods[,i.param])
      
      # dif.drf.rf for estimation = mean
      nom.result <- paste0(param.list[i.param], ".dif.drf.rf.mean")
      col.target = col.target+1
      dif.drf.rf.mean = calculate_dif.drf.rf.mean(true_values = single.param.pods , point_estimates.drf = pred.drf.mean[,i.param], point_estimates.rf = param.pred.rf[[i.param]]$expectation)
      for (i.pods in 1:n.pods) {
        colnames(RESULTS.list.n.pods[[i.pods]])[col.target]=nom.result
        RESULTS.list.n.pods[[i.pods]][k,col.target] = dif.drf.rf.mean[i.pods]
      }
      sum_dif.drf.rf.mean = summary(dif.drf.rf.mean)
      ### Affichage
      summary_str <- paste("Min:", format(sum_dif.drf.rf.mean["Min."], digits = 4, nsmall = 3),
                           "1st Qu.:", format(sum_dif.drf.rf.mean["1st Qu."], digits = 4, nsmall = 3),
                           "Median:", format(sum_dif.drf.rf.mean["Median"], digits = 4, nsmall = 3),
                           "Mean:", format(sum_dif.drf.rf.mean["Mean"], digits = 4, nsmall = 3),
                           "3rd Qu.:", format(sum_dif.drf.rf.mean["3rd Qu."], digits = 4, nsmall = 3),
                           "Max:", format(sum_dif.drf.rf.mean["Max."], digits = 4, nsmall = 3), sep = " ")
      cat("Computation done for n.pods of", nom.result, "-->", summary_str, "\n")
      ###
      
      # Relative dif drf for estimation = mean
      nom.result <- paste0(param.list[i.param], ".dif.drf.mean")
      col.target = col.target+1
      dif.drf.mean = calculate_dif.mean(true_values = single.param.pods , point_estimates = pred.drf.mean[,i.param])
      for (i.pods in 1:n.pods) {
        colnames(RESULTS.list.n.pods[[i.pods]])[col.target]=nom.result
        RESULTS.list.n.pods[[i.pods]][k,col.target] = dif.drf.mean[i.pods]
      }
      sum_dif.drf.mean = summary(dif.drf.mean)
      ### Affichage
      summary_str <- paste("Min:", format(sum_dif.drf.mean["Min."], digits = 4, nsmall = 3),
                           "1st Qu.:", format(sum_dif.drf.mean["1st Qu."], digits = 4, nsmall = 3),
                           "Median:", format(sum_dif.drf.mean["Median"], digits = 4, nsmall = 3),
                           "Mean:", format(sum_dif.drf.mean["Mean"], digits = 4, nsmall = 3),
                           "3rd Qu.:", format(sum_dif.drf.mean["3rd Qu."], digits = 4, nsmall = 3),
                           "Max:", format(sum_dif.drf.mean["Max."], digits = 4, nsmall = 3), sep = " ")
      cat("Computation done for n.pods of", nom.result, "-->", summary_str, "\n")
      ###
      
      # Relative dif rf for estimation = mean
      nom.result <- paste0(param.list[i.param], ".dif.rf.mean")
      col.target = col.target+1
      dif.rf.mean = calculate_dif.mean(true_values = single.param.pods , point_estimates = param.pred.rf[[i.param]]$expectation)
      for (i.pods in 1:n.pods) {
        colnames(RESULTS.list.n.pods[[i.pods]])[col.target]=nom.result
        RESULTS.list.n.pods[[i.pods]][k,col.target] = dif.rf.mean[i.pods]
      }
      sum_dif.rf.mean = summary(dif.rf.mean)
      ### Affichage
      summary_str <- paste("Min:", format(sum_dif.rf.mean["Min."], digits = 4, nsmall = 3),
                           "1st Qu.:", format(sum_dif.rf.mean["1st Qu."], digits = 4, nsmall = 3),
                           "Median:", format(sum_dif.rf.mean["Median"], digits = 4, nsmall = 3),
                           "Mean:", format(sum_dif.rf.mean["Mean"], digits = 4, nsmall = 3),
                           "3rd Qu.:", format(sum_dif.rf.mean["3rd Qu."], digits = 4, nsmall = 3),
                           "Max:", format(sum_dif.rf.mean["Max."], digits = 4, nsmall = 3), sep = " ")
      cat("Computation done for n.pods of", nom.result, "-->", summary_str, "\n")
      ###
      
      # dif.drf.rf for estimation = sd
      nom.result <- paste0(param.list[i.param], ".dif.drf.rf.sd")
      col.target = col.target+1
      dif.drf.rf.sd = calculate_dif.drf.rf.sd(true_values = single.param.pods , point_estimates.drf = pred.drf.sd[,i.param], point_estimates.rf = sqrt(param.pred.rf[[i.param]]$variance.cdf))
      for (i.pods in 1:n.pods) {
        colnames(RESULTS.list.n.pods[[i.pods]])[col.target]=nom.result
        RESULTS.list.n.pods[[i.pods]][k,col.target] = dif.drf.rf.sd[i.pods]
      }
      sum_dif.drf.rf.sd = summary(dif.drf.rf.sd)
      ### Affichage
      summary_str <- paste("Min:", format(sum_dif.drf.rf.sd["Min."], digits = 4, nsmall = 3),
                           "1st Qu.:", format(sum_dif.drf.rf.sd["1st Qu."], digits = 4, nsmall = 3),
                           "Median:", format(sum_dif.drf.rf.sd["Median"], digits = 4, nsmall = 3),
                           "Mean:", format(sum_dif.drf.rf.sd["Mean"], digits = 4, nsmall = 3),
                           "3rd Qu.:", format(sum_dif.drf.rf.sd["3rd Qu."], digits = 4, nsmall = 3),
                           "Max:", format(sum_dif.drf.rf.sd["Max."], digits = 4, nsmall = 3), sep = " ")
      cat("Computation done for n.pods of", nom.result, "-->", summary_str, "\n")
      
      # dif.drf.rf for estimation = lengthCI90
      nom.result <- paste0(param.list[i.param], ".dif.drf.rf.lengthCI90")
      col.target = col.target+1
      dif.drf.rf.lengthCI90 <- calculate_dif.drf.rf.lengthCI90(true_values = single.param.pods, Q5_point_estimates.drf = pred.drf.quant[, i.param],
                                                               Q95_point_estimates.drf = pred.drf.quant[, i.param+length(param.list)],
                                                               Q5_point_estimates.rf = param.pred.rf[[i.param]][["quantiles"]][,1], 
                                                               Q95_point_estimates.rf = param.pred.rf[[i.param]][["quantiles"]][,2])
      for (i.pods in 1:n.pods) {
      colnames(RESULTS.list.n.pods[[i.pods]])[col.target]=nom.result
      RESULTS.list.n.pods[[i.pods]][k,col.target] = dif.drf.rf.lengthCI90[i.pods]
      }
      sum_dif.drf.rf.lengthCI90 = summary(dif.drf.rf.lengthCI90)
      ### Affichage
      summary_str <- paste("Min:", format(sum_dif.drf.rf.lengthCI90["Min."], digits = 4, nsmall = 3),
                           "1st Qu.:", format(sum_dif.drf.rf.lengthCI90["1st Qu."], digits = 4, nsmall = 3),
                           "Median:", format(sum_dif.drf.rf.lengthCI90["Median"], digits = 4, nsmall = 3),
                           "Mean:", format(sum_dif.drf.rf.lengthCI90["Mean"], digits = 4, nsmall = 3),
                           "3rd Qu.:", format(sum_dif.drf.rf.lengthCI90["3rd Qu."], digits = 4, nsmall = 3),
                           "Max:", format(sum_dif.drf.rf.lengthCI90["Max."], digits = 4, nsmall = 3), sep = " ")
      cat("Computation done for n.pods of", nom.result, "-->", summary_str, "\n")
 }

    #colnames(param.pods)
    # "N123"     "N4"       "N5"       "t124"     "Dbn3"     "Nbn3"     "t13"      "t12"      "t15"      "ra"       "t12N123"  "Dbn3Nbn3" "t15N5"   
    # i.pods values for estimation = cor t12.N123 
    nom.result <- "cor.N12.t123"
    col.target = col.target+1
    cor.1.8 = pred.drf.cor$cor.1.8
    for (i.pods in 1:n.pods) {
      colnames(RESULTS.list.n.pods[[i.pods]])[col.target]=nom.result
      RESULTS.list.n.pods[[i.pods]][k,col.target] = cor.1.8[i.pods]
    }
    sum_cor.1.8 = summary(cor.1.8)
    ### Affichage
    summary_str <- paste("Min:", format(sum_cor.1.8["Min."], digits = 4, nsmall = 3),
                         "1st Qu.:", format(sum_cor.1.8["1st Qu."], digits = 4, nsmall = 3),
                         "Median:", format(sum_cor.1.8["Median"], digits = 4, nsmall = 3),
                         "Mean:", format(sum_cor.1.8["Mean"], digits = 4, nsmall = 3),
                         "3rd Qu.:", format(sum_cor.1.8["3rd Qu."], digits = 4, nsmall = 3),
                         "Max:", format(sum_cor.1.8["Max."], digits = 4, nsmall = 3), sep = " ")
    cat("Computation done for n.pods of", nom.result, "-->", summary_str, "\n")
    
    # i.pods values for estimation = cor Dbn3.Nbn3 
    nom.result <- "cor.Dbn3.Nbn3 "
    col.target = col.target+1
    cor.5.6 = pred.drf.cor$cor.5.6
    for (i.pods in 1:n.pods) {
      colnames(RESULTS.list.n.pods[[i.pods]])[col.target]=nom.result
      RESULTS.list.n.pods[[i.pods]][k,col.target] = cor.5.6[i.pods]
    }
    sum_cor.5.6 = summary(cor.5.6)
    ### Affichage
    summary_str <- paste("Min:", format(sum_cor.5.6["Min."], digits = 4, nsmall = 3),
                         "1st Qu.:", format(sum_cor.5.6["1st Qu."], digits = 4, nsmall = 3),
                         "Median:", format(sum_cor.5.6["Median"], digits = 4, nsmall = 3),
                         "Mean:", format(sum_cor.5.6["Mean"], digits = 4, nsmall = 3),
                         "3rd Qu.:", format(sum_cor.5.6["3rd Qu."], digits = 4, nsmall = 3),
                         "Max:", format(sum_cor.5.6["Max."], digits = 4, nsmall = 3), sep = " ")
    cat("Computation done for n.pods of", nom.result, "-->", summary_str, "\n")
    
  # i.pods values for estimation = cor t15.N5 
  nom.result <- "cor.t15.N5"
  col.target = col.target+1
  cor.3.9 = pred.drf.cor$cor.3.9
  for (i.pods in 1:n.pods) {
    colnames(RESULTS.list.n.pods[[i.pods]])[col.target]=nom.result
    RESULTS.list.n.pods[[i.pods]][k,col.target] = cor.3.9[i.pods]
  }
  sum_cor.3.9 = summary(cor.3.9)
  ### Affichage
  summary_str <- paste("Min:", format(sum_cor.3.9["Min."], digits = 4, nsmall = 3),
                       "1st Qu.:", format(sum_cor.3.9["1st Qu."], digits = 4, nsmall = 3),
                       "Median:", format(sum_cor.3.9["Median"], digits = 4, nsmall = 3),
                       "Mean:", format(sum_cor.3.9["Mean"], digits = 4, nsmall = 3),
                       "3rd Qu.:", format(sum_cor.3.9["3rd Qu."], digits = 4, nsmall = 3),
                       "Max:", format(sum_cor.3.9["Max."], digits = 4, nsmall = 3), sep = " ")
  cat("Computation done for n.pods of", nom.result, "-->", summary_str, "\n")

  # i.pods values for estimation = cor N5.ra 
  nom.result <- "cor.N5.ra"
  col.target = col.target+1
  cor.3.10 = pred.drf.cor$cor.3.10
  for (i.pods in 1:n.pods) {
    colnames(RESULTS.list.n.pods[[i.pods]])[col.target]=nom.result
    RESULTS.list.n.pods[[i.pods]][k,col.target] = cor.3.10[i.pods]
  }
  sum_cor.3.10 = summary(cor.3.10)
  ### Affichage
  summary_str <- paste("Min:", format(sum_cor.3.10["Min."], digits = 4, nsmall = 3),
                       "1st Qu.:", format(sum_cor.3.10["1st Qu."], digits = 4, nsmall = 3),
                       "Median:", format(sum_cor.3.10["Median"], digits = 4, nsmall = 3),
                       "Mean:", format(sum_cor.3.10["Mean"], digits = 4, nsmall = 3),
                       "3rd Qu.:", format(sum_cor.3.10["3rd Qu."], digits = 4, nsmall = 3),
                       "Max:", format(sum_cor.3.10["Max."], digits = 4, nsmall = 3), sep = " ")
  cat("Computation done for n.pods of", nom.result, "-->", summary_str, "\n")

####################################################################################################
############## END PRECISION METRICS FOR ind-PODS COMPUTATIONS               #######################
####################################################################################################
    }

# Recurrent saving as a RDS object (each n.tree*n.train pair)
# saveRDS(RESULTS, file = "RESULTS.mean.SCENARIO_5POPS_indSNPs_ADMIXTURE_13parameters_SNP-Param.Pods.PRIOR_various.ntree.ntrain.rds")
# saveRDS(RESULTS.list.n.pods, file = "RESULTS.list.n.pods.SCENARIO_5POPS_indSNPs_ADMIXTURE_13parameters_SNP-Param.Pods.PRIOR_various.ntree.ntrain.rds")
# 
# saveRDS(RESULTS, file = "RESULTS.mean.SCENARIO_5POPS_indSNPs_ADMIXTURE_13parameters_SNP-Param.Pods.PRIOR_various.ntree.SINGLE.ntrain.rds")
# saveRDS(RESULTS.list.n.pods, file = "RESULTS.list.n.pods.SCENARIO_5POPS_indSNPs_ADMIXTURE_13parameters_SNP-Param.Pods.PRIOR_various.ntree.SINGLE.ntrain.rds")

    
# saveRDS(RESULTS, file = "RESULTS.mean.SCENARIO_5POPS_indSNPs_ADMIXTURE_13parameters_SNP-Param.Pods.PRIOR_various.ntree.ntrain.rds")
# saveRDS(RESULTS.list.n.pods, file = "RESULTS.list.n.pods.SCENARIO_5POPS_indSNPs_ADMIXTURE_13parameters_SNP-Param.Pods.PRIOR_various.ntree.ntrain.rds")

saveRDS(RESULTS, file = "RESULTS.mean.TEST.Param.Pods.PRIOR_various.ntree.SINGLE.ntrain.rds")
saveRDS(RESULTS.list.n.pods, file = "RESULTS.list.n.pods.TEST.Param.Pods.PRIOR_various.ntree.SINGLE.ntrain.rds")   
 


###########################################################################################################
####################### ENDING LOOP OVER n.tree AND n.train ###############################################
###########################################################################################################



















if (ANNEX.COMPUTATION.FOR.DRF.RF.EVALUATION.ala.NAF==TRUE){

###########################################################################################################
####################### START CALCULS ANNEXES + FIGURE DESIGN FROM RESULTS.list.n.pods and RESULTS  #######
###########################################################################################################
#RESULTS.ALL <- readRDS("RESULTS.mean.SCENARIO_5POPS_indSNPs_ADMIXTURE_13parameters_SNP-Param.Pods.PRIOR_various.ntree.ntrain.rds")
#RESULTS.list.n.pods <- readRDS("RESULTS.list.n.pods.SCENARIO_5POPS_indSNPs_ADMIXTURE_13parameters_SNP-Param.Pods.PRIOR_various.ntree.ntrain.rds")

RESULTS.list.n.pods <- readRDS("RESULTS.list.n.pods.TEST.Param.Pods.PRIOR_various.ntree.SINGLE.ntrain.rds")

#head(RESULTS.ALL,n=1000)

##################### CALCULS ANNEXES ######################################################## 
# Proportion de valeurs de pods pour lesquelles drf "meilleurs" que rf - ceci pour toutes les metrics "diff" ou drf et rf sont comparees
# Créer une fonction pour filtrer les lignes de RESULTS.list.n.pods sur n.tree et n.train
RESULTS.list.n.pods[[1]]
ntree = 2000
ntrain = 50000
filter_n_tree_n_train <- function(data_frame) {return(data_frame[data_frame$n.tree == ntree & data_frame$n.train == ntrain, ])}
# Appliquer la fonction à chaque élément de RESULTS.list.n.pods
sub.RESULTS.list.n.pods <- lapply(RESULTS.list.n.pods, filter_n_tree_n_train)
sub.RESULTS.list.n.pods[[1]]
neg.prop = as.data.frame(matrix(0,nrow=1,ncol=length(colnames(sub.RESULTS.list.n.pods[[1]]))))
colnames(neg.prop)=colnames(sub.RESULTS.list.n.pods[[1]]) 
neg.proportion = neg.prop[,-c(1,2)]
for (i in 1:length(colnames(neg.proportion))) {
  metric.target <- sapply(sub.RESULTS.list.n.pods, function(df) df[[colnames(neg.proportion)[i]]])
  summary(metric.target)
  proportion = mean(metric.target<0)
  neg.proportion[1,i]=proportion
}

cat("Proportion de valeurs négatives pour les colonnes contenant 'drf.rf'\n")
cat("n.tree =", sub.RESULTS.list.n.pods[[1]]$n.tree, "n.train =", sub.RESULTS.list.n.pods[[1]]$n.train, "\n")
for (i in 1:ncol(neg.proportion)) {
  colname <- names(neg.proportion)[i]
  # Vérifier si le nom de la colonne contient "drf.rf"
  if (grepl("drf.rf", colname)) {
    value <- neg.proportion[1, i]
    cat(colname, "=", value, "\n")
  }
}

##################### FAIRE DES PLOTS A PARTIR DE RESULTS.list.n.pods  ######################################################## 
### Direct calling
colnames(param.pods)
y.target1 = param.pred.rf[[7]]$expectation
y.target = param.pred.rf[[7]]$expectation
y.target2 = pred.drf.mean[,7]
y.target = pred.drf.mean[,7]
head(y.target,n=5)
summary(y.target)
cor(param.pred.rf[[7]]$expectation,pred.drf.mean[,7])
x.target = param.pods$ra
cor(x.target,y.target)
cor(x.target,y.target2)
cor(x.target,y.target1)
cor(x.target1,y.target2)
cor(param.pods$ra,pred.drf.mean[,7])
cor(param.pods$ra,param.pred.rf[[7]]$expectation)
cor(param.pred.rf[[7]]$expectation,pred.drf.mean[,7])
# Indirect calling
y.target1.name = "N12.dif.drf.mean"
y.target2.name = "N12.dif.rf.mean"
y.target.name = "N123.dif.drf.rf.mean"
y.target.name = "ra.dif.drf.rf.mean"
y.target.name="Dbn3Nbn3.dif.drf.rf.lengthCI90"
colnames(param.pods)
#x.target = param.pods$N12
x.target = param.pods$ra
head(x.target, n=5)
summary(x.target,n=5)
library(ggplot2)
library(reshape2)
##############################################
# Créer une fonction pour filtrer les lignes de RESULTS.list.n.pods sur n.tree et n.train
RESULTS.list.n.pods[[1]]
ntree = 2000
ntrain = 50000
filter_n_tree_n_train <- function(data_frame) {return(data_frame[data_frame$n.tree == ntree & data_frame$n.train == ntrain, ])}
# Appliquer la fonction à chaque élément de RESULTS.list.n.pods
sub.RESULTS.list.n.pods <- lapply(RESULTS.list.n.pods, filter_n_tree_n_train)
sub.RESULTS.list.n.pods[[1]]
# Pour extraire les n.pods valeurs de chaque variable col.target 
y.target1 <- sapply(sub.RESULTS.list.n.pods, function(df) df[[y.target1.name]])
y.target2 <- sapply(sub.RESULTS.list.n.pods, function(df) df[[y.target2.name]])
y.target <- sapply(sub.RESULTS.list.n.pods, function(df) df[[y.target.name]])
head(y.target1, n=10)
head(y.target2, n=10)
head(y.target, n=10)
summary(y.target1)
summary(y.target2)
summary(y.target)

# Graphique-plot avec des points issus de une variable
# Création des données
data <- data.frame(x = x.target,y = y.target)
head(data,n=5)
y_min <- 0
y_max <- 1
# Création du graphique
gg <- ggplot(data, aes(x = x, y = y)) +
  geom_point(shape = 1, color = "black") +  # Points en noir avec fond vide
  geom_smooth(method = "lm", formula = y ~ poly(x, degree=3), se = FALSE, color = "red") +  # Ajustement polynomiale
  ylim(y_min, y_max) +  # Limites de l'axe y
  xlab("ra pod value") +  # Légende de l'axe des x
  ylab("ra.dif.drf.rf") +  # Légende de l'axe des y
  theme_minimal()  # Utilisation d'un thème minimaliste
print(gg)

# Graphique-plot avec des points issus de deux variables
# Création des données
data <- data.frame(x = c(x.target, x.target), 
                   y = c(y.target1, y.target2), 
                   type = rep(c("drf", "rf"), each=length(x.target)))
y_min <- 0.0 # Valeur minimale pour l'axe y
y_max <- 1.0 # Valeur maximale pour l'axe y
# Création du graphique
gg <- ggplot(data, aes(x = x, y = y, color = type)) +
  ylim(y_min, y_max) + # Définir les limites de l'axe des y
  geom_point(shape = 1) +  # shape = 1 pour des points avec un fond vide
  scale_color_manual(values = c("drf" = "red", "rf" = "blue")) +
  xlab("N12 pod value") +  # Légende de l'axe des ordonnées (x)
  ylab("N12.dif") +  # Légende de l'axe des abscisses (y)
  theme_minimal() +
  labs(color = "y-Targets")
print(gg)

##################### FAIRE DES BOXPLOT A PARTIR DE RESULTS.list.n.pods  ########################################################
col.target = "N12.dif.drf.rf.mean"
#col.target = "ra.dif.drf.rf.lengthCI90"
col.target = "N12.dif.drf.mean"
library(ggplot2)
library(reshape2)
#
ntree = 1000
filter_n_tree <- function(data_frame) {return(data_frame[data_frame$n.tree == ntree, ])}
# Appliquer la fonction à chaque élément de RESULTS.list.n.pods
sub.RESULTS.list.n.pods <- lapply(RESULTS.list.n.pods, filter_n_tree)
# Pour extraire les n.pods valeurs de chaque variable col.target (cf boxplot):
values_col.target <- sapply(sub.RESULTS.list.n.pods, function(df) df[[col.target]])
values_col.target = t(values_col.target)
df.box.plot = values_col.target
colnames(df.box.plot) = vector.n.train
dim(values_col.target)
head(values_col.target,n=5)
dim(df.box.plot)
head(df.box.plot,n=5)
summary(df.box.plot)

# Remodeler les données en format long
df.box.plot_long <- melt(df.box.plot)
head(df.box.plot_long, n=5)
df.box.plot_long = df.box.plot_long[-1]
colnames(df.box.plot_long)=c("n.train","pod.value")
head(df.box.plot_long, n=5)

# Calculer les quantiles à 5% et 95%
quantiles <- tapply(df.box.plot_long$pod.value, df.box.plot_long$n.train, function(x) quantile(x, c(0.05, 0.95)))
quantiles
# Récupérer les niveaux réels de n.train
levels_n.train <- levels(factor(df.box.plot_long$n.train))
levels_n.train

# Créer le boxplot avec ggplot2
y_min <- -0.5 # Valeur minimale pour l'axe y
y_max <- 1.5 # Valeur maximale pour l'axe y
ggplot(df.box.plot_long, aes(x = factor(n.train), y = pod.value)) +
  geom_boxplot(notch = TRUE, notchidth = 0.95, fill = NA, color = "black", outlier.shape = TRUE, alpha = 0.5, position = "dodge", varwidth = TRUE) +
  labs(title = "n.tree = 2000 - Boxplot for various n.train values",
       x = "n.train",
       y = col.target) +
  scale_x_discrete(labels = levels_n.train) +
  scale_y_continuous(limits = c(y_min, y_max)) +
  # Ajouter la valeur moyenne
  stat_summary(fun = mean, geom = "point", shape = 16, size = 3, fill = "red", color = "red", position = position_dodge(width = 0.75))+
  # Ajouter les quantiles à 5% et 95% pour les different boxplot avec les facteurs n.train
  geom_point(data = subset(df.box.plot_long, n.train == "500"), aes(x = factor(n.train), y = quantiles[["500"]][1]), shape = 16, size = 3, color = "blue", position = position_dodge(width = 0.75)) +
  geom_point(data = subset(df.box.plot_long, n.train == "500"), aes(x = factor(n.train), y = quantiles[["500"]][2]), shape = 16, size = 3, color = "blue", position = position_dodge(width = 0.75)) +
  geom_point(data = subset(df.box.plot_long, n.train == "1000"), aes(x = factor(n.train), y = quantiles[["1000"]][1]), shape = 16, size = 3, color = "blue", position = position_dodge(width = 0.75)) +
  geom_point(data = subset(df.box.plot_long, n.train == "1000"), aes(x = factor(n.train), y = quantiles[["1000"]][2]), shape = 16, size = 3, color = "blue", position = position_dodge(width = 0.75)) +
  geom_point(data = subset(df.box.plot_long, n.train == "2000"), aes(x = factor(n.train), y = quantiles[["2000"]][1]), shape = 16, size = 3, color = "blue", position = position_dodge(width = 0.75)) +
  geom_point(data = subset(df.box.plot_long, n.train == "2000"), aes(x = factor(n.train), y = quantiles[["2000"]][2]), shape = 16, size = 3, color = "blue", position = position_dodge(width = 0.75)) 
  # geom_point(data = subset(df.box.plot_long, n.train == "5000"), aes(x = factor(n.train), y = quantiles[["5000"]][1]), shape = 16, size = 3, color = "blue", position = position_dodge(width = 0.75)) +
  # geom_point(data = subset(df.box.plot_long, n.train == "5000"), aes(x = factor(n.train), y = quantiles[["5000"]][2]), shape = 16, size = 3, color = "blue", position = position_dodge(width = 0.75))


####################################### Figures for mean estimations from RESULTS #############################
# library(ggplot2)
# RESULTS <- readRDS("SCENARIO_5POPS_indSNPs_ADMIXTURE_7parameters_SNP-Param.Pods.PRIOR_various.ntree.ntrain.rds")
# toto = readRDS("SCENARIO_5POPS_indSNPs_ADMIXTURE_7parameters_SNP-Param.Pods.PRIOR_various.ntree.ntrain.rds")
# head(toto[,],n=1000)
# print(toto)
# # Extraire les lignes où n.train=XXXX ou n.tree=XXX
# #sub.RESULTS <- subset(RESULTS, n.train == 20000)
# #sub.RESULTS <- subset(RESULTS, n.train == 5000)
# sub.RESULTS <- subset(RESULTS, n.tree == 1000)
# head(sub.RESULTS[,], n=100)
# colnames(sub.RESULTS)

#################################  NMAE ####################################
# Création de graphiques avec des légendes personnalisées (NMAE = f(n.train | n.tree fixed)) - t124.
x_min <- 500  # Valeur minimale pour l'axe x
x_max <- 50000  # Valeur maximale pour l'axe x
y_min <- 0.10 # Valeur minimale pour l'axe y
y_max <- 0.32 # Valeur maximale pour l'axe y
gg <- ggplot(data = sub.RESULTS, aes(x = n.train, y = t124.NMAE.mean.mean.drf)) +
  geom_point(size = 3, aes(color = "drf")) + 
  geom_point(data = sub.RESULTS, aes(x = n.train, y = t124.NMAE.mean.mean.rf, color = "rf"), size = 3) +  # Points bleus
  xlab("n.train") +  # Légende de l'axe des ordonnées (x)
  ylab("NMAE.mean.mean") +  # Légende de l'axe des abscisses (y)
  ggtitle("PARAM = t124. 
  n.tree=2000 - NMAE.mean.mean=f(n.train) - 7 variable parameters (SNPs)") +  # Titre du graphique
  scale_x_continuous(limits = c(x_min, x_max)) +
  scale_y_continuous(limits = c(y_min, y_max)) + # Réglez l'axe des ordonnées à l'échelle et le range
  scale_x_continuous(limits = c(x_min, x_max), minor_breaks = seq(0, x_max, by = 1000)) +  # Spécifier les marques de l'axe x tous les 1000
  scale_y_continuous(limits = c(y_min, y_max)) + # Réglez l'axe des ordonnées à l'échelle et le range
  scale_color_manual(values = c("rf" = "blue", "drf" = "red")) +
  labs(color = "Method") +
  guides(color = guide_legend(title = "Method")) +
  theme(plot.title = element_text(hjust = 0.5))  # Centrer le titre
print(gg)

#################################  mean Sd ####################################
# Création de graphiques avec des légendes personnalisées (NMAE = f(n.train | n.tree fixed)) - ra.
x_min <- 500  # Valeur minimale pour l'axe x
x_max <- 50000  # Valeur maximale pour l'axe x
y_min <- 0.05 # Valeur minimale pour l'axe y
y_max <- 0.25 # Valeur maximale pour l'axe y
gg <- ggplot(data = sub.RESULTS, aes(x = n.train, y = ra.sd.mean.drf)) +
  geom_point(size = 3, aes(color = "drf")) + 
  geom_point(data = sub.RESULTS, aes(x = n.train, y = ra.sd.mean.rf, color = "rf"), size = 3) +  # Points bleus
  xlab("n.train") +  # Légende de l'axe des ordonnées (x)
  ylab("sd.mean") +  # Légende de l'axe des abscisses (y)
  ggtitle("PARAM = ra. 
  n.tree=2000 - sd.mean=f(n.train) - 7 variable parameters (SNPs)") +  # Titre du graphique
  # scale_x_continuous(limits = c(x_min, x_max)) +
  # scale_y_continuous(limits = c(y_min, y_max)) + # Réglez l'axe des ordonnées à l'échelle et le range
  scale_x_continuous(limits = c(x_min, x_max), minor_breaks = seq(0, x_max, by = 1000)) +  # Spécifier les marques de l'axe x tous les 1000
  scale_y_continuous(limits = c(y_min, y_max)) + # Réglez l'axe des ordonnées à l'échelle et le range
  scale_color_manual(values = c("rf" = "blue", "drf" = "red")) +
  labs(color = "Method") +
  guides(color = guide_legend(title = "Method")) +
  theme(plot.title = element_text(hjust = 0.5))  # Centrer le titre
print(gg)

################### 90% CI ##################################################
# Création de graphiques avec des légendes personnalisées (NMAE = f(n.train | n.tree fixed)) - t124.
x_min <- 500  # Valeur minimale pour l'axe x
x_max <- 50000  # Valeur maximale pour l'axe x
y_min <- 0.85 # Valeur minimale pour l'axe y
y_max <- 1.0 # Valeur maximale pour l'axe y
gg <- ggplot(data = sub.RESULTS, aes(x = n.train, y =  t124.CI90.drf)) +
  geom_point(size = 3, aes(color = "drf")) + 
  geom_point(data = sub.RESULTS, aes(x = n.train, y =  t124.CI90.rf, color = "rf"), size = 3) +  # Points bleus
  xlab("n.train") +  # Légende de l'axe des ordonnées (x)
  ylab("CI90") +  # Légende de l'axe des abscisses (y)
  ggtitle("PARAM = t124.
  n.tree=2000 - CI90% - 7 variable parameters (SNPs)") +  # Titre du graphique
  # scale_x_continuous(limits = c(x_min, x_max)) +
  # scale_y_continuous(limits = c(y_min, y_max)) + # Réglez l'axe des ordonnées à l'échelle et le range
  scale_x_continuous(limits = c(x_min, x_max), minor_breaks = seq(0, x_max, by = 1000)) +  # Spécifier les marques de l'axe x tous les 1000
  scale_y_continuous(limits = c(y_min, y_max)) + # Réglez l'axe des ordonnées à l'échelle et le range
  scale_color_manual(values = c("rf" = "blue", "drf" = "red")) +
  labs(color = "Method") +
  guides(color = guide_legend(title = "Method")) +
  theme(plot.title = element_text(hjust = 0.5))  # Centrer le titre
print(gg)

head(sub.RESULTS[,], n=100)

################ Intersection CI90 + correlations mean et sd ############################
# Création de graphiques avec des légendes personnalisées (correlations mean et sd= f(n.train | n.tree fixed))
gg <- ggplot(data = sub.RESULTS, aes(x = n.train)) +
  geom_point(aes(y = t12.N12.mean.cor, color = "t12.N12.mean.cor", shape = "t12.N12.mean.cor"), size = 3) +
  geom_point(aes(y = t12.N12.sd.cor, color = "t12.N12.sd.cor", shape = "t12.N12.sd.cor"), size = 3) +
  geom_point(aes(y = ra.N3.mean.cor, color = "ra.N3.mean.cor", shape = "ra.N3.mean.cor"), size = 3) +
  geom_point(aes(y = ra.N3.sd.cor, color = "ra.N3.sd.cor", shape = "ra.N3.sd.cor"), size = 3) +
  scale_shape_manual(values = c("t12.N12.mean.cor" = 16, "t12.N12.sd.cor" = 1, 
                                "ra.N3.mean.cor" = 16, "ra.N3.sd.cor" = 1)) +
  scale_color_manual(values = c("t12.N12.mean.cor" = "red", "t12.N12.sd.cor" = "red",
                                "ra.N3.mean.cor" = "green", "ra.N3.sd.cor" = "green")) +
  xlab("n.train") +
  ylab("Cor-mean.or.Cor-sd") +
  ggtitle("n.tree=2000 - Cor-mean.or.Cor-sd=f(n.train) 
          - 7 variable parameters (SNPs)") +
  scale_x_continuous(limits = c(x_min, x_max), minor_breaks = seq(0, x_max, by = 1000)) +
  scale_y_continuous(limits = c(y_min, y_max), minor_breaks = seq(0, y_max, by = 0.1)) +
  labs(color = "", shape = "") +  # Retirer les légendes pour couleur et forme
  guides(color = guide_legend(title = ""), shape = guide_legend(title = "")) +
  theme(plot.title = element_text(hjust = 0.5))
print(gg)

###########################################################################################################
####################### END CALCULS ANNEXES + FIGURE DESIGN FROM RESULTS.list.n.pods and RESULTS  ##########
###########################################################################################################




###########################################################################################################
#######################   SOME MORE CODE                 ##################################################
###########################################################################################################

# About correlations

# Charger les données (remplacer 'mon_dataframe' par votre dataframe)
mon_dataframe <- pred.drf.mean
# Calculer la matrice de corrélation
matrice_correlation <- cor(mon_dataframe)
matrice_correlation
# Afficher un graphique en paires
colnames(param.pods)
pairs(mon_dataframe, 
      panel = function(x, y) {
        points(x, y)
        usr <- par("usr"); on.exit(par(usr))
        par(usr = c(0, 1, 0, 1))
        r <- cor(x, y)
        text(0.5, 0.5, format(r, digits = 2))
      })

plot(param.pods$N12,pred.drf.cor$cor.1.6)
cor(param.pods$N12,pred.drf.cor$cor.1.6)
cor.test(param.pods$N12,pred.drf.cor$cor.1.6)
plot(param.pods$t12,pred.drf.cor$cor.1.6)
cor(param.pods$t12,pred.drf.cor$cor.1.6)
cor.test(param.pods$t12,pred.drf.cor$cor.1.6)

plot(param.pods$t12*1/param.pods$N12,pred.drf.cor$cor.1.6)
cor(param.pods$t12*1/param.pods$N12,pred.drf.cor$cor.1.6)
cor.test(param.pods$t12*1/param.pods$N12,pred.drf.cor$cor.1.6)

plot(param.pods$ra,pred.drf.cor$cor.2.7)
cor(param.pods$ra,pred.drf.cor$cor.2.7)
cor.test(param.pods$ra,pred.drf.cor$cor.2.7)
plot(param.pods$N3,pred.drf.cor$cor.2.7)
cor(param.pods$N3,pred.drf.cor$cor.2.7)
cor.test(param.pods$N3,pred.drf.cor$cor.2.7)


# Some more code
wwAE <- get_sample_weights(model.drf, newdata = stats.pods, num.threads =nbre.threads.used.for.computation)

param = as.data.frame(PARAMS$ra) # only one parameter = ra
param=as.vector(param)
dim(param)
head(param,n=10)

# VI = variableImportance(
#   model.drf,
#   h = NULL,
#   response.scaling = TRUE,
#   type = "raw"
# )

# Importance Variables
Imp.Var = as.data.frame(variable_importance(model.drf))
rownames(Imp.Var)=colnames(stats.pods)
colnames(Imp.Var)="Imp.Var"
head(Imp.Var,n=10)

indices <- order(Imp.Var$Imp.Var, decreasing = TRUE)
Sorted.Imp.Var <- Imp.Var[indices, , drop = FALSE]
head(Sorted.Imp.Var,n=30)





####################################################################################################
####################################################################################################
####################################################################################################
############## NO INFORMATION = START PRECISION METRICS FOR NMAE means metrics COMPUTATIONS  #######
####################################################################################################
####################################################################################################
####################################################################################################
RESULTS.NO.INFO.means <- as.data.frame(matrix(0, nrow = 1 , ncol = dim(param.pods)[2]*2))
dim(RESULTS.NO.INFO.means)
head(RESULTS.NO.INFO.means, n=5)

###########################
cat("NO INFO METRICS - means","\n")
col.target = 0
for (i.param in 1:length(param.list)) 
{
  single.param.pods = as.vector(param.pods[,i.param])
  head(single.param.pods, n=5)
  # Randomization
  single.param.pods.randomized1 <- sample(single.param.pods)
  single.param.pods.randomized2 <- sample(single.param.pods)
  single.param.pods.randomized3 <- sample(single.param.pods)
  length(single.param.pods.randomized3)
  # head(single.param.pods.randomized1, n=5)
  # head(single.param.pods.randomized2, n=5)
  # head(single.param.pods.randomized3, n=5)
  # NMAE normalized by mean
  
  nom.result <- paste0(param.list[i.param], ".NMAE.mean.mean.rf")
  col.target = col.target+1
  NMAE.mean.mean.rf = calculate_NMAE_mean(true_values = single.param.pods.randomized1 , point_estimates = single.param.pods.randomized2)
  assign(nom.result, NMAE.mean.mean.rf)
  cat(nom.result, "=",NMAE.mean.mean.rf ,"\n")
  RESULTS.NO.INFO.means[col.target]=NMAE.mean.mean.rf
  colnames(RESULTS.NO.INFO.means)[col.target] <- nom.result
  
  nom.result <- paste0(param.list[i.param], ".NMAE.mean.mean.drf")
  col.target = col.target+1
  NMAE.mean.mean.drf = calculate_NMAE_mean(true_values = single.param.pods.randomized1 , single.param.pods.randomized3)
  assign(nom.result, NMAE.mean.mean.drf)
  cat(nom.result, "=",NMAE.mean.mean.drf ,"\n")
  RESULTS.NO.INFO.means[col.target]=NMAE.mean.mean.drf
  colnames(RESULTS.NO.INFO.means)[col.target] <- nom.result
  
}

saveRDS(RESULTS.NO.INFO.means, file = "RESULTS.NO.INFO.means.SCENARIO_5POPS_indSNPs_ADMIXTURE_13parameters_SNP-Param.Pods.PRIOR_various.ntree.ntrain.rds")

####################################################################################################
############## NO INFORMATION = END PRECISION METRICS FOR NMAE means metrics COMPUTATIONS    #######
####################################################################################################

####################################################################################################
############## NO INFORMATION = START PRECISION METRICS FOR ind-PODS COMPUTATIONS  #################
####################################################################################################
RESULTS.NO.INFO.npods <- as.data.frame(matrix(0, nrow = dim(param.pods)[1] , ncol = dim(param.pods)[2]*10))
colnames(RESULTS.NO.INFO.npods) = colnames(param.pods)
dim(RESULTS.NO.INFO.npods)
head(RESULTS.NO.INFO.npods, n=5)

###########################
cat("NO INFO METRICS - ipods","\n")
col.target = 0
for (i.param in 1:length(param.list)) 
{
  single.param.pods = as.vector(param.pods[,i.param])
  head(single.param.pods, n=5)
  # Randomization
  single.param.pods.randomized1 <- sample(single.param.pods)
  single.param.pods.randomized2 <- sample(single.param.pods)
  single.param.pods.randomized3 <- sample(single.param.pods)
  length(single.param.pods.randomized3)
  head(single.param.pods.randomized1, n=5)
  head(single.param.pods.randomized2, n=5)
  head(single.param.pods.randomized3, n=5)
  
  # dif.drf.rf for estimation = mean
  nom.result <- paste0(param.list[i.param], ".dif.drf.rf.mean")
  col.target = col.target+1
  dif.drf.rf.mean = calculate_dif.drf.rf.mean(true_values = single.param.pods.randomized1 , point_estimates.drf = single.param.pods.randomized2, point_estimates.rf = single.param.pods.randomized3)
  for (i.pods in 1:n.pods) {
    colnames(RESULTS.NO.INFO.npods)[col.target]=nom.result
    RESULTS.NO.INFO.npods[i.pods,col.target] = dif.drf.rf.mean[i.pods]
  }
  sum_dif.drf.rf.mean = summary(dif.drf.rf.mean)
  ### Affichage
  summary_str <- paste("Min:", format(sum_dif.drf.rf.mean["Min."], digits = 4, nsmall = 3),
                       "1st Qu.:", format(sum_dif.drf.rf.mean["1st Qu."], digits = 4, nsmall = 3),
                       "Median:", format(sum_dif.drf.rf.mean["Median"], digits = 4, nsmall = 3),
                       "Mean:", format(sum_dif.drf.rf.mean["Mean"], digits = 4, nsmall = 3),
                       "3rd Qu.:", format(sum_dif.drf.rf.mean["3rd Qu."], digits = 4, nsmall = 3),
                       "Max:", format(sum_dif.drf.rf.mean["Max."], digits = 4, nsmall = 3), sep = " ")
  cat("Computation done for n.pods of", nom.result, "-->", summary_str, "\n")
  ###
  
  # Relative dif drf for estimation = mean
  nom.result <- paste0(param.list[i.param], ".dif.drf.mean")
  col.target = col.target+1
  dif.drf.mean = calculate_dif.mean(true_values = single.param.pods.randomized1 , point_estimates = single.param.pods.randomized2)
  for (i.pods in 1:n.pods) {
    colnames(RESULTS.NO.INFO.npods)[col.target]=nom.result
    RESULTS.NO.INFO.npods[i.pods,col.target] = dif.drf.mean[i.pods]
  }
  sum_dif.drf.mean = summary(dif.drf.mean)
  ### Affichage
  summary_str <- paste("Min:", format(sum_dif.drf.mean["Min."], digits = 4, nsmall = 3),
                       "1st Qu.:", format(sum_dif.drf.mean["1st Qu."], digits = 4, nsmall = 3),
                       "Median:", format(sum_dif.drf.mean["Median"], digits = 4, nsmall = 3),
                       "Mean:", format(sum_dif.drf.mean["Mean"], digits = 4, nsmall = 3),
                       "3rd Qu.:", format(sum_dif.drf.mean["3rd Qu."], digits = 4, nsmall = 3),
                       "Max:", format(sum_dif.drf.mean["Max."], digits = 4, nsmall = 3), sep = " ")
  ##########################################cat("Computation done for n.pods of", nom.result, "-->", summary_str, "\n")
  ###
  
  # Relative dif rf for estimation = mean
  nom.result <- paste0(param.list[i.param], ".dif.rf.mean")
  col.target = col.target+1
  dif.rf.mean = calculate_dif.mean(true_values = single.param.pods.randomized1 , point_estimates = single.param.pods.randomized2)
  for (i.pods in 1:n.pods) {
    colnames(RESULTS.NO.INFO.npods)[col.target]=nom.result
    RESULTS.NO.INFO.npods[i.pods,col.target] = dif.rf.mean[i.pods]
  }
  sum_dif.rf.mean = summary(dif.rf.mean)
  ### Affichage
  summary_str <- paste("Min:", format(sum_dif.rf.mean["Min."], digits = 4, nsmall = 3),
                       "1st Qu.:", format(sum_dif.rf.mean["1st Qu."], digits = 4, nsmall = 3),
                       "Median:", format(sum_dif.rf.mean["Median"], digits = 4, nsmall = 3),
                       "Mean:", format(sum_dif.rf.mean["Mean"], digits = 4, nsmall = 3),
                       "3rd Qu.:", format(sum_dif.rf.mean["3rd Qu."], digits = 4, nsmall = 3),
                       "Max:", format(sum_dif.rf.mean["Max."], digits = 4, nsmall = 3), sep = " ")
  ##########################################cat("Computation done for n.pods of", nom.result, "-->", summary_str, "\n")
  ###
  
  # dif.drf.rf for estimation = sd
  nom.result <- paste0(param.list[i.param], ".dif.drf.rf.sd")
  col.target = col.target+1
  dif.drf.rf.sd = calculate_dif.drf.rf.sd(true_values = single.param.pods.randomized1 , point_estimates.drf = single.param.pods.randomized2, point_estimates.rf = single.param.pods.randomized3)
  for (i.pods in 1:n.pods) {
    colnames(RESULTS.NO.INFO.npods)[col.target]=nom.result
    RESULTS.NO.INFO.npods[i.pods,col.target] = dif.drf.rf.sd[i.pods]
  }
  sum_dif.drf.rf.sd = summary(dif.drf.rf.sd)
  ### Affichage
  summary_str <- paste("Min:", format(sum_dif.drf.rf.sd["Min."], digits = 4, nsmall = 3),
                       "1st Qu.:", format(sum_dif.drf.rf.sd["1st Qu."], digits = 4, nsmall = 3),
                       "Median:", format(sum_dif.drf.rf.sd["Median"], digits = 4, nsmall = 3),
                       "Mean:", format(sum_dif.drf.rf.sd["Mean"], digits = 4, nsmall = 3),
                       "3rd Qu.:", format(sum_dif.drf.rf.sd["3rd Qu."], digits = 4, nsmall = 3),
                       "Max:", format(sum_dif.drf.rf.sd["Max."], digits = 4, nsmall = 3), sep = " ")
  cat("Computation done for n.pods of", nom.result, "-->", summary_str, "\n")
  
 
}

saveRDS(RESULTS.NO.INFO.npods, file = "RESULTS.NO.INFO.npods.SCENARIO_5POPS_indSNPs_ADMIXTURE_13parameters_SNP-Param.Pods.PRIOR_various.ntree.ntrain.rds")

####################################################################################################
############## NO INFO = END PRECISION METRICS FOR ind-PODS COMPUTATIONS     #######################
####################################################################################################


RESULTS.NO.INFO.npods.RANDOM <- readRDS("RESULTS.NO.INFO.npods.SCENARIO_5POPS_indSNPs_ADMIXTURE_13parameters_SNP-Param.Pods.PRIOR_various.ntree.ntrain.rds")
head(RESULTS.NO.INFO.npods.RANDOM,n=5)
RESULTS.NO.INFO.means.RANDOM <- readRDS("RESULTS.NO.INFO.means.SCENARIO_5POPS_indSNPs_ADMIXTURE_13parameters_SNP-Param.Pods.PRIOR_various.ntree.ntrain.rds")
head(RESULTS.NO.INFO.means.RANDOM,n=5)


}














