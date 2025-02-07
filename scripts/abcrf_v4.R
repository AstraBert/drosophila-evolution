###### RANDOM FOREST MODEL CHOICE - R-SCRIPT FOR DIYABC ############
# Authors: Louis raynal, Jean-Michel Marin, Arnaud Estoup
# Date: 19/06/2018 + 17/08/2024 + 11-10-2024 + 08/11/2024
# Version provided to Astra BERTELLI 08-11-2024
####################################################################

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
#####################################################

###### Loading of the library abcrf and other useful libraries
library(abcrf) # RF
help(package = "abcrf")
library(MASS) # LDA
library(ranger) # CSRF
library(parallel)

#### Various options
options(max.print=10000) # To print tables until 10000 lines
# Number of cores available for your computer
n_cores <- detectCores()
# Nbre of CPU cores I want to use for parallel computation (if not define then all available CPU cores will be used)
how.many.cores.used.for.computation = n_cores-2

###### output text file that will include various numerical outputs
output.file = "diyabc_4_test.txt"
#output.file = "nanne_dataset_pour_graph_lda_OK.txt"
n.run = 1 # Most people do a single RF run (n.run=1)...I prefere to do several runs (e.g. n.run=10) to better appreciate the robustness of my conclusions

###### Key parameters to inform in order to run RF  
# Number of simulations taken from the training dataset (reftable) that will be used to built RF trees
N.train <- 20000
# Number of trees in the forest
ntree <- 1000
# Threshold results (a gadget to use RF output in a specific way = pruning when many scenarios are compared = to be explained later)
results_threshold_fraction_trees = 0.05
# Specific data treatment for large number of statobsRF leading to different winning scenarios
large_number_of_statobsRF_leading_to_different_winning_scenarios = FALSE

###### Give here the names of three key files produced by diyabcRF: (i) the file corresponding to the bin format reference table (with all simulated data summarized with various summary statistics), 
# (ii) the file of the header (with various informations on the simulations), 
# and (iii) the file corresponding54 the observed dataset (summarized with the same set of summary statistics). 
name.reference.table <- "reftableRF.bin"
name.header.file <- "headerRF.txt"
#name.observed.dataset <- "statobsRF_all_statobs_vectors.txt" # Note that this observed dataset may includes several lines (i.e. several vectors of observed data/sumstats) 
name.observed.dataset <- "statobsRF.txt"

###### GROUPING OR NOT GROUPING SCENARIOS IN THE ANALYSIS ?
# Analysis of each scenarios independently ---> grouping.scenarios="NO")
# Analysis of scenarios by groups ---> grouping.scenarios="YES"
grouping.scenarios <- "NO"

if (grouping.scenarios=="YES")
{
  # Definition of the scenarios groups (if grouping = YES)
  group.list = list("1", "2", "11", "12", "13", "14")
  # group.list = list("4","5")
  # group.list = list(c("1", "2", "3", "4", "5"), c("6", "7", "8", "9", "10")) 
  #group.list = list("1","2","3","4")
  n.scen =length(group.list)
}

# Loading data from key files produced by diyabcRF  using the specific fonction of abcrf (reference.table)
# N = total numbre of simuations one wants to load from the reference table
reference.table <- readRefTable(filename = name.reference.table, header=name.header.file)
cat("ciao!")
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
# Check the stat.obs object
cat("Partial display of the statobs file:","\n")
head(stat.obs[,1:5], n=16)

# ####### Exploring the data in the reftable
# # Extraire les valeurs de la colonne "FST_1_4.5"
# fst_values <- reference.table$stats[, "FST_1_4.5"]
# # Vérifier que les données sont numériques pour les analyses
# fst_values <- as.numeric(fst_values)
# # Tracer les valeurs
# plot(fst_values, main = "Plot des valeurs FST_1_4.5", xlab = "Index", ylab = "Valeurs")
# # Calculer et tracer la densité
# fst_density <- density(fst_values)
# plot(fst_density, main = "Densité des valeurs FST_1_4.5", xlab = "Valeurs", ylab = "Densité")
# # Calculer les quantiles
# fst_quantiles <- quantile(fst_values, probs = c(0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95))
# # Afficher les quantiles
# cat("Quantiles des valeurs FST_1_4.5 :\n")
# print(fst_quantiles)
# # Facultatif : Sauvegarder la densité et les quantiles dans des fichiers
# #write.table(data.frame(fst_density), "FST_1_4.5_density.txt", row.names = FALSE)
# write.table(data.frame(Quantiles = fst_quantiles), "FST_1_4.5_quantiles.txt", row.names = TRUE)
# summary(fst_values)
# 
# # Extraire les valeurs de la colonne "NAL_1_4"
# fst_values <- reference.table$stats[, "NAL_1_4"]
# # Vérifier que les données sont numériques pour les analyses
# fst_values <- as.numeric(fst_values)
# # Tracer les valeurs
# plot(fst_values, main = "Plot des valeurs NAL_1_4", xlab = "Index", ylab = "Valeurs")
# # Calculer et tracer la densité
# fst_density <- density(fst_values)
# plot(fst_density, main = "Densité des valeurs NAL_1_4", xlab = "Valeurs", ylab = "Densité")
# # Calculer les quantiles
# fst_quantiles <- quantile(fst_values, probs = c(0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95))
# # Afficher les quantiles
# cat("Quantiles des valeurs NAL_1_4 :\n")
# print(fst_quantiles)
# # Facultatif : Sauvegarder la densité et les quantiles dans des fichiers
# #write.table(data.frame(fst_density), "NAL_1_4_density.txt", row.names = FALSE)
# write.table(data.frame(Quantiles = fst_quantiles), "NAL_1_4_quantiles.txt", row.names = TRUE)
# summary(fst_values)
# 
# # Extraire les valeurs de la colonne "NAL_1_5"
# fst_values <- reference.table$stats[, "NAL_1_5"]
# # Vérifier que les données sont numériques pour les analyses
# fst_values <- as.numeric(fst_values)
# # Tracer les valeurs
# plot(fst_values, main = "Plot des valeurs NAL_1_5", xlab = "Index", ylab = "Valeurs")
# # Calculer et tracer la densité
# fst_density <- density(fst_values)
# plot(fst_density, main = "Densité des valeurs NAL_1_5", xlab = "Valeurs", ylab = "Densité")
# # Calculer les quantiles
# fst_quantiles <- quantile(fst_values, probs = c(0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95))
# # Afficher les quantiles
# cat("Quantiles des valeurs NAL_1_5 :\n")
# print(fst_quantiles)
# # Facultatif : Sauvegarder la densité et les quantiles dans des fichiers
# #write.table(data.frame(fst_density), "NAL_1_5_density.txt", row.names = FALSE)
# write.table(data.frame(Quantiles = fst_quantiles), "NAL_1_5_quantiles.txt", row.names = TRUE)
# summary(fst_values)
# 
# # Extraire les valeurs de la colonne "HET_1_4"
# fst_values <- reference.table$stats[, "HET_1_4"]
# # Vérifier que les données sont numériques pour les analyses
# fst_values <- as.numeric(fst_values)
# # Tracer les valeurs
# plot(fst_values, main = "Plot des valeurs HET_1_4", xlab = "Index", ylab = "Valeurs")
# # Calculer et tracer la densité
# fst_density <- density(fst_values)
# plot(fst_density, main = "Densité des valeurs HET_1_4", xlab = "Valeurs", ylab = "Densité")
# # Calculer les quantiles
# fst_quantiles <- quantile(fst_values, probs = c(0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95))
# # Afficher les quantiles
# cat("Quantiles des valeurs HET_1_4 :\n")
# print(fst_quantiles)
# # Facultatif : Sauvegarder la densité et les quantiles dans des fichiers
# #write.table(data.frame(fst_density), "HET_1_4_density.txt", row.names = FALSE)
# write.table(data.frame(Quantiles = fst_quantiles), "HET_1_4_quantiles.txt", row.names = TRUE)
# summary(fst_values)
# 
# # Extraire les valeurs de la colonne "HET_1_5"
# fst_values <- reference.table$stats[, "HET_1_5"]
# # Vérifier que les données sont numériques pour les analyses
# fst_values <- as.numeric(fst_values)
# # Tracer les valeurs
# plot(fst_values, main = "Plot des valeurs HET_1_5", xlab = "Index", ylab = "Valeurs")
# # Calculer et tracer la densité
# fst_density <- density(fst_values)
# plot(fst_density, main = "Densité des valeurs HET_1_5", xlab = "Valeurs", ylab = "Densité")
# # Calculer les quantiles
# fst_quantiles <- quantile(fst_values, probs = c(0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95))
# # Afficher les quantiles
# cat("Quantiles des valeurs HET_1_5 :\n")
# print(fst_quantiles)
# # Facultatif : Sauvegarder la densité et les quantiles dans des fichiers
# #write.table(data.frame(fst_density), "HET_1_5_density.txt", row.names = FALSE)
# write.table(data.frame(Quantiles = fst_quantiles), "HET_1_5_quantiles.txt", row.names = TRUE)
# summary(fst_values)

sink(file = output.file, split = TRUE)
cat(" ############## MODEL CHOICE analysis using random forest ############### ","\n")
cat("\n")
cat("Name of the statobs file:", name.observed.dataset,"\n")
cat("Number of statobs in the file =", nrow(stat.obs),"\n")
cat("Name of the reference table file:", name.reference.table,"\n")
cat("Name of the headerfile file:", name.header.file,"\n")
cat("Number of simulations loaded from the reference table =",reference.table$nrec, "\n")
cat("Number of scenarios (i.e. models) in the reference table =",reference.table$nscen, "\n")
cat("Number of simulations available for each scenario from the loaded reference table =",reference.table$nrecscen,"\n")
cat("Number of parameters in the reference table =",reference.table$nparam, "\n")
cat("Number of summary statistics in the reference table (without LDA axes) =",ncol(reference.table$stats), "\n")
cat("Number of simulations in the TRAINING DATASET used to built rf trees =", N.train, "\n")
cat("Number of trees in the forest =",ntree, "\n")
cat("Number of cores available = ",n_cores, "\n")
cat("Number of cores used for computation =",how.many.cores.used.for.computation, "\n")

if (grouping.scenarios == "YES") {
  cat("Grouping or selection list of scenarios = ", toString(group.list), "\n")
  cat("Number of ANALYSED scenarios (i.e. models) using RF =",n.scen, "\n")
} else {
  n.scen = reference.table$nscen
  cat("Analysis of each scenario independently (no grouping or selelection of scenarios)","\n")
  cat("Number of ANALYSED scenarios (i.e. models) using RF =",n.scen, "\n")
}
cat("n.run =",n.run,"\n")
cat("\n") 

### List of DataFrames for compilation of outputs over the n.run of nstatobs
# Number of rows in stat.obs
n.statobs <- nrow(stat.obs)
# Initialization of the list of DataFrames
result.allstatobs.allruns <- vector("list", n.run)
# Create column names for the scenarios and the last 3 columns (Prob, Win, Prior_error)
col_names <- c(paste0("S", 1:n.scen), "Prob", "Win", "Prior_error")
result.allstatobs.i.run <- data.frame(matrix(NA, nrow = n.statobs, ncol = n.scen + 3))
colnames(result.allstatobs.i.run) <- col_names
# Loop to fill the list with empty DataFrames
for (i in 1:n.run) result.allstatobs.allruns[[i]] <- result.allstatobs.i.run
# Dataframe for each statobs
result.by.statobs <- vector("list", n.statobs)
for (i in 1:n.statobs) {
  result.by.statobs[[i]] <- data.frame(matrix(NA, nrow = n.run, ncol = n.scen + 3))
  colnames(result.by.statobs[[i]]) <- col_names
}

###### Construction of a data.frame object including the variables of importance for RF analysis 
# that are the scenario index and the summary statistics
reference.table$scenarios <- as.factor(reference.table$scenarios) 
# Conversion as factor
scenario <- reference.table$scenarios[1:N.train]
sumstat <- reference.table$stats[1:N.train,]
dataTrain <- data.frame(scenario, sumstat)

###### Here we go !!! Forest construction (without or with scenario grouping) #################

# We first add five noise variables (randomly drawned into [0;1] uniform distributions) to the forest training dataset 
# to assess which summary statistics bring information
NOISE.datatrain <- matrix(runif(N.train*5), ncol=5)
colnames(NOISE.datatrain) <- c("NOISE1", "NOISE2", "NOISE3", "NOISE4", "NOISE5")
dataTrain <- cbind(dataTrain, NOISE.datatrain)

for (i.run in 1:n.run) {

cat("#############################################################################################","\n")
cat("--------------------> i.run =",i.run,"\n")
cat("\n")

# We can now train the forest

if (grouping.scenarios=="NO") {
  cat("\n")
  cat("RF TRAINING step using all scenarios separately", "\n")
  model.rf <- abcrf(
    scenario ~ ., 
    dataTrain, 
    lda = TRUE, 
    ntree = ntree, 
    ncores=how.many.cores.used.for.computation,
    paral = TRUE, 
    min.node.size = 1, 
    save.memory = FALSE
  )
}

if (grouping.scenarios == "YES") {
  cat("\n")
  cat("RF training step using groups of scenarios", "\n")
  model.rf <- abcrf(
    scenario ~ ., 
    dataTrain, 
    group = group.list,
    lda = TRUE, 
    ntree = ntree, 
    ncores=how.many.cores.used.for.computation,
    paral = TRUE, 
    min.node.size = 1, 
    save.memory = FALSE
  )
}

# Forest first set of outputs (prior error rates and importance measure of the summary statistics as explanatory variables)
cat("\n")
cat("PRIOR ERROR RATES (using the Out-of-Bag statistical procedure): we provide mean Out-of-Bag prior error rate and matrix of errors with the true scenario number as rows and the chosen scenarios as columns","\n")
print(model.rf)
cat("Numerical values of Importance Measures of the summary statistics (i.e. explanatory variables) are stored in the file named importance_measures_of_summary statistics_SORTED.txt","\n")
cat("\n")
# Forest second set of outputs
output <- model.rf$model.rf$variable.importance
decreasing.order <- order(output,decreasing=TRUE)
var.imp.sorted <- data.frame(Importance=output[decreasing.order])
write.table(var.imp.sorted,file="IM_of_SS_SORTED.txt",quote=FALSE,sep="\t\t",row.names=TRUE, col.names="Stat_name           Importance_measure_(as_in_SupMat_of_Raynal_et_al_2018_p17)")

###### Prediction on an observed dataset (cf. stat.obs) ##################
#Note that this observed dataset may includes several lines (i.e. several vectors of observed data/sumstats) 
# We add five noise variables (randomly drawned into [0;1] uniform distributions) to the observed dataset
NOISE.obs <- matrix(runif(5), ncol=5)
colnames(NOISE.obs) <- c("NOISE1", "NOISE2", "NOISE3", "NOISE4", "NOISE5")
stat.obs <- cbind(stat.obs,  NOISE.obs)
# We run the RF prediction on the stat.obs dataset
pred.dataobs <- predict(model.rf, stat.obs, dataTrain, ntree, ncores=how.many.cores.used.for.computation, min.node.size=1, paral.predict = TRUE)

cat("RF MODEL CHOICE PREDICTION on the observed dataset(s): we provide the index of the selected model (i.e. scenario), the votes for each model (over the total number of trees in the forest) and the posterior probability for the selected model","\n")
print(pred.dataobs, row.names = FALSE)
cat("\n")
cat("\n")
#### Pour compilation sur les i.run et les i.statobs
  result.allstatobs.i.run[1:n.scen] <- pred.dataobs$vote
  result.allstatobs.i.run[n.scen + 1] <- pred.dataobs$post.prob
  result.allstatobs.i.run[n.scen + 2] <- pred.dataobs$allocation
  result.allstatobs.i.run[n.scen + 3] <- model.rf$prior.err
  result.allstatobs.allruns[[i.run]] <- result.allstatobs.i.run
}

#sink(file = "output_file_new.txt", split = TRUE) # for testing or independent saving
#############################################################################################
#############################################################################################
cat("\n")
cat("#############################################################################################","\n")
cat("#############################################################################################","\n")
cat("ABSTRACT over the n.run for each vector of STATOBS","\n")
cat("#############################################################################################","\n")
cat("#############################################################################################","\n")
cat("\n")

cat("Recalling what are the different statobs lines", "\n")
fichier <- name.observed.dataset  
lines <- readLines(fichier)
selected.lines <- lines[1:n.statobs]
statobs.selected.lines <- paste0("Statobs_", 1:n.statobs, " = ", selected.lines)
# Writing each statobs line on a different line + without quotes
cat(statobs.selected.lines, sep = "\n")

# Fill up dataframe for all runs of each statobs
for (i.statobs in 1:n.statobs) {
  for (i.run in 1:n.run) result.by.statobs[[i.statobs]][i.run,] =result.allstatobs.allruns[[i.run]][i.statobs,]
}

# Display a summary for each dataframe in the list
cat("\n")
for (i.statobs in 1:n.statobs) {
  cat("Summary for Statobs =", i.statobs, "\n")
  print(summary(result.by.statobs[[i.statobs]]))
  cat("\n")
}
cat("\n")

cat("\n")
# Global results
# Initialize a dataframe to store consolidated results
all.different.winning.scenarios.over.all.statobs <- data.frame(
  Win.Scen.ID = integer(),
  occur = integer(),
  stringsAsFactors = FALSE
)
for (i.statobs in 1:n.statobs) {
  cat("Winner scenarios for Statobs =", i.statobs, "\n")
  print(result.by.statobs[[i.statobs]]$Win, row.names = FALSE)
  cat("\n")
  # Retrieve unique elements and count their occurrences for the current i.statobs
  win_values <- result.by.statobs[[i.statobs]]$Win
  win_counts <- as.data.frame(table(win_values))
  # Rename the columns of the win_counts dataframe
  colnames(win_counts) <- c("Win.Scen.ID", "occur")
  # Convert the Win.Scen.ID column to numeric
  win_counts$Win.Scen.ID <- as.numeric(as.character(win_counts$Win.Scen.ID))
  # Add the results to the global dataframe with aggregation of occurrences
  all.different.winning.scenarios.over.all.statobs <- rbind(all.different.winning.scenarios.over.all.statobs, win_counts)
}
cat("\n")
cat("\n")
# Aggregate the results to have a single row per Win.Scen.ID with the sum of occurrences
all.different.winning.scenarios.over.all.statobs <- aggregate(
  occur ~ Win.Scen.ID, 
  data = all.different.winning.scenarios.over.all.statobs, 
  sum
)
# Display the final consolidated dataframe
print(all.different.winning.scenarios.over.all.statobs)

# RESULTS WITH A THRESHOLD
# Initialize a list to store results by statobs
list.result.threshold <- vector("list", n.statobs)
# Loop through the statobs
for (i.statobs in 1:n.statobs) {
  # Initialize an empty dataframe for each statobs
  result.threshold <- data.frame(
    Scenario = integer(),
    Mean.Votes = numeric(),
    stringsAsFactors = FALSE
  )
  # Loop through the scenarios
  for (scen in 1:n.scen) {
    # Calculate the mean of votes over n.run values for the current scenario
    mean.votes <- mean(result.by.statobs[[i.statobs]][, scen], na.rm = TRUE)
    result.threshold <- rbind(result.threshold, data.frame(Scenario = scen, Mean.Votes = mean.votes))
    # # Check if the mean exceeds the threshold = if sorting on a threshold is needed here... unnecessary since handled below
    # if (mean.votes > threshold) {
    #   # Add the scenario and mean votes to the dataframe
    #   result.threshold <- rbind(result.threshold, data.frame(Scenario = scen, Mean.Votes = mean.votes))
    # }
  }
  # Store the dataframe of each statobs in the list
  list.result.threshold[[i.statobs]] <- result.threshold
}
# Display the list of results
# Iterate through each element of list.result.threshold and display it with a statobs number
cat("\n")
cat("Mean vote numbers for each statobs", "\n")
for (i.statobs in 1:length(list.result.threshold)) {
  cat("\n")
  cat("Statobs #", i.statobs, "\n")
  print(list.result.threshold[[i.statobs]])
}

# FINAL OUTPUT
cat("\n")
# Initialize an empty dataframe to store final results
final.df.over.all.statobs.threshold <- data.frame(
  Scenario = integer(),
  Global.Mean.Votes = numeric(),
  stringsAsFactors = FALSE
)
# Iterate through each i.statobs in list.result.threshold
for (i.statobs in 1:n.statobs) {
  # Add the results of each statobs to the final dataframe
  final.df.over.all.statobs.threshold <- rbind(final.df.over.all.statobs.threshold, list.result.threshold[[i.statobs]])
}

# Calculate the mean of Mean.Votes for each unique scenario
final.df.over.all.statobs.threshold <- aggregate(
  Mean.Votes ~ Scenario,
  data = final.df.over.all.statobs.threshold,
  mean
)
# Affichage des résultats
cat("\nFinal results without threshold:\n")
cat("Global mean number of votes over all statobs:\n")
print(final.df.over.all.statobs.threshold)
# Ajout d'un message final clair
cat("\n### End of Results ###\n")
# Suppression des dépendances à RStudio-specific environments
if (exists(".rs.WorkingDataEnv")) {
  rm(list = ls(envir = .rs.WorkingDataEnv), envir = .rs.WorkingDataEnv)
}
if (exists(".rs.CachedDataEnv")) {
  rm(list = ls(envir = .rs.CachedDataEnv), envir = .rs.CachedDataEnv)
}

cat("#############################################################################################","\n")
cat("#############################################################################################","\n")
cat("########### OTHER PRESENTATION OF FINAL RESULTS WITHOUT ANY THRESHOLD ######################", "\n")
cat("#############################################################################################","\n")
cat("#############################################################################################","\n")

# Initialize two dataframes for mean number of votes and mean fraction of votes
col_names <- c("Scenario", paste0("Statobs", 1:n.statobs))
dataframe.all.results.in.mean.nbre.of.votes <- data.frame(matrix(ncol = n.statobs + 1, nrow = n.scen + 1))
dataframe.all.results.in.mean.fraction.of.votes <- data.frame(matrix(ncol = n.statobs + 1, nrow = n.scen + 1))
colnames(dataframe.all.results.in.mean.nbre.of.votes) <- col_names
colnames(dataframe.all.results.in.mean.fraction.of.votes) <- col_names
# Fill in the "Scenario" column
dataframe.all.results.in.mean.nbre.of.votes$Scenario <- c(1:n.scen, "TotalTrees")
dataframe.all.results.in.mean.fraction.of.votes$Scenario <- c(1:n.scen, "TotalFraction")

# Populate the remaining columns with Mean.Votes data for each statobs
for (i.statobs in 1:n.statobs) {
  # Insert Mean.Votes for the first 28 rows
  dataframe.all.results.in.mean.nbre.of.votes[1:n.scen, paste0("Statobs", i.statobs)] <- list.result.threshold[[i.statobs]]$Mean.Votes
  # Calculate the total for each Statobs column in the mean number of votes dataframe
  dataframe.all.results.in.mean.nbre.of.votes[n.scen + 1, paste0("Statobs", i.statobs)] <- sum(list.result.threshold[[i.statobs]]$Mean.Votes)
  # Insert Mean.Votes as a percentage for the first 28 rows
  dataframe.all.results.in.mean.fraction.of.votes[1:n.scen, paste0("Statobs", i.statobs)] <- list.result.threshold[[i.statobs]]$Mean.Votes / ntree
  # Calculate the total for each Statobs column in the mean fraction of votes dataframe
  dataframe.all.results.in.mean.fraction.of.votes[n.scen + 1, paste0("Statobs", i.statobs)] <- sum(list.result.threshold[[i.statobs]]$Mean.Votes / ntree)
}
# Display the final data frames
cat("\n")
# Afficher les résultats
cat("DataFrame for mean number of votes with TotalTrees row","\n")
print(dataframe.all.results.in.mean.nbre.of.votes)
cat("nDataFrame for mean fraction of votes with Total% row","\n")
print(dataframe.all.results.in.mean.fraction.of.votes)
# Ajout d'un message final clair
cat("\n### End of Results ###\n")
# Vérifications et suppression des environnements spécifiques (optionnel)
if (exists(".rs.WorkingDataEnv") && length(ls(envir = .rs.WorkingDataEnv)) > 0) {
  cat("\nCleaning .rs.WorkingDataEnv...\n")
  rm(list = ls(envir = .rs.WorkingDataEnv), envir = .rs.WorkingDataEnv)
}
if (exists(".rs.CachedDataEnv") && length(ls(envir = .rs.CachedDataEnv)) > 0) {
  cat("\nCleaning .rs.CachedDataEnv...\n")
  rm(list = ls(envir = .rs.CachedDataEnv), envir = .rs.CachedDataEnv)
}

######## Output With threshold
cat("#############################################################################################","\n")
cat("#############################################################################################","\n")
cat("########### OTHER PRESENTATION OF FINAL RESULTS USING A DEFINED THRESHOLD          ##########", "\n")
cat("#############################################################################################","\n")
cat("#############################################################################################","\n")
threshold <- ntree * results_threshold_fraction_trees
cat("ntree =",ntree,"\n")
cat("results_threshold_fraction_trees=",results_threshold_fraction_trees,"\n")
cat("threshold in min number of trees =",threshold,"\n")
cat("\n")
# Filter sub-frames based on thresholds
sub_dataframe_mean_votes <- dataframe.all.results.in.mean.nbre.of.votes[
  apply(dataframe.all.results.in.mean.nbre.of.votes[, -1, drop = FALSE], 1, function(row) any(row >= threshold)), 
]
sub_dataframe_mean_fraction <- dataframe.all.results.in.mean.fraction.of.votes[
  apply(dataframe.all.results.in.mean.fraction.of.votes[, -1, drop = FALSE], 1, function(row) any(row >= results_threshold_fraction_trees)), 
]

# Display results
# Remove last line row (Total)
sub_dataframe_mean_votes <- sub_dataframe_mean_votes[-nrow(sub_dataframe_mean_votes), ]
sub_dataframe_mean_fraction <- sub_dataframe_mean_fraction[-nrow(sub_dataframe_mean_fraction), ]
cat("Sub-dataframe for the scenarios with mean number of votes > ", threshold, "\n")
print(sub_dataframe_mean_votes)
cat("\n")
cat("Sub-dataframe for the scenarios with mean fraction of votes > ", results_threshold_fraction_trees, "\n")
print(sub_dataframe_mean_fraction)


# cat("Sum of the nbre of votes / the fraction of votes for scenarios > threshold", "\n")
cat( "\n")
# Check and select the columns after ‘Scenario’.
cols_to_sum <- sub_dataframe_mean_votes[, -1, drop = FALSE]
# Compute the sum of each column (statobs)
col_sums <- colSums(cols_to_sum)
# Dataframe for results
result <- data.frame(Column = names(col_sums), Sum = col_sums)
cat("Sum of the nbre of votes for scenario relaining after threshold", "\n")
print(result)
cat( "\n")
# Check and select the columns after ‘Scenario’.
cols_to_sum <- sub_dataframe_mean_fraction[, -1, drop = FALSE]
# Compute the sum of each column (statobs)
col_sums <- colSums(cols_to_sum)
# Dataframe for results
result <- data.frame(Column = names(col_sums), Sum = col_sums)
cat("Sum of the nbre of votes for scenario relaining after threshold", "\n")
print(result)

if (large_number_of_statobsRF_leading_to_different_winning_scenarios==TRUE) {
############## Some additional output for specific questions ############################

result.major.class.by.statobs <- numeric(n.statobs)
for (i.statobs in 1:n.statobs) {
  # Récupérer les valeurs de la colonne $Win
  win_values <- result.by.statobs[[i.statobs]]$Win
  # Identifier les fréquences des classes
  freq_table <- table(win_values)
  # Trouver la classe majoritaire (la plus fréquente, ou la plus grande en cas d'ex-æquo)
  major.class <- as.integer(names(freq_table[freq_table == max(freq_table)]))
  major.class <- max(major.class)  # Prendre la plus grande en cas d'ex-æquo
  result.major.class.by.statobs[i.statobs]= major.class
}

for (i.statobs in 1:5) cat("Statobs #",i.statobs," Major class=",result.major.class.by.statobs[i.statobs],"\n")

cat("Recalling what are the different statobs lines in a convenient format", "\n")
fichier <- name.observed.dataset  
lines <- readLines(fichier)
selected.lines <- lines[1:n.statobs]
# Extraire la partie entre ".dat_" et ".txt" dans selected.lines.short
selected.lines.short <- sub(".*\\.dat_(.*?)\\.txt", "\\1", selected.lines)
# Vérifier le résultat
#str(selected.lines.short)

# Initialiser le vecteur n.statobs
n.statobs <- length(result.by.statobs)  # Supposons que n.statobs est basé sur la longueur de result.by.statobs

# Créer un vecteur pour regrouper les informations
all_statobs_major_wins <- vector("character", n.statobs)
# Remplir l'objet all_statobs_major_wins
for (i.statobs in seq_len(n.statobs)) {
  # Récupérer les informations nécessaires
  selected_line <- selected.lines.short[i.statobs]
  major_class <- as.numeric(result.major.class.by.statobs[i.statobs])  # Convertir explicitement en numérique
  win_values <- paste(result.by.statobs[[i.statobs]]$Win, collapse = " ")
  # Combiner les informations sur une même ligne avec 6 espaces entre les éléments
  all_statobs_major_wins[i.statobs] <- paste0(
    selected_line, "    MajorS =", 
    major_class, "     All WinS =",
    win_values
  )
}

for (i.statobs in 1:5) cat("NOT SORTED Statobs ID",i.statobs," Win_pattern:    ",all_statobs_major_wins[i.statobs],"\n")

# Extraire les valeurs de MajorS à partir de chaque chaîne
numeric_values <- sapply(all_statobs_major_wins, function(x) {
  matches <- regmatches(x, regexpr("MajorS =([0-9]+)", x))  # Trouver "MajorS =N"
  as.numeric(sub("MajorS =", "", matches))  # Extraire la valeur numérique
})
# Trier all_statobs_major_wins en fonction des valeurs extraites de MajorS
sorted_indices <- order(numeric_values)  # Indices triés
all_statobs_major_wins_sorted <- all_statobs_major_wins[sorted_indices]

cat("\n")
cat(" ############## Output Winning S ############### ","\n")
cat("\n")
for (i.statobs in 1:n.statobs) cat("SORTED Statobs ID",i.statobs," Win_pattern:    ",all_statobs_major_wins_sorted[i.statobs],"\n")
cat("\n")
#head(all_statobs_major_wins_sorted)

############ Analyses des facteurs impliqués dans Win

# Extraire la première partie (avant "MajorS =") et la valeur de MajorS
data <- data.frame(
  Factors = sapply(all_statobs_major_wins, function(x) {
    # Extraire la partie avant "MajorS ="
    strsplit(trimws(x), "\\s{4,}")[[1]][1]
  }),
  MajorS = as.numeric(sapply(all_statobs_major_wins, function(x) {
    # Extraire la valeur numérique après "MajorS ="
    matches <- regmatches(x, regexpr("MajorS =([0-9]+)", x))
    as.numeric(sub("MajorS =", "", matches))
  }))
)
# Décomposer la colonne Factors en éléments distincts en fonction du symbole "_"
factor_split <- strsplit(as.character(data$Factors), "_")
# Transformer la liste en dataframe avec des colonnes pour chaque facteur
factor_df <- do.call(rbind, factor_split)
factor_df <- as.data.frame(factor_df, stringsAsFactors = TRUE)
# Nommer les colonnes pour les facteurs
colnames(factor_df) <- paste0("Factor", seq_len(ncol(factor_df)))
# Combiner avec MajorS
analysis_data <- cbind(factor_df, MajorS = data$MajorS)
# Vérifier la structure du dataframe final
str(analysis_data)
# Afficher les premières lignes pour vérification
# head(analysis_data)

# Analyse 1 : Régression logistique (si MajorS est catégorique)
# Transformer MajorS en variable binaire : 0 pour 6, 1 pour différent de 6
analysis_data$MajorS_binary <- as.factor(ifelse(analysis_data$MajorS == 6, 0, 1))
# Vérifier les niveaux des facteurs
sapply(analysis_data[, c("Factor1", "Factor2", "Factor3", "Factor4")], nlevels)
# Filtrer les colonnes ayant au moins deux niveaux
valid_factors <- sapply(analysis_data[, c("Factor1", "Factor2", "Factor3", "Factor4")], nlevels) > 1
valid_factor_names <- names(valid_factors[valid_factors])
# Créer une formule basée sur les facteurs valides
formula <- as.formula(paste("MajorS_binary ~", paste(valid_factor_names, collapse = " + ")))
# Modèle de régression logistique avec les facteurs valides
logistic_model <- glm(formula, data = analysis_data, family = binomial)
# Résumé du modèle
summary(logistic_model)

# Analyse 2 - RandomForest
# Charger le package randomForest
library(randomForest)
# Préparer les données
analysis_data$MajorS_factor <- as.factor(analysis_data$MajorS)  # Convertir MajorS en facteur si elle est discrète
# Ajuster un modèle de forêt aléatoire
rf_model <- randomForest(MajorS_factor ~ Factor1 + Factor2 + Factor3 + Factor4,
                         data = analysis_data,
                         ntree = 500,              # Nombre d'arbres
                         mtry = 2,                 # Nombre de variables testées à chaque division
                         importance = TRUE)        # Calcul de l'importance des variables
# Résumé du modèle
print(rf_model)
# Importance des variables
importance_values <- importance(rf_model)

# Ouvrir une nouvelle fenêtre graphique plus grande (uniquement pour RStudio ou un environnement GUI)
dev.new(width = 10, height = 6)  # Ajuster les dimensions si nécessaire
# Graphique de l'importance des facteurs
#varImpPlot(rf_model, main = "Importance des facteurs")
# Spécifier le fichier png et ses dimensions
png(paste0(output.file,"_imporFactRF.png"), width = 800, height = 600)
# Générer le graphique
varImpPlot(rf_model, main = "Importance des facteurs")
# Fermer le périphérique graphique
dev.off()

# Après avoir construit une Random Forest et analysé les effets des facteurs sur 
# une variable comme MajorS, les Partial Dependence Plots (PDP) sont très utiles pour
# comprendre l'impact moyen des facteurs. Cependant, d'autres analyses complémentaires 
# peuvent enrichir votre compréhension, notamment :

# SHAP Values (SHapley Additive exPlanations) :
# Quantifie la contribution de chaque facteur pour chaque prédiction.
# Fournit une explication localisée et globale des prédiction

library(iml)
# Préparer l'explainer
explainer <- Predictor$new(rf_model, data = analysis_data, y = analysis_data$MajorS)
# Calculer les valeurs Shapley pour une observation
shapley <- Shapley$new(explainer, x.interest = analysis_data[1, ])
# Sauvegarder le graphique dans un fichier PNG
png(filename = paste0(output.file, "_Shapley.png"), width = 800, height = 600)
# Générer le graphique
plot(shapley)
# Fermer le périphérique graphique pour enregistrer l'image
dev.off()

# MAIS AUSSI: Résumé des compléments aux PDP :
#   Analyse	Méthode principale	Résultat attendu
# Importance des variables	varImpPlot, SHAP	Identifie les facteurs les plus influents.
# Interaction entre facteurs	PDP (2 facteurs), Interaction Depth	Montre les interactions significatives.
# Effets locaux (par observation)	LIME	Explique pourquoi une prédiction spécifique a été faite.
# Effet conditionnel individuel (ICE)	ICE plots	Examine l’impact des facteurs pour des observations individuelles.
# Analyse des erreurs	Confusion matrix, résidus	Identifie les classes mal prédictes et les erreurs du modèle.
# Importance par sous-groupes	Réentraîner des modèles sur différents groupes	Vérifie si certains facteurs sont plus influents dans des sous-populations spécifiques.
# En combinant ces analyses, vous obtiendrez une compréhension beaucoup plus approfondie des effets des facteurs et des décisions du modèle



# Les **Partial Dependence Plots (PDP)** permettent d’interpréter comment les valeurs d’un facteur (ou d’une paire de facteurs) influencent les prédictions moyennes (`yhat`) du modèle. Voici comment interpréter les valeurs de `yhat` :
#   
#   ---
#   
#   ### **1. Analyse pour un seul facteur**
#   
#   #### **Définition :**
#   - Les valeurs de `yhat` représentent la **prédiction moyenne** du modèle lorsque le facteur étudié est fixé à une valeur donnée, et que toutes les autres variables du modèle sont marginalisées (fixées à leur distribution observée dans les données).
# 
# #### **Interprétation :**
# - **Variation de `yhat`** :
#   - Si `yhat` varie fortement en fonction des valeurs du facteur, cela indique que ce facteur a une influence significative sur la variable cible.
# - Une `yhat` stable signifie que le facteur n’a que peu d’impact.
# - **Direction de l'effet** :
#   - Une augmentation de `yhat` suggère que la valeur du facteur favorise une probabilité ou une réponse plus élevée pour la classe cible (dans le cas de classification) ou pour la variable continue (dans le cas de régression).
#   - Une diminution indique l’effet inverse.
# 
# #### **Exemple d’interprétation :**
# Si le PDP pour `Factor1` montre :
#   - `yhat` = 0.75 pour `Factor1 = "A"`.
#   - `yhat` = 0.50 pour `Factor1 = "B"`.
#   Cela signifie que le modèle prédit en moyenne une probabilité de 75 % pour la classe cible lorsque `Factor1 = "A"`, et de 50 % lorsque `Factor1 = "B"`.
# 
# ---
# 
# ### **2. Analyse pour une paire de facteurs**
# 
# #### **Définition :**
# - Les valeurs de `yhat` représentent la **prédiction moyenne** lorsque les deux facteurs sont fixés à une combinaison donnée de leurs niveaux, tandis que toutes les autres variables du modèle sont marginalisées.
# 
# #### **Interprétation :**
# - **Interaction entre les deux facteurs** :
#   - Si `yhat` varie en fonction des combinaisons des deux facteurs, cela indique une **interaction** entre eux. Une interaction forte signifie que l’effet d’un facteur dépend de la valeur de l’autre.
#   - Si les lignes ou les régions sont parallèles (dans une heatmap ou un graphique 3D), il y a peu ou pas d’interaction.
# - **Synergie ou antagonisme** :
#   - Une augmentation de `yhat` dans certaines combinaisons indique que ces valeurs conjointes des facteurs favorisent une probabilité plus élevée ou une réponse plus forte.
#   - Une diminution indique l’inverse.
# 
# #### **Exemple d’interprétation :**
# Si le PDP pour `Factor1` et `Factor2` montre :
#   - `yhat` = 0.85 pour `Factor1 = "A"` et `Factor2 = "X"`.
#   - `yhat` = 0.60 pour `Factor1 = "A"` et `Factor2 = "Y"`.
#   Cela signifie que le modèle prédit une probabilité moyenne de 85 % pour la classe cible lorsque `Factor1 = "A"` et `Factor2 = "X"`. Cependant, la probabilité diminue à 60 % lorsque `Factor2` passe à `Y`, indiquant une interaction.
# 
# ---
# 
# ### **Conseils pour interpréter les PDP :**
# 
# 1. **Graphique pour un facteur** :
#    - Le graphique typique montre `yhat` en ordonnée (axe y) et les valeurs du facteur en abscisse (axe x).
#    - Analysez la forme (montée, descente, plate) pour comprendre l’effet du facteur.
# 
# 2. **Graphique pour une paire de facteurs** :
#    - **Heatmap** :
#      - Les couleurs indiquent les valeurs de `yhat` (plus clair = plus élevé).
#      - Regardez les combinaisons où `yhat` est le plus élevé ou le plus bas.
#    - **Surface 3D** :
#      - Examinez les pics ou les creux pour identifier les interactions.
# 
# ---
# 
# ### **Résumé :**
# 
# - **Un facteur** :
#   - Variation de `yhat` : Impact direct du facteur.
#   - Forme du graphique : Direction et amplitude de l’effet.
# - **Deux facteurs** :
#   - Changements de `yhat` en fonction des combinaisons : Interaction entre facteurs.
#   - Patterns (heatmap ou surface) : Synergie ou antagonisme.
# 
# Ces analyses vous aident à comprendre le comportement du modèle et à identifier les facteurs clés influençant la variable cible. 😊

# library(ggplot2)
# library(pdp)
# # Boucle pour générer un PDP pour chaque facteur dans factor_df
# for (factor_name in colnames(factor_df)) {
#   
#   # Créer une grille de prédiction pour le facteur actuel
#   pred_grid <- data.frame(
#     factor = levels(analysis_data[[factor_name]])
#   )
#   colnames(pred_grid) <- factor_name  # Nommer correctement la colonne
#   
#   # Générer le PDP pour le facteur actuel
#   pdp_result <- partial(
#     object = rf_model,
#     pred.var = factor_name,
#     pred.grid = pred_grid,
#     train = analysis_data,
#     plot = TRUE,
#     plot.engine = "ggplot2"
#   )
#   
#   # Ajouter un titre spécifique au graphique
#   print(paste("PDP pour", factor_name))
#   print(pdp_result)
# }
# 
# # Boucle pour générer un PDP pour chaque paire de facteurs
# factor_names <- colnames(factor_df)
# n_factors <- length(factor_names)
# 
# for (i in 1:(n_factors - 1)) {
#   for (j in (i + 1):n_factors) {
#     factor1 <- factor_names[i]
#     factor2 <- factor_names[j]
#     
#     # Créer une grille de prédiction pour les deux facteurs
#     interaction_grid <- expand.grid(
#       Factor1 = levels(analysis_data[[factor1]]),
#       Factor2 = levels(analysis_data[[factor2]])
#     )
#     colnames(interaction_grid) <- c(factor1, factor2)  # Nommer correctement les colonnes
#     
#     # Générer le PDP pour les deux facteurs
#     pdp_result <- partial(
#       object = rf_model,
#       pred.var = c(factor1, factor2),
#       pred.grid = interaction_grid,
#       train = analysis_data,
#       plot = TRUE,
#       plot.engine = "ggplot2"
#     )
#     
#     # Ajouter un titre spécifique au graphique
#     print(paste("PDP pour l'interaction entre", factor1, "et", factor2))
#     print(pdp_result)
#   }
# }

library(ggplot2)
library(pdp)

# Boucle pour générer un PDP pour chaque facteur dans factor_df
for (factor_name in colnames(factor_df)) {
  # Créer une grille de prédiction pour le facteur actuel
  pred_grid <- data.frame(
    factor = levels(analysis_data[[factor_name]])
  )
  colnames(pred_grid) <- factor_name  # Nommer correctement la colonne
  
  # Générer le PDP pour le facteur actuel
  pdp_result <- partial(
    object = rf_model,
    pred.var = factor_name,
    pred.grid = pred_grid,
    train = analysis_data,
    plot = TRUE,
    plot.engine = "ggplot2"
  )
  file_name <- paste0(output.file, "_pdp_", factor_name, ".png")
  ggsave(file_name, plot = pdp_result, width = 8, height = 6)
  print(paste("PDP sauvegardé pour", factor_name, "dans", file_name))
}

# Boucle pour générer un PDP pour chaque paire de facteurs
factor_names <- colnames(factor_df)
n_factors <- length(factor_names)
for (i in 1:(n_factors - 1)) {
  for (j in (i + 1):n_factors) {
    factor1 <- factor_names[i]
    factor2 <- factor_names[j]
    # Créer une grille de prédiction pour les deux facteurs
    interaction_grid <- expand.grid(
      Factor1 = levels(analysis_data[[factor1]]),
      Factor2 = levels(analysis_data[[factor2]])
    )
    colnames(interaction_grid) <- c(factor1, factor2)  # Nommer correctement les colonnes
    # Générer le PDP pour les deux facteurs
    pdp_result <- partial(
      object = rf_model,
      pred.var = c(factor1, factor2),
      pred.grid = interaction_grid,
      train = analysis_data,
      plot = TRUE,
      plot.engine = "ggplot2"
    )
    # Sauvegarder la figure dans un fichier
    file_name <- paste0(output.file, "_pdp_", factor1, "_", factor2, ".png")
    ggsave(file_name, plot = pdp_result, width = 8, height = 6)
    print(paste("PDP sauvegardé pour l'interaction entre", factor1, "et", factor2, "dans", file_name))
  }
}

# Fermer proprement les sink(s) ouverts
if (sink.number() > 0) {
  sink()
}
# Nettoyage des environnements spécifiques à RStudio
cat("Attempting cleanup of RStudio-specific environments...\n")
tryCatch({
  if (exists(".rs.WorkingDataEnv") && is.environment(.rs.WorkingDataEnv)) {
    if (length(ls(envir = .rs.WorkingDataEnv)) > 0) {
      rm(list = ls(envir = .rs.WorkingDataEnv), envir = .rs.WorkingDataEnv)
    }
  }
}, error = function(e) {
  cat("Warning: Failed to clean .rs.WorkingDataEnv. Continuing execution...\n")
})
tryCatch({
  if (exists(".rs.CachedDataEnv") && is.environment(.rs.CachedDataEnv)) {
    if (length(ls(envir = .rs.CachedDataEnv)) > 0) {
      rm(list = ls(envir = .rs.CachedDataEnv), envir = .rs.CachedDataEnv)
    }
  }
}, error = function(e) {
  cat("Warning: Failed to clean .rs.CachedDataEnv. Continuing execution...\n")
})
# Message final clair
cat("\n### End of Results ###\n")

}





#########################################  ################################################
# Figures error et graph_lda ##############################################################
#########################################  ################################################
#########################################  ################################################
cat("#############################################################################################","\n")
cat("Drawing figures for the last of the i.run =",i.run,"\n")
cat("#############################################################################################","\n")

if (grouping.scenarios=="NO")
{
  ########### FIGURE OUTPUT ###################################################################################
  #TO BE USED WHEN Analysis of each scenarios independently ---> grouping.scenarios="NO")
  #############################################################################################################
  ## Graphic providing prior error rates for forests with different number of trees (pdf file name = error_vs_ntree.pdf")
  ## computed using an Out-of-Bag procedure, e.g. Fig. 3 in Pudlo et al. 2016
  # Adjust graphical margins and disable the "Press <ENTER> to Continue..." prompt
  par(mar = c(3, 2, 2, 1) + 0.1, ask = FALSE)
  # First, produce two figures (LDA and VarImp) in the Plots window and automatically save them as PDF files
  # Disable the "Press <ENTER> to Continue" prompt
  # par(ask = FALSE)
  original_readline <- readline
  readline <- function(prompt = "") {
    if (grepl("Press <ENTER> to Continue", prompt)) {
      return("")
    } else {
      original_readline(prompt)
    }
  }
  plot(model.rf, dataTrain, obs = stat.obs, pdf = TRUE, n.var = 30)
  
  ########## If "margin problem", the graphs are still produced; THEN PROCEED BY RUNNING THIS PART OF THE CODE !!!!!
  # Save the error = f(n.tree) plot in a PNG file
  png("err_oob_plot.png", width = 800, height = 600)
  # Compute err.oob
  err.oob <- err.abcrf(model.rf, dataTrain, paral = TRUE)
  # Produce the plot
  plot(err.oob)
  # Close the PNG file to save the plot
  dev.off()  # Important to close the device so it returns to the Plots window
  cat("THREE ILLUSTRATIVE GRAPHICS have been produced and saved in three different files:", quote=FALSE)
  cat("File named prior_errors_vs_number_of_trees.png: Graphic providing prior error rates for forests with different number of trees and computed using an Out-of-Bag procedure, e.g. Fig. 3 in Pudlo et al. 2016", "\n")
  cat("File named graph_lda.pdf = LDA projections of the reference table for the different scenarios plus the observed dataset cf. black star in the figure", "\n")
  cat("File name graph_varImpPlot.pdf = the contributions of the 30 most important statistics to the RF (e.g. Fig. S6 and Fig. S7 in Pudlo et al. 2016)", "\n")
  # #####################################################################################
  
  
  # ############################### START RENAME OUTPUT FILES ######################################################
  # Give a specific name or use the one you defined for your analysis
  #output.file.save <- "abcrf_LUIS_USA__GAUT_vs_FRAI_S1&S2_tok_sap_ok_t1000_s20000_pS_10runs.txt"
  output.file.save <- output.file
  
  # Extract the name without extension and the extension of output.file
  nse_output.file.save <- tools::file_path_sans_ext(output.file)
  ext_output.file.save <- tools::file_ext(output.file)
  
  # List of output files to rename
  files_to_rename <- c("graph_lda.pdf",
                       "graph_varImpPlot.pdf",
                       "err_oob_plot.png",
                       "IM_of_SS_SORTED.txt")
  # Rename each output file
  for (file in files_to_rename) {
    # Extract the name without extension and the extension of the output file
    nse4f <- tools::file_path_sans_ext(file)
    ext4f <- tools::file_ext(file)
    # Create the new file name
    new_name <- paste0(nse4f, "_", nse_output.file.save, ".", ext4f)
    # Rename the file
    file.rename(file, new_name)
  }
  # ############################### END RENAME OUTPUT FILES ######################################################
}

if (grouping.scenarios=="YES")
{
  ########### FIGURE OUTPUT BIS ###############################################################################
  #TO BE USED WHEN Analysis groups of scenarios ---> grouping.scenarios="YES")
  #############################################################################################################
  ## Graphic providing prior error rates for forests with different number of trees (pdf file name = error_vs_ntree.pdf")
  ## computed using an Out-of-Bag procedure, e.g. Fig. 3 in Pudlo et al. 2016
  # Adjust graphical margins and disable the "Press <ENTER> to Continue..." prompt
  par(mar = c(3, 2, 2, 1) + 0.1, ask = FALSE)
  # First, produce two figures (LDA and VarImp) in the Plots window and automatically save them as PDF files
  # Disable the "Press <ENTER> to Continue..." prompt
  # dataTrain.SubScen <- dataTrain[dataTrain$scenario %in% c(1,2), ]
  dataTrain.SubScen <- dataTrain[dataTrain$scenario %in% group.list, ]
  # Renumber the row indices
  rownames(dataTrain.SubScen) <- NULL
  # Verify the result
  # head(dataTrain.SubScen)
  # First, produce two figures (LDA and VarImp) in the Plots window and automatically save them as PDF files
  # Disable the "Press <ENTER> to Continue" prompt
  # par(ask = FALSE)
  original_readline <- readline
  readline <- function(prompt = "") {
    if (grepl("Press <ENTER> to Continue", prompt)) {
      return("")
    } else {
      original_readline(prompt)
    }
  }
  plot(model.rf, dataTrain.SubScen, obs = stat.obs, pdf = TRUE, n.var = 30)
  
  ########## WARNING: If "margin problem", the graphs are still produced; TO FINISH THE PROCESS THEN JUST RUN THIS BELOW SECTION OF THE CODE !!!!!
  # Save the error = f(n.tree) plot in a PNG file
  png("err_oob_plot.png", width = 800, height = 600)
  # Compute err.oob
  err.oob <- err.abcrf(model.rf, dataTrain.SubScen, paral = TRUE)
  # Produce the plot
  plot(err.oob)
  # Close the PNG file to save the plot
  dev.off()  # Important to close the device so it returns to the Plots window
  cat("THREE ILLUSTRATIVE GRAPHICS have been produced and saved in three different files:", quote=FALSE)
  cat("File named prior_errors_vs_number_of_trees.png: Graphic providing prior error rates for forests with different number of trees and computed using an Out-of-Bag procedure, e.g. Fig. 3 in Pudlo et al. 2016", "\n")
  cat("File named graph_lda.pdf = LDA projections of the reference table for the different scenarios plus the observed dataset cf. black star in the figure", "\n")
  cat("File name graph_varImpPlot.pdf = the contributions of the 30 most important statistics to the RF (e.g. Fig. S6 and Fig. S7 in Pudlo et al. 2016)", "\n")
  # #####################################################################################
  # ############################### START RENAME OUTPUT FILES ######################################################
  # Give a specific name or use the one you defined for your analysis
  #output.file.save <- "abcrf_LUIS_USA__GAUT_vs_FRAI_S1&S2_tok_sap_ok_t1000_s20000_pS_10runs.txt"
  output.file.save <- output.file
  
  # Extract the name without extension and the extension of output.file
  nse_output.file.save <- tools::file_path_sans_ext(output.file)
  ext_output.file.save <- tools::file_ext(output.file)
  
  # List of output files to rename
  files_to_rename <- c("graph_lda.pdf",
                       "graph_varImpPlot.pdf",
                       "err_oob_plot.png",
                       "IM_of_SS_SORTED.txt")
  # Rename each output file
  for (file in files_to_rename) {
    # Extract the name without extension and the extension of the output file
    nse4f <- tools::file_path_sans_ext(file)
    ext4f <- tools::file_ext(file)
    # Create the new file name
    new_name <- paste0(nse4f, "_", nse_output.file.save, ".", ext4f)
    # Rename the file
    file.rename(file, new_name)
  }
  # ############################### END RENAME OUTPUT FILES ######################################################
}
