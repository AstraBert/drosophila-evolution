###### RANDOM FOREST MODEL CHOICE - R-SCRIPT FOR DIYABC ############
# Authors: Louis raynal, Jean-Michel Marin, Arnaud Estoup
# Date: 19/06/2018 + 17/08/2024 + 11-10-2024 + 08/11/2024
# Version provided to Astra BERTELLI 08-11-2024

###### Loading of the library abcrf and other useful libraries
library(abcrf) # RF
help(package = "abcrf")
library(MASS) # LDA
library(ranger) # CSRF

# Various options
options(max.print=10000) # Pour imprimer des tables jusque 10000 lignes
# Nbre of CPU cores used for parallel computation (if not define then all available CPU cores will be used)
ncores = 32

###### output text file that will include various numerical outputs
output.file = "abcrf_LUIS_ARG_bestOf_NO_GHOST_versus_two_GHOST_S_3scen_t1000_s20000.txt"
n.run = 10 # Most people do a single RF run (n.run=1)...I prefere to do several runs (e.g. n.run=10) to better appreciate the robustness of my conclusions

###### Key parameters to inform to run RF  
# Number of simulations taken from the training dataset (reftable) that will be used to built RF trees
N.train <- 60000
# Number of trees in the forest (ntree)
ntree <- 1000
# Threshold results (a gadget to use RF output in a specific way = prunning when manay scenarios are compared = to be explained later)
results_threshold_fraction_trees = 0.05

###### Reading of the three key files produced by diyabcRF: (i) the file corresponding to the bin format reference table (with all simulated data summarized with various summary statistics), 
# (ii) the file of the header (with various informations on the simulations), 
# and (iii) the file corresponding the observed dataset (summarized with the same set of summary statistics). 
# N = total numbre of simuations one wants to load from the reference table
name.reference.table <- "reftableRF.bin"
name.header.file <- "headerRF.txt"
name.observed.dataset <- "statobsRF.txt" # Note that this observed dataset may includes several lines (i.e. several vectors of observed data/sumstats) 
#name.observed.dataset <- "statobsRF1234_saplia_sapnin_toknin_toklia.txt"

# Loading data using the scpecific fonction of abcrf
reference.table <- readRefTable(filename = name.reference.table, header=name.header.file) #, N=100000)
stat.obs <- read.table(name.observed.dataset, header=TRUE)
n.scen = reference.table$nscen

# GROUPING OR NOT GROUPING SCENARIOS IN THE ANALYSIS ?
# Analysis of each scenarios independently ---> grouping.scenarios="NO")
# Analysis of scenarios by groups ---> grouping.scenarios="YES"
grouping.scenarios <- "NO"

if (grouping.scenarios=="YES")
{
# Definition of the scenarios groups (if grouping = YES)
group.list = list("8", "10")
#group.list = list("7","9")
#group.list = list(c("1", "2", "3", "4"), c("5", "6", "7", "8")) 
#group.list = list("1","2","3","4")
n.scen =length(group.list)
}

sink(file = output.file, split = TRUE)
cat(" ############## MODEL CHOICE analysis using random forest ############### ","\n")
cat("\n")
cat("Name of the statobs file:", name.observed.dataset,"\n")
cat("Number of statobs in the file =", nrow(stat.obs),"\n")
cat("Name of the reference table file:", name.reference.table,"\n")
cat("Name of the headerfile file:", name.header.file,"\n")
cat("Number of simulations loaded from the reference table =",reference.table$nrec, "\n")
cat("Number of scenarios (i.e. models) in the reference table =",reference.table$nscen, "\n")
cat("Number of ANALYSED scenarios (i.e. models) =",n.scen, "\n")
cat("Number of simulations available for each scenario from the loaded reference table =",reference.table$nrecscen,"\n")
cat("Number of parameters in the reference table =",reference.table$nparam, "\n")
cat("Number of summary statistics in the reference table (without LDA axes) =",ncol(reference.table$stats), "\n")
cat("Number of simulations in the TRAINING DATASET used to built rf trees =", N.train, "\n")
cat("Number of trees in the forest =",ntree, "\n")
if (grouping.scenarios == "YES") {
  cat("Grouping list = ", toString(group.list), "\n")
} else {
  cat("Analysis of each scenario independently (no grouping of scenarios)","\n")
}
cat("n.run =",n.run,"\n")
cat("\n") 

### List of DataFrames for compilation of outputs over the n.run of nstatobs
# Number of rows in stat.obs
n.statobs <- nrow(stat.obs)
# Initialization of the list of DataFrames
result.by.statobs <- vector("list", n.statobs)
# Create column names for the scenarios and the last 3 columns (Prob, Win, Prior_error)
col_names <- c(paste0("S", 1:n.scen), "Prob", "Win", "Prior_error")

# Loop to fill the list with empty DataFrames
for (i in 1:n.statobs) {
  df.result.nrun <- data.frame(matrix(NA, nrow = n.run, ncol = n.scen + 3))
  colnames(df.result.nrun) <- col_names
  result.by.statobs[[i]] <- df.result.nrun
}

###### Construction of a data.frame object including the variables of importance for RF analysis that are the scenario index and the summary statistics
reference.table$scenarios <- as.factor(reference.table$scenarios) 
# Conversion as factor
scenario <- reference.table$scenarios[1:N.train]
sumstat <- reference.table$stats[1:N.train,]
dataTrain <- data.frame(scenario, sumstat)

###### Forest construction (without or with scenario grouping) #################

# We add five noise variables (randomly drawned into [0;1] uniform distributions) to the forest training dataset to assess which summary statistics bring information
NOISE.datatrain <- matrix(runif(N.train*5), ncol=5)
colnames(NOISE.datatrain) <- c("NOISE1", "NOISE2", "NOISE3", "NOISE4", "NOISE5")
dataTrain <- cbind(dataTrain, NOISE.datatrain)

for (i.run in 1:n.run) {

cat("--------------------> i.run =",i.run,"\n")
cat("\n")


# We can now train the forest
if (grouping.scenarios=="NO") {
  cat("\n")
  cat("RF TRAINING step using all scenarios separately", "\n")
  model.rf <- abcrf(scenario~., dataTrain, lda= TRUE, ntree=ntree, paral=TRUE, min.node.size=1, save.memory=FALSE)
}
if (grouping.scenarios == "YES") {
  cat("\n")
  cat("RF training step using groups of scenarios", "\n")
  # model.rf <- abcrf(scenario~., dataTrain, group = list("2", "3"), lda = TRUE, ntree = ntree, paral = TRUE, min.node.size = 1, save.memory = FALSE)
  model.rf <- abcrf(
    scenario ~ ., 
    dataTrain, 
    group = group.list,
    #group = list("7","9"),
    #group = list(c("1", "2", "3", "4"), c("5", "6", "7", "8")), 
    #group = list("1","2","3","4"),
    lda = TRUE, 
    ntree = ntree, 
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
write.table(var.imp.sorted,file="Importance_Measures_of_SUMSTATS_SORTED.txt",quote=FALSE,sep="\t\t",row.names=TRUE, col.names="Statistics_name           Importance_measure_(as_in_SupMat_of_Raynal_et_al_2018_p17)")

###### Prediction on an observed dataset (cf. stat.obs) ##################
#Note that this observed dataset may includes several lines (i.e. several vectors of observed data/sumstats) 
# We add five noise variables (randomly drawned into [0;1] uniform distributions) to the observed dataset
NOISE.obs <- matrix(runif(5), ncol=5)
colnames(NOISE.obs) <- c("NOISE1", "NOISE2", "NOISE3", "NOISE4", "NOISE5")
stat.obs <- cbind(stat.obs,  NOISE.obs)
# We run the RF prediction on the stat.obs dataset
pred.dataobs <- predict(model.rf, stat.obs, dataTrain, ntree, paral = TRUE, min.node.size=1, paral.predict = TRUE)

cat("RF MODEL CHOICE PREDICTION on the observed dataset(s): we provide the index of the selected model (i.e. scenario), the votes for each model (over the total number of trees in the forest) and the posterior probability for the selected model","\n")
print(pred.dataobs, row.names = FALSE)
cat("\n")
cat("#############################################################################################","\n")
cat("\n")

#### Pour compilation sur les i.run et les i.statobs
for (i.statobs in 1:n.statobs) {
  # Récupérer la structure du dataframe pour le statobs courant
  df.result.nrun <- result.by.statobs[[i.statobs]]
  
  # Remplir la ligne i.run du dataframe df.result.nrun avec les valeurs de pred.dataobs$vote
  df.result.nrun[i.run, 1:n.scen] <- pred.dataobs$vote[i.statobs,1:n.scen]
  df.result.nrun[i.run, n.scen + 1] <- pred.dataobs$post.prob[[i.statobs]]
  df.result.nrun[i.run, n.scen + 2] <- pred.dataobs$allocation[[i.statobs]]
  df.result.nrun[i.run, n.scen + 3] <- model.rf$prior.err
  
  # Sauvegarder les modifications dans la liste
  result.by.statobs[[i.statobs]] <- df.result.nrun
  }
}
# sink(file = "output.file.new..txt", split = TRUE) # for testing or independent saving
#############################################################################################
#############################################################################################
cat("\n")
cat("#############################################################################################","\n")
cat("#############################################################################################","\n")
cat("ABSTRACT over the n.run for each vector of STATOBS","\n")
cat("#############################################################################################","\n")
cat("#############################################################################################","\n")
cat("\n")

# Display a summary for each dataframe in the list
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
# Display the final dataframe
cat( "\n")
cat("Global mean number of votes over all statobs","\n")
print(final.df.over.all.statobs.threshold)

########### FINAL RESULTS #################
cat( "\n")
cat( "\n")
######## No threshold

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
dataframe.all.results.in.mean.fraction.of.votes$Scenario <- c(1:n.scen, "Total%")

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
print("DataFrame for mean number of votes with TotalTrees row:")
cat("\n")
print(dataframe.all.results.in.mean.nbre.of.votes)
cat("\n")
print("DataFrame for mean fraction of votes with Total% row:")
cat("\n")
print(dataframe.all.results.in.mean.fraction.of.votes)
cat("\n")
cat("\n")

######## With threshold
cat("#############################################################################################","\n")
cat("#############################################################################################","\n")
cat("########### OTHER PRESENTATION OF FINAL RESULTS WITH inititialy defined THRESHOLD ##########", "\n")
cat("#############################################################################################","\n")
cat("#############################################################################################","\n")

threshold <- ntree * results_threshold_fraction_trees
cat("ntree =",ntree,"\n")
cat("results_threshold_fraction_trees=",results_threshold_fraction_trees,"\n")
cat("threshold in min number of trees =",threshold,"\n")
cat("\n")
# Define the indices of the columns to analyze
start_col <- 2  # The first column to include
end_col <- min(n.statobs + 1, ncol(nbre_df))  # Limit to the last available column
# Apply the filter
if (start_col <= end_col) {
  filtered_nbre_rows <- nbre_df[
    apply(nbre_df[, start_col:end_col, drop = FALSE], 2, function(row) any(row > threshold)), 
  ]
} else {
  warning("No valid columns to filter. Check the value of n.statobs.")
  filtered_nbre_rows <- data.frame()  # Return an empty dataframe if the range is invalid
}
# Display the results
print(filtered_nbre_rows)
# Remove the last row from the original dataframes
nbre_df <- dataframe.all.results.in.mean.nbre.of.votes[1:n.scen, ]
percentage_df <- dataframe.all.results.in.mean.fraction.of.votes[1:n.scen,]

# Define the indices of the columns to analyze
start_col <- 2  # The first column to include
end_col <- min(n.statobs + 1, ncol(percentage_df))  # Limit to the last available column
# For the fraction of votes
if (start_col <= end_col) {
  filtered_percentage_rows <- percentage_df[
    apply(percentage_df[, start_col:end_col, drop = FALSE], 2, function(row) any(row > results_threshold_fraction_trees)), 
  ]
} else {
  warning("No valid columns to filter. Check the value of n.statobs.")
  filtered_percentage_rows <- data.frame()  # Return an empty dataframe if the range is invalid
}
# Display the results
print(filtered_percentage_rows)
# Add total row at the end of each filtered dataframe
if (nrow(filtered_nbre_rows) > 0) {
  # Adjust column selection to avoid invalid dimensions
  start_col <- 2
  end_col <- min(n.statobs + 1, ncol(filtered_nbre_rows))
  
  # Calculate column sums only if there are valid columns
  if (start_col <= end_col) {
    total_nbre <- colSums(filtered_nbre_rows[, start_col:end_col, drop = FALSE])
    filtered_nbre_rows <- rbind(filtered_nbre_rows, c("Total", total_nbre))
  } else {
    warning("No valid columns for summation in filtered_nbre_rows.")
  }
}
if (nrow(filtered_percentage_rows) > 0) {
  # Adjust column selection to avoid invalid dimensions
  start_col <- 2
  end_col <- min(n.statobs + 1, ncol(filtered_percentage_rows))
  
  # Calculate column sums only if there are valid columns
  if (start_col <= end_col) {
    total_percentage <- colSums(filtered_percentage_rows[, start_col:end_col, drop = FALSE])
    filtered_percentage_rows <- rbind(filtered_percentage_rows, c("Total", total_percentage))
  } else {
    warning("No valid columns for summation in filtered_percentage_rows.")
  }
}
# Set row names sequentially
if (nrow(filtered_nbre_rows) > 0) {
  rownames(filtered_nbre_rows) <- 1:nrow(filtered_nbre_rows)
}

if (nrow(filtered_percentage_rows) > 0) {
  rownames(filtered_percentage_rows) <- 1:nrow(filtered_percentage_rows)
}

# Create final dataframes
dataframe.all.results.in.mean.nbre.of.votes.threshold <- filtered_nbre_rows
dataframe.all.results.in.mean.fraction.of.votes.threshold <- filtered_percentage_rows

# Display the final dataframes
cat("\n")
print("Filtered DataFrame for mean number of votes with threshold and totals:")
cat("\n")
print(dataframe.all.results.in.mean.nbre.of.votes.threshold)
cat("\n")
print("Filtered DataFrame for mean fraction of votes with threshold and totals:")
cat("\n")
print(dataframe.all.results.in.mean.fraction.of.votes.threshold)

# Clean escaped from the output file
sink()

#########################################  ################################################
#########################################  ################################################
# Figures error et graph_lda ##############################################################
#########################################  ################################################
#########################################  ################################################
cat("#############################################################################################","\n")
cat("Drawing figures for the last of the i.run =",i.run,"\n")
cat("#############################################################################################","\n")

########### FIGURE OUTPUT ###################################################################################
#TO BE USED WHEN Analysis of each scenarios independently ---> grouping.scenarios="NO")
#############################################################################################################
## Graphic providing prior error rates for forests with different number of trees (pdf file name = error_vs_ntree.pdf")
## computed using an Out-of-Bag procedure, e.g. Fig. 3 in Pudlo et al. 2016
# Adjust graphical margins and disable the "Press <ENTER> to Continue..." prompt
par(mar = c(3, 2, 2, 1) + 0.1, ask = FALSE)
# First, produce two figures (LDA and VarImp) in the Plots window and automatically save them as PDF files
# Disable the "Press <ENTER> to Continue..." prompt
plot(model.rf, dataTrain, obs = stat.obs, pdf = TRUE, n.var = 30)

########## If "margin problem", the graphs are still produced; proceed with
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
                     "Importance_Measures_of_SUMSTATS_SORTED.txt")
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


# ########### FIGURE OUTPUT BIS ###############################################################################
# #TO BE USED WHEN Analysis groups of scenarios ---> grouping.scenarios="YES")
# #############################################################################################################
# ## Graphic providing prior error rates for forests with different number of trees (pdf file name = error_vs_ntree.pdf")
# ## computed using an Out-of-Bag procedure, e.g. Fig. 3 in Pudlo et al. 2016
# # Adjust graphical margins and disable the "Press <ENTER> to Continue..." prompt
# par(mar = c(3, 2, 2, 1) + 0.1, ask = FALSE)
# # First, produce two figures (LDA and VarImp) in the Plots window and automatically save them as PDF files
# # Disable the "Press <ENTER> to Continue..." prompt
# # dataTrain.SubScen <- dataTrain[dataTrain$scenario %in% c(1,2), ]
# dataTrain.SubScen <- dataTrain[dataTrain$scenario %in% group.list, ]
# # Renumber the row indices
# rownames(dataTrain.SubScen) <- NULL
# # Verify the result
# # head(dataTrain.SubScen)
# # First, produce two figures (LDA and VarImp) in the Plots window and automatically save them as PDF files
# # Disable the "Press <ENTER> to Continue..." prompt
# # par(ask = FALSE)
# plot(model.rf, dataTrain.SubScen, obs = stat.obs, pdf = TRUE, n.var = 30)
# 
# ########## If "margin problem", the graphs are still produced; proceed with
# # Save the error = f(n.tree) plot in a PNG file
# png("err_oob_plot.png", width = 800, height = 600)
# # Compute err.oob
# err.oob <- err.abcrf(model.rf, dataTrain.SubScen, paral = TRUE)
# # Produce the plot
# plot(err.oob)
# # Close the PNG file to save the plot
# dev.off()  # Important to close the device so it returns to the Plots window
# cat("THREE ILLUSTRATIVE GRAPHICS have been produced and saved in three different files:", quote=FALSE)
# cat("File named prior_errors_vs_number_of_trees.png: Graphic providing prior error rates for forests with different number of trees and computed using an Out-of-Bag procedure, e.g. Fig. 3 in Pudlo et al. 2016", "\n")
# cat("File named graph_lda.pdf = LDA projections of the reference table for the different scenarios plus the observed dataset cf. black star in the figure", "\n")
# cat("File name graph_varImpPlot.pdf = the contributions of the 30 most important statistics to the RF (e.g. Fig. S6 and Fig. S7 in Pudlo et al. 2016)", "\n")
# # #####################################################################################
# # ############################### START RENAME OUTPUT FILES ######################################################
# # Give a specific name or use the one you defined for your analysis
# #output.file.save <- "abcrf_LUIS_USA__GAUT_vs_FRAI_S1&S2_tok_sap_ok_t1000_s20000_pS_10runs.txt"
# output.file.save <- output.file
# 
# # Extract the name without extension and the extension of output.file
# nse_output.file.save <- tools::file_path_sans_ext(output.file)
# ext_output.file.save <- tools::file_ext(output.file)
# 
# # List of output files to rename
# files_to_rename <- c("graph_lda.pdf", 
#                      "graph_varImpPlot.pdf", 
#                      "err_oob_plot.png", 
#                      "Importance_Measures_of_SUMSTATS_SORTED.txt")
# # Rename each output file
# for (file in files_to_rename) {
#   # Extract the name without extension and the extension of the output file
#   nse4f <- tools::file_path_sans_ext(file)
#   ext4f <- tools::file_ext(file)
#   # Create the new file name
#   new_name <- paste0(nse4f, "_", nse_output.file.save, ".", ext4f)
#   # Rename the file
#   file.rename(file, new_name)
# }
# # ############################### END RENAME OUTPUT FILES ######################################################
