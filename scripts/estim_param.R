# Estimation RF parametres avec ou sans PLS 31/05/2018 et 03/2020
# inclu densityPlot {abcrf}	= Plot the posterior density given a new set of summary statistic
# last modifications AE 30/12/2024

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

################# Loading useful libraries #########
library(abcrf) # RF
library(MASS) # LDA
library(ranger) # CSRF
library(rpart) # CART
library(doParallel)
library(pls)
help(package = "abcrf")

#### Various options
options(max.print=10000) # To print tables until 10000 lines
# Number of cores available for your computer
n_cores <- detectCores()
# Nbre of CPU cores I want to use for parallel computation (if not define then all available CPU cores will be used)
how.many.cores.used.for.computation = n_cores-3

# Saving results
save <- TRUE  # Set save to TRUE or FALSE
if (save == TRUE) {
output_file_name = "estimated_params.txt"
# output_file_name = "test.txt"
}

# Key file names
name.reference.table <- "reftableRF.bin"
name.header.file <- "headerRF.txt"
#name.observed.dataset <- "statobsRF_all_statobs_vectors.txt" # Note that this observed dataset may includes several lines (i.e. several vectors of observed data/sumstats) 
#name.observed.dataset <- "statobsRF_target_dataset.mss.dataNANNE.dat_US-Wat_US-Sok_US-Haw_JP-Sap_CN-Lan_BR-Poa.txt"
name.observed.dataset <- "statobsRF.txt"

### RF parameter
selected.model=1
N.train <- 4000 # Je vais entrainer des forêts sur N.train données
Ntest <- 800  # Tests datasets to compute MNAE, MSE, etc sera fait sur Ntest data of the reference.table not used for training
# WARNING: PLS: done on ncomp_total = 0.1 * nstats puis selection of n comp PLS so that 0.99 variance totals des ncomp_total
N.load.reftable = 4000
ntree=1000 #### nbre of trees tree to do the RF

reference.table <- readRefTable(filename = name.reference.table, header=name.header.file, N=N.load.reftable)

# Computation options
param_original = TRUE
param_compound = FALSE
PLS=FALSE
Graphical_density_representation=FALSE

# Parameters to estimate
### Original parameters
list_original_parameters = "NONE"
if (param_original==TRUE) {
  # Liste (globale) obtenue automatiquement a partir du reftable
  
  # Liste des colonnes initiales
  colnames(reference.table$params)
  print(reference.table$params)
  reference.table.dim <- dim(reference.table$params)
  print(reference.table.dim)
  # Créer la liste à partir des noms de colonnes de reference.table$params
  list_original_parameters <- as.list(colnames(reference.table$params))
  # Combiner le contenu de la liste en un vecteur
  list_original_parameters <- unlist(list_original_parameters, use.names = FALSE)
  list_original_parameters <- c(list_original_parameters)
  # Sélection des éléments 1 à 63
  list_original_parameters <- list_original_parameters[1:63]
  #print(list_original_parameters)
  
  # Liste addoc
  
  #list_original_parameters = c("raan1","raan2", "raawat", "raasd", "raaas1", "raaas2")
  list_original_parameters <- c("t34","t23","t12","t32","t24","t10","t31","t21","t14","N1","N2","N3","N4","N24","N23","N12","N32","N24","N10","N31","N21","N14")
}

# Loading key files
# N = total number of simuations one wants to load from the reference table
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

### Préambule et partage des données = focussing on the selected model

reference.table$scenarios <- reference.table$scenarios[indexesModel]
reference.table$stats <- reference.table$stats[indexesModel,]
reference.table$params <- reference.table$params[indexesModel,]
# dim(reference.table$params)
# head(reference.table$params, n=3)
# colnames(reference.table$params)

# Supprimer les colonnes contenant uniquement des NA
reference.table$params <- reference.table$params[, colSums(is.na(reference.table$params)) < nrow(reference.table$params)]
# dim(reference.table$params)
# head(reference.table$params, n=3)
# colnames(reference.table$params)

# Combien on a de données au final
N <- length(reference.table$scenarios) 
print(N)
# Combien de stats
N.stats <- ncol(reference.table$stats)
# Combien de param
N.param <- ncol(reference.table$params)
# Nombre de données potentielles pour faire des tests
N.test_possible <- N - N.train # datatests = potentiellement toutes les donnees n'ayant pas servis pour entrainement

# set.seed(1) # If one want some controle
indicesTrain <- sample(1:N, N.train, replace=FALSE) # indices d'entrainement
indicesTest <- c(1:N)[-indicesTrain] # indices de test

################################################################################
### Compound parameters  #######################################################
################################################################################
# list_compound_parameters <- "None"  # Compound parameteres ra
# if (param_compound == TRUE) {
#   # param1 <- reference.table$params[, "raaas1"]
#   # param2 <- reference.table$params[, "raaas2"]
#   
#   param1 <- reference.table$params[, "DBwat"]
#   param2 <- reference.table$params[, "NBwat"]
#   
#   # Modifier les valeurs de param1 égales à 0.000 en 1/2000
#   param1[param1 == 0.000] <- 1/2000
#   
#   # # Voir si Valeurs nulles ? cf. si le cas alors pose pbl pour calcul 
#   # # de post.NMAE.median post.NMAE.mean prior.NMAE.median prior.NMAE.mean
#   # # Trouver les indices où la valeur est 0.000
#   # indices_zero <- which(param1 == 0.000)
#   # # Compter le nombre de lignes avec la valeur 0.000
#   # nombre_zero <- length(indices_zero)
#   # # Afficher les indices et le nombre
#   # cat("Nombre de lignes avec une valeur de 0.000 :", nombre_zero, "\n")
#   # cat("Indices des lignes correspondantes :", indices_zero, "\n")
#   # # Afficher les lignes correspondantes (si applicable)
#   # if (nombre_zero > 0) {
#   #   print(reference.table$params[indices_zero, ])
#   # } else {
#   #   cat("Aucune ligne avec une valeur de 0.000 trouvée.\n")
#   # }
#   
#   # Définir les objets composés
#   n.compound.param = 4
#   ra.Wat.AS <- 1 - param1
#   ra.Gan2.SA <- param2 * param1
#   ra.nc.AS <- (1 - param2) * param1
#   ra.watGan2.AS <- ra.Wat.AS + ra.Gan2.SA
#   
#   param_compound_dataframe = t(rbind(ra.Wat.AS, ra.Gan2.SA, ra.nc.AS, ra.watGan2.AS))
#   param_compound_dataframe = as.data.frame(param_compound_dataframe)
#   # dim(param_compound_dataframe)
#   # head(param_compound_dataframe, n=3)
#   # colnames(param_compound_dataframe)
#   # str(param_compound_dataframe)
#   # List of compound parameters
#   list_compound_parameters <- colnames(param_compound_dataframe)
#   # Add compound parameters to reference.table$params
#   reference.table$params <- cbind(reference.table$params, param_compound_dataframe)
#   # dim(reference.table$params)
#   # head(reference.table$params, n=3)
#   # colnames(reference.table$params)
#   # str(reference.table$params)
# }

########## Compound parameters for BI = Bottleneck Intensity #################
list_compound_parameters <- "None"
# Activer la génération des variables composites seulement si param_compound est TRUE
if (param_compound == TRUE) {
  #Définir la liste des paramètres
  list_original_parameters <- c("t34","t23","t12","t32","t24","t10","t31","t21","t14","N1","N2","N3","N4","N24","N23","N12","N32","N24","N10","N31","N21","N14")


  # Nombre total de binômes (chaque binôme correspond à deux paramètres consécutifs)
  num_BI_binomes <- length(list_original_parameters) / 2
  # Initialiser un dataframe vide pour stocker les résultats
  BI_binomes_dataframe <- data.frame(matrix(nrow = nrow(reference.table$params), ncol = num_BI_binomes))
  colnames(BI_binomes_dataframe) <- rep("", num_BI_binomes)  # Initialiser les noms de colonnes

  # Générer dynamiquement les binômes
col_index <- 1
  for (i in seq(1, length(list_original_parameters), by = 2)) {
    db_param <- list_original_parameters[i]
    nb_param <- list_original_parameters[i + 1]
   bi_name <- paste0("BI", sub("DB", "", db_param)) # For computing Bottleneck Intensity
   
    # Vérifier si les colonnes existent
    if (db_param %in% colnames(reference.table$params) && nb_param %in% colnames(reference.table$params)) {
      # Calculer le binôme en accédant aux colonnes de la matrice
      binome_col <- reference.table$params[, db_param] / reference.table$params[, nb_param]
      #binome_col <- reference.table$params[, db_param] * reference.table$params[, nb_param]

      # Ajouter le binôme au dataframe
      BI_binomes_dataframe[, col_index] <- binome_col
      colnames(BI_binomes_dataframe)[col_index] <- bi_name  # Nommer la colonne
      col_index <- col_index + 1
    } else {
      # Signaler les colonnes manquantes
      missing_cols <- setdiff(c(db_param, nb_param), colnames(reference.table$params))
      cat("Colonnes manquantes :", paste(missing_cols, collapse = ", "), "\n")
    }
  }

  # Modifier les valeurs de BI égales à 0.000 (si DB = 0) en 1/2000
  # En effet si BI=0 alors pose pbl pour calcul 
  # de post.NMAE.median post.NMAE.mean prior.NMAE.median prior.NMAE.mean
  BI_binomes_dataframe[BI_binomes_dataframe == 0.000] <- 1/2000 # car DB=1 au lieu de 0 et NBmax = 2000

  # Supprimer les colonnes vides dans BI_binomes_dataframe
  BI_binomes_dataframe <- BI_binomes_dataframe[, colnames(BI_binomes_dataframe) != ""]
  # Ajouter les nouvelles colonnes au dataframe original
  reference.table$params <- cbind(reference.table$params, BI_binomes_dataframe)
  list_compound_parameters <- colnames(BI_binomes_dataframe)
  param.list <- colnames(reference.table$params)

  # Afficher les dimensions et les premières lignes pour vérification
  # dim(BI_binomes_dataframe)
  # head(BI_binomes_dataframe)
  # colnames(BI_binomes_dataframe)
  # dim(reference.table$params)
  # head(reference.table$params, n=3)
  # colnames(reference.table$params)
}

################################################################################
################################################################################
sink(file = output_file_name, split = TRUE)

######################################################################################################################
##### TREATMENT FOR A GIVEN PARAMETER   ##############################################################################
######################################################################################################################

# Combiner les deux listes original.param et compound.param en une seule si necessaire
if ((param_original == TRUE)&(param_compound == TRUE)) combined_list <- c(list_original_parameters, list_compound_parameters)
if ((param_original == TRUE)&(param_compound == FALSE)) combined_list <- c(list_original_parameters)
if ((param_original == FALSE)&(param_compound == TRUE)) combined_list <- c(list_compound_parameters)

cat(" ############## PARAMETER ESTIMATION using random forest ############### ","\n")
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
cat("Number of simulations in the TRAINING DATASET used to built rf trees =", N.train, "\n")
cat("Number of trees in the forest =",ntree, "\n")
cat("Number of test datasets that will be used to compute NMAE=",Ntest,"\n")
cat("Nbre of test dataset that can be possibly used to compute NMAE=",N.test_possible,"\n")
cat("Number of cores available = ",n_cores, "\n")
cat("Number of cores used for computation =",how.many.cores.used.for.computation, "\n")
cat("Selected model=",selected.model,"\n")
cat("N train=",N.train,"\n")
cat("N tree: ",ntree,"\n")
cat("N test_calcul_NMAE=",Ntest,"\n")
if (param_original == TRUE) cat("Original parameters estimated = ",list_original_parameters, "\n")
if (param_compound == TRUE) cat("Compound parameters estimated = ",list_compound_parameters, "\n")
cat("Parameters in fine estimated = ",combined_list, "\n")
cat("\n")


# Buiding result table for all parameters
#################################################################
# Dimensions du dataframe
n_rows <- length(combined_list)
col_names <- c("Parameter", "Mean", "Median", "Q5", "Q95", 
               "post.NMAE.median", "post.NMAE.mean", 
               "prior.NMAE.median", "prior.NMAE.mean", "Coverage.OOB")

# Initialisation du dataframe avec NA
RESULT_ESTIM_PARAM <- data.frame(matrix(NA, nrow = n_rows, ncol = length(col_names)))
colnames(RESULT_ESTIM_PARAM) <- col_names
# Remplir la colonne "Parameter" avec les noms des paramètres
RESULT_ESTIM_PARAM$Parameter <- combined_list
##################################################

# PROCESSING ESTIMATION OF PARAMETERS = Parcourir tous les éléments de la liste combinée
param.nbr=0
for (parameter_of_interest in combined_list) {
  param.nbr=param.nbr+1
  cat("\n")  
  cat("\n") 
  cat("######################################################################################################################","\n")
  cat("################################### PROCESSING ESTIMATION OF PARAMETER :", parameter_of_interest,"########################################","\n")
  cat("######################################################################################################################","\n")
  cat("\n") 
  cat("RF analysis - parameter ",parameter_of_interest, "\n")
  cat("Parameter #",param.nbr," over a total of",length(combined_list), "parameters to estimate", "\n")
  y <- reference.table$params[indicesTrain, parameter_of_interest]
  ytest <- reference.table$params[indicesTest, parameter_of_interest]
  cat("LENGTH OF YTEST")
  cat(length(ytest))
  # # Cas particulier avec les param compound ra
  # if (param_original== TRUE) {
  #   if (parameter_of_interest %in% list_original_parameters) {
  #   y <- reference.table$params[indicesTrain, parameter_of_interest]
  #   ytest <- reference.table$params[indicesTest, parameter_of_interest]
  #   }
  # }
  # if (param_compound == TRUE) {
  #   if (parameter_of_interest %in% list_compound_parameters) {
  #   compound_parameter_of_interest <- get(parameter_of_interest, envir = .GlobalEnv)
  #   y <- compound_parameter_of_interest[indicesTrain]
  #   ytest <- compound_parameter_of_interest[indicesTest]
  #   }
  # }
 
  x <- as.data.frame(reference.table$stats[indicesTrain,])
  xtest <- as.data.frame(reference.table$stats[indicesTest,])
  train <- data.frame(r = y, x) # On crée le data.frame pour RF sans PLS
  # On regarde un sous ensemble Ntest defini en amont de données test
  indiceTest <- sample(1:length(ytest), Ntest, replace=FALSE) # Je tire au hasard des données de test
  xtest <- xtest[indiceTest,] # Je réduis mes données de test
  ytest <- ytest[indiceTest]  # pour la réponse ou les covariables
  colnames(xtest) <- colnames(train)[-1]


if (PLS==FALSE) {
###########################   TREATMENTS SANS PLS mais avec NOISE #######################
# AJOUTER 5 variable de NOISE d
monNOISE <- matrix(runif(N.train*5), ncol=5)
colnames(monNOISE) <- c("NOISE1", "NOISE2", "NOISE3", "NOISE4", "NOISE5")
x.NOISE <- cbind(x,monNOISE)
train.NOISE <- data.frame(r = y, x.NOISE) # On cree le data.frame pour les forets avec le NOISE
colnames(x.NOISE) <- colnames(train.NOISE)[-1] # On met les memes noms partout

# RF
rf.abcrf.NOISE <- regAbcrf(r ~ ., data = train.NOISE, ntree = ntree, min.node.size = 5, paral = TRUE, ncores = how.many.cores.used.for.computation)
# Pour visualiser l'effet du nbre de trees sur erreur d'estimation de la moyenne
# Ouvrir un fichier PNG pour l'enregistrement
errorOOB.NOISE <- err.regAbcrf(object = rf.abcrf.NOISE, training = train.NOISE, paral = TRUE)
# Pour visualiser l'effet du nbre de trees sur erreur d'estimation de la moyenne
plot(errorOOB.NOISE)
png(file = paste0(output_file_name,"_error_vs_ntrees_",parameter_of_interest,".png"))
# Créer la figure (assurez-vous d'avoir un appel à la fonction plot() ou similaire)
plot(errorOOB.NOISE)
# Fermer l'appareil graphique pour sauvegarder l'image
dev.off()
# Pour visualiser les variables importances via figure
plot(rf.abcrf.NOISE)
# Pour avoir un graphique ou on peut jouer sur la taille de la police
plotVarImp = variableImpPlot(rf.abcrf.NOISE, n.var=30)
# Ouvrir un fichier PNG pour l'enregistrement
png(file = paste0(output_file_name, "_VarImp.png_",parameter_of_interest,".png"))
# Créer la figure avec variableImpPlot
variableImpPlot(rf.abcrf.NOISE, n.var = 30)
# Fermer l'appareil graphique pour sauvegarder l'image
dev.off()

cat("\n")
cat("ERROR METRICS","\n")
rf.abcrf.NOISE
print(rf.abcrf.NOISE)

####### Pour extraire TOUTES les VarImp du graphique VarImp(stat) et ecrire dans fichier txt ##############
# Pour visualiser toutes les variables importances en txt
#rf.abcrf$model.rf$variable.importance # Pour visualiser les variables importances
#sort(rf.abcrf$model.rf$variable.importance) # Ordre croissant
#sort(rf.abcrf$model.rf$variable.importance, decreasing = TRUE) # Ordre d?croissant
output <- rf.abcrf.NOISE$model.rf$variable.importance
output = data.frame(output)
write.table(output, paste0(output_file_name,"_var_imp_estim_param_NOT_sorted_sans_PLS.txt"), quote=FALSE, row.names=TRUE, col.names="stat_name var_imp")
table_trie = read.table(paste0(output_file_name,"_var_imp_estim_param_NOT_sorted_sans_PLS.txt"),header=TRUE)
table_trie$var_imp = round(table_trie$var_imp, 6)
output_after_trie = table_trie[order(table_trie$var_imp, decreasing=TRUE),]
write.table(output_after_trie, paste0(output_file_name,"_var_imp_estim_param_SORTED_sans_PLS.txt"), quote=FALSE, row.names=FALSE)

### PREDICTIONS on a given DATA_OBS
obs.poi <- as.data.frame(stat.obs)
colnames(obs.poi) <- colnames(xtest)
# On ajoute du NOISE (5 variables)
monNOISEObs <- matrix(runif(5), ncol=5)
colnames(monNOISEObs) <- c("NOISE1", "NOISE2", "NOISE3", "NOISE4", "NOISE5")
obs.poi.NOISE <- cbind(obs.poi, monNOISEObs) # On ajoute les colonnes de scores PLS et NOISEs
colnames(obs.poi.NOISE) <- colnames(train.NOISE)[-1] # Homogenisation noms des colonnes
statobs_predClassic.NOISE <- predict(rf.abcrf.NOISE, obs=obs.poi.NOISE, quantiles = c(0.05,0.95), paral = TRUE, 
                                     paral.predict=TRUE, prior.err.med = TRUE, post.err.med = TRUE, 
                                     ncores=how.many.cores.used.for.computation, min.node.size = 5, training=train.NOISE)

# AFFICHAGE RESULTS + NOISE (NO PLS)
# Résultats avec bruit (sans PLS)
cat("\n")
cat("RESULTS + NOISE (NO PLS)\n")

### Dynamic Building of result Table with all RF results for each parameter
# Extraire les informations de pref.rf for storage
output <- capture.output(statobs_predClassic.NOISE) # Capture les sorties affichées en console
# Identifier et extraire les lignes contenant des informations particulières
NMAE_prior_mean <- as.numeric(sub(".*mean: ", "", grep("Prior out-of-bag normalized mean absolute error computed with mean", output, value = TRUE)))
NMAE_prior_median <- as.numeric(sub(".*median: ", "", grep("Prior out-of-bag normalized mean absolute error computed with median", output, value = TRUE)))
prior_OOB_coverage <- as.numeric(sub(".*coverage: ", "", grep("Prior out-of-bag credible interval coverage", output, value = TRUE)))
# Verify
# cat("Prior out-of-bag normalized mean absolute error (mean):", NMAE_prior_mean, "\n")
# cat("Prior out-of-bag normalized mean absolute error (median):", NMAE_prior_median, "\n")
# cat("Prior out-of-bag credible interval coverage:", prior_OOB_coverage, "\n")

statobs_predClassic.NOISE.key.results <- c(statobs_predClassic.NOISE$expectation,
                         statobs_predClassic.NOISE$med,
                         statobs_predClassic.NOISE$quantiles[1],
                         statobs_predClassic.NOISE$quantiles[2],
                         statobs_predClassic.NOISE$post.NMAE.med,
                         statobs_predClassic.NOISE$post.NMAE.mean,
                         NMAE_prior_median,
                         NMAE_prior_mean, 
                         prior_OOB_coverage)
### Dynamic Building of result Table with the RF key results for each parameter
RESULT_ESTIM_PARAM[param.nbr,2:ncol(RESULT_ESTIM_PARAM)] = statobs_predClassic.NOISE.key.results
print(RESULT_ESTIM_PARAM[param.nbr,])

if (Graphical_density_representation==TRUE) {
  # GRAPHICAL REPRESENTATION OF THE APPROXIMATE POSTERIOR DENSITY
  # Enregistrer le résultat de density plot dans un fichier image au format JPEG avec le nom dérivé de output_file_name
  #densityPlot(object = rf.abcrf.NOISE, obs = obs.poi.NOISE, training = train.NOISE, paral = TRUE, bw = 1000, ylim = c(0, 0.0002), xlim = c(10,16000), lty = 3, main = "Posterior distribution (noPLS)", xlab = "Parameter N1",ylab= "Density")
  # Définir le nom du fichier de sortie au format PNG
  output_file_png <- gsub(".txt", ".png", output_file_name)
  output_file_png = paste0(parameter_of_interest,"_", output_file_png)
  
  # Enregistrer la figure dans un fichier PNG
  png(filename = paste0("Parameter_density_", output_file_png))
  # Générer le density plot
  densityPlot(object = rf.abcrf.NOISE, obs = obs.poi.NOISE, training = train.NOISE, 
              paral = TRUE, bw = 0.05, ylim = c(0, 6), xlim = c(0, 1), 
              lty = 3, main = "Posterior distribution (noPLS)", 
              xlab = "Admixture rate North Carolina", ylab = "Density")
  # Fermer l'appareil de sortie graphique
  dev.off()
}
}

if (PLS==TRUE) {
  ############################## SECTION POUR TRAITEMENT AVEC PLS ################
  ############# Calculs, selection et ajouts de composantes PLS = AE #############
  ################################################################################
  
  ### On laisse tomber le "coude" de la variance expliquee pour la reponse
  ### on laisse tomber les histoire de validation croisee
  ## On va tout simplement selectionner x% de la variance expliquee par les axes PLS
  
  ncomp_total=N.stats # nbre total de composantes PLS = nbre de stats - WARNING: il faut plus de simues que de stats !!!
  ncomp_total # Pour se rappeler !
  ncomp_total = 1.0 * ncomp_total ############ TESTS POUR EVITER DE TROP GROS CALCULS ou alors ncomp_total = 0.1 * ncomp_total
  ncomp_total
  # Analyse PLS sensus stricto
  pls.fit = plsr(r ~ ., data=train, scale=TRUE, ncomp=ncomp_total)
  summary(pls.fit)
  # Pour recuperer le pourcentage de variance de reponse expliquee
  percentYVar <- 100 * drop(R2(pls.fit, estimate = "train",intercept = FALSE)$val)
  plot(1:length(percentYVar), percentYVar, type="b", cex=1.5, cex.axis=1.5, cex.lab=1.5,  xlab="Number of components", ylab="Percent of Variance explained by PLS axes")
  # A veut maintenant selectioner un sous-nbre d'axes PLS (nComposante_sel) afin d'obtenir une proportion p_threshold_PLS = 0.99
  # par exemple de la variance expliquee par PLS: objectif = virer les axes qui amenent rien en terme info
  p_threshold_PLS = 0.95
  p_threshold_PLS
  percentYVar_max = percentYVar[ncomp_total]
  percentYVar_max
  p_var_PLS = p_threshold_PLS * percentYVar_max
  p_var_PLS
  nComposante_sel <- max(which(percentYVar <= p_var_PLS))
  nComposante_sel
  # nComposante_sel = ncomp_total ########  &&&&&&&&&&&&&&&&&&& juste pour voir si tous les axes !!!!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  abline(h=p_var_PLS, col="green") # Pour voir visuellement la ou on va couper avec les nComposantes selectionnees
  poids.PLS = loading.weights(pls.fit)[,1] ### donne les poids de chaque résumé pour la première composante
  poids.PLS
  
  ################### Ajouter dans le training set les scores PLS comme nouvelles covariables AINSI que 5 variables de bruit ############
  # je recupere les scores predits precedemment
  scoresPLS_sel <- scores(pls.fit)[,1:nComposante_sel]
  # Je rajoute les scores et le NOISE (monNOISE = 5 variables UNIF:0-1) comme nouvelles colonne
  monNOISE <- matrix(runif(N.train*5), ncol=5)
  colnames(monNOISE) <- c("NOISE1", "NOISE2", "NOISE3", "NOISE4", "NOISE5")
  x.pls.NOISE <- cbind(x, scoresPLS_sel, monNOISE)
  train.pls.NOISE <- data.frame(r = y, x.pls.NOISE) # On cree le data.frame pour les forets avec PLS et NOISE
  colnames(x.pls.NOISE) <- colnames(train.pls.NOISE)[-1] # On met les memes noms partout
  
  
  ###########################   TREATMENTS AVEC PLS et NOISE #######################
  # RF
  rf.abcrf.pls.NOISE <- regAbcrf(r ~ ., data = train.pls.NOISE, ntree = ntree, min.node.size = 5, paral = TRUE, ncores = how.many.cores.used.for.computation)
  # Pour visualiser l'effet du nbre de trees sur erreur d'estimation de la moyenne
  png(file="param_prior_errors_vs_number_of_trees_PLS_NOISE.png")
  errorOOB_pls_NOISE <- err.regAbcrf(object = rf.abcrf.pls.NOISE, training = train.pls.NOISE, paral = TRUE)
  dev.off()
  # Pour visualiser les variables importances via figure
  plot(rf.abcrf.pls.NOISE, n.var=30)
  # Pour avoir un graphique ou on peut jouer sur la taille de la police
  # variableImpPlot(rf.abcrf.pls.NOISE, cex=0.5, cex.axis=2.0, cex.lab=2.0, n.var=50)
  variableImpPlot(rf.abcrf.pls.NOISE, n.var=50)
  
  # computed using an Out-of-Bag procedure NMAE ?????????? FAUX DEMANDER MPC
  #png(file="param_prior_errors_vs_number_of_trees.png")
  #err.oob <- err.abcrf(rf.abcrf.pls.NOISE, train.pls.NOISE, paral=TRUE)
  #dev.off()
  #plot(err.oob)
  
  ####### Pour extraire TOUTES les VarImp du graphique VarImp(stat) et ecrire dans fichier txt ###### AE #######
  # Pour visualiser toutes les variables importances en txt
  #rf.abcrf.pls.NOISE$model.rf$variable.importance # Pour visualiser les variables importances en txt
  #sort(rf.abcrf.pls.NOISE$model.rf$variable.importance) # Ordre croissant  # Ordre croissant = affichage text
  #sort(rf.abcrf.pls.NOISE$model.rf$variable.importance, decreasing = TRUE) # Ordre decroissant = affichage text
  output <- rf.abcrf.pls.NOISE$model.rf$variable.importance
  output = data.frame(output)
  write.table(output, "var_imp_estim_param_NOT_sorted_AVEC_PLS.txt", quote=FALSE, row.names=TRUE, col.names="stat_name var_imp")
  table_trie = read.table("var_imp_estim_param_NOT_sorted_AVEC_PLS.txt",header=TRUE)
  table_trie$var_imp = round(table_trie$var_imp, 6)
  output_after_trie = table_trie[order(table_trie$var_imp, decreasing=TRUE),]
  write.table(output_after_trie, "var_imp_estim_param_SORTED_AVEC_PLS.txt", quote=FALSE, row.names=FALSE)
  
  #### PREDICTIONS on a given DATA_OBS
  # If obs.poi is not a dataframe or the column names do not match,
  # you can use the following lines:
  obs.poi <- as.data.frame(obs.poi)
  colnames(obs.poi) <- colnames(xtest)
  # On récupère les scores PLS pour les données statobs
  scores_statobs.pls <- predict(pls.fit, obs.poi, type = "scores")
  scores_statobs.pls.sel <- scores_statobs.pls[,1:nComposante_sel, drop = FALSE]
  # cf. evite inverser row et col cf. R fait chier quand un seul dataset !!!!
  # scores_statobs.pls.sel <- scores_statobs.pls[,1:nComposante_sel]
  # scores_statobs.pls.sel <- as.data.frame(t(scores_statobs.pls.sel)) ###WARNING: on doit inverser row et col !!!!
  # On ajoute du NOISE (5 variables)
  monNOISEObs <- matrix(runif(5), ncol=5)
  colnames(monNOISEObs) <- c("NOISE1", "NOISE2", "NOISE3", "NOISE4", "NOISE5")
  obs.poi.pls.NOISE <- cbind(obs.poi, scores_statobs.pls.sel, monNOISEObs) # On ajoute les colonnes de scores PLS et NOISEs
  colnames(obs.poi.pls.NOISE) <- colnames(train.pls.NOISE)[-1] # Homogenisation noms des colonnes
  
  statobs_predClassic.pls.NOISE <- predict(rf.abcrf.pls.NOISE, obs=obs.poi.pls.NOISE, quantiles = c(0.5,0.05,0.95), paral = TRUE,
                                           paral.predict=TRUE, post.err.med = TRUE, ncores=how.many.cores.used.for.computation, min.node.size = 5, training=train.pls.NOISE)
  
  
  df.pls.noise<-c(statobs_predClassic.pls.NOISE$expectation,
                  statobs_predClassic.pls.NOISE$med,
                  statobs_predClassic.pls.NOISE$quantiles[2],
                  statobs_predClassic.pls.NOISE$quantiles[3],
                  statobs_predClassic.pls.NOISE$post.NMAE.med,
                  statobs_predClassic.pls.NOISE$post.NMAE.mean,
                  statobs_predClassic.pls.NOISE$prior.coverage,
                  statobs_predClassic.pls.NOISE$prior.NMAE.med,
                  statobs_predClassic.pls.NOISE$prior.NMAE.mean,
                  statobs_predClassic.pls.NOISE$model.rf$NMAE)
  
  statobs_predClassic.pls.NOISE$med
  statobs_predClassic.pls.NOISE$quantiles[2]
  statobs_predClassic.pls.NOISE$quantiles[3]
  statobs_predClassic.pls.NOISE$variance
  statobs_predClassic.pls.NOISE$post.NMAE.med
  statobs_predClassic.pls.NOISE$post.NMAE.mean
  statobs_predClassic.pls.NOISE$prior.coverage
  statobs_predClassic.pls.NOISE$prior.NMAE.med
  statobs_predClassic.pls.NOISE$prior.NMAE.mean
  statobs_predClassic.pls.NOISE$model.rf$NMAE
  
  # RESULTS + PLS + NOISE
  var.names<-c("mean","median","q5","q95","post.NMAE.median","post.NMAE.mean","prior.NMAE.median","prior.NMAE.mean")
  #  statobs_predClassic.NOISE$prior.coverage NE MARCHE PAS !!!
  var.names
  df.pls.noise
  
  # A GRAPHICAL REPRESENTATION OF THE APPROXIMATE POSTERIOR DENSITY
  #densityPlot(object = rf.abcrf.NOISE, obs = obs.poi.NOISE, training = train.NOISE, paral = TRUE, bw = 0.05, ylim = c(0, 3), xlim = c(0, 1), lty = 3, main = "Posterior distribution without PLS", xlab = "Parameter r_wat = rGHAS1*rGHAS2",ylab= "Density")
  densityPlot(object = rf.abcrf.pls.NOISE, obs = obs.poi.pls.NOISE, training = train.pls.NOISE, paral = TRUE, bw = 0.05, ylim = c(0, 4), xlim = c(0, 1), lty = 3, main = "Posterior distribution WITH PLS", xlab = "Parameter r_wat = rGHAS1*rGHAS2",ylab= "Density")
  # Pour VARIABLE temps ??? --> marche pas !!!
  # densityPlot(object = rf.abcrf.pls.NOISE, obs = obs.poi.pls.NOISE, training = train.pls.NOISE, paral = TRUE, bw = 100, ylim = c(0, 6),from = 100, to = 600, lty = 3, main = "Posterior distribution WITH PLS", xlab = parameter_of_interest)
  
  
  # computed using an Out-of-Bag procedure NMAE ?????????? FAUX DEMANDER MPC
  png(file="param_prior_errors_vs_number_of_trees.png")
  err.oob <- err.abcrf(rf.abcrf.pls.NOISE, train.pls.NOISE, paral=TRUE)
  dev.off()
  plot(err.oob)
  
  ################################################################################
  ############# Calculs, selection et ajouts de composantes PLS = AE #############
  ################################################################################
  
  ### On laisse tomber le "coude" de la variance expliquee pour la reponse
  ### on laisse tomber les histoire de validation croisee
  ## On va tout simplement selectionner x% de la variance expliquee par les axes PLS
  
  ncomp_total=N.stats # nbre total de composantes PLS = nbre de stats - WARNING: il faut plus de simues que de stats !!!
  ncomp_total # Pour se rappeler !
  ncomp_total = 1.0 * ncomp_total ############ TESTS POUR EVITER DE TROP GROS CALCULS ou alors ncomp_total = 0.1 * ncomp_total
  ncomp_total
  # Analyse PLS sensus stricto
  pls.fit = plsr(r ~ ., data=train, scale=TRUE, ncomp=ncomp_total)
  summary(pls.fit)
  # Pour recuperer le pourcentage de variance de reponse expliquee
  percentYVar <- 100 * drop(R2(pls.fit, estimate = "train",intercept = FALSE)$val)
  plot(1:length(percentYVar), percentYVar, type="b", cex=1.5, cex.axis=1.5, cex.lab=1.5,  xlab="Number of components", ylab="Percent of Variance explained by PLS axes")
  # A veut maintenant selectioner un sous-nbre d'axes PLS (nComposante_sel) afin d'obtenir une proportion p_threshold_PLS = 0.99
  # par exemple de la variance expliquee par PLS: objectif = virer les axes qui amenent rien en terme info
  p_threshold_PLS = 0.95
  p_threshold_PLS
  percentYVar_max = percentYVar[ncomp_total]
  percentYVar_max
  p_var_PLS = p_threshold_PLS * percentYVar_max
  p_var_PLS
  nComposante_sel <- max(which(percentYVar <= p_var_PLS))
  nComposante_sel
  # nComposante_sel = ncomp_total ########  &&&&&&&&&&&&&&&&&&& juste pour voir si tous les axes !!!!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  abline(h=p_var_PLS, col="green") # Pour voir visuellement la ou on va couper avec les nComposantes selectionnees
  poids.PLS = loading.weights(pls.fit)[,1] ### donne les poids de chaque résumé pour la première composante
  poids.PLS
  
  ################### Ajouter dans le training set les scores PLS comme nouvelles covariables AINSI que 5 variables de bruit ############
  # je recupere les scores predits precedemment
  scoresPLS_sel <- scores(pls.fit)[,1:nComposante_sel]
  # Je rajoute les scores et le NOISE (monNOISE = 5 variables UNIF:0-1) comme nouvelles colonne
  monNOISE <- matrix(runif(N.train*5), ncol=5)
  colnames(monNOISE) <- c("NOISE1", "NOISE2", "NOISE3", "NOISE4", "NOISE5")
  x.pls.NOISE <- cbind(x, scoresPLS_sel, monNOISE)
  train.pls.NOISE <- data.frame(r = y, x.pls.NOISE) # On cree le data.frame pour les forets avec PLS et NOISE
  colnames(x.pls.NOISE) <- colnames(train.pls.NOISE)[-1] # On met les memes noms partout
  
  
  ###########################   TREATMENTS AVEC PLS et NOISE #######################
  # RF
  rf.abcrf.pls.NOISE <- regAbcrf(r ~ ., data = train.pls.NOISE, ntree = ntree, min.node.size = 5, paral = TRUE, ncores = how.many.cores.used.for.computation)
  # Pour visualiser l'effet du nbre de trees sur erreur d'estimation de la moyenne
  png(file="param_prior_errors_vs_number_of_trees_PLS_NOISE.png")
  errorOOB_pls_NOISE <- err.regAbcrf(object = rf.abcrf.pls.NOISE, training = train.pls.NOISE, paral = TRUE)
  dev.off()
  # Pour visualiser les variables importances via figure
  plot(rf.abcrf.pls.NOISE, n.var=30)
  # Pour avoir un graphique ou on peut jouer sur la taille de la police
  variableImpPlot(rf.abcrf.pls.NOISE, cex=0.5, cex.axis=2.0, cex.lab=2.0, n.var=50)
  
  # computed using an Out-of-Bag procedure NMAE ?????????? FAUX DEMANDER MPC
  #png(file="param_prior_errors_vs_number_of_trees.png")
  #err.oob <- err.abcrf(rf.abcrf.pls.NOISE, train.pls.NOISE, paral=TRUE)
  #dev.off()
  #plot(err.oob)
  
  ####### Pour extraire TOUTES les VarImp du graphique VarImp(stat) et ecrire dans fichier txt ###### AE #######
  # Pour visualiser toutes les variables importances en txt
  #rf.abcrf.pls.NOISE$model.rf$variable.importance # Pour visualiser les variables importances en txt
  #sort(rf.abcrf.pls.NOISE$model.rf$variable.importance) # Ordre croissant  # Ordre croissant = affichage text
  #sort(rf.abcrf.pls.NOISE$model.rf$variable.importance, decreasing = TRUE) # Ordre decroissant = affichage text
  output <- rf.abcrf.pls.NOISE$model.rf$variable.importance
  output = data.frame(output)
  write.table(output, "var_imp_estim_param_NOT_sorted_AVEC_PLS.txt", quote=FALSE, row.names=TRUE, col.names="stat_name var_imp")
  table_trie = read.table("var_imp_estim_param_NOT_sorted_AVEC_PLS.txt",header=TRUE)
  table_trie$var_imp = round(table_trie$var_imp, 6)
  output_after_trie = table_trie[order(table_trie$var_imp, decreasing=TRUE),]
  write.table(output_after_trie, "var_imp_estim_param_SORTED_AVEC_PLS.txt", quote=FALSE, row.names=FALSE)
  
  #### PREDICTIONS on a given DATA_OBS
  # If obs.poi is not a dataframe or the column names do not match,
  # you can use the following lines:
  obs.poi <- as.data.frame(obs.poi)
  colnames(obs.poi) <- colnames(xtest)
  # On récupère les scores PLS pour les données statobs
  scores_statobs.pls <- predict(pls.fit, obs.poi, type = "scores")
  scores_statobs.pls.sel <- scores_statobs.pls[,1:nComposante_sel, drop = FALSE]
  # cf. evite inverser row et col cf. R fait chier quand un seul dataset !!!!
  # scores_statobs.pls.sel <- scores_statobs.pls[,1:nComposante_sel]
  # scores_statobs.pls.sel <- as.data.frame(t(scores_statobs.pls.sel)) ###WARNING: on doit inverser row et col !!!!
  # On ajoute du NOISE (5 variables)
  monNOISEObs <- matrix(runif(5), ncol=5)
  colnames(monNOISEObs) <- c("NOISE1", "NOISE2", "NOISE3", "NOISE4", "NOISE5")
  obs.poi.pls.NOISE <- cbind(obs.poi, scores_statobs.pls.sel, monNOISEObs) # On ajoute les colonnes de scores PLS et NOISEs
  colnames(obs.poi.pls.NOISE) <- colnames(train.pls.NOISE)[-1] # Homogenisation noms des colonnes
  
  statobs_predClassic.pls.NOISE <- predict(rf.abcrf.pls.NOISE, obs=obs.poi.pls.NOISE, quantiles = c(0.5,0.05,0.95), paral = TRUE,
                                           paral.predict=TRUE, post.err.med = TRUE, ncores=how.many.cores.used.for.computation, min.node.size = 5, training=train.pls.NOISE)
  
  
  df.pls.noise<-c(statobs_predClassic.pls.NOISE$expectation,
                  statobs_predClassic.pls.NOISE$med,
                  statobs_predClassic.pls.NOISE$quantiles[2],
                  statobs_predClassic.pls.NOISE$quantiles[3],
                  statobs_predClassic.pls.NOISE$post.NMAE.med,
                  statobs_predClassic.pls.NOISE$post.NMAE.mean,
                  statobs_predClassic.pls.NOISE$prior.coverage,
                  statobs_predClassic.pls.NOISE$prior.NMAE.med,
                  statobs_predClassic.pls.NOISE$prior.NMAE.mean,
                  statobs_predClassic.pls.NOISE$model.rf$NMAE)
  
  statobs_predClassic.pls.NOISE$med
  statobs_predClassic.pls.NOISE$quantiles[2]
  statobs_predClassic.pls.NOISE$quantiles[3]
  statobs_predClassic.pls.NOISE$variance
  statobs_predClassic.pls.NOISE$post.NMAE.med
  statobs_predClassic.pls.NOISE$post.NMAE.mean
  statobs_predClassic.pls.NOISE$prior.coverage
  statobs_predClassic.pls.NOISE$prior.NMAE.med
  statobs_predClassic.pls.NOISE$prior.NMAE.mean
  statobs_predClassic.pls.NOISE$model.rf$NMAE
  
  
}
}

cat("\n")
cat("######################################################################################################################","\n")
cat("################################### ALL PARAMETER ESTIMATION RESULTS + NOISE (NO PLS)             ####################","\n")
cat("######################################################################################################################","\n")
cat("\n")
print(RESULT_ESTIM_PARAM, row.names=FALSE)

sink()




