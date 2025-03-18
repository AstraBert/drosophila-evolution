# Copyright (C) {2024} {GLM, PB, JMM, AE}
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

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

#################################################################################################################
### Delayed start of the script Rscript GOF.PRIOR.POSTERIOR.18-01-2025.R for 7200 seconds (2 hours) #############
#Sys.sleep(7200) # Temporise pendant 7200 secondes (2 heures)
Sys.sleep(20) # Temporise pendant 20 secondes
#################################################################################################################

################################################################################
current_directory <- getwd() # Récupère le répertoire de travail
print(current_directory)     # Affiche le répertoire
################################################################################

################################################################################
library(abcgof)
library(abcrf)
################################################################################
# Un brin de nettoyage plus spécifique: répertoire sim à supprimer
directory.path.to.remove <- "./sim"
# Supprimer le répertoire et son contenu (cf. recursive = TRUE)
unlink(directory.path.to.remove, recursive = TRUE)
# Fonction pour convertir des secondes en heures, minutes et secondes
convert_time <- function(time_in_seconds) {
  hours <- floor(time_in_seconds / 3600)
  minutes <- floor((time_in_seconds %% 3600) / 60)
  seconds <- (time_in_seconds %% 3600) %% 60
  return(sprintf("%02d:%02d:%05.2f", hours, minutes, seconds))
}
################################################################################
## location of data, simulation and results
swd = getwd() # swd = source work directory
folder <- swd
################################################################################
# For controling starting point and bootsrapping
myseed <- 1162
set.seed(myseed)
#################### Nbre of relicates cf. n.boot and nbre of cores for parallelization #########################################################
n.boot = 600
# pour ncores.prior et ncores.posterior
ncores.prior = 11
# From Paul: Quand il essaie de paralléliser, il doit avoir la fonction simulate_swd sur tous les coeurs, ce qui n'est pas le cas par défaut.
# Comme la parallelisation a déjà lieu sur les simus, ce n'est pas très utile de paralléliser sur le reste, j'ai mis 1 sur toutes les simus de mon côté.
if (n.boot==1) ncores.posterior = 1  # Mettre ncores.posterior = 1 si n.boot=1 !!!
if (n.boot>1) ncores.posterior = 31 # Mettre ncores.posterior = 1 ou plus si n.boot>1 !!!
ncores.posterior = 11 # Mettre ncores.posterior = 1 si n.boot=1 !!!
ncores_sim = 11
################################################################################
# Define here the type of analysis ?
REFTABLE.BIN=TRUE # In this case we use a diyabc reftableRF.bin file as starting point
GOF.PRIOR = FALSE
GOF.POSTERIOR = TRUE
post.change.parameter.order.for.simulation = FALSE
post.rej = TRUE
post.rej.loclin = FALSE
post.rej.ridge = FALSE
################################################################################
# Parameters
#n.ref.prior.all <- c(500, 1000, 2000, 3000, 4000, 5000) # total number of points in the reference dataset for prior gof
n.ref.prior.all <- c(500, 1000, 2000) # total number of points in the reference dataset for prior gof

# n.ref.hp.all <- c(100*1000) # total number of points in the reference dataset for the post and freq gof
# n.post.all <- c(1000) # number of points kept in the posterior # On ne parle plus de epsilon mais equivalent
n.ref.hp.all <- c(60*1000) # total number of points in the reference dataset for the post and freq gof
n.post.all <- c(500) # number of points kept in the posterior # On ne parle plus de epsilon mais equivalent
# n.ref.hp.all <- c(50*1000, 100*1000) # total number of points in the reference dataset for the post and freq gof
# n.post.all <- c(500, 1000, 2000) # number of points kept in the posterior # On ne parle plus de epsilon mais equivalent

split <- 0.5 # split of the posterior for the frequentist test
k_all <- 1:20 # k vector
k_range <- c(5, 20) #### range effectif continue de k values

################################################################################
# Output
datestamp_day <- "2024-12-23" #format(Sys.time(), "%Y-%m-%d")#
#res.file.name.prior <- paste0("GOF_HAfinal_estim_param_S17NGS_no_moustach_2017_7_2-1_NZ_",datestamp_day,"_PRIOR.txt")
res.file.name.prior <- "priors.txt"
#res.file.name.posterior = paste0("GOF_HAfinal_estim_param_S17NGS_NOmous_2017_7_2-1_NZ_s50000_e1p100",datestamp_day,"_POSTERIOR_rej_500boot.txt")
res.file.name.posterior <- "posteriors.txt"
################################################################################

if (REFTABLE.BIN==FALSE) {
#####################################################################################
########## WHEN STARTING FROM txt files (as for the human SNP case) #################
#####################################################################################
 ## Example = Human real data files
 file.dataobs = "statobsRF_SNP_1_2000.txt" # For "GOF posterior" only
 #file.dataobs.replica = "statobsRF_snp_2001_4000.txt"
 file.dataobs.replica = "statobsRF_snp_1_2000.txt"
 data.obs <- as.matrix(read.table(file.dataobs,header=TRUE), drop = FALSE)
 data.obs.replica <- as.matrix(read.table(file.dataobs.replica,header=TRUE), drop = FALSE)
 
 # For GOF PRIOR analysis
 n.scenario.prior.gof = 6
 prior.gof.file.ref = "yellow_points_from_human_12000snp_maf_hudson_reftableRF_allS_11000PerScen_param_S"
 #prior.gof.file.ref = "green_points_from_human_12000snp_maf_hudson_reftableRF_allS_11000PerScen_param_S"
 
 # For any GOF POSTERIOR analysis
 scen.ref <- 2 # ID number of the analysed scenario !!!!! Important pour la transformation des parametres
 scenario <- paste0(scen.ref)
 posterior.gof.file.ref = "yellow_points_from_human_12000snp_maf_hudson_reftableRF_S2S3_110000PerScen_param_S2.txt"
 #posterior.gof.file.ref = "green_points_from_human_12000snp_maf_hudson_reftableRF_S2S3_110000PerScen_param_S2.txt"
 posterior.gof.file.param <- "human_12000snp_maf_hudson_reftableRF_S2S3_110000PerScen_param_S2.txt" # Pas optimal
 param.ref.all <- as.matrix(read.table(file.path(folder, posterior.gof.file.param), header = TRUE))[1:100000, -1]
# #############################################################################################
}

#####################################################################################
########## WHEN STARTING FROM reftableRF.BIN ET statobsRF.TXT (as for SWD) ##########
#####################################################################################
if (REFTABLE.BIN==TRUE) {
###### Reading of the three key files produced by diyabcRF: 
# (i) the file corresponding to the bin format reference table (with all simulated data summarized with various summary statistics), 
# (ii) the file of the header (with various informations on the simulations), 
# and (iii) the file corresponding the observed dataset (summarized with the same set of summary statistics). 
# N = numbre of simuations one wants to load from the reference table
name.refTable <- "reftableRF.bin" # i.e. "reftableRF_microsats_13_24.bin"
name.header <- "headerRF.txt"
# name.observed.dataset <- "statobsRF_all_statobs_vectors_without_comment_lines.txt"
# name.observed.dataset.replica = "statobsRF_all_statobs_vectors_without_comment_lines.txt"
#name.observed.dataset <- "statobsRF_snp_1_2_000.txt"
name.observed.dataset <- "statobsRF.txt"
#name.observed.dataset.replica <- "statobsRF_snp_2001_4000.txt"
name.observed.dataset.replica <- "statobsRF_replica.txt"
# Warning: # Nouvel ordre scenario issu de headerRF.txt voir Step6-Q-chatGPT de INSTRUCTION-CHATGPT-constructions de scenarios_05-12-2024.txt )

scen.ref <- 1 # ID number of the analysed scenario when computing posterior GOF !!!!! 
# Important pour la transformation des parametres (étape specifique pour chaque reftable)
scenario <- paste0(scen.ref)
n.multiple.obs.datasets = 1

refTable <- readRefTable(filename = name.refTable, header=name.header)
stat.obs <- read.table(name.observed.dataset, header=TRUE)
stat.obs.replica <- read.table(name.observed.dataset.replica, header=TRUE)
data.obs = as.matrix(stat.obs)
data.obs.replica = as.matrix(stat.obs.replica)
data.obs = data.obs[1:n.multiple.obs.datasets,]
data.obs.replica = data.obs.replica[1:n.multiple.obs.datasets,] 

###### Write key variables of the refTable and RF treatments in the output file
cat("Name of the reference table:", name.refTable, "\n")
cat("Name of the observed dataset:", name.observed.dataset, "\n")
cat("Name of the replica observed dataset:", name.observed.dataset.replica, "\n")
cat("Number of simulations loaded from the reference table:", refTable$nrec, "\n")
cat("Number of scenarios (i.e. models) in the reference table:", refTable$nscen, "\n")
cat("Number of simulations available for each scenario from the loaded reference table:", refTable$nrecscen, "\n")
cat("Number of parameters recovered from the reference table:", refTable$nparam, "\n")
cat("Number of summary statistics in the reference table:", ncol(refTable$stats), "\n")

# Appliquer format pour convertir toutes les valeurs de refTable$params en écriture non scientifique
# Désactiver l'écriture scientifique
# options(scipen = 999)
# # Toutes les valeurs deviennent numeriques
# refTable$params <- apply(refTable$params, 2, as.numeric)
# head(refTable$params, n=5)
}

################################################################################
################################################################################
# START OF MAIN COMPUTATION SECTION
################################################################################
################################################################################

################################################################################
################################################################################
# Prior GOF
################################################################################
################################################################################
 if (GOF.PRIOR==TRUE) {
   
    sink(file = res.file.name.prior, split = TRUE)   
   
    if (REFTABLE.BIN==TRUE) {
    cat("#################### PRIOR GOF #################","\n")
    cat("Work directory =",current_directory,"\n")
    cat("myseed =",myseed,"\n")
    cat("k_all =",k_all,"\n")
    cat("k_range =",k_range,"\n")
  	cat("nboot =",n.boot,"\n")
    cat("\n")
    cat("Name of the reference table:", name.refTable, "\n")
    cat("Name of the observed dataset:", name.observed.dataset, "\n")
    cat("Number of simulations loaded from the reference table:", refTable$nrec, "\n")
    cat("Number of scenarios (i.e. models) in the reference table:", refTable$nscen, "\n")
    cat("Number of simulations available for each scenario from the loaded reference table:", refTable$nrecscen, "\n")
    cat("Number of parameters recovered from the reference table:", refTable$nparam, "\n")
    cat("Number of summary statistics in the reference table:", ncol(refTable$stats), "\n")
    cat("\n")
    n.scenario.prior.gof=refTable$nscen
    }
    
   for (i.scen in 1:n.scenario.prior.gof) {    
   
    ##### ALA HUMAN reftable txt JAS #######
      if (REFTABLE.BIN==FALSE) {
	  cat("#################### PRIOR GOF #################","\n")
      cat("Work directory =",current_directory,"\n")
      cat("myseed =",myseed,"\n")
      cat("k_all =",k_all,"\n")
      cat("k_range =",k_range,"\n")
	    cat("nboot =",n.boot,"\n")
      cat("\n")
      prior.gof.file.ref.i.scen = paste0(prior.gof.file.ref,i.scen,".txt")
      #prior.gof.file.ref.i.scen.replica = paste0(prior.gof.file.ref.replica,i.scen,".txt")
      #prior.gof.file.param.i.scen = paste0(prior.gof.file.param,i.scen,".txt")
      data.ref <- as.matrix(read.table(file.path(folder, prior.gof.file.ref.i.scen), header = FALSE))[1:10000, 2:131]
      #data.ref.replica <- as.matrix(read.table(file.path(folder, prior.gof.file.ref.i.scen.replica), header = FALSE))[1:10000, 2:131]
	  cat("name name of reference table file =",prior.gof.file.ref.i.scen,"\n")
	  cat("name of dataobs file =",file.dataobs,"\n")
      }
        
      #### AVEC reftableRF.bin #########
      if (REFTABLE.BIN==TRUE) {
      # Afficher la structure de refTable pour vérifier les données disponibles
      # str(refTable)
      # Convertir i.scen en caractère pour la comparaison
      i.scen.char <- as.character(i.scen)
      # Identifie les indices de refTable correspondant au scénario i.scen
      indices.scenario.i.scen <- which(refTable$scenarios == i.scen.char)
      # Extraire les statistiques pour les indices trouvés
      stats.scenario.i.scen <- as.data.frame(refTable$stats[indices.scenario.i.scen,])
      data.ref = as.matrix(stats.scenario.i.scen)
      colnames(data.ref) = colnames(data.obs)
      }
       
      cat("\n")
      cat("##########################################################","\n")
      cat("PRIOR GOF analyzed scenario =",i.scen,"\n")
      cat("\n")
      
      for (n.ref.prior in n.ref.prior.all) {
      #res_file <-  here("results", paste0(results_name_nref, "_prior.rds"))
      data.ref.n.ref.prior <- data.ref[1:n.ref.prior,]
      #param.ref = param.ref.all[1:n.ref.prior,]

      ## GOF computation
      # Capturer le temps de début
      start_time <- proc.time()

        gof.prior <- gfit(target=data.obs,
                          sumstat=data.ref.n.ref.prior,
                          nb.replicate = n.ref.prior*split,
                          score = c("lof", "kNN"),
                          k = k_all,
                          k_range = k_range,
                          norm = sd,
                          ncores = ncores.prior,
                          nboot = n.boot)

        # Capturer le temps de fin
        end_time <- proc.time()
        # Calculer la durée d'exécution
        execution_time <- end_time - start_time

        # Convertir le temps écoulé en heures, minutes et secondes
        execution_time_converted <- convert_time(execution_time["elapsed"])

        summary.lof = summary(gof.prior, score = "lof", k = "max", level = 0.95)
        summary.kNN = summary(gof.prior, score = "kNN", k = 1, level = 0.95)

        # Capture des résumés sous forme de vecteurs de chaînes de caractères
        summary.lof_text <- capture.output(print(summary.lof))
        summary.kNN_text <- capture.output(print(summary.kNN))

        cat("--------------> n.sim =", n.ref.prior,"- n.calib =", n.ref.prior*split,"\n")
        print(execution_time_converted)
        cat("\n")
        print(summary.lof)
        cat("\n")
        print(summary.kNN)
        cat("\n")
        }
     }
 }  

    if (GOF.POSTERIOR ==TRUE) {
    ################################################################################
    ################################################################################
    # Posterior GOF (rej, loclin et ridge)
    ################################################################################
    ################################################################################
      
      sink(file = res.file.name.posterior, split = TRUE)
      
      ### Préambule et partage des données 
      indexesModel <- which(refTable$scenarios == scen.ref) #### Choix du scenario pour lequel on veut faire du posterior GOF
      refTable$scenarios <- refTable$scenarios[indexesModel]
      refTable$stats <- refTable$stats[indexesModel,]
      refTable$params <- refTable$params[indexesModel,]
      
      #ELIMINER SI NECESSAIRE LES COLONNES DE PARAMETRES FIXES
      #Pour HA ala NGS = "t4m", "t3m", "t2m"
      # dim(refTable$params)
      # head(refTable$params, n=3)
      # refTable$params <- refTable$params[, !(colnames(refTable$params) %in% c("t4m", "t3m", "t2m"))]
      # dim(refTable$params)
      # head(refTable$params, n=3)
      # dim(refTable$stats)
      # head(refTable$stats[1:3,1:3], n=3)
      
      if (post.change.parameter.order.for.simulation == TRUE) {

      #### Changing param order if necessay for simulation of replicates
      refTable$params_initial_order_with_NA <- refTable$params
      # Remove columns with no data (NA columns)
      refTable$params <- refTable$params[, colSums(!is.na(refTable$params)) > 0]
      
      # dim(refTable$param)
      # head(refTable$param, n=3)
      
      # Définir le nouvel ordre des colonnes (sans NA)
      # from headerRF.txt and Step6-Q-chatGPT )
      nouvel_ordre <- c(
        "Nwat", "Nsok", "Nhw", "Nsap", "Nlia", "Nsd", "Nnc", "Nwis", "Ngen", "Ncol", "Nbra", "NGan2", "NGan1", "NGhw", "NGan3", "NAC",
        "tbra", "DBbra", "NBbra", "raabra", "tnc", "DBnc", "NBnc", "twis", "DBwis", "NBwis", "tgen", "DBgen", "NBgen", "tcol", "DBcol",
        "NBcol", "tan3", "DBan3", "NBan3", "tsd", "DBsd", "NBsd", "raasd", "twat", "DBwat", "NBwat", "raawat", "tsok", "DBsok", "NBsok",
        "tan2", "DBan2", "NBan2", "raan2", "tan1", "DBan1", "NBan1", "raan1", "th", "DBh", "NBh", "tGhw", "DBGhw", "NBGhw", "tj", "tc",
        "µmic_1", "pmic_1", "snimic_1")
      # Nouvel ordre S2 = scenario fraimout complet bra = wat+sd
        # nouvel_ordre <- c(
          # "Nwat", "Nsok", "Nhw", "Nsap", "Nlia", "Nsd", "Nnc", "Nwis", "Ngen", "Ncol",
          # "Nbra", "Narg", "Ncl1", "Ncl2", "Ncl3", "NAC", "tc3", "DBc3", "NBc3", "tc2",
          # "DBc2", "NBc2", "tc1", "DBc1", "NBc1", "tfarg", "DBarg", "NBarg", "raaarg", "tfbra",
          # "DBbra", "NBbra", "raabra", "tfcol", "DBcol", "NBcol", "tfgen", "DBgen", "NBgen",
          # "tfwis", "DBwis", "NBwis", "tfnc", "DBnc", "NBnc", "tfsd", "DBsd", "NBsd", "tfsok",
          # "DBsok", "NBsok", "raasok", "tfwat", "DBwat", "NBwat", "raawat", "tfh", "DBh", "NBh",
          # "tfj", "tfc", "µmic_1", "pmic_1", "snimic_1")
          
      # Vérifier que toutes les colonnes du nouvel ordre existent dans le dataframe d'origine
      if (!all(nouvel_ordre %in% colnames(refTable$params))) stop("Certaines colonnes du nouvel ordre sont manquantes dans le dataframe d'origine.")
       
      # ReCréer le dataframe refTable$params avec le nouvel ordre spécifié et en enlevant les colonnes vides (avec des NA)
      refTable$params <- refTable$params[, nouvel_ordre]
      }

      #head(refTable$params, n=3)
      
      if (REFTABLE.BIN==TRUE) {
      cat("#################### POSTERIOR GOF - analysed scenario ",scen.ref,"#################","\n")
      cat("Work directory =",current_directory,"\n")
      cat("myseed =",myseed,"\n")
      cat("k_all =",k_all,"\n")
      cat("k_range =",k_range,"\n")
      cat("\n")
      cat("Name of the reference table:", name.refTable, "\n")
      cat("Name of the observed dataset:", name.observed.dataset, "\n")
      cat("Name of the replica observed dataset:", name.observed.dataset.replica, "\n")
      cat("Number of simulations loaded from the reference table:", refTable$nrec, "\n")
      cat("Number of scenarios (i.e. models) in the reference table:", refTable$nscen, "\n")
      cat("Number of simulations available for each scenario from the loaded reference table:", refTable$nrecscen, "\n")
      cat("Number of parameters recovered from the reference table:", refTable$nparam, "\n")
      cat("Number of summary statistics in the reference table:", ncol(refTable$stats), "\n")
      cat("\n")
      data.ref <- as.matrix(refTable$stats)
      #data.ref <- as.matrix(refTable$stats)[1:100000,]
      #data.ref.replica <- as.matrix(read.table(file.path(folder, posterior.gof.file.ref.i.scen.replica), header = FALSE))[1:10000, 2:131]
      colnames(data.ref) = colnames(data.obs)
      #colnames(data.ref.replica) = colnames(data.obs)
      param.ref.all <- as.matrix(refTable$params)
      #param.ref.all <- as.matrix(refTable$params)[1:100000,]
      
      # dim(param.ref.all)
      # head(param.ref.all, n=3)
      # colnames(param.ref.all)
      # 
      # dim(data.ref)
      # head(data.ref[1:3,1:3], n=3)
      # colnames(data.ref)
      
      
      # Vérifier si des colonnes contiennent des NA
      cols_with_na <- colnames(refTable$params)[colSums(is.na(refTable$params)) > 0]
      
      # Afficher les colonnes avec des NA
      if (length(cols_with_na) > 0) {
        print(paste("Les colonnes avec des NA sont :", paste(cols_with_na, collapse = ", ")))
       } 
      }
      
      if (REFTABLE.BIN==FALSE) {
      cat("#################### POSTERIOR GOF - analysed scenario ",scen.ref,"#################","\n")
      cat("Work directory =",current_directory,"\n")
      cat("myseed =",myseed,"\n")
      cat("k_all =",k_all,"\n")
      cat("k_range =",k_range,"\n")
      cat("\n")
      
      data.ref <- as.matrix(read.table(file.path(folder, posterior.gof.file.ref), header = FALSE))[1:100000, 2:131]
      #data.ref.replica <- as.matrix(read.table(file.path(folder, posterior.gof.file.ref.i.scen.replica), header = FALSE))[1:10000, 2:131]
      colnames(data.ref) = colnames(data.obs)
      #colnames(data.ref.replica) = colnames(data.obs)
      param.ref.all <- as.matrix(read.table(file.path(folder, posterior.gof.file.param.i.scen), header = TRUE))[1:100000, -1]
      }
      
 
        ################################################################################
        ################################################################################
        # DEBUT PREPARATION OF SIMULATION
        ################################################################################
        ################################################################################
        # Prepare simulation function
        #source("sim_human_AE-01-05-2024.R")
        source("../../../scripts/sim_swd.R")
        path_to_headers <- swd

        seed <- myseed
        ncores_sim <- ncores_sim
        datestamp_simu <- datestamp_day
        #if (!test_in_alternative) datestamp_simu <- paste0(datestamp_simu, "_H0")
        path_to_sim <- setup_sim(path_to_headers, datestamp_simu, myseed, ncores_sim)

        # Fonction originale Paul
        # sim.fun.ref <- function(params, ncores_sim, path_to_sim) {
        #   params[, !grepl("ra", colnames(params))] <- round(params)[, !grepl("ra", colnames(params))]
        #   params <- cbind(scen.ref, params)
        #   colnames(params)[1] <- "scenario"
        #   new_sumstats <- simulate_human(params, scen.ref, path_to_sim, ncores_sim)
        #   return(as.matrix(new_sumstats)[, -1])
        #   }

        #  Ajouter scen.ref et rep.run comme paramètres pour la fonction, afin que leur utilisation soit explicite et contrôlée de manière externe.
        sim.fun.ref <- function(params, ncores_sim, path_to_sim, scen.ref) {
            # # Liste des motifs à exclure pour arrondir (colnames qui contiennent "raa", "µmic_1", "pmic_1", ou "snimic_1")
            # cols_to_exclude <- c("raa", "µmic_1", "pmic_1", "snimic_1")
            # Liste des motifs à exclure pour arrondir (colnames qui contiennent "raa", "µmic_1", "pmic_1", ou "snimic_1")
            cols_to_exclude <- c("ra", "µmic_1", "pmic_1", "snimic_1")
            # Appliquer !grepl pour exclure les colonnes correspondant aux motifs ci-dessus
            cols_to_round <- !grepl(paste(cols_to_exclude, collapse = "|"), colnames(params))
            params[, cols_to_round] <- round(params[, cols_to_round])
            # Ajouter une colonne de scénario référence au début des paramètres
            params <- cbind(scenario = scen.ref, params)
            # Simuler les statistiques sommaires
            #new_sumstats <- simulate_human(params, scen.ref, path_to_sim, ncores_sim)
            new_sumstats <- simulate_swd(params, scen.ref, path_to_sim, ncores_sim)
            # Retourner les statistiques sommaires en excluant la première colonne
            return(as.matrix(new_sumstats)[, -1])
        }
        
        # ########### Transformation of parameters (si regression et si human!!!) ###############
        # # All param: any order
        # N1_bounds <- c(1000.0,100000.0)
        # N2_bounds <- c(1000.0,100000.0)
        # N3_bounds <- c(1000.0,100000.0)
        # N4_bounds <- c(1000.0,100000.0)
        # t1_bounds <- c(1.0,30.0)
        # t2_bounds <- c(100.0,10000.0)
        # d3_bounds <- c(0.0,50.0)
        # Nbn3_bounds <- c(5,500)
        # d4_bounds <- c(0.0,50.0)
        # Nbn4_bounds <- c(5,500)
        # N34_bounds <- c(1000.0,100000.0)
        # t3_bounds <- c(100.0,10000.0)
        # d34_bounds <- c(0.0,50.0)
        # Nbn34_bounds <- c(5,500)
        # t4_bounds <- c(100.0,10000.0)
        # Na_bounds <- c(100.0,10000.0)
        # ra_bounds <- c(0.05,0.95)
        # t11_bounds <- c(1.0,30.0)
        # t22_bounds <- c(100.0,10000.0)
        # t33_bounds <- c(100.0,10000.0)
        # t44_bounds <- c(100.0,10000.0)
        # 
        # # Param dans "le bon ordre"
        # param_bounds <- cbind(N1_bounds, N2_bounds, N3_bounds, N4_bounds,
        #                       t1_bounds, ra_bounds, t2_bounds,
        #                       d3_bounds, Nbn3_bounds,
        #                       d4_bounds, Nbn4_bounds,
        #                       N34_bounds, t3_bounds,
        #                       d34_bounds, Nbn34_bounds,
        #                       t4_bounds, Na_bounds)
        # 
        # param_lower_bound <- param_bounds[1, ]
        # param_upper_bound <- param_bounds[2, ]
        # 
        # ind_times_order <- c(7, 13, 16) # necessaire car t4<t3<t2
        # param_transform <- rep("logit", 17)
        # param_transform[ind_times_order] <- "none"
        # 
        # names(param_transform) <- names(param_lower_bound) <- names(param_upper_bound) <- colnames(param.ref.all) # Tous le meme nom
        # 
        # # Pour les parametres entier (pas du type ra) = un peu avant ou apres pour pouvoir avoir inf et sup apres logit
        # param_lower_bound[!grepl("ra", names(param_lower_bound))] <- param_lower_bound[!grepl("ra", names(param_lower_bound))] - 0.49
        # param_upper_bound[!grepl("ra", names(param_upper_bound))] <- param_upper_bound[!grepl("ra", names(param_upper_bound))] + 0.49
        # 
        # trans_and_back_no_order <- gofabcpkg:::check_param_transform(param.ref.all, param_transform, param_lower_bound, param_upper_bound)
        # 
        # trans_order <- function(x, ind, alpha = 0.25) {
        #   y <- x
        #   for (ii in 2:length(ind)) {
        #     y[, ind[ii]] <- x[, ind[ii]] - 1 - x[, ind[ii-1]]
        #     y[y[, ind[ii]] == 0, ind[ii]] <- y[y[, ind[ii]] == 0, ind[ii]] + alpha
        #   }
        #   z <- y
        #   z[, ind[1]] <- gofabcpkg:::logit(y[, ind[1]],
        #                                    param_lower_bound[ind[1]],
        #                                    param_upper_bound[ind[1]] - 2)
        #   for (ii in 2:length(ind)) {
        #     z[, ind[ii]] <- gofabcpkg:::logitVec(y[, ind[ii]],
        #                                          0,
        #                                          param_upper_bound[ind[ii]] - (length(ind) - ii + 1) - x[, ind[ii-1]])
        #   }
        #   return(z)
        # }
        # trans_back <- function(z, ind, alpha = 0.25) {
        #   x <- y <- z
        #   y[, ind[1]] <- gofabcpkg:::logistic(z[, ind[1]],
        #                                       param_lower_bound[ind[1]],
        #                                       param_upper_bound[ind[1]] - 2)
        #   x[, ind[1]] <- y[, ind[1]]
        #   for (ii in 2:length(ind)) {
        #     y[, ind[ii]] <- gofabcpkg:::logisticVec(z[, ind[ii]],
        #                                             0,
        #                                             param_upper_bound[ind[ii]] - (length(ind) - ii + 1) - x[, ind[ii-1]])
        #     x[, ind[ii]] <- y[, ind[ii]] + x[, ind[ii-1]] + 1
        #   }
        #   # rounding
        #   x[, ind] <- round(x[, ind])
        #   return(x)
        # }
        # colnames(param.ref.all)[ind_times_order]
        # trans_par_order <- function(x) {
        #   y <- trans_and_back_no_order$transform(x)
        #   y <- trans_order(y, ind_times_order)
        #   colnames(y) <- colnames(x)
        #   return(y)
        # }
        # backtrans_par_order <- function(y) {
        #   x <- trans_back(y, ind_times_order)
        #   x <- trans_and_back_no_order$back_transform(x)
        #   colnames(x) <- colnames(y)
        #   return(x)
        # }
        # 
        # # tt <- trans_par_order(param.ref.all[1:3, ])
        # # ttt <- backtrans_par_order(tt)
        # # all.equal(ttt, param.ref.all[1:3, ])
        # # #sim.fun.ref(ttt, 6, path_to_sim)
        # # pp <- param.ref.all[1:2, ]
        # # #pp[, c(7, 13, 16)] <- rbind(c(100, 101, 102), c(9998, 9999, 10000))
        # # tt <- trans_par_order(pp)
        # # ttt <- backtrans_par_order(tt)
        # # all.equal(ttt, pp)
        # 
        # trans_and_back <- list(transform = trans_par_order,
        #                        back_transform = backtrans_par_order)

 
  ##################################
        
      
  ################################################################################
  ################################################################################
  # FIN PREPARATION OF SIMULATION
  ################################################################################
  ################################################################################
        
        ################################################################################
        ################################################################################
        # OUTPUT
        ################################################################################
        ################################################################################

        if (REFTABLE.BIN == FALSE) {
        cat("#####################################","\n")
        cat("Analysed scenario = S",scen.ref,"\n")
        cat("Reftable file name =",posterior.gof.file.ref,"\n")
        #cat("refscen.replica =",posterior.gof.file.ref.i.scen.replica,"\n")
        cat("dataobs =",file.dataobs,"\n")
        cat("dataobs.replica =",file.dataobs.replica,"\n")
        cat("myseed =",myseed,"\n")
        cat("k_all =",k_all,"\n")
        cat("k_range =",k_range,"\n")
        cat("\n")
        
        cat("\n")
        cat("#################### posterior GOF #################","\n")
        cat("refscen =",posterior.gof.file.ref,"\n")
        #cat("refscen.replica =",posterior.gof.file.ref.i.scen.replica,"\n")
        cat("dataobs =",file.dataobs,"\n")
        cat("dataobs.replica =",file.dataobs.replica,"\n")
        cat("\n")
        }
        
        if (REFTABLE.BIN == TRUE) {
          cat("\n")
          cat("#################### posterior GOF #################","\n")
          cat("Analysed scenario = S",scen.ref,"\n")
          cat("Reftable file name =",name.refTable,"\n")
          #cat("refscen.replica =",name.refTable.replica,"\n")
          cat("dataobs =",name.observed.dataset,"\n")
          cat("dataobs.replica =",name.observed.dataset.replica,"\n")
          cat("myseed =",myseed,"\n")
          cat("k_all =",k_all,"\n")
          cat("k_range =",k_range,"\n")
          cat("\n")
          }
          
    ################################################################################
    # Loop over n.ref
    for (n.ref in n.ref.hp.all) {

      data.ref.n.ref <- data.ref[1:n.ref,]
      param.ref.n.ref = param.ref.all[1:n.ref,]

      ####################################################################
      ## loop over n.post and hence eps

      for (n.post in n.post.all) {
        eps <- n.post / n.ref

      if (post.rej == TRUE) {
        ########## Rejection ##############################################
        # Capturer le temps de début
        start_time <- proc.time()

        gof.hp.rej <- hpgfit(target = data.obs,
                             target.replica = data.obs.replica ,
                             param = param.ref.n.ref,
                             sumstat = data.ref.n.ref,
                             sim.fun= sim.fun.ref,
                             method = c("rejection"),
                             kernel = c("epanechnikov"),
                             lambda = c(0.0001, 0.001, 0.01),
                             #param_transform = "none",
                             # trans.fun = trans_and_back$transform,
                             # back.trans.fun = trans_and_back$back_transform,
                             # param_lower_bound = param_lower_bound,
                             # param_upper_bound = param_upper_bound,
                             score = c("lof", "kNN"),
                             k = k_all,
                             k_range = k_range,
                             eps = eps,
                             split = 0.5,
                             norm = sd,
                             ncores = ncores.posterior,
                             nboot = n.boot,
                             ncores_sim = ncores_sim,
                             path_to_sim = path_to_sim,
                             scen.ref = scen.ref)

        # Capturer le temps de fin
        end_time <- proc.time()
        # Calculer la durée d'exécution
        execution_time <- end_time - start_time

        # Convertir le temps écoulé en heures, minutes et secondes
        execution_time_converted <- convert_time(execution_time["elapsed"])

        summary.lof.95 = summary(gof.hp.rej, score = "lof", k = "max", level = 0.95)
        summary.kNN.95 = summary(gof.hp.rej, score = "kNN", k = 1, level = 0.95)
        #summary.lof.90 = summary(gof.hp.rej, score = "lof", k = "max", level = 0.90)
        #summary.kNN.90 = summary(gof.hp.rej, score = "kNN", k = 1, level = 0.90)

        cat("--------------> REJ: n.sim =", n.ref,"  n.post = ",n.post, "  eps = ",eps, "  n.calib =", n.post*split,"\n")
        print(execution_time_converted)
        cat("\n")
        print(summary.lof.95)
        cat("\n")
        print(summary.kNN.95)
        cat("\n")

        if (post.rej.loclin == TRUE) {
          ########## rej.loclin ##############################################
          start_time <- proc.time()

          gof.hp.loclin  <- hpgfit(target = data.obs,
                               target.replica = data.obs.replica ,
                               param = param.ref.n.ref,
                               sumstat = data.ref.n.ref,
                               sim.fun= sim.fun.ref,
                               method = "loclinear",
                               kernel = c("epanechnikov"),
                               lambda = c(0.0001, 0.001, 0.01),
                               #param_transform = "none",
                               trans.fun = trans_and_back$transform,
                               back.trans.fun = trans_and_back$back_transform,
                               # param_lower_bound = param_lower_bound,
                               # param_upper_bound = param_upper_bound,
                               score = c("lof", "kNN"),
                               k = k_all,
                               k_range = k_range,
                               eps = eps,
                               split = 0.5,
                               norm = sd,
                               ncores = ncores.posterior,
                               nboot = n.boot,
                               ncores_sim = ncores_sim,
                               path_to_sim = path_to_sim,
                               scen.ref = scen.ref)
#							 ,
#                             n.post = n.post) ### AE?

          # Capturer le temps de fin
          end_time <- proc.time()
          # Calculer la durée d'exécution
          execution_time <- end_time - start_time

          # Convertir le temps écoulé en heures, minutes et secondes
          execution_time_converted <- convert_time(execution_time["elapsed"])

          summary.lof.95 = summary(gof.hp.loclin, score = "lof", k = "max", level = 0.95)
          summary.kNN.95 = summary(gof.hp.loclin, score = "kNN", k = 1, level = 0.95)
          #summary.lof.90 = summary(gof.hp.loclin, score = "lof", k = "max", level = 0.90)
          #summary.kNN.90 = summary(gof.hp.loclin, score = "kNN", k = 1, level = 0.90)

          cat("--------------> LOCLIN: n.sim =", n.ref,"  n.post = ",n.post, "  eps = ",eps, "  n.calib =", n.post*split,"\n")
          print(execution_time_converted)
          cat("\n")
          print(summary.lof.95)
          cat("\n")
          print(summary.kNN.95)
          cat("\n")
          # cat("\n")
          # print(summary.lof.90)
          # cat("\n")
          # print(summary.kNN.90)
          # cat("\n")

          # Écrire le contenu combiné dans un fichier txt
          cat("\n")
          cat("--------------> LOCLIN: n.sim =", n.ref,"  n.post = ",n.post, "  eps = ",eps, "  n.calib =", n.post*split,"\n")
          print(execution_time_converted)
          cat("\n")
          print(summary.lof.95)
          cat("\n")
          print(summary.kNN.95)
          cat("\n")
          # cat("\n")
          # print(summary.lof.90)
          # cat("\n")
          # print(summary.kNN.90)
          # cat("\n")
        }

        if (post.rej.ridge == TRUE) {
          ########## rej.ridge ##############################################
          start_time <- proc.time()

          gof.hp.ridge  <- hpgfit(target = data.obs,
                                   target.replica = data.obs.replica ,
                                   param = param.ref.n.ref,
                                   sumstat = data.ref.n.ref,
                                   sim.fun= sim.fun.ref,
                                   method = "ridge",
                                   kernel = c("epanechnikov"),
                                   lambda = c(0.0001, 0.001, 0.01),
                                   #param_transform = "none",
                                   trans.fun = trans_and_back$transform,
                                   back.trans.fun = trans_and_back$back_transform,
                                   # param_lower_bound = param_lower_bound,
                                   # param_upper_bound = param_upper_bound,
                                   score = c("lof", "kNN"),
                                   k = k_all,
                                   k_range = k_range,
                                   eps = eps,
                                   split = 0.5,
                                   norm = sd,
                                   ncores = ncores.posterior,
                                   nboot = n.boot,
                                   ncores_sim = ncores_sim,
                                   path_to_sim = path_to_sim,
                                   scen.ref = scen.ref)

          # Capturer le temps de fin
          end_time <- proc.time()
          # Calculer la durée d'exécution
          execution_time <- end_time - start_time

          # Convertir le temps écoulé en heures, minutes et secondes
          execution_time_converted <- convert_time(execution_time["elapsed"])

          summary.lof.95 = summary(gof.hp.ridge, score = "lof", k = "max", level = 0.95)
          summary.kNN.95 = summary(gof.hp.ridge, score = "kNN", k = 1, level = 0.95)
          # summary.lof.90 = summary(gof.hp.ridge, score = "lof", k = "max", level = 0.90)
          # summary.kNN.90 = summary(gof.hp.ridge, score = "kNN", k = 1, level = 0.90)

          cat("--------------> RIDGE: n.sim =", n.ref,"  n.post = ",n.post, "  eps = ",eps, "  n.calib =", n.post*split,"\n")
          print(execution_time_converted)
          cat("\n")
          print(summary.lof.95)
          cat("\n")
          print(summary.kNN.95)
          cat("\n")
          # cat("\n")
          # print(summary.lof.90)
          # cat("\n")
          # print(summary.kNN.90)
          # cat("\n")

          # Écrire le contenu combiné dans un fichier txt
          cat("\n")
          cat("--------------> RIDGE: n.sim =", n.ref,"  n.post = ",n.post, "  eps = ",eps, "  n.calib =", n.post*split,"\n")
          print(execution_time_converted)
          cat("\n")
          print(summary.lof.95)
          cat("\n")
          print(summary.kNN.95)
          cat("\n")
          # cat("\n")
          # print(summary.lof.90)
          # cat("\n")
          # print(summary.kNN.90)
          # cat("\n")
          }
        }
      }
    }
  }

sink()


