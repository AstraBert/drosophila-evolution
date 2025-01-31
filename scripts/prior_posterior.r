# prior_posterior check_with_nbest_datasets_from_euclidian_distance: AE+JMM 13-01-2024
# WARNING: put the file distance.cpp in the same directory !!!

#### Initial cleaning and general options #################
# Remove all objects in the global environment
rm(list = ls())
# Free memory used by deleted objects
gc()
# Clear all open graphics (if using RStudio)
if (!is.null(dev.list())) dev.off()
# Reset options
options(default = TRUE)
#options(prompt = FALSE)
options(ask = FALSE)
# Clear history (optional)
cat("", file = ".Rhistory")

################
################
library(abcrf)
library(doRNG)
library(Rcpp)
library(dplyr)
sourceCpp("../../scripts/distances.cpp") # Warning: fast computation of Euclidean distances
options(max.print=100000)
SEED = TRUE
if (SEED) # If we want to set the seed for the random number generator: usually 1967
{
  seed = 1967
  set.seed(seed)  
}

# File names
outputfile <- "prior_posterior_checking_.txt"
# Reference file
file.reference.sim = "reftableRF.bin"
# Header file
file.header = "headerRF.txt"
# Observed datafile
file.obs.poi = "statobsRF.txt"

# Parameter values to implement
nparticules.reftable <- 60000
compute.pval.stats.prior = TRUE
compute.pval.stats.posterior = TRUE  # If we want to compute the p-values of the summary statistics for each scenario
if (compute.pval.stats.posterior) Threshold.SumStats.Distrib = 0.05 # Threshold to determine the number of (best) particles SELECTED for computing SumStats p-values (equal to 1 for prior checking)

# Import the reference table and observed data
# reftable <- readRefTable(filename = file.reference.sim, header=file.header)
# N.rec <- reftable$nrec
# N.rec
###################
reftable <- readRefTable(filename = file.reference.sim, header=file.header,N=nparticules.reftable)
obs.poi <- read.table(file.obs.poi,header=TRUE)
PARAM.S <- as.data.frame(reftable$params) # the parameter values
STAT.S <- as.data.frame(reftable$stats) # the summary statistics
#index_scen <- as.vector(reftable$scenarios) # the model indexes
index_scen <- reftable$scenarios # the model indexes
N.scen <- reftable$nscen # the number of scenarios
N.rec <- reftable$nrec # the number of simulations
N.param <- ncol(reftable$param) # the number of parameters
N.stat <- ncol(reftable$stats) # the number of summary statistics
n.param.scen = matrix(0,N.scen)
pvalue = Dobs = mean.dsim = matrix(0,N.scen)
n.Q=9 # Number of quantiles
Q.Dsim = array(0,dim = c(N.scen,n.Q))
obs.poi.mad = as.matrix(obs.poi)
scen.stat.pval = matrix(0,N.stat,N.scen)
PRIOR.scen.stat.pval = matrix(0,N.stat,N.scen)
stat.names = matrix(colnames(STAT.S))
scen.names = paste("Scen", 1:N.scen, sep = "") 
colnames(scen.stat.pval) = colnames(PRIOR.scen.stat.pval) = scen.names
rownames(scen.stat.pval) = rownames(PRIOR.scen.stat.pval) = stat.names
SUM.outliers.1p1000=SUM.outliers.1p100=SUM.outliers.5p100=SUM.outliers.95p100=SUM.outliers.99p100=SUM.outliers.999p1000= matrix(0,N.scen)
scen.stat.pval.iter = matrix(0,N.stat,N.scen)

sink(file = outputfile, split = TRUE)
### Some details on my data (to run first to adjust analysis design)
cat("Output file name:",outputfile,"\n")
cat("Reference file name:",file.reference.sim,"\n")
cat("Header file name:",file.header,"\n")
cat("StatObs file name:",file.obs.poi,"\n")
cat("\n")
cat(sprintf("Number of simulations in the reference table:"),N.rec,"\n")
cat(sprintf("Number of scenarios (i.e. models) in the reference table:"),N.scen,"\n")
cat(sprintf("Number of parameters recovered from the reference table:"),N.param,"\n")
cat(sprintf("Summary statistics for the analysis:"),N.stat,"\n")
cat("Details about the data to be analyzed","\n")
for (j in 1:N.scen) cat("Data scenario #",j, "Number of particles = ",nrow(STAT.S[index_scen==j,]),"\n")
cat(sprintf("Number of summary statistics for the analysis:"),N.stat,"\n")
cat(sprintf("Posterior data Threshold:"),Threshold.SumStats.Distrib,"\n")
cat("\n")

############# Computation of P-values for the summary statistics (from PRIOR distributions cf. threshold = 1.0) ###################################################
if (compute.pval.stats.prior)
{ 
  cat("######################################,\n")
  cat("########### PRIOR ANALYSIS ###########","\n")
  cat("######################################,\n")
  
  N.sel.Sum.Stat.Row.Scen.Prior = matrix(0,N.scen)  
  for (j in 1:N.scen) 
  {
    N.sel.Sum.Stat.Row.Scen.Prior[j] = round(nrow(STAT.S[index_scen==j,])*1.0) # Computation on priors and hence 1.0
    sumstat.scen = matrix(0,nrow(STAT.S[index_scen==j,]),ncol(STAT.S))
    sumstat.scen.mad = matrix(0,nrow(STAT.S[index_scen==j,]),ncol(STAT.S))
    ordered.sumstat.scen.mad = matrix(0,N.sel.Sum.Stat.Row.Scen.Prior[j],ncol(STAT.S))
    sumstat.scen = STAT.S[index_scen==j,]
    sumstat.scen.mad = sumstat.scen 
    # For the TESTED scenario (PRIOR PODS) calculate the MAD, normalize using sumstats from the reftable and sumstats from the pods
    for (i in 1:N.stat)
    {
      mado.stat = mad(data.matrix(STAT.S[index_scen==j,i]))
      if (mado.stat > 0)
      {
        sumstat.scen.mad[,i] = sumstat.scen[,i]/mado.stat
        obs.poi.mad[,i] = obs.poi[,i]/mado.stat
      }
    }  
    sumstat.scen.mad = as.matrix(sumstat.scen.mad)
    sumstat.scen = as.matrix(sumstat.scen)
    
    for (n in 1:nrow(obs.poi))
    {
      # For the TESTED scenario, select the n best pods 
      target.pod = obs.poi.mad[n,]
      dist.obs.sim= calDIST(unlist(target.pod),unlist(sumstat.scen.mad))
      list.order.dist.pods = order(dist.obs.sim)[1:N.sel.Sum.Stat.Row.Scen.Prior[j]]
      ordered.sumstat.scen.mad = sumstat.scen.mad[list.order.dist.pods,]
      # Computation of the probability of each summary stat
      for (i in 1:N.stat) 
      {
        PRIOR.scen.stat.pval[i,j] <- mean(obs.poi.mad[n,i] >= ordered.sumstat.scen.mad[,i])
      }  
    }
  }  
  
  # Displaying SumStats distributions and P-values in the output file (WARNING: computation using PRIOR DISTRIBUTIONS)
  cat("\n")
  cat("Location of the observed SumStats within the distributions of simulated SumStats (WARNING: computation using PRIOR DISTRIBUTIONS cf. idem Threshold = 1.0)","\n")
  cat("\n")
  for (j in 1:N.scen)
  {  
    cat("Data scenario #",j,"\n")
    cat("Number of particles AVAILABLE for computing SumStats p-values = ",nrow(STAT.S[index_scen==j,]),"\n")
    cat("Number of particles SELECTED for computing SumStats p-values = ",N.sel.Sum.Stat.Row.Scen.Prior[j],"\n")
  }
  cat("\n") 
  cat("Probability that the SumStat value of the SIMULATED datasets is INFERIOR or EQUAL to those of the OBSERVED dataset (PRIORS)","\n") 
  print(PRIOR.scen.stat.pval)
  cat("\n")
  cat(sprintf("Mean value and Quantiles of SumStat p-values (0.05,0.10,0.50,0.90,0.95) for each scenario"),"\n")
  for (j in 1:N.scen) 
  {
    cat(sprintf("Scenario "),j,sprintf(" Mean: "),round(mean(PRIOR.scen.stat.pval[,j],prob=c(0.05,0.1,0.5,0.90,0.95)),digits=3)
        ,sprintf(" Q: "),round(quantile(PRIOR.scen.stat.pval[,j],prob=c(0.05,0.1,0.5,0.90,0.95)),digits=3),"\n")
  }
  cat("\n")
  cat("Number of outlying SumStats over a total of",N.stat,"summary statistics","\n")
  for (j in 1:N.scen) 
  {
    SUM.outliers.1p1000[j] <- sum(PRIOR.scen.stat.pval[,j] <= 0.001) 
    SUM.outliers.1p100[j] <- sum(PRIOR.scen.stat.pval[,j] <= 0.01)
    SUM.outliers.5p100[j] <- sum(PRIOR.scen.stat.pval[,j] <= 0.05)
    SUM.outliers.95p100[j] <- sum(PRIOR.scen.stat.pval[,j] >= 0.95) 
    SUM.outliers.99p100[j] <- sum(PRIOR.scen.stat.pval[,j] >= 0.99)
    SUM.outliers.999p1000[j] <- sum(PRIOR.scen.stat.pval[,j] >= 0.999)
    cat("Scenario #",j,"\n")
    cat("SUM.outliers.1p1000 = ",SUM.outliers.1p1000[j]," ( Expectation = ",round(N.stat*1/1000),")","\n")
    cat("SUM.outliers.1p100 = ",SUM.outliers.1p100[j]," ( Expectation = ",round(N.stat*1/100),")","\n")
    cat("SUM.outliers.5p100 = ",SUM.outliers.5p100[j]," ( Expectation = ",round(N.stat*5/100),")","\n")
    cat("SUM.outliers.95p100 = ",SUM.outliers.95p100[j]," ( Expectation = ",round(N.stat*5/100),")","\n")
    cat("SUM.outliers.99p100 = ",SUM.outliers.99p100[j]," ( Expectation = ",round(N.stat*1/100),")","\n")
    cat("SUM.outliers.999p1000 = ",SUM.outliers.999p1000[j]," (Expectation = ",round(N.stat*1/1000),")","\n")
    cat("\n")
  }
}

############# Generate posterior reftables with threshold = Threshold.SumStats.Distrib ###################################################
if (compute.pval.stats.posterior)
  cat("\n")
cat("\n")
cat("######################################,\n")
cat("########### POSTERIOR ANALYSIS #######","\n")
cat("######################################,\n")
{ 
  N.sel.Sum.Stat.Row.Scen = matrix(0,N.scen)  
  for (j in 1:N.scen) 
  {
    N.sel.Sum.Stat.Row.Scen[j] = round(nrow(STAT.S[index_scen==j,])*Threshold.SumStats.Distrib)
    sumstat.scen = matrix(0,nrow(STAT.S[index_scen==j,]),ncol(STAT.S))
    PARAM.S.scen = matrix(0,nrow(PARAM.S[index_scen==j,]),ncol(PARAM.S))
    sumstat.scen.mad = matrix(0,nrow(STAT.S[index_scen==j,]),ncol(STAT.S))
    ordered.sumstat.scen.mad = matrix(0,N.sel.Sum.Stat.Row.Scen[j],ncol(STAT.S))
    sumstat.scen = STAT.S[index_scen==j,]
    PARAM.S.scen = PARAM.S[index_scen==j,]
    sumstat.scen.mad = sumstat.scen 
    
    for (i in 1:N.stat)
    {
      mado.stat = mad(data.matrix(STAT.S[index_scen==j,i]))
      if (mado.stat > 0)
      {
        sumstat.scen.mad[,i] = sumstat.scen[,i]/mado.stat
        obs.poi.mad[,i] = obs.poi[,i]/mado.stat
      }
    }  
    sumstat.scen.mad = as.matrix(sumstat.scen.mad)
    sumstat.scen = as.matrix(sumstat.scen)
    
    for (n in 1:nrow(obs.poi))
    {
      # For the TESTED scenario, select the n best pods 
      target.pod = obs.poi.mad[n,]
      dist.obs.sim= calDIST(unlist(target.pod),unlist(sumstat.scen.mad))
      list.order.dist.pods = order(dist.obs.sim)[1:N.sel.Sum.Stat.Row.Scen[j]]
      ordered.sumstat.scen.mad = sumstat.scen.mad[list.order.dist.pods,]
      # Computation of the probability of each summary stat
      for (i in 1:N.stat) 
      {
        scen.stat.pval[i,j] <- mean(obs.poi.mad[n,i] >= ordered.sumstat.scen.mad[,i])
      }  
    }
  }    
}

############# Computation of P-values for the summary statistics 
##### (using POSTERIOR distributions)###################################################

# Displaying SumStats distributions and P-values
cat("\n")
cat("Location of the observed SumStats within the distributions of simulated SumStats: POSTERIOR DISTRIBUTIONS with Threshold = ",Threshold.SumStats.Distrib,"\n")
cat("\n")
for (j in 1:N.scen)
{  
  cat("Data scenario #",j,"\n")
  cat("Number of particles AVAILABLE for computing SumStats p-values = ",nrow(STAT.S[index_scen==j,]),"\n")
  cat("Number of particles SELECTED for computing SumStats p-values = ",N.sel.Sum.Stat.Row.Scen[j],"(Threshold for computing the distributions of SumStats = ",Threshold.SumStats.Distrib,")","\n")
}
cat("\n") 
cat("Probability that the SumStat value of the SIMULATED datasets is INFERIOR or EQUAL to those of the OBSERVED dataset","\n") 
cat("\n") 
print(scen.stat.pval)
cat("\n")
cat(sprintf("Mean value and Quantiles of SumStat p-values (0.05,0.10,0.50,0.90,0.95) for each scenario"),"\n")
for (j in 1:N.scen) 
{
  cat(sprintf("Scenario "),j,sprintf(" Mean: "),round(mean(scen.stat.pval[,j],prob=c(0.05,0.1,0.5,0.90,0.95)),digits=3)
      ,sprintf(" Q: "),round(quantile(scen.stat.pval[,j],prob=c(0.05,0.1,0.5,0.90,0.95)),digits=3),"\n")
}
cat("\n")
cat("Number of outlying SumStats over a total of",N.stat,"summary statistics","\n")
for (j in 1:N.scen) 
{
  SUM.outliers.1p1000[j] <- sum(scen.stat.pval[,j] <= 0.001) 
  SUM.outliers.1p100[j] <- sum(scen.stat.pval[,j] <= 0.01)
  SUM.outliers.5p100[j] <- sum(scen.stat.pval[,j] <= 0.05)
  SUM.outliers.95p100[j] <- sum(scen.stat.pval[,j] >= 0.95) 
  SUM.outliers.99p100[j] <- sum(scen.stat.pval[,j] >= 0.99)
  SUM.outliers.999p1000[j] <- sum(scen.stat.pval[,j] >= 0.999)
  cat("Scenario #",j,"\n")
  cat("SUM.outliers.1p1000 = ",SUM.outliers.1p1000[j]," ( Expectation = ",round(N.stat*1/1000),")","\n")
  cat("SUM.outliers.1p100 = ",SUM.outliers.1p100[j]," ( Expectation = ",round(N.stat*1/100),")","\n")
  cat("SUM.outliers.5p100 = ",SUM.outliers.5p100[j]," ( Expectation = ",round(N.stat*5/100),")","\n")
  cat("SUM.outliers.95p100 = ",SUM.outliers.95p100[j]," ( Expectation = ",round(N.stat*5/100),")","\n")
  cat("SUM.outliers.99p100 = ",SUM.outliers.99p100[j]," ( Expectation = ",round(N.stat*1/100),")","\n")
  cat("SUM.outliers.999p1000 = ",SUM.outliers.999p1000[j]," (Expectation = ",round(N.stat*1/1000),")","\n")
  cat("\n")
}

cat("#####################################################################,\n")
cat("############ SOME USEFULL INSIGHTS FOR MODEL IMPROVEMENT ############","\n")
cat("#####################################################################,\n")
cat("\n")
for (j in 1:N.scen) 
{
  cat("Scenario #",j,"\n")
  cat("\n")
  for (i in 1: N.stat)
    if (scen.stat.pval[i,j]<=0.05) cat("Posterior.Outlier.5%",stat.names[i],"OBS-value = ",obs.poi[1,i],"POSTERIOR P-value=",scen.stat.pval[i,j],"PRIOR P-value=",PRIOR.scen.stat.pval[i,j],"[P_Posterior - P_Prior] =",scen.stat.pval[i,j]-PRIOR.scen.stat.pval[i,j],"\n")
  cat("\n")
  for (i in 1: N.stat)
    if (scen.stat.pval[i,j]>=0.95) cat("Posterior.Outlier.95%",stat.names[i],"OBS-value = ",obs.poi[1,i],"POSTERIOR P-value=",scen.stat.pval[i,j],"PRIOR P-value=",PRIOR.scen.stat.pval[i,j],"[P_Posterior - P_Prior] =",scen.stat.pval[i,j]-PRIOR.scen.stat.pval[i,j],"\n")
  cat("\n")
  for (i in 1: N.stat)
    if (scen.stat.pval[i,j]<=0.01) cat("Posterior.Outlier.1%",stat.names[i],"OBS-value = ",obs.poi[1,i],"POSTERIOR P-value=",scen.stat.pval[i,j],"PRIOR P-value=",PRIOR.scen.stat.pval[i,j],"[P_Posterior - P_Prior] =",scen.stat.pval[i,j]-PRIOR.scen.stat.pval[i,j],"\n")
  cat("\n")
  for (i in 1: N.stat)
    if (scen.stat.pval[i,j]>=0.99) cat("Posterior.Outlier.99%",stat.names[i],"OBS-value = ",obs.poi[1,i],"POSTERIOR P-value=",scen.stat.pval[i,j],"PRIOR P-value=",PRIOR.scen.stat.pval[i,j],"[P_Posterior - P_Prior] =",scen.stat.pval[i,j]-PRIOR.scen.stat.pval[i,j],"\n")
  cat("\n")
}
sink()


