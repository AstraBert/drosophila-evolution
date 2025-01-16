require(poolfstat)
# Loading required package: poolfstat
drosevol.df <- read.csv("drosevol_readcount.csv")
drosevol.dgn <- drosevol.df[,1]
drosevol.cnxj <- drosevol.df[,2]
drosevol.cnother <- drosevol.df[,3]
drosevol.cnqtp <- drosevol.df[,4]
drosevol.isr <- drosevol.df[,5]
drosevol.drossim <- drosevol.df[,6]
load("eu_pops_data.RData")
eu.selected.df <- read.csv("selected_dest_eu_snps.csv")
eu.selected.snps <- eu.selected.df[,1]
eu.pops.count <- eu.pops.data@refallele.readcount
dim(eu.pops.count)
# [1] 1217825     347
length(eu.selected.snps)
# [1] 1137891
eu.selected.snps <- eu.selected.snps + 1
eu.pops.count.filtered <- eu.pops.count[eu.selected.snps, ]
dim(eu.pops.count.filtered)  # Should return 1137891, 347
# [1] 1137891     347
eu.dgn <- cbind(eu.pops.count.filtered, drosevol.dgn)
eu.dgn.cnxj <- cbind(eu.pops.count.filtered, drosevol.cnxj)
eu.dgn.cnxj.cnother <- cbind(eu.pops.count.filtered, drosevol.cnother)
eu.dgn <- cbind(eu.pops.count.filtered, drosevol.dgn)
eu.dgn.cnxj <- cbind(eu.dgn, drosevol.cnxj)
eu.dgn.cnxj.cnother <- cbind(eu.dgn.cnxj, drosevol.cnother)
eu.dgn.cnxj.cnother.cnqtp <- cbind(eu.dgn.cnxj.cnother, drosevol.cnqtp)
eu.dgn.cnxj.cnother.cnqtp.isr <- cbind(eu.dgn.cnxj.cnother.cnqtp, drosevol.isr)
eu.dgn.cnxj.cnother.cnqtp.isr.drossim <- cbind(eu.dgn.cnxj.cnother.cnqtp.isr, drosevol.drossim)
eu.pops.count <- eu.dgn.cnxj.cnother.cnqtp.isr.drossim
colnames(eu.pops.count) <- NULL
eu.pops.data@refallele.readcount <- eu.pops.count
eu.pops.cov <- eu.pops.data@readcoverage
eu.pops.cov.filtered <- eu.pops.cov[eu.selected.snps, ]
rceu.dgn <- cbin(eu.pops.cov.filtered, drosevol.dgn)
rceu.dgn <- cbind(eu.pops.cov.filtered, drosevol.dgn)
rceu.dgn.cnxj <- cbind(rceu.dgn, drosevol.cnxj)
rceu.dgn.cnxj.cnother <- cbind(rceu.dgn.cnxj, drosevol.cnother)
rceu.dgn.cnxj.cnother.cnqtp <- cbind(rceu.dgn.cnxj.cnother, drosevol.cnqtp)
rceu.dgn.cnxj.cnother.cnqtp.isr <- cbind(rceu.dgn.cnxj.cnother.cnqtp, drosevol.isr)
rceu.dgn.cnxj.cnother.cnqtp.isr.drossim <- cbind(rceu.dgn.cnxj.cnother.cnqtp.isr, drosevol.drossim)
eu.pops.cov <- rceu.dgn.cnxj.cnother.cnqtp.isr.drossim
eu.pops.data@readcoverage <- eu.pops.cov
eu.pops.poolsizes <- c(eu.pops.data@poolsizes, 40, 25, 50, 50, 32, 100)
length(eu.pops.poolsizes)
# [1] 353
eu.pops.data@poolsizes <- eu.pops.poolsizes
eu.pops.poolnames <- eu.pops.data@poolnames
eu.pops.poolnames <- c(eu.pops.poolnames, "DGN", "CNXJ", "CnOther", "CnQTP", "ISR", "DrosSim")
eu.pops.data@npools <- 353
all.pops.fstata <- compute.fstats(eu.pops.data,  nsnp.per.bjack.block = 10000, computeDstat = FALSE,verbose=TRUE)
# Block-Jackknife specification
# 120 Jackknife blocks identified with 1200000 SNPs (out of 1217825 ).
# SNPs map to 5 different chrom/scaffolds
# Average (min-max) Block Sizes: 1.022 ( 0.49 - 6.97 ) Mb
# Estimating Q1
# 0%   10   20   30   40   50   60   70   80   90   100%
# [----|----|----|----|----|----|----|----|----|----|
# **************************************************|
# Estimating Q2
# 0%   10   20   30   40   50   60   70   80   90   100%
# [----|----|----|----|----|----|----|----|----|----|
# **************************************************|
# Estimating within-population heterozygosities
# --Error in `.rowNamesDF<-`(x, value = value) : invalid 'row.names' length

