require(poolfstat)

load("results/all_pops_data_wosim.RData") # -> all.pops.data

dros.sim.df <- read.csv("results/drossim_snps.csv")
dros.sim.column <- dros.sim.df[,1]
all.pops.count <- cbind(all.pops.data@refallele.readcount)
colnames(all.pops.count) <- NULL
all.pops.data@refallele.readcount <- all.pops.count
all.pops.data@poolsizes <- rep(1000,21)  
all.pops.data@poolnames <- c(all.pops.data@poolnames, "Pool21")
all.pops.data@npools <- 21
all.pops.cov <- all.pops.data@readcoverage
all.pops.cov <- cbin(all.pops.cov, dros.sim.column)
colnames(all.pops.cov) <- NULL
all.pops.data@readcoverage <- all.pops.cov
all.pops.fstats <- compute.fstats(all.pops.data, nsnp.per.bjack.block = 1000, computeDstat = TRUE,verbose=TRUE)
save(all.pops.fstats, file="results/all_pops_fstats_wsim.RData")    