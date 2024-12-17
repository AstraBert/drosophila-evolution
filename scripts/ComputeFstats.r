all.pops <- vcf2pooldata(vcf.file="/gatk_modified/userdata/abertelli/drosophila-evolution/results/fake_pools.vcf", min.maf="0.05", poolsizes = rep(100000,20))
all.pops.fstats<-compute.fstats(all.pops,nsnp.per.bjack.block = 1000, computeDstat = TRUE,verbose=FALSE)
all.pops.fst<-computeFST(all.pops, nsnp.per.bjack.block = 100, verbose=FALSE)
