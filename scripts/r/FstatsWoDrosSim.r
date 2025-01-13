require(poolfstat)

all.pops.data <- vcf2pooldata(vcf.file="/gatk_modified/userdata/abertelli/drosophila-evolution/results/drosophila_evolution.bcftools_fakepools_wholegen.vcf.gz", min.maf="0.05", poolsizes=rep(100000,20))

save(all.pops.data, file="/gatk_modified/userdata/abertelli/drosophila-evolution/results/all_pops_data.RData")
