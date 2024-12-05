library("poolfstat")

allpops.readcount <-vcf2pooldata(vcf.file="/gatk_modified/userdata/abertelli/drosophila-evolution/results/drosophila_evolution.bcftools_2R.vcf.gz",poolsizes=rep(80,335),min.maf="0.01")
