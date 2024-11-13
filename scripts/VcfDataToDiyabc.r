library(poolfstat)

poolvec <- c(392, 80, 80, 80)

dest4pop.readcount <-vcf2pooldata(vcf.file="/gatk_modified/drosophila-project/data/pop_data/test_diyabc/DEST.four_pops.neutral.vcf.gz",poolsizes=poolvec,min.maf="0.01")

# [1] OUTPUT: 
# Reading Header lines
# VarScan like format detected for allele count data:
#  the AD field contains allele depth
# for the alternate allele and RD field for the reference allele
# (N.B., positions with more than one alternate allele will be ignored)
# Parsing allele counts
# 1e+05  lines processed in 0 h  0 m  2 s : 45921 SNPs found
# Data consists of 45921 SNPs for 4 Pools

pooldata2diyabc(
    dest4pop.readcount,
    writing.dir = "/gatk_modified/drosophila-project/data/pop_data/test_diyabc/",
    prefix = "DEST.four_pops.neutral",
    diyabc.mrc = 1,
)

