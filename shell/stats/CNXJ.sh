source activate python_deps
bcftools stats /gatk_modified/userdata/abertelli/drosophila-evolution/results/groups/CNXJ/CNXJ_all.bcftools.vcf.gz > /gatk_modified/userdata/abertelli/drosophila-evolution/results/groups/CNXJ/CNXJ_all.vchk
plot-vcfstats -p /gatk_modified/userdata/abertelli/drosophila-evolution/results/groups/CNXJ/plots /gatk_modified/userdata/abertelli/drosophila-evolution/results/groups/CNXJ/CNXJ_all.vchk
conda deactivate
