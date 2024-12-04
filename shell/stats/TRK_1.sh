source activate python_deps
bcftools stats /gatk_modified/userdata/abertelli/drosophila-evolution/results/groups/TRK/TRK_1/TRK_1_all.bcftools.vcf.gz > /gatk_modified/userdata/abertelli/drosophila-evolution/results/groups/TRK/TRK_1/TRK_1_all.vchk
plot-vcfstats -p /gatk_modified/userdata/abertelli/drosophila-evolution/results/groups/TRK/TRK_1/plots /gatk_modified/userdata/abertelli/drosophila-evolution/results/groups/TRK/TRK_1/TRK_1_all.vchk
conda deactivate
