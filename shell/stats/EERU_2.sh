source activate python_deps
bcftools stats /gatk_modified/userdata/abertelli/drosophila-evolution/results/groups/EE/EERU_2/EERU_2_all.bcftools.vcf.gz > /gatk_modified/userdata/abertelli/drosophila-evolution/results/groups/EE/EERU_2/EERU_2_all.vchk
plot-vcfstats -p /gatk_modified/userdata/abertelli/drosophila-evolution/results/groups/EE/EERU_2/plots /gatk_modified/userdata/abertelli/drosophila-evolution/results/groups/EE/EERU_2/EERU_2_all.vchk
conda deactivate
