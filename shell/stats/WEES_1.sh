source activate python_deps
bcftools stats /gatk_modified/userdata/abertelli/drosophila-evolution/results/groups/WE/WEES_1/WEES_1_all.bcftools.vcf.gz > /gatk_modified/userdata/abertelli/drosophila-evolution/results/groups/WE/WEES_1/WEES_1_all.vchk
plot-vcfstats -p /gatk_modified/userdata/abertelli/drosophila-evolution/results/groups/WE/WEES_1/plots /gatk_modified/userdata/abertelli/drosophila-evolution/results/groups/WE/WEES_1/WEES_1_all.vchk
conda deactivate
