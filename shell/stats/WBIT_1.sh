source activate python_deps
bcftools stats /gatk_modified/userdata/abertelli/drosophila-evolution/results/groups/WB/WBIT_1/WBIT_1_all.bcftools.vcf.gz > /gatk_modified/userdata/abertelli/drosophila-evolution/results/groups/WB/WBIT_1/WBIT_1_all.vchk
plot-vcfstats -p /gatk_modified/userdata/abertelli/drosophila-evolution/results/groups/WB/WBIT_1/plots /gatk_modified/userdata/abertelli/drosophila-evolution/results/groups/WB/WBIT_1/WBIT_1_all.vchk
conda deactivate
