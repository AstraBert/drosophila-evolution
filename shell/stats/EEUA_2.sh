source activate python_deps
bcftools stats /gatk_modified/userdata/abertelli/drosophila-evolution/results/groups/EE/EEUA_2/EEUA_2_all.bcftools.vcf.gz > /gatk_modified/userdata/abertelli/drosophila-evolution/results/groups/EE/EEUA_2/EEUA_2_all.vchk
plot-vcfstats -p /gatk_modified/userdata/abertelli/drosophila-evolution/results/groups/EE/EEUA_2/plots /gatk_modified/userdata/abertelli/drosophila-evolution/results/groups/EE/EEUA_2/EEUA_2_all.vchk
conda deactivate
