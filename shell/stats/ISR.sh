source activate python_deps
bcftools stats /gatk_modified/userdata/abertelli/drosophila-evolution/results/groups/ISR/ISR_all.bcftools.vcf.gz > /gatk_modified/userdata/abertelli/drosophila-evolution/results/groups/ISR/ISR_all.vchk
plot-vcfstats -p /gatk_modified/userdata/abertelli/drosophila-evolution/results/groups/ISR/plots /gatk_modified/userdata/abertelli/drosophila-evolution/results/groups/ISR/ISR_all.vchk
conda deactivate
