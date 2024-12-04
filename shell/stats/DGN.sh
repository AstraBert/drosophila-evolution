source activate python_deps
bcftools stats /gatk_modified/userdata/abertelli/drosophila-evolution/results/groups/DGN/DGN_all.bcftools.vcf.gz > /gatk_modified/userdata/abertelli/drosophila-evolution/results/groups/DGN/DGN_all.vchk
plot-vcfstats -p /gatk_modified/userdata/abertelli/drosophila-evolution/results/groups/DGN/plots /gatk_modified/userdata/abertelli/drosophila-evolution/results/groups/DGN/DGN_all.vchk
conda deactivate
