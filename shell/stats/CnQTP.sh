source activate python_deps
bcftools stats /gatk_modified/userdata/abertelli/drosophila-evolution/results/groups/CnQTP/CnQTP_all.bcftools.vcf.gz > /gatk_modified/userdata/abertelli/drosophila-evolution/results/groups/CnQTP/CnQTP_all.vchk
plot-vcfstats -p /gatk_modified/userdata/abertelli/drosophila-evolution/results/groups/CnQTP/plots /gatk_modified/userdata/abertelli/drosophila-evolution/results/groups/CnQTP/CnQTP_all.vchk
conda deactivate
