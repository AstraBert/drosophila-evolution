source activate python_deps
bcftools stats /gatk_modified/userdata/abertelli/drosophila-evolution/results/groups/CnOther/CnOther_all.bcftools.vcf.gz > /gatk_modified/userdata/abertelli/drosophila-evolution/results/groups/CnOther/CnOther_all.vchk
plot-vcfstats -p /gatk_modified/userdata/abertelli/drosophila-evolution/results/groups/CnOther/plots /gatk_modified/userdata/abertelli/drosophila-evolution/results/groups/CnOther/CnOther_all.vchk
conda deactivate
