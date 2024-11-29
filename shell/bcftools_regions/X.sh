wd=/gatk_modified/userdata/abertelli/drosophila-evolution

source activate gatk_modified

bcftools mpileup \
    -C 50 \
    -Ou \
    -f $wd/data/reference/dmel-6.59.fa \
    -b $wd/data/freebayes_inputs/bamfiles.txt \
    -q 20 \
    -Q 20 \
    -r "X" \
    -a DP,AD | bcftools call -mv -Oz  --format-fields GQ,GP > $wd/results/drosophila_evolution.bcftools_X.vcf.gz

conda deactivate