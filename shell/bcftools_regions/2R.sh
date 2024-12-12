wd=/gatk_modified/userdata/abertelli/drosophila-evolution

source activate gatk_modified

bcftools mpileup \
    -C 50 \
    -Ou \
    -f $wd/data/reference/dmel-6.59.fa \
    -b $wd/data/freebayes_inputs/bamfiles_DGN_subs.txt \
    -Q 20 \
    -q 0 \
    -r "2R" \
    -a DP,AD | bcftools call -mv -Oz  --format-fields GQ,GP | gunzip -c | head -n 100000

conda deactivate