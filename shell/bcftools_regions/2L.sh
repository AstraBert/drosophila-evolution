eval "$(conda shell.bash hook)"

wd=/home/abert/media1/projects/drosophila-evolution

source activate gatk_modified

bcftools mpileup \
    -Ou \
    -f $wd/data/reference/dmel-6.62.fa \
    -b $wd/data/bamfiles/bamfiles.txt \
    -q 20 \
    -Q 20 \
    -r "2L" \
    --threads 120 \
    -a DP,AD | bcftools call -mv -Oz  --format-fields GQ,GP > $wd/data/vcf/2L.vcf.gz

conda deactivate
