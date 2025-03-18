wd="/home/abert/media1/projects/drosophila-evolution/data/vcf"

eval "$(conda shell.bash hook)"

conda activate gatk_modified

bcftools filter -O z \
    -o $wd/wholegenome.noinv.vcf.gz \
    -t ^"2L:1725744-13654180,2R:14891154-20776334,3L:2673046-16808841,3R:20906917-29531297,3R:15932209-25244010,3R:19596867-32579331,3R:11250567-26640370,X:13125736-20079328,X:17328912-20093711" \
    $wd/wholegenome.vcf.gz \
    --threads 120
    
conda deactivate

