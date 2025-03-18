wd="/home/abert/media1/projects/drosophila-evolution/data/vcf"

eval "$(conda shell.bash hook)"

conda activate gatk_modified

bcftools filter -O z \
    -o $wd/wholegenome.noinv.nohet.vcf.gz \
    -t ^"2L:0-5000000,2R:0-5000000,3L:0-5000000,3R:0-5000000,2L:18513712-23513712,2R:20286936-25286936,3L:23110227-28110227,3R:27079331-32079331,X:0-5000000,X:18542271-23542271" \
    $wd/wholegenome.noinv.vcf.gz \
    --threads 120
    
conda deactivate



