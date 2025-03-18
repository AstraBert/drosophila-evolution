wd="/home/abert/media1/projects/drosophila-evolution"

bash $wd/shell/wget_orchestration.sh # get DEST BAM files
bash $wd/shell/download_fastqfiles.sh # get DGN FASTQ files and map them

bash $wd/shell/variant_call.sh  # call variants
cd $wd/data/vcf/

conda activate gatk_modified
bcftools concat -O z -o wholegenome.vcf.gz 2R.vcf.gz 2L.vcf.gz 3R.vcf.gz 3L.vcf.gz --threads 120
conda deactivate

# filtering inversions and heterochromatic regions
bash $wd/shell/filter_inversions.sh
bash $wd/shell/filter_het.sh

# Pseudo Poolification
conda activate python_deps
python3 $wd/scripts/python/PseudoPoolify.py
conda deactivate
