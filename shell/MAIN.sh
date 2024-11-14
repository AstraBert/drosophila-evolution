#!/bin/bash

## GET FASTQS FOR CHINESE AND ISRAELI SAMPLES
source activate gatk_modified

bash shell/download_pop_fastqs.sh

cp data/*/*/*.dedup.bam* data/

## GET THE REFERENCE GENOME FOR D. melanogaster FROM FlyBase
mkdir -p data/reference/

cd data/reference/

wget -O dmel-6.59.fa.gz https://ftp.flybase.net/genomes/Drosophila_melanogaster/current/fasta/dmel-all-chromosome-r6.59.fasta.gz

## INDEX THE REFERENCE GENOME
bwa-mem2 index dmel-6.59.fa.gz
gunzip -c dmel-6.59.fa.gz >dmel-6.59.fa
samtools faidx dmel-6.59.fa
samtools dict dmel-6.59.fa >dmel-6.59.dict

cd ../../

conda deactivate

## RENAME BAM FILES IN ORDER TO MATCH WITH THEIR ORIGINAL POPULATION
conda activate python_deps

python3 scripts/RenameBamFiles.py

conda deactivate

## DOWNLOAD FOUR EASTERN EUROPEAN AND FOUR WESTERN EUROPEAN SAMPLES FROM DESTv2
cd data/

while IFS= read -r url; do
    wget "$url" 
done < bam_download_link.txt

## RENAME SAMPLES FROM DESTv2 ACCORDING TO WHAT DID BEFORE
source activate python_deps

python3 scripts/RenameDestSamples.py

conda deactivate

## RENAME AND MOVE DGN DATA
bash shell/re_african_data.sh

rm -rf data/DGN/

## MOVE ALL THE BAM FILES TOGETHER
mkdir -p data/bamfiles && mv data/*.bam* data/bamfiles

conda activate gatk_modified

## CREATE THE INPUTS FOR FreeBayes
mkdir data/freebayes_inputs/

### 1. Find all 250000bp chunks
fastalength -f data/reference/dmel-6.59.fa | awk -f scripts/generate_chunk.awk -v chunk_size=250000 - > data/freebayes_inputs/all.chunks

conda deactivate

### 2. Divide them into X and autosomal chunks
source activate python_deps

python3 scripts/AutoVsSexChunks.py

conda deactivate

### 3. Create a file listing all the BAM files
for f in /gatk_modified/drosophila-project/data/bamfiles/*.bam
do
    echo $f >> data/freebayes_inputs/bamfiles.txt
done

source activate freebayes-env

freebayes-parallel \
    <(fasta_generate_regions.py \
        /gatk_modified/drosophila-project/data/reference/dmel-6.59.fa.fai \
        100000) \
    100 \
    -f /gatk_modified/drosophila-project/data/reference/dmel-6.59.fa \
    -L /gatk_modified/drosophila-project/data/freebayes_inputs/bamfiles.txt \
    --pooled-continuous |
    gzip >/gatk_modified/drosophila-project/results/drosophila_evolution.vcf.gz

conda deactivate