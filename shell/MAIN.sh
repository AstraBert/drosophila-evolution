#!/bin/bash

source activate gatk_modified

bash shell/download_pop_fastqs.sh

cp data/*/*/*.dedup.bam* data/

mkdir -p data/reference/

cd data/reference/

## GET THE REFERENCE DNA FOR D. melanogaster FROM FlyBase
wget -O dmel-6.59.fa.gz https://ftp.flybase.net/genomes/Drosophila_melanogaster/current/fasta/dmel-all-chromosome-r6.59.fasta.gz

## INDEX
bwa-mem2 index dmel-6.59.fa.gz
gunzip -c dmel-6.59.fa.gz >dmel-6.59.fa
samtools faidx dmel-6.59.fa
samtools dict dmel-6.59.fa >dmel-6.59.dict

cd ../../

conda deactivate

conda activate python_deps

python3 scripts/RenameBamFiles.py

conda deactivate

cd data/

while IFS= read -r url; do
    wget "$url" 
done < bam_download_link.txt

source activate python_deps

python3 scripts/RenameDestSamples.py

conda deactivate

bash shell/re_african_data.sh

rm -rf data/DGN/

mkdir -p data/bamfiles && mv data/*.bam* data/bamfiles

conda activate gatk_modified

mkdir data/freebayes_inputs/

fastalength -f data/reference/dmel-6.59.fa | awk -f scripts/generate_chunk.awk -v chunk_size=250000 - > data/freebayes_inputs/all.chunks

conda deactivate

source activate python_deps

python3 scripts/AutoVsSexChunks.py

conda deactivate

for f in /gatk_modified/drosophila-project/data/bamfiles/*.bam
do
    echo $f >> data/freebayes_inputs/bamfiles.txt
done

