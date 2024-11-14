#!/bin/bash

wd=/gatk_modified/drosophila-project/data/pop_data

for sg in ISR_1 ISR_2
do
    input_file=${wd}/download_${sg}.txt 

    while IFS= read -r accession; do
        outp=${wd}/data/${sg}/${accession} 
        mkdir -p $outp/logs
        
        bash /gatk_modified/drosophila-project/data/pop_data/shell/mapping_pipeline.sh \
            -r /gatk_modified/drosophila-project/data/drosophila_starting_test/dmel-6.59.fa.gz \
            -fq1 ${wd}/data/${sg}/${accession}_1.fastq.gz \
            -fq2 ${wd}/data/${sg}/${accession}_2.fastq.gz \
            -o $outp \
            -t 32
    done < "$input_file"
done