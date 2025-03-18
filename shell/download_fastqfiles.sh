wd="/home/abert/media1/projects/drosophila-evolution/data"

for sg in isr
do
    outp=$wd/interim/${sg}
    mkdir -p $outp
    cd $outp
    input_file=$wd/download_india_${sg}.txt
    while IFS= read -r accession; do
        echo "Downloading and processing: $accession"
        
        # Download the SRA file
        prefetch "$accession"
        
        mkdir -p $wd/mapping/${sg}/${accession}

        # Convert the SRA file to FASTQ format
        fasterq-dump -e 120 --split-files "$accession"
        pigz -p 120 $wd/interim/${sg}/${accession}_1.fastq
        pigz -p 120 $wd/interim/${sg}/${accession}_2.fastq
        rm -rf $wd/interim/${sg}/${accession}
        bash $wd/../shell/mapping_pipeline.sh \
            -fq1 $wd/interim/${sg}/${accession}_1.fastq.gz \
            -fq2 $wd/interim/${sg}/${accession}_2.fastq.gz \
            -r $wd/reference/dmel-6.62.fa.gz \
            -o $wd/mapping/${sg}/${accession} \
            -t 120
    done < "$input_file"
done
