# Unraveling the shared evolutionary history of European humans and Drosophila melanogaster
## Bachelor's thesis project by Astra Bertelli

**Supervisors**: Martin Kapun, Lino Ometto

**Collaborators**: Alan Bergland, Arnoud Estoup, Stefan Eichart, Mathieu Gautier, Joaquin Nunez

This repository collects the data, the code and the analysis workflow for the project.

## Environment

To reproduce the analysis, we strongly advise that you use the Docker image that we purposedly built for this project. You should get Docker for your platform following [this tutorial](https://dev.to/astrabert/1mindocker-2-get-docker-kh), adjust the path to the data you want to mount inside your Docker image in [`env`](./.env) under `ÙSERDATA_PATH` and then simply run:

```bash
docker compose up -d
docker exec -it $(sudo docker ps -qf "name=drosophila_project_env") /bin/bash
```

Or, if you are on Linux/macOS:

```bash
bash compose.sh
```

You'll be put inside the container (a semi-isolated virtual machine) in which our Docker image is up and running, and you'll find all your mounted data under `/gatk_modified/drosophila-project`.

If you want to check new releases for the image as well as past versions, you can do that on [Docker Hub](https://hub.docker.com/repository/docker/astrabert/silly-gat-kay/general) or on the [GitHub repository](https://github.com/AstraBert/silly-gat-kay).

## Analysis workflow

Find all the steps of the workflow in [MAIN.sh](./shell/MAIN.sh).

**Getting the data**: We stored all the download links or SRA accessions in specific files that we subsequently used to get the data either through `prefetch`+`fasterq-dump` or through `wget`. 

```bash
wd=/gatk_modified/drosophila-project/data

## GET FASTQS FOR CHINESE AND ISRAELI SAMPLES
for sg in CN_XJ CnOther CnQTP ISR_1 ISR_2
do
    outp=${wd}/data_files/${sg}
    
    mkdir -p $outp
    
    cd $outp

    input_file=${wd}/download_${sg}.txt 

    while IFS= read -r accession; do
        echo "Downloading and processing: $accession"
        
        prefetch "$accession"
        
        fasterq-dump -e 50 --split-files "$accession"
        pigz -p 50 ${wd}/data_files/${sg}/${accession}_1.fastq
        pigz -p 50 ${wd}/data_files/${sg}/${accession}_2.fastq
        rm -rf ${wd}/data_files/${sg}/${accession} 
    done < "$input_file"
done

cp $wd/*/*/*.dedup.bam* $wd/

cd $wd

## DOWNLOAD FOUR EASTERN EUROPEAN AND FOUR WESTERN EUROPEAN SAMPLES FROM DESTv2
while IFS= read -r url; do
    wget "$url" 
done < bam_download_link.txt
```

**Preprocessing the data**: We preprocessed the data modifying the naming in order for them to match with the population from which they're coming from. 

```bash
## RENAME BAM FILES IN ORDER TO MATCH WITH THEIR ORIGINAL POPULATION
conda activate python_deps

python3 scripts/RenameBamFiles.py
python3 scripts/RenameDestSamples.py

conda deactivate

## RENAME AND MOVE DGN DATA
bash shell/re_african_data.sh
rm -rf data/DGN/
```

And then we moved all the BAM files together:

```bash
mkdir -p data/bamfiles && mv data/*.bam* data/bamfiles
```

**Preparing the input for FreeBayes**: We got the reference genome and indexed it:

```bash
source activate gatk_modified

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
```

Then we extracted all the 250.000bp-long chunks from the reference genome with and AWK script, chunks that we divided into autosomal and X-linked (excluding chromosome 4 and chromosome Y).

```bash
mkdir data/freebayes_inputs/

### 1. Find all 250000bp chunks
fastalength -f data/reference/dmel-6.59.fa | awk -f scripts/generate_chunk.awk -v chunk_size=250000 - > data/freebayes_inputs/all.chunks

conda deactivate

### 2. Divide them into X and autosomal chunks
source activate python_deps

python3 scripts/AutoVsSexChunks.py

conda deactivate
```

We then created a file with the list of all the BAM files we have:

```bash
### 3. Create a file listing all the BAM files
for f in /gatk_modified/drosophila-project/data/bamfiles/*.bam
do
    echo $f >> data/freebayes_inputs/bamfiles.txt
done
```

**Variant calling with FreeBayes**: TBD