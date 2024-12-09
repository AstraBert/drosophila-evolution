# Unraveling the shared evolutionary history of European humans and Drosophila melanogaster
## Bachelor's thesis project by Astra Bertelli

**Supervisors**: Martin Kapun, Lino Ometto

**Collaborators**: Alan Bergland, Arnoud Estoup, Stefan Eichert, Mathieu Gautier, Joaquin Nunez

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

**Set your working directory**: Set the working repository so that it matches your local file system:

```bash
wd="/gatk_modified/userdata/abertelli/drosophila-evolution/"
```

### Data getting workflow

<div align="center">
    <img src="./imgs/data_getting.jpg" alt="Pipeline for getting data">
    <p><i>Workflow for getting needed data</i></p>
    <sbr>
</div>

**Getting the reference genome and index it**: Download the latest release of _Drosophila melanogaster_'s genome from [FlyBase]() and index it with `samtools`+`bwa-mem2`:

```bash
## GET THE REFERENCE GENOME FOR D. melanogaster FROM FlyBase
mkdir -p $wd/data/reference/

cd $wd/data/reference/

wget -O dmel-6.59.fa.gz https://ftp.flybase.net/genomes/Drosophila_melanogaster/current/fasta/dmel-all-chromosome-r6.59.fasta.gz

## INDEX THE REFERENCE GENOME
bwa-mem2 $wd/data/reference/index dmel-6.59.fa.gz
gunzip -c $wd/data/reference/dmel-6.59.fa.gz >$wd/data/reference/dmel-6.59.fa
samtools faidx $wd/data/reference/dmel-6.59.fa
samtools dict $wd/data/reference/dmel-6.59.fa >$wd/data/reference/dmel-6.59.dict
```

**Get already mapped data files**: The DGN data were available from local storage. DESTv2 data can be instead downloaded from the database using `wget`:

```bash
## DOWNLOAD FOUR EASTERN EUROPEAN AND FOUR WESTERN EUROPEAN SAMPLES FROM DESTv2
cd $wd/data/bamfiles/

while IFS= read -r url; do
    wget "$url" 
done < $wd/data/bam_download_link.txt

## DOWNLOAD ADDITIONAL EASTERN EUROPEAN AND WESTERN EUROPEAN SAMPLES FROM DESTv1 AND DESTv2
cnt=0
mkdir $wd/shell/wget_data/

while IFS= read -r url; do
    ((cnt++))
    echo "wget "$url"" > $wd/shell/wget_data/${cnt}.sh 
done < $wd/data/additional_dest_bams.txt

source activate freebayes-env
echo "will cite" | parallel --citation >/dev/null 2>&1
parallel --bar -j 100 bash ::: $wd/shell/wget_data/*.sh
conda deactivate
```

### Mapping data pipeline

<div align="center">
    <img src="./imgs/data_mapping.jpg" alt="Pipeline for mapping data">
    <p><i>Raw sequencing data mapping pipeline</i></p>
    <sbr>
</div>

**Get the raw sequencing data and map them**: Store all the SRA accessions in specific files and subsequently use them to get the data either through `prefetch`+`fasterq-dump`. Use the [mapping_pipeline.sh](./shell/mapping_pipeline.sh) script to map them on the fly. 

```bash
for sg in CN_XJ CnOther CnQTP ISR_1 ISR_2
do
    outp=${wd}/data_files/${sg}
    
    mkdir -p $outp
    
    cd $outp

    input_file=${wd}/download_${sg}.txt 

    while IFS= read -r accession; do
        echo "Downloading and processing: $accession"
        
        # Download the SRA file
        prefetch "$accession"
        
        mkdir -p ${wd}/data/mapping/${sg}/${accession}

        # Convert the SRA file to FASTQ format
        fasterq-dump -e 50 --split-files "$accession"
        pigz -p 50 ${wd}/data_files/${sg}/${accession}_1.fastq
        pigz -p 50 ${wd}/data_files/${sg}/${accession}_2.fastq
        rm -rf ${wd}/data_files/${sg}/${accession} 
        bash $wd/shell/mapping_pipeline.sh \
            -fq1 ${wd}/data_files/${sg}/${accession}_1.fastq \
            -fq2 ${wd}/data_files/${sg}/${accession}_2.fastq \
            -r ${wd}/data/reference/dmel-6.59.fa.gz \
            -o ${wd}/data/mapping/${sg}/${accession} \
            -t 50
    done < "$input_file"
done

## RETAIN ONLY THE BAM FILES AND MOVE THEM
rm -rf ${wd}/data_files/
mkdir ${wd}/data/bamfiles
mv ${wd}/data/mapping/*/*.bam* ${wd}/data/bamfiles
```

### Processing the data

**Preprocessing the data**: We preprocessed the data modifying the naming in order for them to match with the population from which they're coming from. 

```bash
conda activate python_deps

## RENAME BAM FILES IN ORDER TO MATCH WITH THEIR ORIGINAL POPULATION
python3 $wd/scripts/RenameBamFiles.py

## RENAME SAMPLES FROM DESTv2 ACCORDING TO WHAT DID BEFORE
python3 $wd/scripts/RenameDestSamples.py

## RENAME ADDITIONAL SAMPLES FROM DESTv1/DESTv2 ACCORDING TO WHAT DID BEFORE
python3 $wd/scripts/RenameAdditionalDest.py

conda deactivate

counter=0
for f in $wd/data/DGN/*.bam
do
    ((counter++))
    mv $f $wd/data/bamfiles/DGN_${counter}.bam
done
```

**Adjusting the RG line in DGN BAM files and indexing them**: The `@RG` line in the DGN BAM files either reports `sample` as `SM` tag. In order to avoid confusion in the variant calling process, substitute the `@RG` line using `samtools` and custom python scripts:

```bash
for n in {1..165} 
do
    f="$wd/data/bamfiles/DGN_${n}.bam"
    fres="$wd/data/dgn_renamed/DGN_${n}.bam"

    source activate python_deps

    rgline=$(python3 $wd/scripts/ExtractRgLineText.py -i $f)
    rgid=$(python3 $wd/scripts/ExtractIdFromRgLine.py -rgl "$rgline")

    conda deactivate

    source activate gatk_modified

    samtools addreplacerg --threads 100 -r "$rgline" -w -o $fres $f
    samtools index -@ 100 $fres

    conda deactivate
done

## REMOVE THE UNUSED DIRECTORIES
rm -rf $wd/data/bamfiles/DGN_*.bam
mv $wd/data/dgn_renamed/*.bam* $wd/data/bamfiles/
rm -rf $wd/data/dgn_renamed/
```

**Adjusting the RG line in DEST BAM files and indexing them**: The `@RG` line in the DEST BAM files reports `sample_name` as `SM` tag. In order to avoid confusion in the variant calling process, substitute the `@RG` line using `samtools` and custom python scripts:

```bash
mkdir -p $wd/data/dest_renamed/

for sample in EERU EEUA WEES WEFR
do
    for n in 1 2
    do 
        f="$wd/data/bamfiles/${sample}_${n}.bam"
        fres="$wd/data/dest_renamed/${sample}_${n}.bam"

        source activate python_deps

        rgline=$(python3 $wd/scripts/RgLineForDest.py -i $f)

        conda deactivate

        source activate gatk_modified

        samtools addreplacerg --threads 100 -r "$rgline" -w -o $fres $f
        samtools index -@ 100 $fres

        conda deactivate
    done
done

mkdir -p $wd/data/add_dest_renamed/


for sample in EBAT WBIT EBPL WBDE EBHU CYP TRK
do
    for n in 1
    do 
        f="$wd/data/bamfiles/${sample}_${n}.bam"
        fres="$wd/data/add_dest_renamed/${sample}_${n}.bam"

        source activate python_deps

        rgline=$(python3 $wd/scripts/RgLineForDest.py -i $f)

        conda deactivate

        echo $rgline
        source activate gatk_modified

        samtools addreplacerg --threads 100 -r "$rgline" -w -o $fres $f
        samtools index -@ 100 $fres

        conda deactivate
    done
done

rm -rf $wd/data/bamfiles/EE??_?.bam*
rm -rf $wd/data/bamfiles/WE??_?.bam*
rm -rf $wd/data/bamfiles/WB??_?.bam*
rm -rf $wd/data/bamfiles/EB??_?.bam*
rm -rf $wd/data/bamfiles/CYP_?.bam*
rm -rf $wd/data/bamfiles/TRK_?.bam*
mv $wd/data/dest_renamed/*.bam* $wd/data/bamfiles/
mv $wd/data/add_dest_renamed/*.bam* $wd/data/bamfiles/
rm -rf $wd/data/dest_renamed/
rm -rf $wd/data/add_dest_renamed/
``` 


**Preparing the input for FreeBayes and BCFTools**: Create a file with the list of all the BAM files we have:

```bash
for f in $wd/data/bamfiles/*.bam
do
    echo $f >> $wd/data/freebayes_inputs/bamfiles.txt
done
```

### Variant calling

**Variant calling with FreeBayes**: Use multi-threaded FreeBayes to call the variants from BAM files, storing them into a gzipped VCF file. 

```bash
source activate freebayes-env

export TMPDIR=/gatk_modified/userdata/tmp/
freebayes-parallel \
    <(fasta_generate_regions.py \
        $wd/data/reference/dmel-6.59.fa.fai \
        100000) \
    80 \
    -f $wd/data/reference/dmel-6.59.fa \
    -L $wd/data/freebayes_inputs/bamfiles.txt \
    -C 1 \
    -F 0.02 \
    -G 5 \
    --limit-coverage 250 \
    --use-best-n-alleles 4 \
    --strict-vcf \
    --pooled-continuous |
    gzip >$wd/results/drosophila_evolution.freebayes.vcf.gz

conda deactivate
```

**Variant calling with BCFtools**: Use `bcftools` to perform variant calling on the data (useful for comparison between the two methods, but also faster then FreeBayes and so easier to get the data). Parallelize the variant calling duividing it on different chromosomes and spawn it on different threads with `GNU parallel`:

```bash
source activate freebayes-env
echo "will cite" | parallel --citation >/dev/null 2>&1
parallel --bar -j 40 bash ::: $wd/shell/bcftools_regions/*.sh
conda deactivate
```

This is a code snippet for the variant calling command:

```bash
wd=/gatk_modified/userdata/abertelli/drosophila-evolution

source activate gatk_modified
bcftools mpileup \
    -C 50 \
    -Ou \
    -f $wd/data/reference/dmel-6.59.fa \
    -b $wd/data/freebayes_inputs/bamfiles.txt \
    -q 20 \
    -Q 20 \
    -r "2L" \
    -a DP,AD | bcftools call -mv -Oz --format-fields GQ,GP > $wd/results/drosophila_evolution.bcftools_2L.vcf.gz
conda deactivate
```

### VCF Post-Processing

From `bcftools` you will obtain a set of five VCF files, each representing one of the chromosomes of interest (2R, 2L, 3R, 3L and X).

**Concatenate VCF files**: You can use `bcftools concat` to concatenate the five VCF files obtained into one:

```bash
bcftools concat \
    -O z \
    --threads 16 \
    -o $wd/results/drosophila_evolution.bcftools_all.vcf.gz \
    $wd/results/drosophila_evolution.bcftools_2R.vcf.gz \
    $wd/results/drosophila_evolution.bcftools_2L.vcf.gz \
    $wd/results/drosophila_evolution.bcftools_3R.vcf.gz \
    $wd/results/drosophila_evolution.bcftools_3L.vcf.gz \
    $wd/results/drosophila_evolution.bcftools_X.vcf.gz  
```

**Calculate simple VCF stats**: You can calculate simple VCF stats with `bcftools stats` and plot them with `plot-vcfstats`.

```bash
## EXTRACT STATS FROM VCF FILES

source activate python_deps
bcftools stats $wd/results/drosophila_evolution.bcftools_all.vcf.gz > $wd/results/drosophila_evolution.bcftools_all.vchk
plot-vcfstats -p $wd/results/bcftools_plots/ $wd/results/drosophila_evolution.bcftools_all.vchk
conda deactivate
```