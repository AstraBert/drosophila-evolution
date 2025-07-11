#!/bin/bash

wd="/gatk_modified/userdata/abertelli/drosophila-evolution/"

conda activate gatk_modified

## GET THE REFERENCE GENOME FOR D. melanogaster FROM FlyBase
mkdir -p $wd/data/reference/

cd $wd/data/reference/

wget -O dmel-6.59.fa.gz https://ftp.flybase.net/genomes/Drosophila_melanogaster/current/fasta/dmel-all-chromosome-r6.59.fasta.gz

## INDEX THE REFERENCE GENOME
bwa-mem2 $wd/data/reference/index dmel-6.59.fa.gz
gunzip -c $wd/data/reference/dmel-6.59.fa.gz >$wd/data/reference/dmel-6.59.fa
samtools faidx $wd/data/reference/dmel-6.59.fa
samtools dict $wd/data/reference/dmel-6.59.fa >$wd/data/reference/dmel-6.59.dict

## GET FASTQS FOR CHINESE AND ISRAELI SAMPLES
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

conda activate python_deps

python3 $wd/scripts/DgnSelect.py

conda deactivate

for sg in dgn
do
    outp=${wd}/data_files/${sg}
    
    mkdir -p $outp
    
    cd $outp

    input_file=${wd}/data/download_${sg}.txt 

    while IFS= read -r accession; do
        echo "Downloading and processing: $accession"
        
        # Download the SRA file
        prefetch "$accession"
        
        mkdir -p ${wd}/data/mapping/${sg}/${accession}

        # Convert the SRA file to FASTQ format
        fasterq-dump -e 100 --split-files "$accession"
        pigz -p 100 ${wd}/data_files/${sg}/${accession}_1.fastq
        pigz -p 100 ${wd}/data_files/${sg}/${accession}_2.fastq
        rm -rf ${wd}/data_files/${sg}/${accession} 
        bash $wd/shell/mapping_pipeline.sh \
            -fq1 ${wd}/data_files/${sg}/${accession}_1.fastq.gz \
            -fq2 ${wd}/data_files/${sg}/${accession}_2.fastq.gz \
            -r ${wd}/data/reference/dmel-6.59.fa.gz \
            -o ${wd}/data/mapping/${sg}/${accession} \
            -t 100
    done < "$input_file"
done

conda activate python_deps

python3 $wd/scripts/RenameDGN.py

conda deactivate

rm -rf ${wd}/data_files/
mkdir ${wd}/data/bamfiles
mv ${wd}/data/mapping/*/*.bam* ${wd}/data/bamfiles

conda deactivate

## DOWNLOAD FOUR EASTERN EUROPEAN AND FOUR WESTERN EUROPEAN SAMPLES FROM DESTv2
cd $wd/data/bamfiles/

while IFS= read -r url; do
    wget "$url" 
done < $wd/data/bam_download_link.txt

## DOWNLOAD AADITIONAL EASTERN EUROPEAN AND WESTERN EUROPEAN SAMPLES FROM DESTv1 AND DESTv2
cd $wd/data/bamfiles/

cnt=0
mkdir $wd/shell/wget_data/

while IFS= read -r url; do
    ((cnt++))
    echo "wget "$url"" > $wd/shell/wget_data/${cnt}.sh 
done < $wd/data/additional_dest_bams.txt

conda activate freebayes-env
echo "will cite" | parallel --citation >/dev/null 2>&1
parallel --bar -j 100 bash ::: $wd/shell/wget_data/*.sh
conda deactivate

conda activate python_deps

## RENAME BAM FILES IN ORDER TO MATCH WITH THEIR ORIGINAL POPULATION
python3 $wd/scripts/RenameBamFiles.py

## RENAME SAMPLES FROM DESTv2 ACCORDING TO WHAT DID BEFORE
python3 $wd/scripts/RenameDestSamples.py

## RENAME ADDITIONAL SAMPLES FROM DESTv1/DESTv2 ACCORDING TO WHAT DID BEFORE
python3 $wd/scripts/RenameAdditionalDest.py

conda deactivate

mkdir -p $wd/data/dest_renamed/ 


for sample in EERU EEUA WEES WEFR
do
    for n in 1 2
    do 
        f="$wd/data/bamfiles/${sample}_${n}.bam"
        fres="$wd/data/dest_renamed/${sample}_${n}.bam"

        conda activate python_deps

        rgline=$(python3 $wd/scripts/RgLineForDest.py -i $f)

        conda deactivate

        conda activate gatk_modified

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

        conda activate python_deps

        rgline=$(python3 $wd/scripts/RgLineForDest.py -i $f)

        conda deactivate

        echo $rgline
        conda activate gatk_modified

        samtools addreplacerg --threads 100 -r "$rgline" -w -o $fres $f
        samtools index -@ 100 $fres

        conda deactivate
    done
done


rm -rf $wd/data/bamfiles/DGN_*.bam
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

## CREATE THE INPUTS FOR FreeBayes
mkdir $wd/data/freebayes_inputs/

### 3. Create a file listing all the BAM files
for f in $wd/data/bamfiles/*.bam
do
    echo $f >> $wd/data/freebayes_inputs/bamfiles.txt
done

## VARIANT CALLING WITH BCFTOOLS

conda activate freebayes-env
echo "will cite" | parallel --citation >/dev/null 2>&1
parallel --bar -j 100 bash ::: $wd/shell/bcftools_regions/*.sh
conda deactivate

## CONCATENATE VCF FILE

bcftools concat \
    -O z \
    --threads 100 \
    -o $wd/results/drosophila_evolution.bcftools_all.vcf.gz \
    $wd/results/drosophila_evolution.bcftools_2R.vcf.gz \
    $wd/results/drosophila_evolution.bcftools_2L.vcf.gz \
    $wd/results/drosophila_evolution.bcftools_3R.vcf.gz \
    $wd/results/drosophila_evolution.bcftools_3L.vcf.gz \
    $wd/results/drosophila_evolution.bcftools_X.vcf.gz  

## PSEUDO-POOLIFICATION

conda activate python_deps

python3 $wd/scripts/python/RandomAlleleFromVcf.py

conda deactivate

## CONVERSION TO VCF

head -n 1901 > $wd/results/drosophila_evolution.bcftools_fakepools_wholegen.vcf
cat $wd/results/fake_pools_all.tsv.gz >> $wd/results/drosophila_evolution.bcftools_fakepools_wholegen.vcf

conda activate gatk_modified

bcftools view -O z -o $wd/results/drosophila_evolution.bcftools_fakepools_wholegen.vcf.gz $wd/results/drosophila_evolution.bcftools_fakepools_wholegen.vcf

conda deactivate

## CALCULATION OF F-STATS

conda activate R

### 1. Without Drosophila simulans

Rscript $wd/scripts/r/FstatsWoDrosSim.r


### 2. With Drosophila simulans

Rscript $wd/scripts/r/FstatsWithDrosSim.r

## RANDOM ALLELE PCA

Rscript $wd/scripts/r/RandomAllelePca.r

conda deactivate

## KRIGING INTERPOLATION

conda activate python_deps

python3 $wd/scripts/python/KrigingWithMap.py

conda deactivate