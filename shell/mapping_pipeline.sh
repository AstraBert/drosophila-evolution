#!/bin/bash

# Usage function
usage() {
    echo "Usage: mapping_pipeline.sh -r,--reference REFERENCE -fq,--fastq FASTQ_FILE -o,--outputdir OUTPUTDIR [-t, --threads THREADS] 

    REQUIRED ARGUMENTS:
        -r, --reference REFERENCE: Compressed reference FASTA file (.fa.gz, .fasta.gz and other extensions.gz)
        -fq1, --fastq_file1 FASTQ_FILE1: Compressed or uncompressed input FASTQ file with raw paired-end reads
        -fq2, --fastq_file2 FASTQ_FILE2: Compressed or uncompressed input FASTQ file with raw paired-end reads
        -o, --outputdir OUTPUTDIR: Output directory where the results will be stored
    OPTIONAL ARGUMENTS:
        -t, --threads THREADS: integer number of threads (default: 200)

    Input mapping_pipeline.sh -h to show this message again.
    
    DOCS: https://github.com/AstraBert/drosophila-evolution/blob/main/docs/mapping_pipeline.md"
    exit 1
}

eval "$(conda shell.bash hook)"

reference_file=""
read_1=""
read_2=""
outputdir=""
threads=200

while [ $# -gt 0 ]; do
    case "$1" in
        -h|--help)
            usage
            ;;
        -r|--reference)
            reference_file="$2"
            shift 2
            ;;
        -t|--threads)
            threads="$2"
            shift 2
            ;;
        -fq1|--fastq1)
            read_1="$2"
            shift 2
            ;;
        -fq2|--fastq2)
            read_2="$2"
            shift 2
            ;;
        -o|--outputdir)
            outputdir="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done 

if [ -z $reference_file ] ||  [ -z $read_1 ] || [ -z $read_2 ]  ||  [ -z $outputdir ]; then
    echo "Missing required arguments!"
    usage
fi

if [ ! -f $reference_file ] ||  [ ! -f $read_1 ] || [ ! -f $read_2 ]; then
    echo "One or more files you provided seem not to exist: please check them and retry"
    usage
fi

smpl=$(basename $read_1)
sample="${smpl%_*}"
rfrnc="${reference_file%.*.*}"

mkdir -p $outputdir
mkdir -p $outputdir/logs

echo """ARGUMENTS FOR THE ANALYSIS:
reference_file = $reference_file
fastq_file1 = $read_1
fastq_file2 = $read_2
outputdir = $outputdir
threads = $threads
""" > $outputdir/logs/analysis_args.txt

echo """ARGUMENTS FOR THE ANALYSIS:
reference_file = $reference_file
fastq_file1 = $read_1
fastq_file2 = $read_2
outputdir = $outputdir
threads = $threads
"""

echo ""
echo ""

echo "STARTING THE ANALYSIS"
echo ""
echo "(You will find the logs of each step at: $outputdir/logs)"
echo ""
echo "0%>                         100%"

conda activate gatk_modified

echo "Trimming of reads"
echo "0%---------->                 100%"
cutadapt \
    -q 18 \
    --minimum-length 25 \
    -o $outputdir/$sample.trimmed1.fq.gz \
    -p $outputdir/$sample.trimmed2.fq.gz \
    -b ACACTCTTTCCCTACACGACGCTCTTCCGATC \
    -B CAAGCAGAAGACGGCATACGAGAT \
    -O 15 \
    -n 3 \
    --cores=$threads \
    $read_1 $read_2 > $outputdir/logs/cutadapt.log 2>&1

echo "Mapping and sorting"
echo "0%--------------->              100%"
bwa-mem2 mem \
    -t 100 \
    -M \
    -R "@RG\tID:$sample\tSM:$sample\tPL:illumina\tLB:lib1" \
    $reference_file \
    $outputdir/$sample.trimmed1.fq.gz \
    $outputdir/$sample.trimmed2.fq.gz |
    samtools view -@ $threads -Sbh -q 20 -F 0x100 - |
    samtools sort -@ $threads - >$outputdir/$sample.sorted_merged.bam

echo "Deduplicating"
echo "0%-------------------->         100%"

conda deactivate

conda activate picard

picard MarkDuplicates \
    REMOVE_DUPLICATES=true \
    I=$outputdir/$sample.sorted_merged.bam \
    O=$outputdir/$sample.dedup.bam \
    M=$outputdir/$sample.mark_duplicates_report.txt \
    VALIDATION_STRINGENCY=SILENT > $outputdir/logs/dedup.log 2>&1

conda deactivate

conda activate gatk_modified

samtools index $outputdir/$sample.dedup.bam > $outputdir/logs/samtools_ind.log 2>&1

conda deactivate

echo "END OF THE ANALYSIS"
echo "0%------------------------------>100%"
