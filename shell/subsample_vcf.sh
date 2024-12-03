#!/bin/bash

# Usage function
usage() {
    echo "Usage: subsample_vcf.sh -v VCF_FILE -s SAMPLES -o OUTPUTDIR

    REQUIRED ARGUMENTS:
        -v, --vcf VCF_FILE: Compressed VCF file to subsample (.vcf.gz or .extension.gz)
        -s, --samples SAMPLES: text file containing a list of bamfiles (more in the docs)
        -o, --outputdir OUTPUTDIR: Output directory where the results will be stored
    OPTIONAL ARGUMENTS:
        -t, --threads THREADS: integer number of threads (default: 1)
        -sgn, --samples_group_name SAMPLES_GROUP_NAME: name of the samples group to append at the name of the resulting VCF (more in the docs, default: 'samplegroup')

    Input subsample_vcf.sh -h to show this message again."
    exit 1
}

vcf_file=""
samples=""
outputdir=""
samplegroupname="samplegroup"
threads=1

while [ $# -gt 0 ]; do
    case "$1" in
        -h|--help)
            usage
            ;;
        -v|--vcf)
            vcf_file="$2"
            shift 2
            ;;
        -s|--samples)
            samples="$2"
            shift 2
            ;;
        -o|--outputdir)
            outputdir="$2"
            shift 2
            ;;
        -t|--threads)
            threads="$2"
            shift 2
            ;;
        -sgn|--samples_group_name)
            samplegroupname="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done 

if [ -z $vcf_file ] ||  [ -z $samples ]  ||  [ -z $outputdir ]; then
    echo "Missing required arguments!"
    usage
fi

if [ ! -f $vcf_file ] ||  [ ! -f $samples ]; then
    echo "One or more files you provided seem not to exist: please check them and retry"
    usage
fi

smpl=$(basename $vcf_file)
namespec="${smpl%.*.*}"

mkdir -p $outputdir
mkdir -p $outputdir/logs

echo """ARGUMENTS FOR THE ANALYSIS:
vcf_file = $vcf_file
fastq_file1 = $samples
samples = $outputdir
threads = $threads
""" > $outputdir/logs/analysis_args.txt

echo """ARGUMENTS FOR THE ANALYSIS:
vcf_file = $vcf_file
fastq_file1 = $samples
samples = $outputdir
threads = $threads
"""

echo ""
echo ""

echo "STARTING THE ANALYSIS"
echo ""
echo "(You will find the logs of each step at: $outputdir/logs)"
echo ""

# Capture start time
start_time=$(date +%s)

source activate python_deps

mkdir -p $outputdir/data/

outsmpl=$(basename $samples)
outsampls="${outsmpl%.*}"
outsamples=${outsmpl}.samples.txt

scriptdir=$(dirname $0)

python3 $scriptdir/../scripts/SamplesFromBamFiles.py \
    -i $samples \
    -o $outputdir/data/$outsamples 2> $outputdir/logs/SamplesFromBamFiles.log

conda deactivate

source activate gatk_modified

mkdir -p $outputdir/results/

bcftools view \
    -S $outputdir/data/$outsamples \
    -O z \
    -o $outputdir/results/$namespec.$samplegroupname.vcf.gz \
    $vcf_file 2> $outputdir/logs/bcftools.log

conda deactivate

end_time=$(date +%s)

execution_time=$((end_time - start_time))

execution_time_hr=$((execution_time / 3600))
execution_time_min=$(((execution_time % 3600) / 60))
execution_time_sec=$((execution_time % 60))

echo "YOUR ANALYSIS IS COMPLETE IN $execution_time_hr hours, $execution_time_min minutes, and $execution_time_sec seconds"
