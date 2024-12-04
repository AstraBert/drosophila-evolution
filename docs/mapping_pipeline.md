# Mapping Pipeline

<div align="center">
    <img src="../imgs/data_mapping.jpg" alt="Pipeline for mapping data">
    <p><i>Raw sequencing data mapping pipeline</i></p>
    <sbr>
</div>

[mapping_pipeline.sh](../shell/mapping_pipeline.sh) is a custom shellscript that combines the capabilities of several well-established bioinformatics tools to perform mapping on raw sequencing, paired-end read data deriving from Illumina technologies.

## Workflow

As represented by the flowchart at the beginning of this page, the workflow of the pipeline is straightforward: 

- The two input FASTQ files are trimmed with `cutadapt`, producing one joint FASTQ file
- The resulting FASTQ file is mapped against the indexed reference genome with `bwa-mem2`
- The resulting BAM file is deduplicated with `picard MarkDuplicates`
- The deduplicated BAM file is indexed with `samtools index`

## Inputs

As input, the pipeline takes:

- Two FASTQ files representing raw, paired-end reads
- The already-indexed reference genome in compressed FASTA format

## Options

### REQUIRED ARGUMENTS:
- `-r`, `--reference` REFERENCE: Compressed reference FASTA file (.fa.gz, .fasta.gz and other extensions.gz)
- `-fq1`, `--fastq_file1` FASTQ_FILE1: Compressed or uncompressed input FASTQ file with raw paired-end reads
- `-fq2`, `--fastq_file2` FASTQ_FILE2: Compressed or uncompressed input FASTQ file with raw paired-end reads
- `-o`, `--outputdir` OUTPUTDIR: Output directory where the results will be stored

### OPTIONAL ARGUMENTS:
- `-t`, `--threads` THREADS: integer number of threads (default: 200)

Input `mapping_pipeline.sh -h` to show the help message

## Outputs

In the output directory you will find:

- The combined FASTQ file resulting from trimming
- A `.bam` file (the original BAM file)
- A `.dedup.bam` file (the deduplicated BAM file)
- A `.dedup.bam.bai` file (the index for the deduplicated BAM file)

## Example usage

One could use the pipeline with the following example command:

```bash
# Example usage
bash shell/mapping_pipeline.sh \
    -r ./data/reference.fa.gz \
    -fq1 ./raw_reads/sample_1.fastq.gz \
    -fq2 ./raw_reads/sample_2.fastq.gz \
    -o ./results_directory/ \
    -t 100
```

