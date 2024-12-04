# VCF Subsampling Pipeline

[subsample_vcf.sh](../shell/subsample_vcf.sh) is a custom shellscript that combines the capabilities of several well-established bioinformatics tools to perform subsampling of VCF files.

## Workflow

The workflow of this pipeline is simple:

- The custom python script [`SamplesFromBamFiles.py`](../scripts/SamplesFromBamFiles.py) is employed to create a sample file that should serve as input for `bcftools view`, extracting the sample name from the `@RG` line of the BAM files listed in the input text file with the BAM files list
- `bcftools view` carries the actual subsampling starting from the samples file obtained in the previous step

## Inputs

As input, the pipeline takes:

- A compressed VCF file
- An text file with the list of BAM files that were used as source for the variant calling and that refer to the samples contained in the VCF

## Options

### REQUIRED ARGUMENTS:

- `-v`, `--vcf` VCF_FILE: Compressed VCF file to subsample (.vcf.gz or .extension.gz)
- `-s`, `--samples` SAMPLES: text file containing a list of bamfiles (more in the docs)
- `-o`, `--outputdir` OUTPUTDIR: Output directory where the results will be stored

### OPTIONAL ARGUMENTS:

- `-t`, `--threads` THREADS: integer number of threads (default: 1)
- `-sgn`, `--samples_group_name` SAMPLES_GROUP_NAME: name of the samples group to append at the name of the resulting VCF (more in the docs, default: 'samplegroup')

Input `subsample_vcf.sh -h` to show the help message.

> *The `--samples_group_name` option allows the user to specify a name for their samples group (for example 'European'): this will be added to the resulting VCF file, which will be then named BASENAME.European.vcf.gz (BASENAME is the basename of the input VCF file)*

## Outputs

In the results directory, you will find:

- `logs`: logs of the two steps are stored here
- `data`: samples files with the BAM list as input for the subsampling pipeline are stored here
- `results`: samples-specific VCFs are stored here

## Example usage

One could use the pipeline with the following example command:

```bash
# Example usage
bash shell/subsample_vcfs.sh \
    -s data/bamfiles.txt \
    -v data/all_samples.vcf.gz \
    -o results/subsamples/ \
    -sgn European \
    -t 100
```

