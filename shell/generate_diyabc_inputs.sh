wd="/home/abert/media1/projects/drosophila-evolution"
vcf_file="$wd/results/drosophila_evolution.bcftools_fakepools_wholegen.vcf.gz"

mkdir -p $wd/data/simulation_vcfs/
mkdir -p $wd/data/simulation_data/logs/

for i in {1..4}
do
	conda activate gatk_modified
    bcftools view \
        -S $wd/data/simulation_data/subsample_${i}.txt \
        -O z \
        -o $wd/data/simulation_vcfs/subsample_${i}.vcf.gz \
        $vcf_file 2> $wd/data/simulation_data/logs/bcftools.log
	conda deactivate

    mkdir -p $wd/data/diyabc_inputs_${i}/

    conda activate R

    Rscript $wd/scripts/r/Vcf2DiyabcGen.r \
        $wd/data/simulation_vcfs/subsample_${i}.vcf.gz \
        $wd/data/diyabc_inputs_${i}/
    
    conda deactivate

    mv $wd/data/diyabc_inputs_${i}/PopData.diyabc $wd/data/diyabc_inputs_${i}/POOL_PopData.snp
done
