wd="/home/abert/media1/projects/drosophila-evolution"
vcf_file="$wd/results/drosophila_evolution.bcftools_fakepools_wholegen.vcf.gz"

mkdir -p $wd/data/simulation_data/logs/

for i in {1..4}
do

   	input_vcf="$wd/data/simulation_vcfs/subsample_${i}.vcf.gz"
    mkdir -p $wd/data/diyabc_inputs_${i}/

    source activate R

    Rscript $wd/scripts/r/Vcf2DiyabcGen.r \
        $input_vcf \
        $wd/data/diyabc_inputs_${i}/
    
    conda deactivate

    mv $wd/data/diyabc_inputs_${i}/PopData.diyabc $wd/data/diyabc_inputs_${i}/POOL_PopData.snp
done
