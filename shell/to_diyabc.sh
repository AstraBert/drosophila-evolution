eval "$(conda shell.bash hook)"

wd="/home/abert/media1/projects/drosophila-evolution"

mkdir -p $wd/data/simulation_vcfs/

for i in {1..4}
do
    mkdir -p $wd/data/diyabc_inputs_${i}/

    conda activate R

    Rscript $wd/scripts/r/Vcf2DiyabcGen.r \
        $wd/data/simulation_vcfs/subsample_${i}_noinv.vcf.gz \
        $wd/data/diyabc_inputs_${i}/
    
    conda deactivate
done
