eval "$(conda shell.bash hook)"
wd="/home/abert/media1/projects/drosophila-evolution"
conda activate freebayes-env
echo "will cite" | parallel --citation >/dev/null 2>&1
parallel --bar -j 120 bash ::: $wd/shell/bcftools_regions/*.sh
conda deactivate
