wd="/home/abert/media1/projects/drosophila-evolution/data/diyabc_inputs_2/"

cd $wd

source activate R
Rscript ../../scripts/diyabc_problem_analysis/Rscrip_to_select_SNPs_for_poolseq_analyses_AE_15-01-2025.R
mv POOL_PopData_df_poolseq_readmin_20.txt POOL_PopData_1.snp
conda deactivate

cat POOL_PopData.snp | head -n 2 > POOL_PopData_2.snp
cat POOL_PopData_1.snp >> POOL_PopData_2.snp
micro POOL_PopData_2.snp

cat ../sample_header_diyabc.txt > header.txt
micro header.txt

sleep 15

cd $wd/../diyabc_executable/
./diyabc-RF-linux-v1.1.54 -p $wd -n "t:100;c:1;s:2011;f:f"
./diyabc-RF-linux-v1.1.54 -p $wd -R "ALL" -t 100 -g 120 -r 12000 -m

cd $wd
source activate R
Rscript ../../scripts/abcrf_v4.R
conda deactivate

