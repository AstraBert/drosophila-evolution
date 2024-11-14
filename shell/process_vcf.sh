source activate python_deps

python3 /gatk_modified/drosophila-project/scripts/RegionsToYaml.py \
    -i /gatk_modified/drosophila-project/data/pop_data/test_diyabc/regs_w_invs.txt \
    -o /gatk_modified/drosophila-project/data/pop_data/test_diyabc/regs_w_invs

pigz -d -p 100 /gatk_modified/drosophila-project/data/pop_data/test_diyabc/DEST_but_is_four_pops.vcf.gz

# only allow -im if ou want the intermediate files to be saved
python3 /gatk_modified/drosophila-project/scripts/FilterVcf.py \
    -i /gatk_modified/drosophila-project/data/pop_data/test_diyabc/DEST_but_is_four_pops.vcf \
    -ff /gatk_modified/drosophila-project/data/pop_data/test_diyabc/regs_w_invs.yaml \
    -hd "CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,MA_Tan_Lar_1_2021-09-06,FR_Cot_Vie_1_2019-10-19,UA_Don_Mar_1_2021-08-23,US_Vir_Cha_1_2016-10-03" \
    -o  /gatk_modified/drosophila-project/data/pop_data/test_diyabc/DEST_but_is_four_filtered_pops \
    -im   

pigz -p 100 /gatk_modified/drosophila-project/data/pop_data/test_diyabc/DEST_but_is_four_pops.vcf
pigz -p 100 /gatk_modified/drosophila-project/data/pop_data/test_diyabc/DEST_but_is_four_filtered_pops*

conda deactivate

