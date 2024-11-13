source activate gatk_modified

bcftools view -s 'ZM_Sou_Sia_1_2010-07-16,FR_Cot_Vie_1_2019-10-19,UA_Don_Mar_1_2021-08-23,US_Vir_Cha_1_2016-10-03' /gatk_modified/drosophila-project/data/pop_data/test_diyabc/DEST.spns.neutral.vcf.gz -O z -o /gatk_modified/drosophila-project/data/pop_data/test_diyabc/DEST.four_pops.neutral.vcf.gz

conda deactivate