source activate gatk_modified

wd="/gatk_modified/drosophila-project/data/pop_data/test_diyabc"

bcftools filter -O z -o ${wd}/DEST_but_is_four_pops.bcffilt.vcf.gz -t "^2L:1725744-13654180,2R:14891154-20776334,3L:2673046-16808841,3R:20906917-29531297,3R:15932209-25244010,3R:19596867-32579331,3R:11250567-26640370,X,Y,4" ${wd}/DEST_but_is_four_pops.vcf.gz

conda deactivate



