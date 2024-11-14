source activate gatk_modified
counter=0
for f in /gatk_modified/drosophila-project/data/DGN/*.bam
do
    ((counter++))
    samtools index $f
    mv $f /gatk_modified/drosophila-project/data/DGN_${counter}.bam
    mv ${f}.bai /gatk_modified/drosophila-project/data/DGN_${counter}.bam.bai
done
conda deactivate