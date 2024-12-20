sudo docker run -it -v /home/minionuser/miniondata1/:/gatk_modified/userdata/ astrabert/silly-gat-kay:latest

df1 = pl.DataFrame({"CHROM": ["chr1", "chr2", "chr3", "chr4"],"POS": [100, 200, 300, 400],"VALUE": [10, 20, 30, 40]})
df2 = pl.DataFrame({"CHROM": ["chr1", "chr3"],"POS": [100, 300]})
