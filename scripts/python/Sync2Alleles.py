import polars as pl

snp_file = input("SNP file: ")

csv = pl.read_csv(snp_file)
sync = pl.read_csv("dros_sim.sync.gz", separator="\t")

csv = csv.select(pl.col("Chromosome"), pl.col("Position"), pl.col("RefAllele"), pl.col("AltAllele"))
sync = sync.rename({"CHROM": "Chromosome", "POS": "Position"})
sync = sync.drop("COUNT")
drossim = csv.join(sync, on=["Chromosome", "Position"], how="inner")

# Add the comparison column
drossim = drossim.with_columns([
    pl.when(pl.col("ALL") == pl.col("RefAllele"))
    .then(100)
    .otherwise(0)
    .alias("RefAlleleMatch")
])

counts = drossim["RefAlleleMatch"].to_list()
counts = [str(count)+"\n" for count in counts]
f = open("drossim_snps_finalsims.csv", "w")
counts = ["DrosSim\n"]+counts
f.writelines(counts)
f.close()
