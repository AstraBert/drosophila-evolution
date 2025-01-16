import polars as pl

def read_vcf(vcffile: str) -> pl.DataFrame:
    df = pl.read_csv(vcffile, separator="\t", n_threads=22, skip_rows=1909)
    return df

if __name__ == "__main__":
    vcf = read_vcf("drosevol.noneurope.vcf.gz")
    vcf.rename({"#CHROM": "Chromosome", "POS": "Position", "REF": "RefAllele", "ALT": "AltAllele"})
    vcf.drop(["ID", "QUAL", "FILTER", "INFO", "FORMAT"])
    df = pl.read_csv("dest_eu_snps.tsv", separator="\t", columns=["Chromosome", "Position", "RefAllele", "AltAllele"])
    dest_dros = vcf.join(df, on=["Chromosome", "Position", "RefAllele", "AltAllele"], how="inner")
    m = dest_dros["DGN"].to_list()
    dest_dros["DGN"] = pl.Series([40 for i in range(len(m))])
    dest_dros["CNXJ"] = pl.Series([25 for i in range(len(m))])
    dest_dros["CnOther"] = pl.Series([50 for i in range(len(m))])
    dest_dros["CnQTP"] = pl.Series([50 for i in range(len(m))])
    dest_dros["ISR"] = pl.Series([32 for i in range(len(m))])
    dff = dest_dros.select(["DGN", "CNXJ", "CnOther", "CnQTP", "ISR"])
    dff.write_csv("drosevol_readcount.csv")
    
