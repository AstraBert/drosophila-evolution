import polars as pl

def read_vcf(vcffile: str) -> pl.DataFrame:
    df = pl.read_csv(vcffile, separator="\t", n_threads=22, skip_rows=1909)
    return df

if __name__ == "__main__":
    vcf = read_vcf("drosevol.noneurope.vcf.gz")
    vcf = vcf.rename({"#CHROM": "Chromosome", "POS": "Position", "REF": "RefAllele", "ALT": "AltAllele"})
    vcf = vcf.drop(["ID", "QUAL", "FILTER", "INFO", "FORMAT"])
    print(vcf.head())
    df = pl.read_csv("dest_eu_snps.tsv", separator="\t", columns=["Chromosome", "Position", "RefAllele", "AltAllele"])
    print(df.head())
    dest_dros = vcf.join(df, on=["Chromosome", "Position", "RefAllele", "AltAllele"], how="inner")
    dgn = dest_dros["DGN"].to_list()
    dgn_allelic_status = [0 if el.split(":")[0]=="1/1" or el.split(":")[0]=="./." else 40 for el in dgn]
    cnxj = dest_dros["DGN"].to_list()
    cnxj_allelic_status = [0 if el.split(":")[0]=="1/1" or el.split(":")[0]=="./." else 25 for el in cnxj]
    cnother = dest_dros["DGN"].to_list()
    cnother_allelic_status = [0 if el.split(":")[0]=="1/1" or el.split(":")[0]=="./." else 50 for el in cnother]
    cnqtp = dest_dros["DGN"].to_list()
    cnqtp_allelic_status = [0 if el.split(":")[0]=="1/1" or el.split(":")[0]=="./." else 50 for el in cnother]
    isr = dest_dros["DGN"].to_list()
    dest_dros["DGN"] = dgn_allelic_status
    dest_dros["CNXJ"] = cnxj_allelic_status
    dest_dros["CnOther"] = cnother_allelic_status
    dest_dros["CnQTP"] = cnqtp_allelic_status
    dff = dest_dros.select(["DGN", "CNXJ", "CnOther", "CnQTP", "ISR"])
    dff.write_csv("drosevol_readcount.csv")
    
