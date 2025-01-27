import polars as pl

def read_vcf(vcffile: str) -> pl.DataFrame:
    df = pl.read_csv(vcffile, separator="\t", n_threads=22, skip_rows=1911)
    return df

if __name__ == "__main__":
    sync = pl.read_csv("dros_sim.sync.gz", separator="\t")
    sync = sync.rename({"CHROM": "Chromosome", "POS": "Position", "ALL": "DrosSim"})
    sync = sync.drop("COUNT")
    vcf = read_vcf("drosevol.noneurope.noinv.vcf.gz")
    vcf = vcf.rename({"#CHROM": "Chromosome", "POS": "Position", "REF": "RefAllele", "ALT": "AltAllele"})
    vcf = vcf.drop(["ID", "QUAL", "FILTER", "INFO", "FORMAT"])
    print(vcf.head())
    df = pl.read_csv("dest.europe.trk.noinv_snpinfo.csv", separator="\t")
    print(df.head())
    dest_dros = vcf.join(df, on=["Chromosome", "Position", "RefAllele", "AltAllele"], how="inner")
    dest_dros_sim = dest_dros.join(sync, on=["Chromosome", "Position"], how="inner") 
    selected_snps = dest_dros_sim["id"].to_list()
    t = open("selected_dest_eu_snps_noinv.csv", "w")
    t.write("POS\n")
    selected_snps = [str(snp)+"\n" for snp in selected_snps]
    t.writelines(selected_snps)
    t.close()
    drossim = dest_dros_sim["DrosSim"].to_list()
    refalleles = dest_dros_sim["RefAllele"].to_list()
    altalleles = dest_dros_sim["AltAllele"].to_list()
    drossim_allelic_status = [
        0 if drossim[i] == altalleles[i]
        else 100
        for i in range(len(drossim))
    ]
    dgn = dest_dros_sim["DGN"].to_list()
    dgn_allelic_status = [
        0 if el.split(":")[0] == "1/1" or el.split(":")[0] == "./."
        else 20 if el.split(":")[0] == "0/1" or el.split(":")[0] == "1/0"
        else 40
        for el in dgn
    ]
    cnxj = dest_dros_sim["CNXJ"].to_list()
    cnxj_allelic_status = [
        0 if el.split(":")[0] == "1/1" or el.split(":")[0] == "./."
        else 13 if el.split(":")[0] == "0/1" or el.split(":")[0] == "1/0"
        else 25
        for el in cnxj
    ]
    cnother = dest_dros_sim["CnOther"].to_list()
    cnother_allelic_status = [
        0 if el.split(":")[0] == "1/1" or el.split(":")[0] == "./."
        else 25 if el.split(":")[0] == "0/1" or el.split(":")[0] == "1/0"
        else 50
        for el in cnother
    ]
    cnqtp = dest_dros_sim["CnQTP"].to_list()
    cnqtp_allelic_status = [
        0 if el.split(":")[0] == "1/1" or el.split(":")[0] == "./."
        else 25 if el.split(":")[0] == "0/1" or el.split(":")[0] == "1/0"
        else 50
        for el in cnqtp
    ]
    isr = dest_dros_sim["ISR"].to_list()
    isr_allelic_status = [
        0 if el.split(":")[0] == "1/1" or el.split(":")[0] == "./."
        else 16 if el.split(":")[0] == "0/1" or el.split(":")[0] == "1/0"
        else 32
        for el in cnqtp
    ]
    f = open("drosevol_readcount_noinv.csv", "w")
    f.write("DGN,CNXJ,CnOther,CnQTP,ISR,DrosSim\n")
    for el in range(len(isr_allelic_status)):
        f.write(f"{dgn_allelic_status[el]},{cnxj_allelic_status[el]},{cnother_allelic_status[el]},{cnqtp_allelic_status[el]},{isr_allelic_status[el]},{drossim_allelic_status[el]}\n")
    f.close()
    
