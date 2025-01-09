import polars as pl


def read_vcf(vcffile: str) -> pl.DataFrame:
    df = pl.read_csv(vcffile, separator="\t", n_threads=22, skip_rows=1907)
    return df

def read_sync(syncfile: str) -> pl.DataFrame:
    df = pl.read_csv(syncfile, separator="\t", n_threads=22)
    return df

def add_dros_sim(vcfdf: pl.DataFrame, syncdf: pl.DataFrame):
    syncdf = syncdf.rename({"CHROM": "Chromosome"})
    syncdf = syncdf.rename({"POS": "Position"})
    findf = syncdf.join(vcfdf, on=["Chromosome", "Position"], how="inner")
    print("Perfomerd join")
    finaldf = findf.drop(findf.columns[4:])
    ls = []
    positionsvcf = vcfdf["Position"].to_list()
    chromsvcf = vcfdf["Chromosome"].to_list() 
    refallelesvcf = vcfdf["RefAllele"].to_list()
    altallelesvcf = vcfdf["AltAllele"].to_list() 
    vcf = {f"{chromsvcf[i]}.{positionsvcf[i]}": [refallelesvcf[i], altallelesvcf[i]] for i in range(len(positionsvcf))}  
    syncchrom = finaldf["Chromosome"].to_list() 
    syncpos = finaldf["Position"].to_list() 
    syncalleles = finaldf["ALL"].to_list() 
    sync = {f"{syncchrom[i]}.{syncpos[i]}": syncalleles[i] for i in range(len(syncpos))} 
    print("Created dicts")
    for k in vcf:
        if sync[k] == vcf[k][0]:
            ls.append("100")
        elif sync[k] == vcf[k][1]:
            ls.append("0")
        else:
            ls.append("0") 
    print("Finished loop")
    snps_counts = ["Pool21\n"]+[f"{l}\n" for l in ls]
    return snps_counts 


if __name__ == "__main__":
    vcfdf = pl.read_csv("/gatk_modified/userdata/abertelli/drosophila-evolution/results/all_pops_snps.csv", n_threads=20)
    print(vcfdf.head())
    print(vcfdf.height)
    syncdf = read_sync("/gatk_modified/userdata/abertelli/drosophila-evolution/data/dros_sim/dros_sim.sync.gz")
    print(syncdf.head())
    print(syncdf.height)
    snps_counts = add_dros_sim(vcfdf, syncdf)
    f = open("/gatk_modified/userdata/abertelli/drosophila-evolution/results/drossim_snps.csv", "w")
    f.writelines(snps_counts)
    f.close()

    



