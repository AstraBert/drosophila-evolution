import polars as pl


def read_vcf(vcffile: str) -> pl.DataFrame:
    df = pl.read_csv(vcffile, separator="\t", n_threads=22, skip_rows=1907)
    return df

def read_sync(syncfile: str) -> pl.DataFrame:
    df = pl.read_csv(syncfile, separator="\t", n_threads=22)
    return df

def add_dros_sim(vcfdf: pl.DataFrame, syncdf: pl.DataFrame):
    syncdf = syncdf.rename({"CHROM": "#CHROM"})
    findf = syncdf.join(vcfdf, on=["#CHROM", "POS"], how="inner")
    print("Perfomerd join")
    finaldf = findf.drop(findf.columns[4:])
    ls = []
    positionsvcf = vcfdf["POS"].to_list()
    chromsvcf = vcfdf["#CHROM"].to_list() 
    refallelesvcf = vcfdf["REF"].to_list()
    altallelesvcf = vcfdf["ALT"].to_list() 
    vcf = {f"{chromsvcf[i]}.{positionsvcf[i]}": [refallelesvcf[i], altallelesvcf[i]] for i in range(len(positionsvcf))}  
    syncchrom = finaldf["#CHROM"].to_list() 
    syncpos = finaldf["POS"].to_list() 
    syncalleles = finaldf["ALL"].to_list() 
    sync = {f"{syncchrom[i]}.{syncpos[i]}": syncalleles[i] for i in range(len(syncpos))} 
    print("Created dicts")
    for k in vcf:
        if sync[k] == vcf[k][0]:
            ls.append("0/0:1,1,1:100:100,0:1,0,0:10")
        elif sync[k] == vcf[k][1]:
            ls.append("1/1:1,1,1:100:0,100:1,0,0:10")
        else:
            ls.append("./.:0,0,0:0:0,0:0,0,0:0") 
    print("Finished loop")
    drossim = pl.Series("DrosSim", ls)
    updateddf = vcfdf.insert_column(len(vcfdf.columns), drossim)
    print("Updated DF")
    return updateddf


if __name__ == "__main__":
    vcfdf = read_vcf("/gatk_modified/userdata/abertelli/drosophila-evolution/results/drosophila_evolution.bcftools_fakepools_wholegen.vcf.gz")
    print(vcfdf.head())
    print(vcfdf.height)
    syncdf = read_sync("/gatk_modified/userdata/abertelli/drosophila-evolution/data/dros_sim/dros_sim.sync.gz")
    print(syncdf.head())
    print(syncdf.height)
    updateddf = add_dros_sim(vcfdf, syncdf)
    print(updateddf.head())
    print(updateddf.height)
    updatedpd = updateddf.to_pandas()
    print("Transformed to Pandas")
    updatedpd.to_csv("/gatk_modified/userdata/abertelli/drosophila-evolution/results/drosophila_evolution.bcftools_fakepools_withsim.tsv.gz", sep="\t", index=False)



