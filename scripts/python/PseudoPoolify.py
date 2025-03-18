import polars as pl  
from typing import List, Dict, Tuple
import pandas as pd
import random as r
from collections import Counter
import json

def read_vcf(vcffile: str, skip_rows: int, n_threads: int = 20) -> pl.DataFrame:
    df = pl.read_csv(vcffile, separator="\t", n_threads=n_threads, skip_rows=skip_rows)
    return df

def snps(df: pl.DataFrame) -> pl.DataFrame:
    snpsdf = df.filter((pl.col("REF").is_in(["A","T","C","G"])) & (pl.col("ALT").is_in(["A","T","C","G"])))
    return snpsdf

def randomly_subsample(ref: str, alt: str, infostr: str):
    refad, altad = int(infostr.split(":")[3].split(",")[0]), int(infostr.split(":")[3].split(",")[1])
    refls = [ref for i in range(refad)] 
    altls = [alt for i in range(altad)] 
    total = refls + altls
    if len(total) == 0:
        return (".", 0, 0, "NONE")
    else:
        allele = r.choice(total)
        if allele in refls:
            return (allele, len(refls), "REF")
        return (allele, len(altls), "ALT")

def most_common_el(ls: list):
    counter = Counter(ls)
    element = counter.most_common(1)[0][0] 
    return element 

def is_ref(allele: str, ls: list):
    for i in ls:
        if allele == i[0]:
            if i[2] == "REF":
                return True, "REF"
            return False, i[2]

def build_infostr(ls: list, originfostr: str):
    infostr = ""
    alleles = [l[0] for l in ls]
    alleles_without = [allele for allele in alleles if allele != "."] 
    if len(set(alleles_without)) == 2:
        infostr+="0/1:"
    elif len(set(alleles_without)) == 1:
        allele = alleles_without[0]
        if is_ref(allele,ls)[0]:
            infostr+="0/0:"
        else:
            infostr+="1/1:"
    else:
        infostr="./.:0,0,0:0:0,0:0,0,0:0"
        return infostr
    infostr+=originfostr.split(":")[1]+":"
    infostr+=str(sum([l[1] for l in ls]))+":"
    if len(set(alleles)) >= 2:
        reflen = sum([l[1] for l in ls if l[2]=="REF"])
        altlen = sum([l[1] for l in ls if l[2]=="ALT"])
        infostr+=f"{reflen},{altlen}:"
    else:
        allele = alleles[0]
        if is_ref(allele,ls)[0]:
            infostr+=f"{sum([l[1] for l in ls])},0:"
        else:
            infostr+=f"0,{sum([l[1] for l in ls])}:"
    infostr+=originfostr.split(":")[4]+":"
    infostr+=originfostr.split(":")[5]
    return infostr
            

def snps_map(df: pl.DataFrame, population_to_samples: Dict[str,List[str]]) -> pl.DataFrame:
    population_to_snp_freq = {k: [] for k in population_to_samples}  
    refs = df["REF"].to_list() 
    alts = df["ALT"].to_list() 
    for k in population_to_samples:
        if len(population_to_samples[k])>1:
            pops_freqs = [] 
            origsamplestr = df[population_to_samples[k][0]].to_list()
            for sample in population_to_samples[k]:
                samplelist = df[sample].to_list() 
                freqslist = [randomly_subsample(refs[el], alts[el], samplelist[el]) for el in range(len(samplelist))]
                pops_freqs.append(freqslist)  
            pop_freq = []  
            for i in range(len(pops_freqs[0])):
                l = [] 
                for j in range(len(pops_freqs)):
                    l.append(pops_freqs[j][i]) 
                pop_freq.append(build_infostr(l, origsamplestr[i]))
            population_to_snp_freq[k] = pop_freq
    actual_df = pl.DataFrame(population_to_snp_freq)
    return actual_df 

if __name__ == "__main__":
    nskip = int(input("Number of rows to skip: "))
    n_thr = int(input("Number of threads to use: "))
    js_pop = input("Path to the JSON file with populations: ")
    vcf_pt = input("Path to the VCF file: ")
    opt = input("Output TSV path: ")
    fjs = open(js_pop)
    js = json.load(fjs) 
    fjs.close()   
    pops2samples = js["populations"]
    print("Reading VCF...")
    df = read_vcf(vcf_pt, nskip, n_thr)
    print(df.head())
    print(df.height)
    print("Read VCF!")
    print("Selecting only SNPs and subsampling")
    snps_df = snps(df)
    print(snps_df.head())
    print(snps_df.height)
    print("Selected only SNPs and subsampling!")
    print("Mapping")
    newdf = snps_map(snps_df, pops2samples)
    print(newdf.head())
    print(newdf.height)
    for k in pops2samples:
        snps_df = snps_df.drop(pops2samples[k])
    cols = newdf.columns
    for c in cols:
        snps_df.insert_column(len(snps_df.columns), pl.Series(c, newdf[c].to_list()))
    snps_df.write_csv(opt, separator="\t")
