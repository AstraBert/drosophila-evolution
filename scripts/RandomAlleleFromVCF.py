import polars as pl  
from typing import List, Dict
import pandas as pd
import random as r
from collections import Counter
import pysam

def read_vcf(vcffile: str) -> pl.DataFrame:
    df = pl.read_csv(vcffile, separator="\t", n_threads=16, skip_rows=1901)
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
        return "."
    else:
        allele = r.choice(total)
        return allele

def most_common_el(ls: list):
    counter = Counter(ls)
    element = counter.most_common(1)[0][0] 
    return element 


def snps_map(df: pl.DataFrame, population_to_samples: Dict[str,List[str]]) -> List[List[float]]:
    population_to_snp_freq = {k: [] for k in population_to_samples}  
    refs = df["REF"].to_list() 
    alts = df["ALT"].to_list() 
    poss = df["POS"].to_list() 
    for k in population_to_samples:
        if len(population_to_samples[k])>1:
            pops_freqs = [] 
            for sample in population_to_samples[k]:
                samplelist = df[sample].to_list() 
                freqslist = [randomly_subsample(refs[el], alts[el], samplelist[el]) for el in range(len(samplelist))]
                pops_freqs.append(freqslist)  
            pop_freq = []  
            for i in range(len(pops_freqs[0])):
                l = [] 
                for j in range(len(pops_freqs)):
                    l.append(pops_freqs[j][i]) 
                pop_freq.append(most_common_el(l))
            population_to_snp_freq[k] = pop_freq
    print(population_to_snp_freq["CNXJ"])
    pseudo_df = {"population": list(population_to_snp_freq.keys())}
    for j in range(len(poss)):
        key = f"{poss[j]}"
        values = [population_to_snp_freq[pseudo_df["population"][el]][j] for el in range(len(pseudo_df["population"]))] 
        pseudo_df.update({key: values})
    actual_df = pl.DataFrame(pseudo_df)
    return actual_df 


def extract_sample_from_rgline(bamfile: str) -> str:
    bam = pysam.AlignmentFile(bamfile, "rb")
    head = bam.header
    strhead = str(head)
    lines = strhead.split("\n")
    rgline = [line for line in lines if line.startswith("@RG")]
    rgline = rgline[0]
    rgline_comps = rgline.split("\t")
    sampletag = [comp for comp in rgline_comps if comp.startswith("SM:")]
    sample = sampletag[0].split(":")[1]
    return sample

def find_pops_from_bamlist(bamfile: str, population_list: List[str]) -> Dict[str,List[str]]:
    populations = {el: [] for el in population_list}
    bamfs = open(bamfile, "r")
    bamlist = [bam.replace("\n", "") for bam in bamfs.readlines()]
    for k in populations:
        for bam in bamlist:
            bamname = bam.split("/")[-1].split(".")[0]
            if bamname.startswith(k):
                populations[k].append(bam) 
    population_to_samples = {el: [] for el in population_list}
    for k in populations:
        samples = [extract_sample_from_rgline(el) for el in populations[k]]
        population_to_samples[k] = samples
    return population_to_samples

if __name__ == "__main__":
    bamfile = "/gatk_modified/userdata/abertelli/drosophila-evolution/data/freebayes_inputs/bamfiles.txt"
    pops = ['CNXJ', 'CnOther', 'CnQTP', 'ISR']
    pops2samples = find_pops_from_bamlist(bamfile, pops)
    print("Reading VCF...")
    df = read_vcf("/gatk_modified/userdata/abertelli/drosophila-evolution/results/example.vcf")
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
    pddf = newdf.to_pandas()
    pddf.to_csv("/gatk_modified/userdata/abertelli/drosophila-evolution/results/fake_pools.tsv", sep="\t", index=False)