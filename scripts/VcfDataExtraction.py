import polars as pl  
from typing import List, Dict
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt 
import seaborn as sns
import numpy as np
import pysam

def read_vcf(vcffile: str) -> pl.DataFrame:
    df = pl.read_csv(vcffile, separator="\t", n_threads=16, skip_rows=1901)
    return df

def snps_and_random_subsample(df: pl.DataFrame) -> pl.DataFrame:
    snpsdf = df.filter((pl.col("REF").is_in(["A","T","C","G"])) & (pl.col("ALT").is_in(["A","T","C","G"])))
    reddf = snpsdf.sample(50000, seed=42)
    return reddf

def find_ref_copies(string: str):
    refalt = string.split(":")[0]
    ref, alt = refalt.split("/")[0], refalt.split("/")[1]
    if ref == "." or alt == ".":
        return -1
    return 2 - int(ref) - int(alt)

def snps_map(df: pl.DataFrame, population_to_samples: Dict[str,List[str]]) -> List[List[float]]:
    population_to_snp_freq = {k: [] for k in population_to_samples}  
    for k in population_to_samples:
        if len(population_to_samples[k])>1:
            pops_freqs = [] 
            for sample in population_to_samples[k]:
                samplelist = df[sample].to_list() 
                freqslist = [find_ref_copies(el) for el in samplelist]
                pops_freqs.append(freqslist)  
            pop_freq = []  
            for i in range(len(pops_freqs[0])):
                l = 0
                for j in range(len(pops_freqs)):
                    l += pops_freqs[j][i] 
                pop_freq.append(l / len(pops_freqs))
            population_to_snp_freq[k] = pop_freq
        else:
            sample = population_to_samples[k][0] 
            samplelist = df[sample].to_list() 
            freqslist = [find_ref_copies(el) for el in samplelist]
            population_to_snp_freq[k] = freqslist
    pseudo_df = {"population": list(population_to_snp_freq.keys())}
    for j in range(len(population_to_snp_freq[pseudo_df["population"][0]])):
        key = f"snp_{j}"
        values = [population_to_snp_freq[pseudo_df["population"][el]][j] for el in range(len(pseudo_df["population"]))] 
        pseudo_df.update({key: values})
    actual_df = pl.DataFrame(pseudo_df)
    return actual_df 


def perform_PCA(df: pl.DataFrame):
    pd_df = df.to_pandas()
    X = pd_df.drop(columns=['population'])
    y = pd_df['population']
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X) 
    pca = PCA(n_components=0.98)  
    X_pca = pca.fit_transform(X_scaled)
    selected_component_names = []

    for i in range(pca.n_components_):
        component_loadings = pca.components_[i]
        top_feature_index = np.argmax(np.abs(component_loadings))
        top_feature_name = X.columns[top_feature_index]
        selected_component_names.append(top_feature_name)
    
    newdf = pd_df.loc[:, ["population"]+selected_component_names]
    return newdf

def plot_heatmap(df: pd.DataFrame) -> None:
    snp_values = df.drop(columns=["population"]).T
    annotations = df["population"].values 
    plt.figure(figsize=(10,5))
    sns.heatmap(
        snp_values,
        fmt="s",
        cmap="coolwarm",
        cbar_kws = {"label": "REF Count"},
        linewidths= 0.5 
    )
    plt.xlabel("Populations")
    plt.ylabel("SNPs")
    plt.title("REF Count for SNPs for population")
    plt.xticks(ticks=[i+0.5 for i in range(len(annotations))], labels=annotations, rotation=45)
    plt.show()

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

def create_heatmap_with_annotations(df, annotation_col, output_file=None):
    """
    Creates a heatmap from a DataFrame with annotations and numeric values.

    Parameters:
        df (pd.DataFrame): The DataFrame containing the data.
                          One column should contain annotations, and the others numeric data.
        annotation_col (str): The name of the column containing annotations.
        output_file (str, optional): Path to save the figure. If None, the figure is not saved.

    Returns:
        None
    """
    # Separate annotations and numeric data
    annotations = df[annotation_col].values
    data = df.drop(columns=[annotation_col]).values

    # Create the figure and axis
    fig, ax = plt.subplots(figsize=(10, 6))

    # Generate the heatmap
    cax = ax.imshow(data, cmap="coolwarm", aspect="auto")

    # Add colorbar
    cbar = fig.colorbar(cax, ax=ax, label="Values")

    # Set axis labels and ticks
    ax.set_xticks(range(data.shape[1]))
    ax.set_xticklabels(df.columns.drop(annotation_col), rotation=45)
    ax.set_yticks(range(data.shape[0]))
    ax.set_yticklabels(annotations)
    ax.set_xlabel("Columns")
    ax.set_ylabel("Annotations")
    ax.set_title("Heatmap with Annotations")

    # Save the figure if an output path is provided
    if output_file:
        fig.savefig(output_file, dpi=300, bbox_inches="tight")
        print(f"Figure saved to {output_file}")

    # Show the plot
    plt.show()


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
    pops = ['CNXJ', 'CnOther', 'CnQTP', 'DGN', 'ISR', 'WE', 'EE', 'EB', 'WB', 'TRK', 'CYP']
    pops2samples = find_pops_from_bamlist(bamfile, pops)
    print("Reading VCF...")
    df = read_vcf("/gatk_modified/userdata/abertelli/drosophila-evolution/results/drosophila_evolution.bcftools_2R.vcf.gz")
    print(df.head())
    print(df.height)
    print("Read VCF!")
    print("Selecting only SNPs and subsampling")
    snps_df = snps_and_random_subsample(df)
    print(snps_df.head())
    print(snps_df.height)
    print("Selected only SNPs and subsampling!")
    print("Mapping")
    newdf = snps_map(snps_df, pops2samples)
    print(newdf.head())
    print(newdf.height)
    print("Mapped!")
    print("PCAing...")
    pca_df = perform_PCA(newdf)
    print(pca_df)
    print("PCAed!")
    create_heatmap_with_annotations(pca_df, "population", "heatmap.png")








    