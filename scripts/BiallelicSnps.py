import polars as pl  
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-i","--infile",help="Input VCF file", type=str)
args = parser.parse_args()

inf = args.infile

def read_vcf(vcffile: str) -> pl.DataFrame:
    df = pl.read_csv(vcffile, separator="\t", n_threads=22, skip_rows=1901)
    return df

def snps(df: pl.DataFrame) -> pl.DataFrame:
    snpsdf = df.filter((pl.col("REF").is_in(["A","T","C","G"])) & (pl.col("ALT").is_in(["A","T","C","G"])))
    return snpsdf


if __name__ == "__main__":
    df = read_vcf(inf)
    print(f"Total Variants for: {inf} ", df.height)
    snpsdf = snps(df)
    print(f"Biallelic SNPsfor: {inf} ", snpsdf.height)