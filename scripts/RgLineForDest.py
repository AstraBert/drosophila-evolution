import pysam
import os
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-i", "--inputfile", type=str, required=True, help="Input BAM file")

args = parser.parse_args()
inf = args.inputfile


def extract_rgline(bamfile: str) -> str:
    bam = pysam.AlignmentFile(bamfile, "rb")
    head = bam.header
    strhead = str(head)
    lines = strhead.split("\n")
    rgline = [line for line in lines if line.startswith("@RG")]
    rgline = rgline[0]
    samplename = rgline.split("\t")[2]
    samplename = samplename.split(":")[0] + ":" + bamfile.split("/")[-1].split(".")[0]
    rgline_split = rgline.split("\t")
    rgline_split[2] = samplename
    rgline = ("\t").join(rgline_split)
    return rgline


if __name__ == "__main__":
    rgline = extract_rgline(inf)
    print(rgline)
