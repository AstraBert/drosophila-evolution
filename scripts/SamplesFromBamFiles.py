import pysam
import sys
from argparse import ArgumentParser
import time

parser = ArgumentParser()
parser.add_argument("-i", "--inputfile", type=str, required=True, help="Input TXT file with a list of BAM files")
parser.add_argument("-o", "--outputfile", type=str, required=True, help="Output samples file")

args = parser.parse_args()
inf = args.inputfile
outf = args.outputfile


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


if __name__ == "__main__":
    strt = time.time()
    inff = open(inf, "r")
    lines = inff.readlines()
    bamfiles = [l.replace("\n","") for l in lines]
    outff = open(outf, "w")
    for bamfile in bamfiles: 
        samplename = extract_sample_from_rgline(bamfile)
        outff.write(samplename+"\n")
        print(f"Processed {bamfile}: found {samplename} as sample name", file=sys.stderr)
    outff.close()
    endt = time.time()
    print(f"Execution completed in {endt-strt} s", file=sys.stderr)

