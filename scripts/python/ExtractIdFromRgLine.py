import pysam
import os
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument(
    "-rgl",
    "--rgline",
    type=str,
    required=True,
    help="Input RG Line from whi to extract the ID",
)

args = parser.parse_args()
inf = args.rgline


def extract_id_from_rgline(rgline: str) -> str:
    idx = rgline.split("\t")[1]
    return idx


if __name__ == "__main__":
    idx = extract_id_from_rgline(inf)
    print(idx)
