import os


def reorder_files(files):
    # Define the custom order
    custom_order = ["2R", "2L", "3R", "3L", "X"]
    # Sort the files based on the chromosome suffix using the custom order
    sorted_files = sorted(
        files,
        key=lambda x: custom_order.index(x.split("/")[-1].split(".")[1].split("_")[-1]),
    )
    return sorted_files


if __name__ == "__main__":

    walking = os.walk(
        "/gatk_modified/userdata/abertelli/drosophila-evolution/results/groups/"
    )

    dct = {}

    for root, dirs, files in walking:
        if len(dirs) > 0 and len(files) == 0 and "results" in dirs:
            dct.update({root: []})

    for k in dct:
        vcftoconcat = [
            os.path.join(k, "results", f)
            for f in os.listdir(os.path.join(k, "results"))
            if os.path.isfile(os.path.join(k, "results", f))
        ]
        dct[k] = reorder_files(vcftoconcat)

        shellscript = open(
            f"/gatk_modified/userdata/abertelli/drosophila-evolution/shell/concatenate/{k.split('/')[-1]}.sh",
            "w",
        )
        shellscript.write(
            f"bcftools concat -O z --threads 10 -o {k}/{k.split('/')[-1]}_all.bcftools.vcf.gz {' '.join(dct[k])}\n"
        )
        shellscript.close()
