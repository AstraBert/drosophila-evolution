import os




if __name__ == "__main__":

    walking = os.walk(
        "/gatk_modified/userdata/abertelli/drosophila-evolution/results/groups/"
    )

    dct = {}

    for root, dirs, files in walking:
        if len(dirs) > 0 and len(files) >= 0 and "results" in dirs:
            dct.update({root: []})

    for k in dct:
        vcftoconcat = [
            os.path.join(k, f)
            for f in os.listdir(os.path.join(k))
            if os.path.isfile(os.path.join(k,f))
        ]

        dct[k] = vcftoconcat

        shellscript = open(
            f"/gatk_modified/userdata/abertelli/drosophila-evolution/shell/stats/{k.split('/')[-1]}.sh",
            "w",
        )
        shellscript.write(
            f'source activate python_deps\n'
        )
        shellscript.write(
            f'bcftools stats {dct[k][0]} > {dct[k][0].split(".")[0]}.vchk\n'
        )
        shellscript.write(
            f'plot-vcfstats -p {k}/plots {k}/{dct[k][0].split("/")[-1].split(".")[0]}.vchk\n'
        )
        shellscript.write(
            f'conda deactivate\n'
        )
        shellscript.close()

