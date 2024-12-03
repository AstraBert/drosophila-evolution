import os

if __name__ == "__main__":
    indseq_groups = {"CNXJ": [], "CnOther": [], "CnQTP": [], "DGN": [], "ISR": []} 
    poolseq_groups = {"WE": [], "EE": [], "EB": [], "WB": [], "TRK": [], "CYP": []} 
    regions = ["2L", "3L", "2R", "3R", "X"] 
    bamfiles = "/gatk_modified/userdata/abertelli/drosophila-evolution/data/freebayes_inputs/bamfiles.txt"
    bamfs = open(bamfiles, "r")
    bamlist = [bam.replace("\n","") for bam in bamfs.readlines()] 
    for k in indseq_groups:
        for bam in bamlist:
            bamname = bam.split("/")[-1].split(".")[0]
            if bamname.startswith(k):
                indseq_groups[k].append(bam)
    for k in poolseq_groups:
        for bam in bamlist:
            bamname = bam.split("/")[-1].split(".")[0]
            if bamname.startswith(k):
                poolseq_groups[k].append(bam)
    for k,v in indseq_groups.items():
        fll = f"/gatk_modified/userdata/abertelli/drosophila-evolution/data/bamlists/{k}.txt"
        f = open(fll, "w")
        for el in v:
            f.write(el+"\n")
        f.close()
        shellscript = open(f"/gatk_modified/userdata/abertelli/drosophila-evolution/shell/subsamples/{k}.sh", "w")
        for region in regions:
            shellscript.write(f"bash /gatk_modified/userdata/abertelli/drosophila-evolution/shell/subsample_vcf.sh -v /gatk_modified/userdata/abertelli/drosophila-evolution/results/drosophila_evolution.bcftools_{region}.vcf.gz -sgn {k} -t 10 -s {fll} -o /gatk_modified/userdata/abertelli/drosophila-evolution/results/groups/{k}\n")
        shellscript.close()
    for k,v in poolseq_groups.items():
        os.makedirs(f"/gatk_modified/userdata/abertelli/drosophila-evolution/data/bamlists/{k}")
        bamnames = [] 
        for el in v:
            bamname = el.split("/")[-1].split(".")[0] 
            fll = f"/gatk_modified/userdata/abertelli/drosophila-evolution/data/bamlists/{k}/{bamname}.txt"
            f = open(fll, "w")
            f.write(el+"\n")
            f.close()
            bamnames.append((fll, bamname))
        shellscript = open(f"/gatk_modified/userdata/abertelli/drosophila-evolution/shell/subsamples/{k}.sh", "w")
        for filess in bamnames:
            for region in regions:
                shellscript.write(f"bash /gatk_modified/userdata/abertelli/drosophila-evolution/shell/subsample_vcf.sh -v /gatk_modified/userdata/abertelli/drosophila-evolution/results/drosophila_evolution.bcftools_{region}.vcf.gz -sgn {filess[1]} -t 10 -s {filess[0]} -o /gatk_modified/userdata/abertelli/drosophila-evolution/results/groups/{k}/{filess[1]}\n")
        shellscript.close()