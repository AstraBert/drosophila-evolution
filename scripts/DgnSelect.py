import pandas as pd

csv = pd.read_csv("/gatk_modified/userdata/abertelli/drosophila-evolution/data/DGN.tsv", sep="\t")

idxes = csv["Stock ID"].to_list()
srrs = csv["SRA Accession"].to_list()

c = 0
accs = [] 
for i in range(len(srrs)):
    if idxes[i].startswith("ZI1") or idxes[i].startswith("ZI"):
        accs.append(srrs[i])
        c+=1
        if c==40:
            break

f = open("/gatk_modified/userdata/abertelli/drosophila-evolution/data/download_dgn.txt", "w")
accs = [acc+"\n" for acc in accs]
f.writelines(accs)
f.close() 
