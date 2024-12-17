import pandas as pd
import os
import shutil

f = pd.read_csv("/gatk_modified/userdata/abertelli/drosophila-evolution/data/rename_dgn.csv")
accs = f["ACCESSION"].to_list()
nms = f["NAME"].to_list()
nm2acc = {accs[i]: nms[i] for i in range(len(accs))} 
#fls = [(os.path.join("/gatk_modified/userdata/abertelli/drosophila-evolution/data/bamfiles/", nm2acc[acc] + ".bam"), os.path.join("/gatk_modified/userdata/abertelli/drosophila-evolution/data/mapping/dgn/", acc, acc+".dedup.bam")) for acc in accs if os.path.exists(os.path.join("/gatk_modified/userdata/abertelli/drosophila-evolution/data/mapping/dgn/", acc, acc+".dedup.bam"))] 
flsi = [(os.path.join("/gatk_modified/userdata/abertelli/drosophila-evolution/data/bamfiles/", nm2acc[acc] + ".bam.bai"), os.path.join("/gatk_modified/userdata/abertelli/drosophila-evolution/data/mapping/dgn/", acc, acc+".dedup.bam.bai")) for acc in accs if os.path.exists(os.path.join("/gatk_modified/userdata/abertelli/drosophila-evolution/data/mapping/dgn/", acc, acc+".dedup.bam.bai"))] 

#for fl in fls:
#    shutil.copy(fl[1],fl[0])

for fli in flsi:
    shutil.copy(fli[1],fli[0])