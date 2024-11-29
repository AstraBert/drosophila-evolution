import os
import shutil

ORIGIN_DIRS = [
    "/gatk_modified/userdata/abertelli/drosophila-evolution/data/mapping/CN_XJ",
    "/gatk_modified/userdata/abertelli/drosophila-evolution/data/mapping/CnQTP",
    "/gatk_modified/userdata/abertelli/drosophila-evolution/data/mapping/ISR_1",
    "/gatk_modified/userdata/abertelli/drosophila-evolution/data/mapping/ISR_2",
]

DESTINATION_DIR = "/gatk_modified/userdata/abertelli/drosophila-evolution/data/bamfiles"

for d in ORIGIN_DIRS:
    flnames = [f for f in os.listdir(d) if os.path.isdir(os.path.join(d, f))]
    filess = [
        os.path.join(DESTINATION_DIR, f)
        for f in os.listdir(DESTINATION_DIR)
        if f.split(".")[0] in flnames
    ]
    c = 0
    for j in range(len(filess)):
        if filess[j].endswith(".bai"):
            continue
        else:
            c += 1
            newflnm = d.split("/")[-1].replace("_", "") + "_" + str(c) + ".bam"
            newflnm1 = d.split("/")[-1].replace("_", "") + "_" + str(c) + ".bam.bai"
        newflpath = os.path.join(DESTINATION_DIR, newflnm)
        newflpath1 = os.path.join(DESTINATION_DIR, newflnm1)
        try:
            shutil.move(filess[j], newflpath)
            shutil.move(f"{filess[j]+'.bai'}", newflpath1)
        except Exception as e:
            print(e)
