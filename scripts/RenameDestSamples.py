import os
import shutil

DESTINATION_DIR = "/gatk_modified/drosophila-project/data"
ls = os.listdir(DESTINATION_DIR)
fls = [f for f in os.listdir(DESTINATION_DIR) if os.path.isfile(os.path.join(DESTINATION_DIR,f))]
dictionnaire = {"ES": ["WEES",0], "FR": ["WEFR",0], "RU": ["EERU",0], "UA": ["EEUA",0]}
for fl in fls:
    if fl.split("_")[0] in dictionnaire:
        if fl.endswith(".bai"):
            continue
        else:
            dictionnaire[fl.split("_")[0]][1]+=1
            newflname = dictionnaire[fl.split("_")[0]][0]+"_"+str(dictionnaire[fl.split("_")[0]][1])+".bam"
            newflname1 = dictionnaire[fl.split("_")[0]][0]+"_"+str(dictionnaire[fl.split("_")[0]][1])+".bam.bai"
            newpath = os.path.join(DESTINATION_DIR,newflname)
            newpath1 = os.path.join(DESTINATION_DIR,newflname1)
            shutil.move(os.path.join(DESTINATION_DIR,fl),newpath)
            shutil.move(f"{os.path.join(DESTINATION_DIR,fl)}.bai",newpath1)

