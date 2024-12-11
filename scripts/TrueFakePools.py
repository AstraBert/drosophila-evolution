f = open("data/freebayes_inputs/bamfiles_1.txt")
lines = f.readlines()
bamfiles = [line.replace("\n","") for line in lines]
group_pops = {"DGN": [], "CNXJ": [], "CnOther": [], "CnQTP": [], "ISR": []}
for k in group_pops:
    for i in range(len(bamfiles)):
        if k in bamfiles[i]: 
            group_pops[k].append(str(i+1))    

f.close()

w = open("scripts/FakePools.r", "w")

for k,v in group_pops.items():
    w.write(f"fake.pools.{k} = c({', '.join(v)})\n")
