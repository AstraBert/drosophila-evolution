import pandas as pd

csv = open("divergence_all.tsv", "r+")
lines = csv.readlines()
lines = ["\t".join(line.split(" ")) for line in lines]
lines =  ["Pools\t"+lines[0].replace("bjack\tmean","bjack_mean").replace("bjack\ts.e.","bjack_s.e.")] + lines[1:]
csv.seek(0)
csv.truncate()
csv.writelines(lines)
csv.close()

df = pd.read_csv("divergence_all.tsv", sep="\t")
pools_df = pd.read_csv("pools.csv")
names = pools_df["NAME"].to_list()
pools = pools_df["POOL"].to_list()
pools2names = {pools[i]: names[i] for i in range(len(pools))}
f3pops = df["Pools"].to_list() 
focals = [pools2names[f3pop.split(";")[0]] for f3pop in f3pops]  
pop1s = [pools2names[f3pop.split(";")[1].split(",")[1]] for f3pop in f3pops]  
pop2s = [pools2names[f3pop.split(";")[1].split(",")[0]] for f3pop in f3pops] 
pops1 =[pools2names[f3pop.split(",")[0]] for f3pop in f3pops]  
pops2 =[pools2names[f3pop.split(",")[1]] for f3pop in f3pops]  

df.insert(1, "Pop2", pop2s)
df.insert(1, "Pop1", pop1s)
df.insert(1, "Focal", focals)
df.to_csv("f3stats.tsv", sep="\t", index=False)

