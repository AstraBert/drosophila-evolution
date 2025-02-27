import pandas as pd

f = open("POOL_PopData.txt.id","r")
lines = f.readlines()
idd = [int(l.split(" ")[0]) for l in lines]
f.close()
csv = pd.read_csv("selected_snps_noMAF.csv")
df = csv[csv["id"].isin(idd)]
arms = df["Chromosome"].to_list()
chrs = list(set(arms))
counts = {chr: [arms.count(chr), arms.count(chr)*100/len(arms)] for chr in chrs}
dfp = pd.DataFrame(counts)
print(dfp)
