import pandas as pd

df = pd.read_csv("selected_snps.csv")
arms = df["Chromosome"].to_list()
chrs = list(set(arms))
counts = {chr: [arms.count(chr), arms.count(chr)*100/len(arms)] for chr in chrs}
dfp = pd.DataFrame(counts)
print(dfp)
